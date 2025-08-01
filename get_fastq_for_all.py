#!/usr/bin/env python3
"""
FASTQ File Splitter

This script splits a large FASTQ/FASTQ.gz file into smaller chunks based on 
predefined read ID sets stored in pickle files. It supports three different 
processing modes:

1. Single-threaded mode (--threads 1)
2. Multi-threaded mode (--threads > 1) 
3. seqkit-based mode (--use_seqkit)

The seqkit mode leverages the high-performance seqkit tool for faster processing.
To use seqkit mode, you must have seqkit installed and available in your PATH.

Installation of seqkit:
- Download from: https://github.com/shenwei356/seqkit/releases
- Or install via conda: conda install -c bioconda seqkit
- Or install via go: go install github.com/shenwei356/seqkit@latest

Usage examples:
    # Using built-in single-threaded mode
    python get_fastq_for_all.py -i input.fastq.gz -c chunks.txt -p pickle_dir -o output_dir --pickle_suffix="_read_ids.pkl" --threads 1
    
    # Using built-in multi-threaded mode  
    python get_fastq_for_all.py -i input.fastq.gz -c chunks.txt -p pickle_dir -o output_dir --pickle_suffix="_read_ids.pkl" --threads 8
    
    # Using seqkit mode (recommended for best performance)
    python get_fastq_for_all.py -i input.fastq.gz -c chunks.txt -p pickle_dir -o output_dir --pickle_suffix="_read_ids.pkl" --use_seqkit
"""

import argparse
import gzip
import os
import pickle
from contextlib import ExitStack
import threading
import queue
from collections import defaultdict, deque
import time
import subprocess
import tempfile

def parse_chunk_list(chunk_list_path):
    """
    Parses the chunk list file.

    Args:
        chunk_list_path (str): Path to the chunk list file.
            Expected format: one line per chunk with 'chrom<tab>chunk_idx'.

    Returns:
        list: A list of tuples, where each tuple is (chrom, chunk_idx_str).
              chunk_idx_str is kept as a string for consistent file naming.
    """
    chunks = []
    print(f"Parsing chunk list file: {chunk_list_path}")
    try:
        with open(chunk_list_path, 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith('#'): # Skip empty lines/comments
                    continue
                parts = line.split()
                if len(parts) == 2:
                    chrom, chunk_idx = parts
                    chunks.append((chrom, chunk_idx))
                else:
                    print(f"Warning: Skipping malformed line {i+1} in chunk list: '{line}'")
    except FileNotFoundError:
        print(f"Error: Chunk list file not found at {chunk_list_path}")
        raise
    except Exception as e:
        print(f"Error reading chunk list file {chunk_list_path}: {e}")
        raise
    print(f"Found {len(chunks)} chunks to process.")
    return chunks

def load_read_id_sets(pickle_dir_path, chunks_to_process, pickle_suffix):
    """
    Loads read ID sets from pickle files for the specified chunks
    and creates a reverse map from read_id to chunk identifier.

    Args:
        pickle_dir_path (str): Base directory containing chromosome subdirectories
                               which in turn contain the pickle files.
        chunks_to_process (list): List of (chrom, chunk_idx_str) tuples.
        pickle_suffix (str): The expected suffix for the pickle filenames
                             (e.g., '_read_ids.pkl', '_2_read_ids.pkl').

    Returns:
        dict: A dictionary mapping read IDs (str) to their corresponding
              chunk identifier tuple (chrom, chunk_idx_str).
              Returns None if errors occur during loading.
    """
    read_id_to_chunk_map = {}
    print(f"Loading read ID sets from base directory: {pickle_dir_path}")
    print(f"Expecting pickle file suffix: {pickle_suffix}")
    loaded_count = 0
    skipped_count = 0
    total_read_ids = 0

    for chrom, chunk_idx in chunks_to_process:
        # Construct the filename using the provided suffix
        pickle_filename = f"{chrom}_chunk{chunk_idx}{pickle_suffix}"
        # Construct the full path including the chromosome subdirectory
        pickle_path = os.path.join(pickle_dir_path, chrom, pickle_filename)
        chunk_id = (chrom, chunk_idx) # Use tuple as key/value

        try:
            # Check if the specific pickle file exists before trying to open
            if not os.path.exists(pickle_path):
                 print(f"Warning: Pickle file not found, skipping: {pickle_path}")
                 skipped_count += 1
                 continue

            with open(pickle_path, 'rb') as pf:
                read_ids_set = pickle.load(pf)
                if not isinstance(read_ids_set, set):
                    print(f"Warning: Pickle file {pickle_filename} in {chrom} subdir does not contain a set. Skipping.")
                    skipped_count += 1
                    continue

                for read_id in read_ids_set:
                    # If a read ID is in multiple chunks, this assigns it to the *last* chunk processed containing it.
                    # This behavior is consistent with the original script.
                    read_id_to_chunk_map[read_id] = chunk_id
                loaded_count += 1
                total_read_ids += len(read_ids_set)

        # Removed FileNotFoundError check here as we check with os.path.exists now
        except pickle.UnpicklingError:
            print(f"Error: Could not unpickle file, skipping: {pickle_path}")
            skipped_count += 1
        except Exception as e:
            print(f"Error loading pickle file {pickle_path}: {e}")
            skipped_count += 1

    print(f"Successfully loaded {loaded_count} pickle files.")
    if skipped_count > 0:
        print(f"Skipped {skipped_count} pickle files due to errors or missing files.")
    print(f"Total unique read IDs across loaded sets: {len(read_id_to_chunk_map)}")
    if len(read_id_to_chunk_map) != total_read_ids:
         print(f"Note: Some read IDs might be present in multiple input sets ({total_read_ids} total entries vs {len(read_id_to_chunk_map)} unique). Each read ID is assigned to the last output chunk encountered.")

    if not read_id_to_chunk_map and loaded_count == 0: # Check if *any* file was loaded
        print("Error: No read ID sets were successfully loaded.")
        return None

    return read_id_to_chunk_map

class FastqRecord:
    """Class to represent a FASTQ record"""
    __slots__ = ['header', 'seq', 'plus', 'qual']
    
    def __init__(self, header, seq, plus, qual):
        self.header = header
        self.seq = seq
        self.plus = plus
        self.qual = qual
    
    def get_read_id(self):
        return self.header.split()[0][1:]
    
    def write_to_file(self, file_handle):
        file_handle.write(self.header)
        file_handle.write(self.seq)
        file_handle.write(self.plus)
        file_handle.write(self.qual)

def fastq_reader_worker(input_fastq_path, record_queue, batch_size, stop_event):
    """
    Worker function to read FASTQ records and put them into a queue.
    """
    # Determine if input is gzipped
    is_gzipped = input_fastq_path.lower().endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' # Read in text mode
    
    try:
        with open_func(input_fastq_path, mode) as infile:
            batch = []
            while not stop_event.is_set():
                header = infile.readline()
                if not header: # End of file
                    break
                seq = infile.readline()
                plus = infile.readline()
                qual = infile.readline()

                if not (header.startswith('@') and plus.startswith('+')):
                    continue  # Skip malformed records

                record = FastqRecord(header, seq, plus, qual)
                batch.append(record)
                
                if len(batch) >= batch_size:
                    try:
                        record_queue.put(batch, timeout=1)
                        batch = []
                    except queue.Full:
                        if stop_event.is_set():
                            break
                        continue
                        
            # Put remaining records
            if batch:
                try:
                    record_queue.put(batch, timeout=5)
                except queue.Full:
                    pass
                    
    except Exception as e:
        print(f"Error in reader worker: {e}")
    finally:
        # Signal end of reading
        try:
            record_queue.put(None, timeout=5)
        except queue.Full:
            pass

def processing_worker(record_queue, write_queue, read_id_to_chunk_map, worker_id, stats, stop_event):
    """
    Worker function to process FASTQ records and route them to appropriate chunks.
    """
    processed_count = 0
    
    try:
        while not stop_event.is_set():
            try:
                batch = record_queue.get(timeout=1)
                if batch is None:  # End signal
                    break
                    
                # Group records by chunk
                chunk_records = defaultdict(list)
                
                for record in batch:
                    read_id = record.get_read_id()
                    target_chunk_id = read_id_to_chunk_map.get(read_id)
                    
                    if target_chunk_id:
                        chunk_records[target_chunk_id].append(record)
                        stats['written'] += 1
                    else:
                        stats['ignored'] += 1
                    
                    processed_count += 1
                    stats['processed'] += 1
                
                # Send grouped records to write queue
                for chunk_id, records in chunk_records.items():
                    try:
                        write_queue.put((chunk_id, records), timeout=1)
                    except queue.Full:
                        if stop_event.is_set():
                            break
                        continue
                        
            except queue.Empty:
                continue
            except Exception as e:
                print(f"Error in processing worker {worker_id}: {e}")
                break
                
    except Exception as e:
        print(f"Error in processing worker {worker_id}: {e}")

def writing_worker(write_queue, output_files, buffer_size, stats, stop_event):
    """
    Worker function to write records to output files with buffering.
    """
    # Buffers for each output file
    buffers = {chunk_id: deque() for chunk_id in output_files.keys()}
    
    try:
        while not stop_event.is_set():
            try:
                data = write_queue.get(timeout=1)
                if data is None:  # End signal
                    break
                    
                chunk_id, records = data
                
                if chunk_id in buffers:
                    buffers[chunk_id].extend(records)
                    
                    # Flush buffer if it's getting large
                    if len(buffers[chunk_id]) >= buffer_size:
                        flush_buffer(chunk_id, buffers[chunk_id], output_files, stats)
                        
            except queue.Empty:
                continue
            except Exception as e:
                print(f"Error in writing worker: {e}")
                break
                
        # Flush all remaining buffers
        for chunk_id, buffer in buffers.items():
            if buffer:
                flush_buffer(chunk_id, buffer, output_files, stats)
                
    except Exception as e:
        print(f"Error in writing worker: {e}")

def flush_buffer(chunk_id, buffer, output_files, stats):
    """Flush a buffer to the corresponding output file."""
    if chunk_id in output_files and buffer:
        outfile = output_files[chunk_id]
        for record in buffer:
            record.write_to_file(outfile)
        stats['chunk_counts'][chunk_id] += len(buffer)
        buffer.clear()

def split_fastq_multithreaded(input_fastq_path, read_id_to_chunk_map, output_dir_path, 
                             chunks_to_process, num_workers=4, batch_size=1000, 
                             buffer_size=5000, queue_size=100):
    """
    Reads the input FASTQ file and writes records to corresponding chunk files using multiple threads.

    Args:
        input_fastq_path (str): Path to the input FASTQ or FASTQ.gz file.
        read_id_to_chunk_map (dict): Map from read ID to (chrom, chunk_idx) tuple.
        output_dir_path (str): Directory to write the output chunk FASTQ files.
        chunks_to_process (list): List of (chrom, chunk_idx_str) tuples used
                                  to define output files.
        num_workers (int): Number of processing worker threads.
        batch_size (int): Number of records to read in each batch.
        buffer_size (int): Number of records to buffer before writing.
        queue_size (int): Maximum size of queues.
    """
    print(f"Starting multithreaded FASTQ splitting process for: {input_fastq_path}")
    print(f"Using {num_workers} worker threads, batch size: {batch_size}, buffer size: {buffer_size}")
    
    os.makedirs(output_dir_path, exist_ok=True)
    
    # Thread-safe statistics
    manager = threading.Lock()
    stats = {
        'processed': 0,
        'written': 0,
        'ignored': 0,
        'chunk_counts': defaultdict(int)
    }
    
    # Create queues
    record_queue = queue.Queue(maxsize=queue_size)
    write_queue = queue.Queue(maxsize=queue_size)
    
    # Event to signal workers to stop
    stop_event = threading.Event()
    
    # Use ExitStack to manage multiple output file handles cleanly
    with ExitStack() as stack:
        # Open output files
        output_files = {}
        print(f"Preparing output files in: {output_dir_path}")
        for chunk_id in chunks_to_process:
            chrom, chunk_idx = chunk_id
            output_filename = f"{chrom}_chunk{chunk_idx}.fastq"
            output_path = os.path.join(output_dir_path, output_filename)
            try:
                outfile = stack.enter_context(open(output_path, 'w'))
                output_files[chunk_id] = outfile
            except IOError as e:
                print(f"Error: Could not open output file {output_path} for writing: {e}")
                raise

        if not output_files:
            print("Warning: No output files were opened (perhaps no valid chunks?). Exiting.")
            return

        print(f"{len(output_files)} output files opened successfully.")
        
        # Start worker threads
        threads = []
        
        # Start reader thread
        reader_thread = threading.Thread(
            target=fastq_reader_worker,
            args=(input_fastq_path, record_queue, batch_size, stop_event)
        )
        reader_thread.start()
        threads.append(reader_thread)
        
        # Start processing workers
        for i in range(num_workers):
            worker_thread = threading.Thread(
                target=processing_worker,
                args=(record_queue, write_queue, read_id_to_chunk_map, i, stats, stop_event)
            )
            worker_thread.start()
            threads.append(worker_thread)
        
        # Start writing worker
        writer_thread = threading.Thread(
            target=writing_worker,
            args=(write_queue, output_files, buffer_size, stats, stop_event)
        )
        writer_thread.start()
        threads.append(writer_thread)
        
        # Monitor progress
        start_time = time.time()
        last_count = 0
        
        try:
            while reader_thread.is_alive() or not record_queue.empty() or not write_queue.empty():
                time.sleep(5)  # Report progress every 5 seconds
                current_count = stats['processed']
                if current_count > last_count:
                    elapsed = time.time() - start_time
                    rate = current_count / elapsed if elapsed > 0 else 0
                    print(f"  Processed {current_count//1000}k records ({rate:.0f} records/sec)")
                    last_count = current_count
                    
        except KeyboardInterrupt:
            print("\nInterrupted by user. Stopping workers...")
            stop_event.set()
        
        # Signal end of processing
        try:
            write_queue.put(None, timeout=5)
        except queue.Full:
            pass
        
        # Wait for all threads to finish
        for thread in threads:
            thread.join(timeout=30)
        
        if stop_event.is_set():
            print("Some threads may not have finished cleanly.")

    print("\nFASTQ splitting finished.")
    print(f"Total records processed: {stats['processed']}")
    print(f"Records written to chunks: {stats['written']}")
    print(f"Records ignored (not in any target set): {stats['ignored']}")

    print("\nReads written per chunk:")
    for chunk_id, count in stats['chunk_counts'].items():
        if count > 0: # Only print chunks that received reads
            print(f"  {chunk_id[0]}_chunk{chunk_id[1]}: {count} reads")

# Legacy single-threaded function for compatibility
def split_fastq(input_fastq_path, read_id_to_chunk_map, output_dir_path, chunks_to_process):
    """
    Reads the input FASTQ file and writes records to corresponding chunk files.
    This is the original single-threaded implementation.
    """
    print(f"Starting single-threaded FASTQ splitting process for: {input_fastq_path}")
    os.makedirs(output_dir_path, exist_ok=True)

    # Determine if input is gzipped
    is_gzipped = input_fastq_path.lower().endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' # Read in text mode

    processed_records = 0
    written_records = 0
    ignored_records = 0
    chunk_counts = {chunk_id: 0 for chunk_id in chunks_to_process} # Track counts per chunk

    # Use ExitStack to manage multiple output file handles cleanly
    with ExitStack() as stack:
        # Open input file
        infile = stack.enter_context(open_func(input_fastq_path, mode))
        print("Input FASTQ opened successfully.")

        # Open output files
        output_files = {}
        print(f"Preparing output files in: {output_dir_path}")
        for chunk_id in chunks_to_process:
            chrom, chunk_idx = chunk_id
            output_filename = f"{chrom}_chunk{chunk_idx}.fastq"
            output_path = os.path.join(output_dir_path, output_filename)
            try:
                outfile = stack.enter_context(open(output_path, 'w'))
                output_files[chunk_id] = outfile
            except IOError as e:
                print(f"Error: Could not open output file {output_path} for writing: {e}")
                raise # Reraise error to stop execution

        if not output_files:
             print("Warning: No output files were opened (perhaps no valid chunks?). Exiting.")
             return

        print(f"{len(output_files)} output files opened successfully.")
        print("Processing records...")

        # Iterate through FASTQ records (4 lines at a time)
        while True:
            try:
                header = infile.readline()
                if not header: # End of file
                    break
                seq = infile.readline()
                plus = infile.readline()
                qual = infile.readline()

                if not (header.startswith('@') and plus.startswith('+')):
                     print(f"Warning: Malformed FASTQ record found near record {processed_records + 1}. Attempting to continue.")
                     continue 

                processed_records += 1

                read_id = header.split()[0][1:]

                target_chunk_id = read_id_to_chunk_map.get(read_id)

                if target_chunk_id:
                    if target_chunk_id in output_files:
                        outfile = output_files[target_chunk_id]
                        outfile.write(header)
                        outfile.write(seq)
                        outfile.write(plus)
                        outfile.write(qual)
                        written_records += 1
                        chunk_counts[target_chunk_id] += 1
                    else:
                         print(f"Warning: Read ID '{read_id}' mapped to chunk {target_chunk_id}, but no output file is open for it. Ignoring.")
                         ignored_records += 1
                else:
                    ignored_records += 1

                if processed_records % 100000 == 0:
                    print(f"  Processed {processed_records // 1000}k records...")

            except Exception as e:
                print(f"Error processing record {processed_records + 1}: {e}")
                raise # Reraise to stop execution on error


    print("\nFASTQ splitting finished.")
    print(f"Total records processed: {processed_records}")
    print(f"Records written to chunks: {written_records}")
    print(f"Records ignored (not in any target set): {ignored_records}")

    print("\nReads written per chunk:")
    for chunk_id, count in chunk_counts.items():
         if count > 0: # Only print chunks that received reads
             print(f"  {chunk_id[0]}_chunk{chunk_id[1]}: {count} reads")

def split_fastq_with_seqkit(input_fastq_path, read_id_to_chunk_map, output_dir_path, chunks_to_process):
    """
    Uses seqkit to split the input FASTQ file into chunk files based on read IDs.
    This implementation leverages seqkit's performance for faster processing.
    
    Args:
        input_fastq_path (str): Path to the input FASTQ or FASTQ.gz file.
        read_id_to_chunk_map (dict): Map from read ID to (chrom, chunk_idx) tuple.
        output_dir_path (str): Directory to write the output chunk FASTQ files.
        chunks_to_process (list): List of (chrom, chunk_idx_str) tuples used
                                  to define output files.
    """
    print(f"Starting seqkit-based FASTQ splitting process for: {input_fastq_path}")
    os.makedirs(output_dir_path, exist_ok=True)
    
    # Check if seqkit is available
    try:
        result = subprocess.run(['seqkit', 'version'], capture_output=True, text=True, check=True)
        print(f"Using seqkit version: {result.stdout.strip()}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: seqkit is not installed or not in PATH. Please install seqkit first.")
        print("You can install it from: https://github.com/shenwei356/seqkit")
        raise RuntimeError("seqkit not available")
    
    # Group read IDs by chunk
    chunk_read_ids = defaultdict(set)
    for read_id, chunk_id in read_id_to_chunk_map.items():
        if chunk_id in chunks_to_process:
            chunk_read_ids[chunk_id].add(read_id)
    
    print(f"Preparing to split into {len(chunk_read_ids)} chunks")
    
    # Create temporary directory for read ID files
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_files = {}
        
        # Create read ID files for each chunk
        print("Creating temporary read ID files...")
        for chunk_id in chunks_to_process:
            if chunk_id in chunk_read_ids and chunk_read_ids[chunk_id]:
                chrom, chunk_idx = chunk_id
                temp_filename = f"{chrom}_chunk{chunk_idx}_ids.txt"
                temp_path = os.path.join(temp_dir, temp_filename)
                
                with open(temp_path, 'w') as f:
                    for read_id in chunk_read_ids[chunk_id]:
                        f.write(f"{read_id}\n")
                
                temp_files[chunk_id] = temp_path
                print(f"  Created {temp_filename} with {len(chunk_read_ids[chunk_id])} read IDs")
        
        if not temp_files:
            print("Warning: No read ID files created. No reads to process.")
            return
        
        # Process each chunk with seqkit
        print("\nExtracting reads with seqkit...")
        chunk_counts = {}
        
        for chunk_id, temp_file in temp_files.items():
            chrom, chunk_idx = chunk_id
            output_filename = f"{chrom}_chunk{chunk_idx}.fastq"
            output_path = os.path.join(output_dir_path, output_filename)
            
            print(f"  Processing chunk {chrom}_chunk{chunk_idx}...")
            
            try:
                # Use seqkit grep to extract reads
                # -f: file containing patterns (read IDs)
                # -o: output file
                # --quiet: suppress warnings about missing reads
                cmd = [
                    'seqkit', 'grep',
                    '-f', temp_file,
                    '-o', output_path,
                    '--quiet',
                    input_fastq_path
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                
                # Count reads in output file
                count_cmd = ['seqkit', 'stats', '-T', output_path]
                count_result = subprocess.run(count_cmd, capture_output=True, text=True, check=True)
                
                # Parse the count from seqkit stats output (tab-separated, second column is num_seqs)
                lines = count_result.stdout.strip().split('\n')
                if len(lines) >= 2:
                    count = int(lines[1].split('\t')[3])  # num_seqs is the 4th column (0-indexed: 3)
                    chunk_counts[chunk_id] = count
                    print(f"    Extracted {count} reads to {output_filename}")
                else:
                    chunk_counts[chunk_id] = 0
                    print(f"    Warning: Could not determine read count for {output_filename}")
                
            except subprocess.CalledProcessError as e:
                print(f"    Error processing chunk {chrom}_chunk{chunk_idx}: {e}")
                print(f"    seqkit stderr: {e.stderr}")
                chunk_counts[chunk_id] = 0
            except Exception as e:
                print(f"    Unexpected error processing chunk {chrom}_chunk{chunk_idx}: {e}")
                chunk_counts[chunk_id] = 0
    
    # Summary statistics
    total_written = sum(chunk_counts.values())
    total_expected = sum(len(read_ids) for read_ids in chunk_read_ids.values())
    
    print("\nFASTQ splitting with seqkit finished.")
    print(f"Total reads expected: {total_expected}")
    print(f"Total reads written: {total_written}")
    
    if total_written < total_expected:
        print(f"Note: {total_expected - total_written} reads were not found in the input file")
    
    print("\nReads written per chunk:")
    for chunk_id, count in chunk_counts.items():
        if count > 0:
            print(f"  {chunk_id[0]}_chunk{chunk_id[1]}: {count} reads")

def main():
    parser = argparse.ArgumentParser(
        description="Split a large FASTQ/FASTQ.gz file into smaller chunks based on predefined read ID sets stored in pickle files."
    )
    parser.add_argument(
        "-i", "--input_fastq",
        required=True,
        help="Path to the input FASTQ file (can be .fastq or .fastq.gz)."
    )
    parser.add_argument(
        "-c", "--chunk_list",
        required=True,
        help="Path to the text file listing chunks to process (format: 'chrom<tab>chunk_idx' per line)."
    )
    parser.add_argument(
        "-p", "--pickle_dir",
        required=True,
        help="Base directory containing chromosome subdirectories with pickled read ID sets (e.g., {pickle_dir}/{chrom}/{chrom}_chunk{chunk_idx}_suffix.pkl)."
    )
    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="Directory where the output FASTQ chunk files will be written."
    )
    parser.add_argument(
        "--pickle_suffix",
        required=True,
        help="Suffix for the pickle filenames, including the extension (e.g., '_read_ids.pkl', '_2_read_ids.pkl')."
    )
    
    # Performance optimization parameters
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of worker threads for processing (default: 4). Set to 1 for single-threaded mode."
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=1000,
        help="Number of FASTQ records to process in each batch (default: 1000)."
    )
    parser.add_argument(
        "--buffer_size",
        type=int,
        default=5000,
        help="Number of records to buffer before writing to files (default: 5000)."
    )
    parser.add_argument(
        "--queue_size",
        type=int,
        default=100,
        help="Maximum size of internal queues (default: 100)."
    )
    parser.add_argument(
        "--use_seqkit",
        action="store_true",
        help="Use seqkit for FASTQ splitting instead of the built-in implementation. Requires seqkit to be installed."
    )

    args = parser.parse_args()

    try:
        # 1. Parse chunk list
        chunks_to_process = parse_chunk_list(args.chunk_list)
        if not chunks_to_process:
            print("No valid chunks found in the chunk list. Exiting.")
            return

        # 2. Load read ID sets into a reverse map using the suffix
        read_id_to_chunk_map = load_read_id_sets(args.pickle_dir, chunks_to_process, args.pickle_suffix)
        if read_id_to_chunk_map is None:
             print("Failed to load read ID sets. Exiting.")
             return # Exit if loading failed critically

        # 3. Process FASTQ and split
        if args.use_seqkit:
            time_start = time.time()
            split_fastq_with_seqkit(
                args.input_fastq, 
                read_id_to_chunk_map, 
                args.output_dir, 
                chunks_to_process
            )
            time_end = time.time()
            print(f"Time taken: {time_end - time_start} seconds")
        elif args.threads > 1:
            split_fastq_multithreaded(
                args.input_fastq, 
                read_id_to_chunk_map, 
                args.output_dir, 
                chunks_to_process,
                num_workers=args.threads,
                batch_size=args.batch_size,
                buffer_size=args.buffer_size,
                queue_size=args.queue_size
            )
        else:
            split_fastq(args.input_fastq, read_id_to_chunk_map, args.output_dir, chunks_to_process)

        print("\nScript finished successfully.")

    except FileNotFoundError as e:
        print(f"\nError: Input file not found. {e}")
    except IOError as e:
         print(f"\nError: An I/O error occurred. {e}")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == "__main__":
    main()