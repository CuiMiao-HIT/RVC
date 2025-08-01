import os
import shlex
import argparse
from collections import Counter, defaultdict
from subprocess import Popen, PIPE
from typing import Tuple, Dict, Set # Added Set to typing hints
import pickle
from multiprocessing import Pool
import logging
import re

def subprocess_popen(args,shell=False, stdin=None, stdout=PIPE, bufsize=8388608):
    return Popen(args, shell=shell, stdin=stdin, stdout=stdout, bufsize=bufsize, universal_newlines=True)

def get_chrom_length(bam_file, chrom, samtools_execute_command="samtools"):
    cmd = f"{samtools_execute_command} idxstats {bam_file}"
    process = subprocess_popen(shlex.split(cmd), shell=False)
    for line in process.stdout:
        chr_name, length, *_ = line.strip().split('\t')
        if chr_name == chrom:
            process.terminate()
            return int(length)
    process.stdout.close() # Ensure stdout is closed
    process.wait() # Wait for process to terminate
    return None

# Modified function to return read_info dict and read_id set
def save_read_info(bam_file, output_dir, chrom, chunk_start, chunk_end, samtools_execute_command="samtools") -> Tuple[Dict, Set[str]]:
    os.makedirs(output_dir, exist_ok=True)

    cmd = f"{samtools_execute_command} view --excl-flags 2316 {bam_file} {chrom}:{chunk_start}-{chunk_end}"
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file {bam_file} does not exist")

    process = subprocess_popen(shlex.split(cmd), shell=False)
    # Check process immediately after creation (basic check)
    if process.poll() is not None and process.returncode != 0:
         # Attempt to read stderr if possible, though Popen setup might not capture it easily here
        raise RuntimeError(f"samtools view command failed to start or exited immediately with code {process.returncode}")

    read_info = {}
    read_id_set = set() # Initialize the set for read IDs

    try:
        for line in process.stdout:
            fields = line.strip().split('\t')
            if len(fields) < 11: # Basic check for expected number of fields
                logging.warning(f"Skipping malformed SAM line in chunk {chrom}:{chunk_start}-{chunk_end}: {line.strip()}")
                continue
            read_id = fields[0]
            read_id_set.add(read_id) # Add read_id to the set

            ref_start_pos = int(fields[3])
            seq_length = len(fields[9])
            strand_info = "+" if int(fields[1]) & 16 == 0 else "-"
            cigar_string = fields[5]
            cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')
            cigar_tuples = cigar_pattern.findall(cigar_string)

            has_start_soft_clip = False
            has_end_soft_clip = False
            read_start_pos = 0
            read_end_pos = seq_length - 1 # Initial assumption, might change
            ref_end_pos = ref_start_pos - 1
            end_soft_seq = ""
            start_soft_seq = ""
            read_pos = 0
            ref_pos = ref_start_pos

            for count, operation in cigar_tuples:
                count = int(count)

                if operation == 'S':
                    if read_pos == 0:
                        has_start_soft_clip = True
                        read_start_pos = count # Start read position after soft clipping
                        start_soft_seq = fields[9][:count]
                        # Don't advance ref_pos for initial soft clips
                    # Check if this soft clip is at the end
                    # Need to calculate total read length consumed by non-S/H operations first
                    # Simpler approach: check if *remaining* operations before this S are all S/H
                    elif all(op in 'SH' for _, op in cigar_tuples[cigar_tuples.index((str(count), operation))+1:]):
                         # This assumes the current implementation correctly calculates read_pos up to this point
                        if read_pos + count == seq_length:
                            has_end_soft_clip = True
                            end_soft_seq = fields[9][-count:]
                            # Don't advance ref_pos for terminal soft clips
                    read_pos += count # Always advance read_pos for S
                elif operation == 'M' or operation == '=' or operation == 'X':
                    read_pos += count
                    ref_pos += count
                elif operation == 'I':
                    read_pos += count
                    # Don't advance ref_pos for insertions
                elif operation == 'D' or operation == 'N':
                    # Don't advance read_pos for deletions/skips
                    ref_pos += count
                elif operation == 'H':
                    pass # Hard clipping doesn't consume read bases or ref pos
                elif operation == 'P':
                    pass # Padding doesn't consume read or ref bases

            # Recalculate ref_end_pos based on final ref_pos
            ref_end_pos = ref_pos - 1
            # Store comprehensive read information
            if read_id not in read_info:
                 read_info[read_id] = {
                    'ref_start': ref_start_pos,
                    'ref_end': ref_end_pos,
                    'read_start': read_start_pos, 
                    'read_end': read_end_pos,     
                    'read_length': seq_length,
                    'has_start_soft_clip': has_start_soft_clip,
                    'has_end_soft_clip': has_end_soft_clip,
                    'strand_info': strand_info,
                    'start_soft_seq': start_soft_seq,
                    'end_soft_seq': end_soft_seq,
                    'update_round': 1,
                    'alignment_round': 1
                }
            else:
                # This case should ideally not happen with --excl-flags 2316 unless supplementary alignments
                # are within the same chunk and not filtered. Log it.
                logging.warning(f"Read {read_id} encountered multiple times within chunk {chrom}:{chunk_start}-{chunk_end}. Keeping first encountered record.")
                # print(f"Read {read_id} already exists in read_info")
                # print(read_info[read_id])
                continue # Keep the first one encountered

    finally:
        # Ensure resources are cleaned up
        process.stdout.close()
        ret_code = process.wait()
        if ret_code != 0:
            # Log error even if some lines were processed
             logging.error(f"samtools view command finished with non-zero exit code {ret_code} for chunk {chrom}:{chunk_start}-{chunk_end}")
             # Decide if partial results are acceptable or if an exception should be raised
             # For now, log and return potentially partial results.
             # raise RuntimeError(f"samtools view command failed with code {ret_code}")


    return read_info, read_id_set # Return both dict and set

# Modified function to handle both outputs
def run(bam_file, output_dir, chrom, chunk_num, chunk_idx, chunk_size=None, samtools_execute_command="samtools"):
    os.makedirs(output_dir, exist_ok=True)
    log_file_prefix = f"{chrom}_chunk{chunk_idx}"
    log_file_path = os.path.join(output_dir, f"{log_file_prefix}_processing.log") # Changed log filename slightly

    # Setup file handler for logging
    log_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    log_handler = logging.FileHandler(log_file_path, mode='w') # Overwrite log each run
    log_handler.setFormatter(log_formatter)

    # Get root logger and add handler
    logger = logging.getLogger()
    # Remove existing handlers if any (e.g., from basicConfig in previous runs or imports)
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    logger.addHandler(log_handler)
    logger.setLevel(logging.INFO)


    # Check if file exists
    if not os.path.exists(bam_file):
        error_msg = f"BAM file {bam_file} does not exist"
        logging.error(error_msg)
        raise FileNotFoundError(error_msg)

    if not os.path.exists(f"{bam_file}.bai"):
        error_msg = f"BAM index file {bam_file}.bai does not exist"
        logging.error(error_msg)
        raise FileNotFoundError(error_msg)

    # Get chromosome length and calculate chunk boundaries
    chrom_length = get_chrom_length(bam_file, chrom, samtools_execute_command)
    if not chrom_length:
        error_msg = f"Chromosome {chrom} not found in BAM file {bam_file}"
        logging.error(error_msg)
        raise ValueError(error_msg)

    # Calculate chunk start and end positions and total chunks
    total_chunks = 0
    if chunk_size:
        if chunk_size <= 0:
            raise ValueError("chunk_size must be positive")
        chunk_start = 1 + (chunk_idx - 1) * chunk_size
        chunk_end = min(chunk_idx * chunk_size, chrom_length)
        total_chunks = (chrom_length + chunk_size - 1) // chunk_size
        logging.info(f"Using fixed chunk size: {chunk_size}. Processing chunk {chunk_idx} of {total_chunks} for chromosome {chrom}")
    else:
        if chunk_num is None or chunk_num <= 0:
             raise ValueError("chunk_num must be provided and positive if chunk_size is not specified")
        calculated_chunk_size = (chrom_length + chunk_num - 1) // chunk_num
        chunk_start = 1 + (chunk_idx - 1) * calculated_chunk_size
        chunk_end = min(chunk_idx * calculated_chunk_size, chrom_length)
        total_chunks = chunk_num
        logging.info(f"Dividing into {chunk_num} chunks. Processing chunk {chunk_idx} of {total_chunks} for chromosome {chrom}")

    # Validate chunk index against calculated total chunks
    if chunk_idx < 1 or chunk_idx > total_chunks:
        error_msg = f"Chunk index {chunk_idx} is out of range. Must be between 1 and {total_chunks} for chromosome {chrom}."
        logging.error(error_msg)
        raise ValueError(error_msg)

    # Define output file paths
    read_info_file = os.path.join(output_dir, f"{chrom}_chunk{chunk_idx}_read_info.pkl")
    read_id_set_file = os.path.join(output_dir, f"{chrom}_chunk{chunk_idx}_read_ids.pkl") # New file path for the set


    # Check if the calculated chunk is valid (start <= end)
    if chunk_start > chunk_end:
        logging.warning(f"Chunk {chunk_idx} ({chrom}:{chunk_start}-{chunk_end}) is empty or invalid (start > end). This likely means the chunk index is too high for the given chromosome length and chunking parameters.")
        # Save empty dictionary and empty set for consistency
        read_info = {}
        read_id_set = set()
        with open(read_info_file, 'wb') as f_info:
            pickle.dump(read_info, f_info, protocol=4)
        with open(read_id_set_file, 'wb') as f_set: # Save empty set
            pickle.dump(read_id_set, f_set, protocol=4)
        logging.info(f"Saved empty read_info and read_ids files for empty chunk {chunk_idx}.")
        logging.info(f"Completed processing empty chunk {chunk_idx} of {total_chunks}")
        return 0 # Return 0 reads processed

    logging.info(f"Processing region {chrom}:{chunk_start}-{chunk_end}")

    # Process reads for this chunk, getting both info dict and id set
    try:
        read_info, read_id_set = save_read_info(bam_file, output_dir, chrom, chunk_start, chunk_end, samtools_execute_command)
    except Exception as e:
        logging.error(f"Error processing chunk {chunk_idx} ({chrom}:{chunk_start}-{chunk_end}): {e}")
        # Depending on requirements, you might want to re-raise the exception
        # or return an error indicator, or save empty files.
        # Saving empty files to avoid downstream failures that expect files.
        read_info = {}
        read_id_set = set()
        with open(read_info_file, 'wb') as f_info:
            pickle.dump(read_info, f_info, protocol=4)
        with open(read_id_set_file, 'wb') as f_set:
            pickle.dump(read_id_set, f_set, protocol=4)
        logging.info(f"Saved empty read_info and read_ids files due to error in chunk {chunk_idx}.")
        raise # Re-raise the exception after logging and cleanup attempt

    # Save read information dictionary for this chunk
    with open(read_info_file, 'wb') as f_info:
        pickle.dump(read_info, f_info, protocol=4)
    logging.info(f"Saved read info dictionary to {read_info_file}")

    # Save read ID set for this chunk
    with open(read_id_set_file, 'wb') as f_set:
        pickle.dump(read_id_set, f_set, protocol=4)
    logging.info(f"Saved read ID set to {read_id_set_file}")


    num_reads_in_dict = len(read_info)
    num_unique_ids = len(read_id_set)
    logging.info(f"Completed processing chunk {chunk_idx} of {total_chunks}")
    logging.info(f"Extracted info for {num_reads_in_dict} primary alignments (unique keys in dict).")
    logging.info(f"Found {num_unique_ids} unique read IDs in this chunk.")

    return num_reads_in_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract read information and read IDs from BAM file for a specific chromosome chunk.")
    parser.add_argument("--bam_file", "-b", type=str, required=True, help="The path to the indexed BAM file")
    parser.add_argument("--output_dir", "-o", type=str, required=True, help="The path to the output directory")
    parser.add_argument("--chrom", type=str, required=True, help="The chromosome to be processed")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--chunk_num", type=int, help="Total number of chunks to divide the chromosome into (alternative to --chunk_size)")
    group.add_argument("--chunk_size", type=int, help="Fixed size for each chunk (alternative to --chunk_num)")
    parser.add_argument("--chunk_idx", type=int, required=True, help="The index of the chunk to process (1-based)")
    parser.add_argument("--samtools_path", type=str, default="samtools", help="Path to the samtools executable (default: 'samtools')")


    args = parser.parse_args()

    if args.chunk_idx < 1:
        # No logging is set up yet, print to stderr and exit
        print("Error: Chunk index must be 1-based (>= 1)", file=os.sys.stderr)
        exit(1)
        # raise ValueError("Chunk index must be 1-based (>= 1)") # Can't log yet

    try:
        num_reads = run(
            args.bam_file,
            args.output_dir,
            args.chrom,
            args.chunk_num,
            args.chunk_idx,
            args.chunk_size,
            args.samtools_path # Use the provided samtools path
        )

        # Logging is configured inside 'run', print final messages to stdout
        print(f"Successfully processed chunk {args.chunk_idx} for chromosome {args.chrom}.")
        print(f"Output files saved in: {args.output_dir}")
        # The exact number reported might depend on interpretation (unique IDs vs dict entries)
        # Sticking to dict entries based on run's return value
        print(f"Extracted information for {num_reads} primary alignment records.")

    except FileNotFoundError as e:
        print(f"Error: {e}", file=os.sys.stderr)
        exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=os.sys.stderr)
        exit(1)
    except RuntimeError as e:
         # RuntimeErrors often come from subprocesses, logging inside run should capture details
        print(f"Runtime Error during processing: {e}", file=os.sys.stderr)
        print(f"Check log file in {args.output_dir} for details ({args.chrom}_chunk{args.chunk_idx}_processing.log).", file=os.sys.stderr)
        exit(1)
    except Exception as e:
        # Catch any other unexpected errors
        print(f"An unexpected error occurred: {e}", file=os.sys.stderr)
        # Attempt to guide user to logs if run() was entered
        log_file_path = os.path.join(args.output_dir, f"{args.chrom}_chunk{args.chunk_idx}_processing.log")
        print(f"Check log file ({log_file_path}) if it exists.", file=os.sys.stderr)
        exit(1)