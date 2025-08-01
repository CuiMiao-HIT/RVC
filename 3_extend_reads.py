import os
import argparse
import pickle
import logging
from Bio import SeqIO
from subprocess import Popen, PIPE
import shlex
import  glob
import ksw2_python as ksw2
import numpy as np
import time
import gzip
import gc
from collections import defaultdict
def subprocess_popen(args,shell=False, stdin=None, stdout=PIPE, bufsize=8388608):
    return Popen(args, shell=shell, stdin=stdin, stdout=stdout, bufsize=bufsize, universal_newlines=True)

dna_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
mat = np.array([
    [1, -2, -2, -2],
    [-2, 1, -2, -2],
    [-2, -2, 1, -2],
    [-2, -2, -2, 1]
], dtype=np.int8).flatten()
q = 2
e = 1
flag = ksw2.EZ_EQX
def dna_to_array(seq):
    return np.array([dna_map[c] for c in seq.upper() if c in dna_map], dtype=np.uint8)
def load_chromosome_candidates(candidate_dir: str, chrom: str):
    """
    Load all candidate positions from chunked pickle files for a specific chromosome
    
    Args:
        candidate_dir: Directory containing the pickle files
        chrom: Chromosome name
        
    Returns:
        dict: Combined dictionary of all candidate positions
    """

    pickle_pattern = os.path.join(candidate_dir, f"{chrom}_chunk*_candidate_positions.pkl")
    pickle_files = glob.glob(pickle_pattern)
    
    if not pickle_files:
        raise FileNotFoundError(f"No candidate files found matching pattern: {pickle_pattern}")

    all_candidates = {}
    
    for pickle_file in sorted(pickle_files):  
        try:
            with open(pickle_file, 'rb') as f:
                chunk_candidates = pickle.load(f)
                all_candidates.update(chunk_candidates)
        except Exception as e:
            print(f"Error loading {pickle_file}: {str(e)}")
            continue
    
    print(f"Loaded {len(all_candidates)} candidate positions from {len(pickle_files)} chunks for chromosome {chrom}")
    
    return all_candidates
def load_single_chunk_candidate_positions(candidate_dir: str, chrom: str, chunk_idx: int,round_num: int):
    """
    Load candidate positions from a specific chunk pickle file for a chromosome
    
    Args:
        candidate_dir: Directory containing the pickle files
        chrom: Chromosome name
        chunk_idx: Chunk index number
        
    Returns:
        dict: Dictionary of candidate positions for the specified chunk
    """
    if int(round_num) >=3:
        pickle_file = os.path.join(candidate_dir, f"{chrom}_chunk{chunk_idx}_{int(round_num)-1}_candidate_positions.pkl")
    else:
        pickle_file = os.path.join(candidate_dir, f"{chrom}_chunk{chunk_idx}_candidate_positions.pkl")
    
    try:
        with open(pickle_file, 'rb') as f:
            chunk_candidates = pickle.load(f)
            print(f"Loaded {len(chunk_candidates)} candidate positions from {pickle_file}")
            return chunk_candidates
    except FileNotFoundError:
        print(f"Warning: Pickle file not found: {pickle_file}")
        return {}
    except Exception as e:
        print(f"Error loading {pickle_file}: {str(e)}")
        return {}

def reference_sequence_from(fasta_file_path, regions, samtools_execute_command="samtools"):
    refernce_sequences = []
    region_value_for_faidx = " ".join(regions)

    samtools_faidx_process = subprocess_popen(
        shlex.split("{} faidx {} {}".format(samtools_execute_command, fasta_file_path, region_value_for_faidx))
    )
    while True:
        row = samtools_faidx_process.stdout.readline()
        is_finish_reading_output = row == '' and samtools_faidx_process.poll() is not None
        if is_finish_reading_output:
            break
        if row:
            refernce_sequences.append(row.rstrip())

    # first line is reference name ">xxxx", need to be ignored
    reference_sequence = "".join(refernce_sequences[1:])

    # uppercase for masked sequences
    reference_sequence = reference_sequence.upper()

    samtools_faidx_process.stdout.close()
    samtools_faidx_process.wait()
    if samtools_faidx_process.returncode != 0:
        return None

    return reference_sequence
def batch_read_sequences_from_seqkit(fastq_file_path, read_ids, seqkit_execute_command="seqkit"):
    """
    Batch read sequences from a fastq file using seqkit grep.
    
    Args:
        fastq_file_path: Path to the fastq file
        read_ids: List of read IDs to extract
        seqkit_execute_command: Command to execute seqkit
        
    Returns:
        Dictionary mapping read IDs to their sequences.
    """
    if not read_ids:
        return {}
    
    read_sequences = {}
    
    # Create a temporary file with read IDs
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        for read_id in read_ids:
            temp_file.write(f"{read_id}\n")
        temp_file_path = temp_file.name
    
    try:
        # Use seqkit grep to extract sequences by ID
        seqkit_process = subprocess_popen(
            shlex.split(f"{seqkit_execute_command} grep -f {temp_file_path} {fastq_file_path}")
        )
        
        stdout = seqkit_process.stdout
        
        # Parse the FASTQ output
        current_id = None
        line_count = 0
        
        while True:
            line = stdout.readline()
            # print(line)
            if not line:  # End of stream
                break
                
            line = line.strip()
            if not line:
                continue
            # print(line)
            line_position = line_count % 4
            
            if line_position == 0:  # Header line
                if line.startswith('@'):
                    current_id = line[1:].split()[0]
                else:
                    logging.warning(f"Unexpected line format, expected header: {line}")
                    current_id = None
            elif line_position == 1:  # Sequence line
                if current_id is not None:
                    read_sequences[current_id] = line.upper()
            # line_position 2 is '+' line, line_position 3 is quality line - we skip both
            
            line_count += 1
        
        seqkit_process.stdout.close()
        seqkit_process.wait()
        
        # Check return code
        if seqkit_process.returncode != 0:
            logging.error(f"Seqkit grep process failed with return code {seqkit_process.returncode}")
        
    finally:
        # Clean up temporary file
        if os.path.exists(temp_file_path):
            os.unlink(temp_file_path)
    
    return read_sequences
def read_sequence_from(fastq_file_path, read_id, samtools_execute_command="samtools"):
    region_value_for_faidx = "{}".format(read_id)

    samtools_faidx_process = subprocess_popen(
        shlex.split("{} fqidx {} {}".format(samtools_execute_command, fastq_file_path, region_value_for_faidx))
    )
    
    sequence = None
    line_count = 0
    while True:
        row = samtools_faidx_process.stdout.readline()
        is_finish_reading_output = row == '' and samtools_faidx_process.poll() is not None
        if is_finish_reading_output:
            break
        
        line_count += 1
        if line_count == 2: 
            sequence = row.rstrip().upper()  
    samtools_faidx_process.stdout.close()
    samtools_faidx_process.wait()
    if samtools_faidx_process.returncode != 0:
        return None

    return sequence

def convert_read_id_to_idx(read_info):
    read_id_list = list(read_info.keys())
    read_id_to_idx = {read_id: idx for idx, read_id in enumerate(read_id_list)}
    return read_id_list, read_id_to_idx


EVENT_TYPE = {
    "ORI_READ_START": 0,
    "NEW_READ_START": 1,
    "SNP_START": 2,
    "INSERTION": 3,
    "DELETION_START": 4,
    "DELETION_END": 5,
    "NEW_READ_END": 6,
    "ORI_READ_END": 7
}


#EVENT:(position, type, idx)
def construct_events(candidate_positions, read_info, read_length=50):
    events = []
    read_id_list, read_id_to_idx = convert_read_id_to_idx(read_info)
    print("INFO: Constructing candidate events")
    
    # Create events for candidate positions
    for position, info in candidate_positions.items():
        if info["has_snp"]:
            events.append((position, EVENT_TYPE["SNP_START"], position))
        if info["has_indel"]:
            base_list = info["base_list"]
            max_length = 0
            max_length_type = None
            for base in base_list:
                if len(base[1]) > max_length:
                    max_length = len(base[1])
                    max_length_type = base[1][0]
            if max_length_type == "+":
                events.append((position, EVENT_TYPE["INSERTION"], position))
            else:
                events.append((position, EVENT_TYPE["DELETION_START"], position))
                events.append((position + max_length, EVENT_TYPE["DELETION_END"], position))
    
    candidate_events_num = len(events)
    print(f"INFO: Constructed {candidate_events_num} candidate events")
    print("INFO: Constructing read events")
    
    for read_id, info in read_info.items():
        read_idx = read_id_to_idx[read_id]
        strand_info = info["strand_info"]
        
        if strand_info == "+":  # Positive strand
            if info["has_end_soft_clip"]:
                soft_length = len(info['end_soft_seq'])
                start_pos = info["ref_end"] + 1
                end_pos = start_pos + read_length + soft_length - 1
            else:
                start_pos = info["ref_end"] + 1
                end_pos = start_pos + read_length - 1
        else:  # Negative strand
            if info["has_start_soft_clip"]:
                soft_length = len(info['start_soft_seq'])
                end_pos = info["ref_start"] - 1
                start_pos = end_pos - read_length - soft_length + 1
            else:
                end_pos = info["ref_start"] - 1
                start_pos = end_pos - read_length + 1
        
        # Ensure start_pos is always less than end_pos
        if start_pos > end_pos:
            start_pos, end_pos = end_pos, start_pos
        
        events.append((start_pos, EVENT_TYPE["NEW_READ_START"], read_idx))
        events.append((end_pos, EVENT_TYPE["NEW_READ_END"], read_idx))
    
    read_events_num = len(events) - candidate_events_num
    print(f"INFO: Constructed {read_events_num} read events")
    
    # Sort events by position
    sorted_events = sorted(events, key=lambda x: (x[0], x[1]))
    return sorted_events, read_id_list, read_id_to_idx


def batch_read_sequences_from(fastq_file_path, read_ids, samtools_execute_command="samtools"):
    """
    Batch read sequences from a fastq file using samtools fqidx.
    Handles cases where a read might have no sequence data (represented by 2 lines: header and '+').
    
    Args:
        fastq_file_path: Path to the fastq file
        read_ids: List of read IDs to extract
        samtools_execute_command: Command to execute samtools
        
    Returns:
        Dictionary mapping read IDs to their sequences. Empty string for sequence-less reads.
    """
    if not read_ids:
        return {}
    
    # Join read IDs with spaces for samtools fqidx command
    region_value_for_fqidx = " ".join(read_ids)
    
    read_sequences = {}
    
    samtools_fqidx_process = subprocess_popen(
        shlex.split(f"{samtools_execute_command} fqidx -c {fastq_file_path} {region_value_for_fqidx}")
    )
    
    stdout = samtools_fqidx_process.stdout
    
    while True:
        header_line = stdout.readline()
        if not header_line:  # End of stream
            break
        
        header = header_line.rstrip()
        
        # Expecting a header line starting with '@'
        if not header.startswith('@'):
            logging.warning(f"Unexpected line format, expected header: {header}")
            continue

        current_read_id = header[1:].split()[0]
        
        # Read the next line: this is either the sequence or a '+' indicating no sequence
        line2 = stdout.readline()
        if not line2: 
            logging.warning(f"Premature EOF after header for read ID: {current_read_id}")
            break
            
        line2_stripped = line2.rstrip()
        
        if line2_stripped == '+':
            continue
        else:
            # This is a normal read; line2_stripped is the sequence
            read_sequences[current_read_id] = line2_stripped.upper()
            
            # Read and discard the next two lines: the '+' separator and the quality score
            plus_line = stdout.readline()
            if not plus_line:
                logging.warning(f"Premature EOF after sequence for read ID: {current_read_id}, expected '+' line.")
                break 
            
            quality_line = stdout.readline()
            if not quality_line:
                logging.warning(f"Premature EOF after '+' for read ID: {current_read_id}, expected quality line.")
                break
    
    samtools_fqidx_process.stdout.close()
    samtools_fqidx_process.wait()
    
    # Optionally, check samtools_fqidx_process.returncode for errors
    if samtools_fqidx_process.returncode != 0:
        logging.error(f"Samtools fqidx process failed with return code {samtools_fqidx_process.returncode}")

    return read_sequences          
    
def get_min_max_pos(read_batch, read_info, flanking_size=200):
    """
    Get the minimum and maximum reference positions for a batch of reads with flanking regions.
    
    Args:
        read_batch: List of read IDs to process
        read_info: Dictionary containing read information
        flanking_size: Size of flanking region to add on both sides
        
    Returns:
        tuple: (min_pos, max_pos, position_offset)
        - min_pos: Minimum reference position with flanking
        - max_pos: Maximum reference position with flanking 
        - position_offset: Value to subtract from reference positions for array indexing
    """
    # Find the min and max positions in the batch
    min_pos = None
    max_pos = None
    for read_id in read_batch:
        if read_id not in read_info:
            continue
            
        # Get the reference end position for this read
        ref_end = read_info[read_id]['ref_end']
        
        # Update min and max positions
        if min_pos is None or ref_end < min_pos:
            min_pos = ref_end
        if max_pos is None or ref_end > max_pos:
            max_pos = ref_end
    
    # Add flanking regions
    # For min_pos, we add flanking at the start
    min_pos = max(1, min_pos - flanking_size)  # Ensure position is at least 1
    
    # For max_pos, we extend for alignment and add flanking
    max_pos = max_pos + 200 + flanking_size  # 200 is the default extension for alignment
    
    # The offset is the minimum position (for converting absolute positions to array indices)
    position_offset = min_pos - 1  # Subtract 1 because reference is 1-based
    
    return min_pos, max_pos, position_offset

def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.
    
    Args:
        seq: DNA sequence string
        
    Returns:
        Reverse complemented sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                  'N': 'N', 'n': 'n'}
    return ''.join([complement.get(base, base) for base in reversed(seq)])

def complement_sequence(seq):
    """
    Return the complement of a DNA sequence.
    
    Args:
        seq: DNA sequence string
        
    Returns:
        Complemented sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                  'N': 'N', 'n': 'n'}
    return ''.join([complement.get(base, base) for base in seq])

def complement_base(base):
    """
    Return the complement of a DNA base.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                  'N': 'N', 'n': 'n'}
    return complement.get(base, base)
def align_and_update_candidates_batch(new_fastq, last_fastq, read_batch, read_info, candidate_positions,
                                      ref_fasta, chrom, round_num, middle_fastq_needed_set,
                                      recall_dict, debug_mode=False, identity_threshold=0.95):
    # Batch load sequences from the new fastq
    if debug_mode:
        debug_output = open("./debug_output", "a")
    start_time = time.time()
    new_read_sequences = batch_read_sequences_from_seqkit(new_fastq, read_batch)
    end_time = time.time()
    print(f"INFO: Time taken to batch read {len(read_batch)} new read sequences: {end_time - start_time} seconds")
    
    # Filter to reads that we have sequences for
    middle_fastq_needed = []
    if int(round_num) > 2:
        middle_fastq_needed = [read_id for read_id in read_batch 
                            if read_id in new_read_sequences and 
                            read_id in middle_fastq_needed_set]
    print(f"INFO: Middle fastq needed: {len(middle_fastq_needed)}")
    
    # Batch load sequences from the last fastq if needed
    middle_read_sequences = {}
    if middle_fastq_needed:
        start_time = time.time()
        middle_read_sequences = batch_read_sequences_from_seqkit(last_fastq, middle_fastq_needed)
        end_time = time.time()
        print(f"INFO: Time taken to batch read {len(middle_fastq_needed)} middle fastq sequences: {end_time - start_time} seconds")
    del middle_fastq_needed
    gc.collect()
    
    # Separate reads by strand
    pos_strand_reads = []
    neg_strand_reads = []
    
    for read_id in read_batch:
        if read_id in new_read_sequences:
            if read_info[read_id]['strand_info'] == "+":
                pos_strand_reads.append(read_id)
            else:
                neg_strand_reads.append(read_id)
    
    print(f"INFO: Processing {len(pos_strand_reads)} positive strand reads and {len(neg_strand_reads)} negative strand reads")
    # del read_batch
    gc.collect()
    if pos_strand_reads:
        min_pos, max_pos, position_offset = get_min_max_pos(pos_strand_reads, read_info, flanking_size=200)
        start_time = time.time()
        reference_range = f"{chrom}:{min_pos}-{max_pos}"
        reference_sequence = reference_sequence_from(ref_fasta, [reference_range])
        if reference_sequence is None:
            logging.error(f"Failed to load reference sequence for range {reference_range}")
            return read_info
            
        ref_length = len(reference_sequence)
        end_time = time.time()
        print(f"INFO: Time taken to read reference sequence for positive strand ({reference_range}, length={ref_length}): {end_time - start_time} seconds")
        
        # Process each positive strand read
        start_deletion = 0
        pos_time_start = time.time()
        for read_id in pos_strand_reads:
            # Skip if read sequence wasn't found
            if read_id not in new_read_sequences:
                continue
            delete_flag = False
            new_read_sequence = new_read_sequences[read_id]
            
            # Get middle sequence if needed
            middle_fastq_sequence = ""
            if read_id in middle_read_sequences:
                middle_fastq_sequence = middle_read_sequences[read_id]
            
            # Create the complete read to align
            read_to_align = read_info[read_id]['end_soft_seq'] + middle_fastq_sequence + new_read_sequence
            # Get the reference start position for this read
            reference_start_to_align = read_info[read_id]['ref_end'] + 1
            
            # Calculate the offset within our reference sequence
            local_start_pos = reference_start_to_align - min_pos
            
            # Ensure the local_start_pos is valid
            if local_start_pos < 0:
                local_start_pos = 0
            
            # Calculate end position with extra buffer for alignment
            local_end_pos = min(local_start_pos + len(read_to_align) + 200, ref_length)
            
            # Extract the relevant part of the reference for alignment
            ref_for_alignment = reference_sequence[local_start_pos:local_end_pos]
            
            # Convert to alignment arrays
            read_array = dna_to_array(read_to_align)
            ref_array = dna_to_array(ref_for_alignment)
            
            # Perform alignment
            align_result = ksw2.align(read_array, ref_array, mat, q, e, flag)
            cigar_tuple = align_result['cigar']
            
            read_pos = 0
            ref_pos = 0  # Relative to the extracted reference segment
            soft_clip_flag = False
            
            if cigar_tuple is None:
                continue
                
            # Process CIGAR operations and update candidates
            match_count = 0
            mismatch_count = 0
            pending_updates = []  # Store pending updates: (type, *args)
            
            for idx, (length, op) in enumerate(cigar_tuple):
                if op == ksw2.CIGAR_MATCH:  # Match or mismatch
                    for i in range(length):
                        if read_pos >= len(read_to_align) or ref_pos >= len(ref_for_alignment):
                            break
                            
                        read_base = read_to_align[read_pos]
                        ref_base = ref_for_alignment[ref_pos]
                        
                        # Count matches and mismatches for identity calculation
                        if read_base.upper() == ref_base.upper():
                            match_count += 1
                        else:
                            mismatch_count += 1
                        
                        # Convert back to genome coordinates for candidate lookup
                        genome_pos = local_start_pos + ref_pos + min_pos
                        if genome_pos in candidate_positions:
                            pending_updates.append(('snp', genome_pos, read_id, read_base.upper(), ""))
                        else:
                            if read_base != ref_base:
                                pending_updates.append(('recall', str(genome_pos)+"_"+ref_base+"_"+read_base))
                        read_pos += 1
                        ref_pos += 1
                        
                elif op == ksw2.CIGAR_INS:  # Insertion
                    if read_pos >= len(read_to_align):
                        continue
                    if read_pos + length >= len(read_to_align):
                        # Insertion at the end, treat as soft clip
                        soft_clip_flag = True
                        read_info[read_id]['end_soft_seq'] = read_to_align[read_pos:]
                        read_info[read_id]['has_end_soft_clip'] = True
                        read_info[read_id]['alignment_round'] = int(round_num)
                        read_info[read_id]['update_round'] = int(round_num)
                        read_pos += length
                    else:                # Convert to genome coordinates
                        genome_pos = local_start_pos + ref_pos + min_pos - 1
                        insert_seq = read_to_align[read_pos:read_pos+length].upper()
                        insert_seq = "+" + insert_seq
                        if genome_pos in candidate_positions and read_pos + length <= len(read_to_align):
                            pending_updates.append(('indel', genome_pos, read_id, "", insert_seq.upper()))
                        read_pos += length
                        
                elif op == ksw2.CIGAR_DEL:  # Deletion
                    if debug_mode and idx == 0:
                        debug_output.write(f"INFO: Deletion: read_id: {read_id}, \n length: {length}, \n op: {op}, \n read_pos: {read_pos}, \n ref_pos: {ref_pos}, \n ref_for_alignment: {ref_for_alignment}, \n read_to_align: {read_to_align}, \n middle_fastq_sequence: {middle_fastq_sequence}, \n new_read_sequence: {new_read_sequence}\n")
                    if idx == 0:
                        start_deletion += 1 
                        if length > 10:
                            break
                        
                    if read_pos < 0 or ref_pos + length >= len(ref_for_alignment):
                        break
                    if length > 50:
                        break
                    genome_pos = local_start_pos + ref_pos + min_pos - 1
                    del_end = min(ref_pos + length, len(ref_for_alignment))
                    deletion_seq = "-" + ref_for_alignment[ref_pos:del_end]
                    if genome_pos in candidate_positions:
                        pending_updates.append(('indel', genome_pos, read_id, "", deletion_seq.upper()))
                    ref_pos += length
            
            # Calculate sequence identity and apply updates if threshold is met
            total_aligned = match_count + mismatch_count
            if total_aligned > 0:
                identity = match_count / total_aligned
                if identity >= identity_threshold:
                    # Apply all pending updates
                    for update in pending_updates:
                        if update[0] == 'snp':
                            _, genome_pos, read_id_update, base, indel = update
                            candidate_positions[genome_pos]['read_name_list'].append(read_id_update)
                            candidate_positions[genome_pos]['base_list'].append([base, indel])
                            candidate_positions[genome_pos]['read_name_dict'] = len(
                                candidate_positions[genome_pos]['read_name_list'])
                        elif update[0] == 'indel':
                            _, genome_pos, read_id_update, base, indel = update
                            candidate_positions[genome_pos]['read_name_list'].append(read_id_update)
                            candidate_positions[genome_pos]['base_list'].append([base, indel])
                            candidate_positions[genome_pos]['read_name_dict'] = len(
                                candidate_positions[genome_pos]['read_name_list'])
                        elif update[0] == 'recall':
                            recall_dict[update[1]] += 1
            
            # Update read info with alignment results
            genome_ref_end = local_start_pos + ref_pos + min_pos
            if not delete_flag:
                if soft_clip_flag:
                    read_info[read_id]['ref_end'] = genome_ref_end
                else:
                    if read_pos < len(read_to_align):
                        read_info[read_id]['ref_end'] = genome_ref_end 
                        read_info[read_id]['has_end_soft_clip'] = True
                        read_info[read_id]['end_soft_seq'] = read_to_align[read_pos:]
                        read_info[read_id]['alignment_round'] = int(round_num)
                        read_info[read_id]['update_round'] = int(round_num)
                    else:
                        read_info[read_id]['ref_end'] = genome_ref_end
                        read_info[read_id]['has_end_soft_clip'] = False
                        read_info[read_id]['end_soft_seq'] = ""
                        read_info[read_id]['alignment_round'] = int(round_num)
                        read_info[read_id]['update_round'] = int(round_num)
        print(f"{start_deletion} start deletions")
        pos_time_end = time.time()
        print(f"INFO: Time taken to process positive strand reads: {pos_time_end - pos_time_start} seconds")
    pos_del_time_start = time.time()
    del pos_strand_reads
    gc.collect()
    pos_del_time_end = time.time()
    print(f"INFO: Time taken to delete positive strand reads: {pos_del_time_end - pos_del_time_start} seconds")
    neg_time_start = time.time()
    if neg_strand_reads:
        neg_min_pos = float('inf')
        neg_max_pos = 0
        start_deletion = 0
        for read_id in neg_strand_reads:
            if read_id not in read_info:
                continue
                
            ref_start = read_info[read_id]['ref_start']
            read_length = len(new_read_sequences[read_id]) if read_id in new_read_sequences else 0
            current_min = max(1, ref_start - read_length - 200)
            current_max = ref_start + 200
            
            neg_min_pos = min(neg_min_pos, current_min)
            neg_max_pos = max(neg_max_pos, current_max)
        
        neg_min_pos = int(neg_min_pos)
        neg_max_pos = int(neg_max_pos)
        start_time = time.time()
        neg_reference_range = f"{chrom}:{neg_min_pos}-{neg_max_pos}"
        neg_reference_sequence = reference_sequence_from(ref_fasta, [neg_reference_range])
        if neg_reference_sequence is None:
            logging.error(f"Failed to load reference sequence for range {neg_reference_range}")
            return read_info
            
        neg_ref_length = len(neg_reference_sequence)
        end_time = time.time()
        print(f"INFO: Time taken to read reference sequence for negative strand ({neg_reference_range}, length={neg_ref_length}): {end_time - start_time} seconds")
        for read_id in neg_strand_reads:
            # Skip if read sequence wasn't found
            if read_id not in new_read_sequences:
                continue 
            delete_flag = False
            new_read_sequence = new_read_sequences[read_id]
            middle_fastq_sequence = ""
            if read_id in middle_read_sequences:
                middle_fastq_sequence = middle_read_sequences[read_id] 
            read_to_align = ''.join(reversed(read_info[read_id].get('start_soft_seq', ''))) + middle_fastq_sequence + new_read_sequence
            reference_start_to_align = read_info[read_id]['ref_start'] 
            local_start_pos = reference_start_to_align - neg_min_pos
            if local_start_pos < 0:
                continue
            if local_start_pos >= neg_ref_length:
                continue
            local_start_for_align = local_start_pos
            local_end_for_align = local_start_pos - len(read_to_align) - 200
            if local_end_for_align < 0:
                continue
            if local_start_for_align > neg_ref_length:
                continue
            ref_for_alignment = neg_reference_sequence[local_end_for_align:local_start_for_align]
            ref_for_alignment = reverse_complement(ref_for_alignment)
            read_array = dna_to_array(read_to_align)
            ref_array = dna_to_array(ref_for_alignment)
            align_result = ksw2.align(read_array, ref_array, mat, q, e, flag)
            cigar_tuple = align_result['cigar']
            
            read_pos = 0
            ref_pos = 0  # Relative to the extracted reference segment
            soft_clip_flag = False
            
            if cigar_tuple is None:
                continue
                
            # Process CIGAR operations and update candidates for negative strand
            match_count = 0
            mismatch_count = 0
            pending_updates = []  # Store pending updates: (type, *args)
            
            for idx, (length, op) in enumerate(cigar_tuple):
                if op == ksw2.CIGAR_MATCH:  # Match or mismatch
                    for i in range(length):
                            
                        read_base = read_to_align[read_pos]
                        ref_base = ref_for_alignment[ref_pos]
                        
                        # Count matches and mismatches for identity calculation
                        if read_base.upper() == ref_base.upper():
                            match_count += 1
                        else:
                            mismatch_count += 1
                        
                        # Convert back to genome coordinates for candidate lookup
                        # For negative strand, we need to adjust the coordinates
                        genome_pos = local_start_for_align - ref_pos + neg_min_pos - 1
                        if genome_pos in candidate_positions:

                            pending_updates.append(('snp', genome_pos, read_id, complement_base(read_base.lower()), ""))
                        else:
                            if complement_base(read_base.upper()) != ref_base:
                                pending_updates.append(('recall', str(genome_pos)+"_"+ref_base+"_"+complement_base(read_base.upper())))
                            
                        read_pos += 1
                        ref_pos += 1
                        
                elif op == ksw2.CIGAR_INS:  # Insertion
                    if read_pos >= len(read_to_align):
                        continue
                        
                    if read_pos + length >= len(read_to_align):
                        soft_clip_flag = True
                        read_info[read_id]['start_soft_seq'] = reverse_complement(read_to_align[read_pos:])
                        read_info[read_id]['has_start_soft_clip'] = True
                        read_info[read_id]['alignment_round'] = int(round_num)
                        read_info[read_id]['update_round'] = int(round_num)
                        read_pos += length
                    else:
                        # Convert to genome coordinates
                        genome_pos = local_start_for_align - ref_pos + neg_min_pos - 1 
                        insert_seq = read_to_align[read_pos:read_pos+length].lower()
                        insert_seq = reverse_complement(insert_seq)
                        insert_seq = "+" + insert_seq
                        if genome_pos in candidate_positions and read_pos + length <= len(read_to_align):
                            pending_updates.append(('indel', genome_pos, read_id, "", insert_seq.lower()))
                        read_pos += length
                        
                elif op == ksw2.CIGAR_DEL:  # Deletion
                    if read_pos < 0 or ref_pos + length >= len(ref_for_alignment):
                        # ref_pos += min(length, len(ref_for_alignment) - ref_pos)
                        break
                    if length > 50:
                        break
                        
                    if idx == 0:
                        start_deletion += 1 
                        if length > 10:
                            break
                    genome_pos = local_start_for_align - ref_pos + neg_min_pos-1  - length
                    del_end = min(ref_pos + length, len(ref_for_alignment))
                    deletion_seq = "-" + reverse_complement(ref_for_alignment[ref_pos:del_end])
                    if genome_pos in candidate_positions:
                        pending_updates.append(('indel', genome_pos, read_id, "", deletion_seq.lower()))
                    ref_pos += length
            
            # Calculate sequence identity and apply updates if threshold is met
            total_aligned = match_count + mismatch_count
            if total_aligned > 0:
                identity = match_count / total_aligned
                if identity >= identity_threshold:
                    # Apply all pending updates
                    for update in pending_updates:
                        if update[0] == 'snp':
                            _, genome_pos, read_id_update, base, indel = update
                            candidate_positions[genome_pos]['read_name_list'].append(read_id_update)
                            candidate_positions[genome_pos]['base_list'].append([base, indel])
                            candidate_positions[genome_pos]['read_name_dict'] = len(
                                candidate_positions[genome_pos]['read_name_list'])
                        elif update[0] == 'indel':
                            _, genome_pos, read_id_update, base, indel = update
                            candidate_positions[genome_pos]['read_name_list'].append(read_id_update)
                            candidate_positions[genome_pos]['base_list'].append([base, indel])
                            candidate_positions[genome_pos]['read_name_dict'] = len(
                                candidate_positions[genome_pos]['read_name_list'])
                        elif update[0] == 'recall':
                            recall_dict[update[1]] += 1
            
            # Update read info with alignment results for negative strand
            if not delete_flag:
                genome_ref_start = local_start_for_align - ref_pos + neg_min_pos
                
                if soft_clip_flag:
                    read_info[read_id]['ref_start'] = genome_ref_start
                else:
                    if read_pos < len(read_to_align):
                        read_info[read_id]['ref_start'] = genome_ref_start 
                        read_info[read_id]['has_start_soft_clip'] = True
                        read_info[read_id]['start_soft_seq'] = reverse_complement(read_to_align[read_pos:])
                        read_info[read_id]['alignment_round'] = int(round_num)
                        read_info[read_id]['update_round'] = int(round_num)
                    else:
                        read_info[read_id]['ref_start'] = genome_ref_start
                        read_info[read_id]['has_start_soft_clip'] = False
                        read_info[read_id]['start_soft_seq'] = ""
                        read_info[read_id]['alignment_round'] = int(round_num)
                        read_info[read_id]['update_round'] = int(round_num)
        print(f"{start_deletion} start deletions")
    neg_time_end = time.time()
    print(f"INFO: Time taken to process negative strand reads: {neg_time_end - neg_time_start} seconds")
    neg_del_time_start = time.time()
    del neg_strand_reads
    gc.collect()
    neg_del_time_end = time.time()
    print(f"INFO: Time taken to delete negative strand reads: {neg_del_time_end - neg_del_time_start} seconds")
    if debug_mode:
        print(f"INFO: Start deletion: {start_deletion}")
        debug_output.close()


    return read_info, candidate_positions, recall_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update variant candidates with extended reads")
    parser.add_argument("--new_fastq", type=str, required=True, help="The path to the fastq file")
    parser.add_argument("--last_fastq", type=str, required=True, help="The path to the fastq file")
    parser.add_argument("--round", type=str, required=True, help="The path to the fastq bam file")
    parser.add_argument("--ref_fasta", type=str, required=True, help="The path to the reference fasta file")
    parser.add_argument("--candidate_path", type=str, required=True, help="The path to the candidate positions pickle file")
    parser.add_argument("--read_info_pickle", type=str, required=True, help="The path to the read info pickle file")
    parser.add_argument("--chrom", type=str, required=True, help="The chromosome to be processed")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save output files")
    parser.add_argument("--log_level", type=str, default="INFO", help="Logging level")
    parser.add_argument("--batch_size", type=int, default=500, help="Number of reads to process in each batch")
    parser.add_argument("--tmp_dir", type=str, default="./tmp", help="Directory for temporary files")
    parser.add_argument("--chunk_idx", type=int, default=None, help="Process a specific chunk index (1-based)")
    parser.add_argument("--chunk_num", type=int, default=None, help="Total number of chunks")
    parser.add_argument("--recall_file", type=str, default=None, help="Recall file")
    parser.add_argument("--identity_threshold", type=float, default=0.95, help="Sequence identity threshold (0.95-0.99) for updating candidates")
    parser.add_argument("--read_length", type=int, default=50, help="Read length for constructing events")

    args = parser.parse_args()
    
    # Validate identity threshold
    if not 0.0 <= args.identity_threshold <= 1.0:
        raise ValueError(f"Identity threshold must be between 0.0 and 1.0, got: {args.identity_threshold}")
    
    # Validate read length
    if args.read_length <= 0:
        raise ValueError(f"Read length must be positive, got: {args.read_length}")
    
    # Set up logging
    numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {args.log_level}")
    logging.basicConfig(level=numeric_level, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # Get chromosome from args
    chrom = args.chrom
    
    # Load candidate positions based on whether a specific chunk is requested
    if args.chunk_idx is not None:
        
        logging.info(f"Processing chunk {args.chunk_idx} for chromosome {chrom}")
        candidate_positions = load_single_chunk_candidate_positions(args.candidate_path, chrom, args.chunk_idx, args.round)
        output_suffix = f"_chunk{args.chunk_idx}"
    else:
        # Process all chunks
        logging.info(f"Processing all chunks for chromosome {chrom}")
        candidate_positions = load_chromosome_candidates(args.candidate_path, chrom)
        output_suffix = ""


    if candidate_positions is None or len(candidate_positions) == 0:
        logging.info(f"No candidate positions found for chromosome {chrom}, skipping the rest of the program")
        if args.chunk_idx is not None:
            pickle.dump(candidate_positions, open(os.path.join(args.output_dir, f"{chrom}_chunk{args.chunk_idx}_{args.round}_candidate_positions.pkl"), "wb"))
        else:
            pickle.dump(candidate_positions, open(os.path.join(args.output_dir, f"{chrom}_{args.round}_candidate_positions.pkl"), "wb"))
        exit(0)
    recall_dict = None
    if args.recall_file is not None:
        recall_dict = pickle.load(open(args.recall_file, "rb"))
    else:
        recall_dict = defaultdict(int)
    read_info = pickle.load(open(args.read_info_pickle, "rb"))

    print(f"INFO: Constructing events for chromosome {chrom} round {args.round} with read_length {args.read_length}")
    sorted_events, read_id_list, read_id_to_idx = construct_events(candidate_positions, read_info, read_length=args.read_length)
    compressed_fastq = args.new_fastq
    reads_to_realign = set()
    middle_fastq_needed = set()
    active_new_read_set = set()
    active_candidate_dict = {}
    for (pos, event_type, idx) in sorted_events:
        if event_type == EVENT_TYPE["NEW_READ_START"]:
            active_new_read_set.add(idx)
        elif event_type == EVENT_TYPE["NEW_READ_END"]:
            active_new_read_set.remove(idx)
            if idx not in reads_to_realign:
                read_id = read_id_list[idx]
                # Check strand information to determine how to update read coordinates
                if read_info[read_id]['strand_info'] == "+":
                    # For positive strand, update ref_end as before
                    read_info[read_id]['ref_end'] = pos
                else:
                    read_length = args.read_length
                    # Update ref_start to extend in the upstream direction
                    read_info[read_id]['ref_start'] = max(1, read_info[read_id]['ref_start'] - read_length)
                
                # Update the round information regardless of strand
                read_info[read_id]['update_round'] = int(args.round)

        elif event_type == EVENT_TYPE["SNP_START"]:
            if active_new_read_set:
                for read_idx in active_new_read_set:
                    if read_idx not in reads_to_realign:
                        reads_to_realign.add(read_idx)
                        if read_info[read_id_list[read_idx]]['update_round'] != read_info[read_id_list[read_idx]]['alignment_round']:
                            if int(args.round) > 2:
                                middle_fastq_needed.add(read_id_list[read_idx])

        elif event_type == EVENT_TYPE["INSERTION"]:
            for read_idx in active_new_read_set:
                if read_idx not in reads_to_realign:
                    reads_to_realign.add(read_idx)
                    if read_info[read_id_list[read_idx]]['update_round'] != read_info[read_id_list[read_idx]]['alignment_round']:
                        if int(args.round) > 2:
                            middle_fastq_needed.add(read_id_list[read_idx])

        elif event_type == EVENT_TYPE["DELETION_START"]:
            active_candidate_dict[idx] = {
                "start_covered":[],
                "end_covered":[]
            }
            for read_idx in active_new_read_set:
                active_candidate_dict[idx]["start_covered"].append(read_idx)
        
        elif event_type == EVENT_TYPE["DELETION_END"]:
            for read_idx in active_new_read_set:
                active_candidate_dict[idx]["end_covered"].append(read_idx)
            
            for read_idx in active_candidate_dict[idx]["start_covered"]:
                if read_idx not in reads_to_realign:
                    reads_to_realign.add(read_idx)
                    if read_info[read_id_list[read_idx]]['update_round'] != read_info[read_id_list[read_idx]]['alignment_round']:
                        if int(args.round) > 2:
                            middle_fastq_needed.add(read_id_list[read_idx])
            
            for read_idx in active_candidate_dict[idx]["end_covered"]:
                if read_idx not in reads_to_realign:
                    reads_to_realign.add(read_idx)
                    if read_info[read_id_list[read_idx]]['update_round'] != read_info[read_id_list[read_idx]]['alignment_round']:
                        if int(args.round) > 2:
                            middle_fastq_needed.add(read_id_list[read_idx])
    del active_candidate_dict
    del sorted_events
    del active_new_read_set
    gc.collect()
    
    # Process collected reads in batches
    read_ids_to_realign = [read_id_list[idx] for idx in reads_to_realign]
    total_reads = len(read_ids_to_realign)
    logging.info(f"Total reads to realign: {total_reads}")
    
    for i in range(0, total_reads, args.batch_size):
        batch = read_ids_to_realign[i:i+args.batch_size]
        logging.info(f"Processing batch {i//args.batch_size + 1}/{(total_reads+args.batch_size-1)//args.batch_size}, size: {len(batch)}")
        read_info, candidate_positions, recall_dict = align_and_update_candidates_batch(
            compressed_fastq, 
            args.last_fastq,
            batch, 
            read_info, 
            candidate_positions, 
            args.ref_fasta, 
            chrom, 
            args.round,
            middle_fastq_needed,
            recall_dict,
            identity_threshold=args.identity_threshold
            )
    
    logging.info("Finished updating variant candidates with extended reads")
    #save read_info's keys into a set
    read_info_keys = set(read_info.keys())
    # Save the read_info, recall_dict, and candidate_positions
    if args.chunk_idx is not None:
        pickle.dump(read_info, open(os.path.join(args.output_dir, f"{chrom}_chunk{args.chunk_idx}_{args.round}_read_info.pkl"), "wb"))
        pickle.dump(read_info_keys, open(os.path.join(args.output_dir, f"{chrom}_chunk{args.chunk_idx}_{args.round}_read_ids.pkl"), "wb"))
        pickle.dump(recall_dict, open(os.path.join(args.output_dir, f"{chrom}_chunk{args.chunk_idx}_{args.round}_recall_dict.pkl"), "wb"))
        pickle.dump(candidate_positions, open(os.path.join(args.output_dir, f"{chrom}_chunk{args.chunk_idx}_{args.round}_candidate_positions.pkl"), "wb"))
    else:
        pickle.dump(read_info, open(os.path.join(args.output_dir, f"{chrom}_{args.round}_read_info.pkl"), "wb"))
        pickle.dump(read_info_keys, open(os.path.join(args.output_dir, f"{chrom}_{args.round}_read_ids.pkl"), "wb"))
        pickle.dump(recall_dict, open(os.path.join(args.output_dir, f"{chrom}_{args.round}_recall_dict.pkl"), "wb"))
        pickle.dump(candidate_positions, open(os.path.join(args.output_dir, f"{chrom}_{args.round}_candidate_positions.pkl"), "wb"))
