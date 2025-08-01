import os
import shlex
import argparse
from collections import Counter, defaultdict
from subprocess import Popen, PIPE
from typing import Tuple, List, Dict
import pickle  # Add this import at the top with other imports
from multiprocessing import Pool
import logging
import re  # Add this import at the top of your file if not already present
import random


def subprocess_popen(args,shell=False, stdin=None, stdout=PIPE, bufsize=8388608):
    return Popen(args, shell=shell, stdin=stdin, stdout=stdout, bufsize=bufsize, universal_newlines=True)

def phredscore2raw_score(qual):
    return ord(qual) - 33

def get_chrom_length(bam_file, chrom, samtools_execute_command="samtools"):
    cmd = f"{samtools_execute_command} idxstats {bam_file}"
    process = subprocess_popen(shlex.split(cmd), shell=False)
    for line in process.stdout:
        chr_name, length, *_ = line.strip().split('\t')
        if chr_name == chrom:
            process.terminate()
            return int(length)
    return None

def decode_pileup_bases(pileup_bases, ref_base):
    base_idx = 0
    base_list = []
    while base_idx < len(pileup_bases):
        base = pileup_bases[base_idx]
        if base == '+' or base == '-':
            base_idx += 1
            advance = 0
            while True:
                num = pileup_bases[base_idx]
                if num.isdigit():
                    advance = advance * 10 + int(num)
                    base_idx += 1
                else:
                    break
            base_list[-1][1] = base + pileup_bases[base_idx: base_idx + advance]  # add indel seq
            base_idx += advance - 1

        elif base in "ACGTNacgtn*<>#":
            base_list.append([base, ""])
        elif base in ".,":
            base_list.append([ref_base.upper(), ""]) if base == "." else base_list.append([ref_base.lower(), ""])
        elif base == '^':  # start of read, next base is mq, update mq info
            base_idx += 1
        else:
            base_idx += 1
            continue
        base_idx += 1
    return base_list



def decode_pileup_row(pileup_row: str, min_alt_freq=0.2, min_indel_freq=0.2):
    columns = pileup_row.strip().split('\t')
    contig = columns[0]
    pos = int(columns[1])
    ref_base = columns[2]
    base_list = decode_pileup_bases(columns[4], ref_base)

    read_name_list = columns[7].split(',')
    assert len(base_list)  == len(read_name_list)
    
    # Create read name dictionary
    read_name_dict = {}
    for idx, read_name in enumerate(read_name_list):
        read_name_dict[read_name] = idx
    
    depth = len(base_list)
    if depth == 0:
        return contig, pos, ref_base, read_name_list, base_list, read_name_dict, False, False, depth
    base_to_count_list = []
    for base_char, indel_info in base_list:
        # Create combined representation (e.g., "A", "T+N", "A-GG")
        # Use uppercase base for the combined representation key
        combined_repr_key = base_char.upper()
        if indel_info:
            combined_repr_key += indel_info.upper() # Ensure indel part is also upper for key consistency
        base_to_count_list.append(combined_repr_key)
    base_counter = Counter(base_to_count_list)
    has_snp = False
    has_indel = False
    for base, count in base_counter.items():
        if base == ref_base.upper():
            continue
        alt_freq = count / depth
        length = len(base)
        if length == 1:
            if alt_freq >= min_alt_freq:
                has_snp = True
        elif length > 1:
            if alt_freq >= min_indel_freq:
                has_indel = True
    return contig, pos, ref_base, read_name_list, base_list, read_name_dict, has_snp, has_indel, depth

def process_region(args):
    region, bam_file, ref_fasta, min_depth, min_alt_freq, min_indel_freq, min_mq, min_bq, filter_flag, samtools_execute_command = args
    mq_option = f' --min-MQ {min_mq}'
    bq_option = f' --min-BQ {min_bq}'
    flags_option = f' --excl-flags {filter_flag}'
    read_name_option = "--output-QNAME"
    mp_option = "-s"
    reference_option = f"-f {ref_fasta}"
    region_option = f"-r {region}"
    mpileup_command = f"{samtools_execute_command} mpileup {mq_option} {bq_option} {flags_option} {read_name_option} {mp_option} {reference_option} {region_option} {bam_file}"
    mpileup_process = subprocess_popen(shlex.split(mpileup_command), shell=False)
    
    # snp_candidates = []
    # indel_candidates = []
    candidate_positions = {}
    depth_dict = defaultdict(int)
    for row in mpileup_process.stdout:
        result = decode_pileup_row(row, min_alt_freq, min_indel_freq)
        if result[8] == 0:
            continue
        if result[8] >= min_depth and (result[6] or result[7]):
                candidate_positions[result[1]] = {
                    "contig": result[0],
                    "pos": result[1],
                    "ref_base": result[2],
                    "read_name_list": result[3],
                    "base_list": result[4],
                    "read_name_dict": result[5],
                    "has_snp": result[6],
                    "has_indel": result[7]
                }
        else:
            depth_dict[result[1]] = result[8]
    return candidate_positions, depth_dict


def process_candidate_sites(candidate_positions, output_file, min_depth=8, min_freq=0.15):
    contig = None
    try:
        with open(output_file, 'w') as f_out:
            for pos_key, data in candidate_positions.items():

                contig = data.get("contig")
                position = data.get("pos")
                ref_base = data.get("ref_base")
                base_list = data.get("base_list")
                # print(contig, position, ref_base, base_list)
                if None in [contig, position, ref_base, base_list]:
                    print(f"Warning: Skipping incomplete data for key {pos_key}")
                    continue
                # Process variants (both SNVs and indels)
                process_position_variants(contig, position, ref_base, base_list, f_out, min_depth, min_freq)
    except Exception as e:
        print(f"An unexpected error occurred during processing: {e}")

def process_position_variants(contig, position, ref_base, base_list, f_out, min_depth, min_freq):
    ref_base = ref_base.upper()
    if ref_base == 'N':
        return
    total_depth = len(base_list)
    if total_depth < min_depth:
        return
    raw_observations = []

    for base_char, indel_info in base_list:
        if base_char == "":
            if indel_info.startswith("+") or indel_info.startswith("-"):
                strand_indicator = 1 if indel_info[-1].isupper() else 0
                raw_observations.append((ref_base.upper() + indel_info.upper(), strand_indicator))
        else:
            if base_char == "*":
                #random get strand indicator
                strand_indicator = random.randint(0, 1)
                raw_observations.append((base_char, strand_indicator))
            else:
                strand_indicator = 1 if base_char.isupper() else 0
                raw_observations.append((base_char.upper(), strand_indicator))
                if indel_info.startswith("+") or indel_info.startswith("-"):
                    raw_observations.append((ref_base.upper() + indel_info.upper(), strand_indicator))
    allele_counts = Counter(obs[0] for obs in raw_observations)
    for up_allele_repr, alt_depth in allele_counts.items():
        alt_freq = alt_depth / total_depth if total_depth > 0 else 0
        if alt_freq < min_freq:
            continue
        vcf_ref = ""
        vcf_alt = ""
        is_indel = 0 # 0 for SNV, 1 for INDEL
        if up_allele_repr in 'ACGT': # Potential SNV
            if up_allele_repr != ref_base:
                vcf_ref = ref_base
                vcf_alt = up_allele_repr
                is_indel = 0
            else:
                continue

        elif '+' in up_allele_repr:
            try:
                base_part, ins_seq = up_allele_repr.split('+', 1)
                if base_part != ref_base:
                    print(f"")
                if not ins_seq: continue
                vcf_ref = ref_base
                vcf_alt = ref_base + ins_seq
                is_indel = 1
            except ValueError:
                continue


        elif '-' in up_allele_repr:
            try:
                base_part, del_seq = up_allele_repr.split('-', 1)
                if base_part != ref_base:
                    print(f"")
                if not del_seq: continue
                vcf_ref = ref_base + del_seq
                vcf_alt = ref_base
                is_indel = 1
            except ValueError:
                continue

        else:
            continue

        alt_pos_num = sum(obs[1] for obs in raw_observations if obs[0] == up_allele_repr)
        alt_neg_num = alt_depth - alt_pos_num
        other_depth = total_depth - alt_depth
        
        # Check for strand bias before writing
        if alt_depth >= 4:  # Only check strand bias if we have sufficient depth
            # Calculate strand bias ratio
            if alt_pos_num == 0 or alt_neg_num == 0:
                # Complete strand bias - all reads on one strand
                continue
            

            min_strand_ratio = 0.1
            pos_ratio = alt_pos_num / alt_depth
            neg_ratio = alt_neg_num / alt_depth
            
            if pos_ratio < min_strand_ratio or neg_ratio < min_strand_ratio:
                continue
        
        output_line = (
            f"{contig}_{position}_{is_indel}_{vcf_ref}_{vcf_alt}_"
            f"{total_depth}_{other_depth}_{alt_depth}_"
            f"{alt_pos_num}_{alt_neg_num}"
        )
        f_out.write(output_line + "\n")

def run(bam_file, ref_fasta, output_dir, chrom, min_depth, min_alt_freq, min_indel_freq, min_mq, min_bq, chunk_num, chunk_idx, chunk_size=None, filter_flag="2316", samtools_execute_command="samtools"):
    # Setup logging
    os.makedirs(output_dir, exist_ok=True)
    # Keep the log file name consistent, regardless of chunking method
    log_file_prefix = f"{chrom}_chunk{chunk_idx}"

    logging.basicConfig(
        filename=os.path.join(output_dir, f"{log_file_prefix}_pipeline.log"),
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Check if files exist
    if not os.path.exists(bam_file):
        error_msg = f"BAM file {bam_file} does not exist"
        logging.error(error_msg)
        raise FileNotFoundError(error_msg)

    if not os.path.exists(f"{bam_file}.bai"):
        error_msg = f"BAM index file {bam_file}.bai does not exist"
        logging.error(error_msg)
        raise FileNotFoundError(error_msg)

    if not os.path.exists(ref_fasta):
        error_msg = f"Reference fasta file {ref_fasta} does not exist"
        logging.error(error_msg)
        raise FileNotFoundError(error_msg)

    chrom_length = get_chrom_length(bam_file, chrom, samtools_execute_command)
    if not chrom_length:
        raise ValueError(f"Chromosome {chrom} not found in BAM file")

    # Calculate chunk start and end positions
    if chunk_size:
        # Chunking by fixed size
        chunk_start = 1 + (chunk_idx - 1) * chunk_size
        chunk_end = min(chunk_idx * chunk_size, chrom_length)
        # Calculate total number of chunks based on size for logging/info purposes
        total_chunks = (chrom_length + chunk_size - 1) // chunk_size
        logging.info(f"Using fixed chunk size: {chunk_size}. Processing chunk {chunk_idx} of {total_chunks}")
    else:
        # Original chunking by number of chunks
        if chunk_num <= 0:
             raise ValueError("chunk_num must be positive if chunk_size is not provided")
        calculated_chunk_size = (chrom_length + chunk_num - 1) // chunk_num  # Round up division
        chunk_start = 1 + (chunk_idx - 1) * calculated_chunk_size
        chunk_end = min(chunk_idx * calculated_chunk_size, chrom_length)
        total_chunks = chunk_num
        logging.info(f"Dividing into {chunk_num} chunks. Processing chunk {chunk_idx} of {total_chunks}")

    # Validate chunk index against calculated total chunks
    if chunk_idx < 1 or chunk_idx > total_chunks:
        raise ValueError(f"Chunk index {chunk_idx} is out of range. Must be between 1 and {total_chunks}.")

    # Check if the calculated chunk is valid (start <= end)
    if chunk_start > chunk_end:
        logging.warning(f"Chunk {chunk_idx} starts after it ends ({chunk_start} > {chunk_end}). This likely means the chunk index is too high for the given chromosome length and chunk size/number. No processing needed for this chunk.")
        return 0 # Return 0 candidates as this chunk is empty or invalid

    # Create the region for this chunk
    region = f"{chrom}:{chunk_start}-{chunk_end}"
    logging.info(f"Processing region {region}")

    # Process the specified chunk
    process_args = (region, bam_file, ref_fasta, min_depth, min_alt_freq, min_indel_freq,
                   min_mq, min_bq, filter_flag, samtools_execute_command)

    candidate_positions, depth_dict = process_region(process_args)

    output_first_candidates = os.path.join(output_dir, f"{chrom}_chunk{chunk_idx}_first_candidates")
    process_candidate_sites(candidate_positions, output_first_candidates, min_depth, min_alt_freq)

    # Keep the output pickle file name consistent
    output_filename = f"{chrom}_chunk{chunk_idx}_candidate_positions.pkl"
    with open(os.path.join(output_dir, output_filename), 'wb') as f:
        pickle.dump(candidate_positions, f, protocol=4)
    with open(os.path.join(output_dir, f"{chrom}_chunk{chunk_idx}_depth_dict.pkl"), 'wb') as f:
        pickle.dump(depth_dict, f, protocol=4)


    logging.info(f"Completed processing chunk {chunk_idx} of {total_chunks}")
    logging.info(f"Found {len(candidate_positions)} candidate positions")

    return len(candidate_positions)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get candidate sites from BAM file. Process a specific chunk for parallel execution.")
    parser.add_argument("--bam_file","-b", type=str, required=True, help="The path to the indexed BAM file (with .bai index file)")
    parser.add_argument("--ref_fasta","-r", type=str, required=True, help="The path to the reference fasta file")
    parser.add_argument("--output_dir","-o", type=str, required=True, help="The path to the output directory")

    parser.add_argument("--chrom", type=str, required=True, help="The chromosome to be processed")
    # Chunking arguments - allow specifying either number of chunks or fixed size
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--chunk_num", type=int, help="Total number of chunks to divide the chromosome into (alternative to --chunk_size)")
    group.add_argument("--chunk_size", type=int, help="Fixed size for each chunk (alternative to --chunk_num)")

    parser.add_argument("--chunk_idx", type=int, required=True, help="The index of the chunk to process (1-based)")

    parser.add_argument("--min_depth", type=int, default=10, help="The minimum depth of the candidate site")
    parser.add_argument("--min_alt_freq", type=float, default=0.2, help="The minimum alternative frequency of the candidate site")
    parser.add_argument("--min_indel_freq", type=float, default=0.2, help="The minimum indel frequency of the candidate site")
    parser.add_argument("--min_mq", type=int, default=5, help="The minimum mapping quality of the candidate site")
    parser.add_argument("--min_bq", type=int, default=0, help="The minimum base quality of the candidate site")

    args = parser.parse_args()

    # Validation for chunk_idx >= 1 is implicitly handled by the logic inside run,
    # but we add an explicit check here for clarity.
    if args.chunk_idx < 1:
        raise ValueError("Chunk index must be 1-based (>= 1)")

    # Pass chunk_size to the run function
    num_candidate_positions = run(
        args.bam_file,
        args.ref_fasta,
        args.output_dir,
        args.chrom,
        args.min_depth,
        args.min_alt_freq,
        args.min_indel_freq,
        args.min_mq,
        args.min_bq,
        args.chunk_num, # Pass chunk_num (might be None if chunk_size is used)
        args.chunk_idx,
        args.chunk_size # Pass chunk_size (might be None if chunk_num is used)
    )

    # Determine how to report completion based on which chunking method was used
    if args.chunk_size:
        print(f"Processed chunk {args.chunk_idx} (size {args.chunk_size}) for chromosome {args.chrom}")
    else:
        print(f"Processed chunk {args.chunk_idx} of {args.chunk_num} for chromosome {args.chrom}")

    print(f"Found {num_candidate_positions} candidate positions")
    