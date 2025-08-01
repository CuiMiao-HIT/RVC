import os
import argparse
import pickle
from collections import defaultdict, Counter
import random
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
def process_candidate_sites(candidate_file, output_file, recall_file, depth_file, min_depth=10, min_freq=0.2):
    contig = None
    try:
        with open(candidate_file, 'rb') as f:
            candidate_positions = pickle.load(f)
    except FileNotFoundError:
        return
    except Exception as e:
        return
    if len(candidate_positions) == 0:
        return
    try:
        with open(output_file, 'w') as f_out:
            idx = 0
            for pos_key, data in candidate_positions.items():

                contig = data.get("contig")
                position = data.get("pos")
                ref_base = data.get("ref_base")
                base_list = data.get("base_list")
                if None in [contig, position, ref_base, base_list]:
                    continue

                # Process variants (both SNVs and indels)
                process_position_variants(contig, position, ref_base, base_list, f_out, min_depth, min_freq)
    except IOError as e:
        print(f"")
    except Exception as e:
        print(f"")


def process_position_variants(contig, position, ref_base, base_list, f_out, min_depth, min_freq):
    ref_base = ref_base.upper()
    if ref_base == 'N':
        return
    total_depth = len(base_list)
    if total_depth < min_depth:
        print(f"{contig}\t{position}\t{position}")
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
            print(f"{contig}\t{position}\t{position}")
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
                print(f"{contig}\t{position}\t{position}")
                continue
            

            min_strand_ratio = 0.1
            pos_ratio = alt_pos_num / alt_depth
            neg_ratio = alt_neg_num / alt_depth
            
            if pos_ratio < min_strand_ratio or neg_ratio < min_strand_ratio:
                print(f"{contig}\t{position}\t{position}")
                # print(f"Warning: pos_ratio < min_strand_ratio or neg_ratio < min_strand_ratio at {contig}:{position}. Skipping.")
                continue
        
        output_line = (
            f"{contig}_{position}_{is_indel}_{vcf_ref}_{vcf_alt}_"
            f"{total_depth}_{other_depth}_{alt_depth}_"
            f"{alt_pos_num}_{alt_neg_num}"
        )
        f_out.write(output_line + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Process candidate sites to generate unified variant statistics (SNVs and INDELs)."
    )
    parser.add_argument("--candidate_file", "-c", required=True,
                        help="Path to candidate positions pickle file")
    parser.add_argument("--output_file", "-o", required=True,
                        help="Path to output file for variant statistics")
    parser.add_argument("--recall_file", "-r", required=True,
                        help="Path to recall file")
    parser.add_argument("--depth_file", "-d", required=True,
                        help="Path to depth file")
    parser.add_argument("--min_depth", "-m", type=int, default=10,
                        help="Minimum total depth required to consider a position (default: 10)")
    parser.add_argument("--min_freq", "-f", type=float, default=0.2,
                        help="Minimum allele frequency for SNVs and INDELs (default: 0.2)")

    args = parser.parse_args()
    process_candidate_sites(args.candidate_file, args.output_file, args.recall_file, args.depth_file, args.min_depth, args.min_freq)


if __name__ == "__main__":
    main() 