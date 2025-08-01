#!/bin/bash
set -e 
set -o pipefail 

# ================================================================
# CONFIGURATION PARAMETERS - MODIFY THESE AS NEEDED
# ================================================================

SRC_DIR="/home/user/cuimiao/RVC"
BAM_FILE="/home/user/cuimiao/processed_fastq/scheme2_130_50_50_20/alignments/reads_130_50_50_20_segment1_130bp.bam"
REF_FASTA="/home/user/yuxian/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
# --- FASTQ Files (in processing order) ---
    # First FASTQ is used for initial alignment/processing
    # Subsequent FASTQs are used for iterative read extension
FASTQ_FILES=(
"/home/user/cuimiao/processed_fastq/scheme2_130_50_50_20/reads_130_50_50_20_segment1_130bp.fastq"
"/home/user/cuimiao/processed_fastq/scheme2_130_50_50_20/reads_130_50_50_20_segment2_50bp.fastq"
"/home/user/cuimiao/processed_fastq/scheme2_130_50_50_20/reads_130_50_50_20_segment3_50bp.fastq"
"/home/user/cuimiao/processed_fastq/scheme2_130_50_50_20/reads_130_50_50_20_segment4_20bp.fastq"
    # Add more FASTQ files here if needed for additional extension rounds
)
# --- Read Lengths (corresponding to each FASTQ file) ---
    # First read length is for the original FASTQ file
    # Subsequent read lengths are for each extension round
    # Add more values if you have more FASTQ files
READ_LENGTHS=(130 50 50 20)  
# --- Output Configuration ---
OUTPUT_DIR="${SRC_DIR}/RUN_RESULT"


# --- Processing Parameters ---
    # Chromosomes to Process (default:chr1-chr22)(eg:chr22-->"CHRS=(22)")
CHRS=($(seq 1 22))
THREAD_NUM=22
CHUNK_SIZE=20000000 
BATCH_SIZE=2000000 
IDENTITY_THRESHOLD=0.99
    # Quality Thresholds
MIN_BQ=0
MIN_MQ=0
VI_MIN_DEPTH=8
VI_MIN_FREQ=0.15
IV_QUAL=1

# ================================================================
# END OF CONFIGURATION
# ================================================================

# --- Tool Commands ---
SAMTOOLS_CMD="samtools"
CALLER="${SRC_DIR}/Release/VRT_ caller"
# --- Others ---
CHR_PREFIX=chr
TIMING_LOG="${OUTPUT_DIR}/step_timings.log"

# Validate FASTQ files array
if [ ${#FASTQ_FILES[@]} -lt 2 ]; then
    echo "Error: At least 2 FASTQ files are required (first for initial processing, others for extension rounds)"
    exit 1
fi

# Validate READ_LENGTHS array
if [ ${#READ_LENGTHS[@]} -ne ${#FASTQ_FILES[@]} ]; then
    echo "Error: Number of READ_LENGTHS (${#READ_LENGTHS[@]}) must match number of FASTQ_FILES (${#FASTQ_FILES[@]})"
    exit 1
fi

ORI_FASTQ="${FASTQ_FILES[0]}"
NUM_EXTENSION_ROUNDS=$((${#FASTQ_FILES[@]} - 1))

echo "Configuration Summary:"
echo "- Number of FASTQ files: ${#FASTQ_FILES[@]}"
echo "- Number of extension rounds: ${NUM_EXTENSION_ROUNDS}"
echo "- Original FASTQ: ${ORI_FASTQ}"
for i in $(seq 1 ${NUM_EXTENSION_ROUNDS}); do
    echo "- Extension round ${i} FASTQ: ${FASTQ_FILES[$i]}"
    echo "- Extension round ${i} read length: ${READ_LENGTHS[$i]}"
done
echo "- Output directory: ${OUTPUT_DIR}"
echo "- Identity threshold: ${IDENTITY_THRESHOLD}"

# --- Setup Output and Logging ---
mkdir -p ${OUTPUT_DIR}
JOBLIST_DIR="${OUTPUT_DIR}/joblists"
mkdir -p ${JOBLIST_DIR}

# Create or clear TIMING_LOG
> ${TIMING_LOG}

start_time_total=$(date +%s)
last_step_time=${start_time_total}
echo "Starting Real-time Variant Calling Pipeline (Flexible Extension Mode)..." >> ${TIMING_LOG}
echo "Output Directory: ${OUTPUT_DIR}" >> ${TIMING_LOG}
echo "Chunk Size: ${CHUNK_SIZE}" >> ${TIMING_LOG}
echo "Reference: ${REF_FASTA}" >> ${TIMING_LOG}
echo "BAM File: ${BAM_FILE}" >> ${TIMING_LOG}
echo "Number of FASTQ files: ${#FASTQ_FILES[@]}" >> ${TIMING_LOG}
echo "Number of extension rounds: ${NUM_EXTENSION_ROUNDS}" >> ${TIMING_LOG}
echo "Identity threshold: ${IDENTITY_THRESHOLD}" >> ${TIMING_LOG}
echo "Timing log: ${TIMING_LOG}" >> ${TIMING_LOG}
echo "Joblists will be stored in: ${JOBLIST_DIR}" >> ${TIMING_LOG}
echo "----------------------------------------" >> ${TIMING_LOG}

# --- Tool Dependency Checks ---
if ! command -v "${SAMTOOLS_CMD}" &> /dev/null; then
    echo "Error: ${SAMTOOLS_CMD} command not found. Please install samtools." | tee -a ${TIMING_LOG}
    exit 1
fi

if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel command not found. Please install parallel." | tee -a ${TIMING_LOG}
    exit 1
fi

if ! command -v bgzip &> /dev/null; then
    echo "Error: bgzip command not found. Please install bgzip." | tee -a ${TIMING_LOG}
    exit 1
fi

if ! command -v seqkit &> /dev/null; then
    echo "Error: seqkit command not found. Please install seqkit." | tee -a ${TIMING_LOG}
    exit 1
fi

# --- Pre-calculation: Generate Job Lists ---
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pre-calculating job lists..."

# Get all chromosome lengths from idxstats in one go
IDXSTATS_OUTPUT=$(${SAMTOOLS_CMD} idxstats ${BAM_FILE})
if [ $? -ne 0 ]; then
    echo "Error: samtools idxstats failed on ${BAM_FILE}" | tee -a ${TIMING_LOG}
    exit 1
fi

if [ -z "${IDXSTATS_OUTPUT}" ]; then
    echo "Error: samtools idxstats returned no output for ${BAM_FILE}. Is the BAM file indexed and contains data?" | tee -a ${TIMING_LOG}
    exit 1
fi

# Define job list files
CHUNK_LIST="${OUTPUT_DIR}/chunk_list.txt"
JOBLIST_STEP1A="${JOBLIST_DIR}/joblist_step1A.txt"
JOBLIST_STEP1B="${JOBLIST_DIR}/joblist_step1B.txt"
JOBLIST_STEP2="${JOBLIST_DIR}/joblist_step2.txt"

# Dynamic joblist files for extension rounds
declare -a JOBLIST_GET_FASTQ_FILES
declare -a JOBLIST_EXTEND_FILES
for i in $(seq 1 ${NUM_EXTENSION_ROUNDS}); do
    JOBLIST_GET_FASTQ_FILES[$i]="${JOBLIST_DIR}/joblist_get_fastq_round${i}.txt"
    JOBLIST_EXTEND_FILES[$i]="${JOBLIST_DIR}/joblist_extend_round${i}.txt"
done

JOBLIST_STEP5A="${JOBLIST_DIR}/joblist_step5a_variant_info.txt"
JOBLIST_STEP5B="${JOBLIST_DIR}/joblist_step5b_merge_info.txt"

# Clear existing job lists
> ${JOBLIST_STEP1A}
> ${JOBLIST_STEP1B}
> ${JOBLIST_STEP2}
for i in $(seq 1 ${NUM_EXTENSION_ROUNDS}); do
    > ${JOBLIST_GET_FASTQ_FILES[$i]}
    > ${JOBLIST_EXTEND_FILES[$i]}
done
> ${JOBLIST_STEP5A}
> ${JOBLIST_STEP5B}

total_chunks_calculated=0
for chr_num in ${CHRS[@]}; do
    chr_name="${CHR_PREFIX}${chr_num}"
    chr_length=$(echo "${IDXSTATS_OUTPUT}" | awk -v chr="${chr_name}" '$1 == chr {print $2; exit}')
    if [[ -z "$chr_length" || "$chr_length" -eq 0 ]]; then
        echo "Warning: Chromosome ${chr_name} not found or length is 0 in ${BAM_FILE} idxstats. Skipping." >> ${TIMING_LOG}
        continue
    fi
    num_chunks_for_chr=$(( (chr_length + CHUNK_SIZE - 1) / CHUNK_SIZE ))
    echo "Chromosome ${chr_name}: Length=${chr_length}, Chunks=${num_chunks_for_chr}" >> ${TIMING_LOG}
    total_chunks_calculated=$(( total_chunks_calculated + num_chunks_for_chr ))
    
    # Define step directories
    STEP1_DIR="${OUTPUT_DIR}/1_get_candidate_site"
    STEP2_DIR="${OUTPUT_DIR}/2_get_read_info"
    
    # Dynamic directories for extension rounds
    declare -a GET_FASTQ_DIRS
    declare -a EXTEND_DIRS
    for i in $(seq 1 ${NUM_EXTENSION_ROUNDS}); do
        GET_FASTQ_DIRS[$i]="${OUTPUT_DIR}/$((i+2))_get_fastq_for_all"
        EXTEND_DIRS[$i]="${OUTPUT_DIR}/$((i+2))_extend_reads"
    done
    
    STEP5_DIR="${OUTPUT_DIR}/$((NUM_EXTENSION_ROUNDS+3))_variant_info_dp${VI_MIN_DEPTH}_freq${VI_MIN_FREQ}"

    # Create directories
    mkdir -p ${STEP1_DIR}/${chr_name} ${STEP1_DIR}/logs
    mkdir -p ${STEP2_DIR}/${chr_name} ${STEP2_DIR}/logs
    
    for i in $(seq 1 ${NUM_EXTENSION_ROUNDS}); do
        mkdir -p ${GET_FASTQ_DIRS[$i]}
        mkdir -p ${EXTEND_DIRS[$i]}/${chr_name}/tmp ${EXTEND_DIRS[$i]}/logs
    done
    
    mkdir -p ${STEP5_DIR}/${chr_name} ${STEP5_DIR}/logs

    for (( chunk_idx=1; chunk_idx<=num_chunks_for_chr; chunk_idx++ )); do
        chunk_name="${chr_name} ${chunk_idx}"
        echo "${chunk_name}" >> ${CHUNK_LIST}

        # Step 1A Command
        echo "/usr/bin/time -v python ${SRC_DIR}/1_get_candidate_site.py --bam_file ${BAM_FILE} --ref_fasta ${REF_FASTA} --output_dir ${STEP1_DIR}/${chr_name} --chrom ${chr_name} --min_depth 0 --min_alt_freq 0.1 --min_indel_freq 0.1 --min_mq ${MIN_MQ} --min_bq ${MIN_BQ} --chunk_size ${CHUNK_SIZE} --chunk_idx ${chunk_idx} &> ${STEP1_DIR}/logs/${chr_name}_chunk${chunk_idx}.log" >> ${JOBLIST_STEP1A}

        # Step 2 Command
        echo "/usr/bin/time -v python ${SRC_DIR}/2_get_read_info.py --bam_file ${BAM_FILE} --output_dir ${STEP2_DIR}/${chr_name} --chrom ${chr_name} --chunk_size ${CHUNK_SIZE} --chunk_idx ${chunk_idx} &> ${STEP2_DIR}/logs/${chr_name}_chunk${chunk_idx}.log" >> ${JOBLIST_STEP2}

        # Dynamic extension round commands
        for i in $(seq 1 ${NUM_EXTENSION_ROUNDS}); do
            round_num=$i  # Round numbering starts from 2
            
            # Determine input sources
            if [ $i -eq 1 ]; then
                # First extension round
                read_info_pickle="${STEP2_DIR}/${chr_name}/${chr_name}_chunk${chunk_idx}_read_info.pkl"
                candidate_path="${STEP1_DIR}/${chr_name}/"
                last_fastq="${ORI_FASTQ}"
                recall_file_option=""
            else
                # Subsequent extension rounds
                prev_extend_dir="${EXTEND_DIRS[$((i-1))]}"
                read_info_pickle="${prev_extend_dir}/${chr_name}/${chr_name}_chunk${chunk_idx}_${round_num}_read_info.pkl"
                candidate_path="${prev_extend_dir}/${chr_name}/"
                last_fastq="${GET_FASTQ_DIRS[$((i-1))]}/${chr_name}_chunk${chunk_idx}.fastq"
                recall_file_option="--recall_file ${prev_extend_dir}/${chr_name}/${chr_name}_chunk${chunk_idx}_${round_num}_recall_dict.pkl"
            fi
            
            # Extension command
            current_read_length=${READ_LENGTHS[$i]}
            echo "/usr/bin/time -v python ${SRC_DIR}/3_extend_reads.py --new_fastq ${GET_FASTQ_DIRS[$i]}/${chr_name}_chunk${chunk_idx}.fastq --last_fastq ${last_fastq} --round $((round_num+1)) --ref_fasta ${REF_FASTA} --candidate_path ${candidate_path} --read_info_pickle ${read_info_pickle} --chrom ${chr_name} --output_dir ${EXTEND_DIRS[$i]}/${chr_name} --log_level INFO --batch_size ${BATCH_SIZE} --tmp_dir ${EXTEND_DIRS[$i]}/${chr_name}/tmp --chunk_idx ${chunk_idx} ${recall_file_option} --identity_threshold ${IDENTITY_THRESHOLD} --read_length ${current_read_length} &> ${EXTEND_DIRS[$i]}/logs/${chr_name}_chunk${chunk_idx}.log" >> ${JOBLIST_EXTEND_FILES[$i]}
        done

        # Step 5a Command (Variant Info per chunk) - uses final extension round results
        final_extend_dir="${EXTEND_DIRS[${NUM_EXTENSION_ROUNDS}]}"
        final_round_num=$((NUM_EXTENSION_ROUNDS + 1))
        echo "python ${SRC_DIR}/4_variant_recall.py -c ${final_extend_dir}/${chr_name}/${chr_name}_chunk${chunk_idx}_${final_round_num}_candidate_positions.pkl -o ${STEP5_DIR}/${chr_name}/${chr_name}_chunk${chunk_idx}_variant_info.txt --recall_file ${final_extend_dir}/${chr_name}/${chr_name}_chunk${chunk_idx}_${final_round_num}_recall_dict.pkl --depth_file ${STEP1_DIR}/${chr_name}/${chr_name}_chunk${chunk_idx}_depth_dict.pkl --min_depth ${VI_MIN_DEPTH} --min_freq ${VI_MIN_FREQ} &> ${STEP5_DIR}/logs/${chr_name}_chunk${chunk_idx}_variant_info.log" >> ${JOBLIST_STEP5A}
    done

    # Step 1b Command (Merge first class variant info per chromosome)
    echo "cat ${STEP1_DIR}/${chr_name}/${chr_name}_chunk*_first_candidates > ${STEP1_DIR}/${chr_name}_first_candidates" >> ${JOBLIST_STEP1B}

    # Step 5b Command (Merge variant info per chromosome)
    echo "cat ${STEP5_DIR}/${chr_name}/${chr_name}_chunk*_variant_info.txt > ${STEP5_DIR}/${chr_name}_variant_info.txt" >> ${JOBLIST_STEP5B}

    echo "cat ${STEP5_DIR}/logs/${chr_name}_chunk*_variant_info.log > ${STEP5_DIR}/logs/${chr_name}_variant_info.bed" >> ${JOBLIST_STEP5B}
    echo "cat ${STEP5_DIR}/logs/*_variant_info.bed > ${STEP5_DIR}/logs/tmp_cat.bed" >> ${JOBLIST_STEP5B}
    echo "sort -k1.4n -k2,2n ${STEP5_DIR}/logs/tmp_cat.bed | uniq > ${STEP5_DIR}/logs/tmpcandidates.bed" >> ${JOBLIST_STEP5B}

done

echo "Total chunks to process across all specified chromosomes: ${total_chunks_calculated}" >> ${TIMING_LOG}
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished pre-calculating job lists."
echo "----------------------------------------" | tee -a ${TIMING_LOG}

# --- Step 1a: Get Candidate Sites ---
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting Step 1: Get Candidate Sites"
STEP1_DIR="${OUTPUT_DIR}/1_get_candidate_site"
parallel -j ${THREAD_NUM} --joblog ${STEP1_DIR}/parallel_get_candidate_site.log < ${JOBLIST_STEP1A}

current_time=$(date +%s)
duration=$((current_time - last_step_time))
echo "Step 1 (Get Candidate Sites) completed in: $(date -u -d @${duration} +'%H hours %M minutes %S seconds')" >> ${TIMING_LOG}
last_step_time=${current_time}
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished Step 1."
echo "----------------------------------------" | tee -a ${TIMING_LOG}

# Part 1b: Merge first class variant info per chromosome
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Merging first class variant info per chromosome..."
bash ${JOBLIST_STEP1B}
    # --- Variant Calling ---
STEP1b_DIR="${OUTPUT_DIR}/1_first_variant_calling_dp${VI_MIN_DEPTH}_freq${VI_MIN_FREQ}"
mkdir -p ${STEP1b_DIR}
parallel -j ${THREAD_NUM} --joblog ${STEP1b_DIR}/parallel_variant_calling.log "${CALLER} --chr ${CHR_PREFIX}{1} \
--qual ${IV_QUAL} \
--fin_file ${STEP1_DIR}/${CHR_PREFIX}{1}_first_candidates \
--bed_file ${STEP1b_DIR}/${CHR_PREFIX}{1}_candidatepos.txt \
--fout_vcf ${STEP1b_DIR}/${CHR_PREFIX}{1}.vcf" ::: ${CHRS[@]}

    # --- VCF Merge and Sort ---
MERGE_OUTPUT_DIR1b="${STEP1b_DIR}" 
OUTPUT_FILE1b="${OUTPUT_DIR}/First_part.vcf"
COMPRESSED_FILE1b="${OUTPUT_FILE1b}.gz"
TEMP_DIR=$(mktemp -d -p ${OUTPUT_DIR})

VCF_FILES1b=()
valid_chrs_for_merge1b=()
for chr_num in ${CHRS[@]}; do
    chr_name="${CHR_PREFIX}${chr_num}"
    vcf_file="${MERGE_OUTPUT_DIR1b}/${chr_name}.vcf"
    if [ -f "$vcf_file" ] && [ -s "$vcf_file" ]; then
        VCF_FILES1b+=("$vcf_file")
        valid_chrs_for_merge1b+=("$chr_name")
    else
         echo "Warning: First VCF file ${vcf_file} not found or empty. Skipping from merge." >> ${TIMING_LOG}
    fi
done
    # Create standard VCF header
cat > "${TEMP_DIR}/header.vcf" << EOF
##fileformat=VCFv4.2
##fileDate=$(date +%Y%m%d)
##source=run_realtime_v2_identity_flexible.sh_MergeScript
##reference=${REF_FASTA}
EOF
    # Add contig lines only for chromosomes being merged
for chr_name in ${valid_chrs_for_merge1b[@]}; do
    chr_length=$(echo "${IDXSTATS_OUTPUT}" | awk -v chr="${chr_name}" '$1 == chr {print $2; exit}')
    if [[ -n "$chr_length" && "$chr_length" -gt 0 ]]; then
         echo "##contig=<ID=${chr_name},length=${chr_length}>" >> "${TEMP_DIR}/header.vcf"
    else
         echo "##contig=<ID=${chr_name}>" >> "${TEMP_DIR}/header.vcf"
    fi
done
cat >> "${TEMP_DIR}/header.vcf" << EOF
##INFO=<ID=PASS,Number=0,Type=Flag,Description="Passed variant filter">
##INFO=<ID=LowQual,Number=0,Type=Flag,Description="Low quality variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type (SNP or INDEL)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
EOF
    # Process VCF files for merging
TMP_LIST1b="${TEMP_DIR}/vcf_list.txt"
> ${TMP_LIST1b}
processed_files_count1b=0
for file in "${VCF_FILES1b[@]}"; do
    base_name=$(basename "$file")
    out_file="${TEMP_DIR}/${base_name}.processed"
    if grep -q "^##fileformat=VCF" "$file"; then
         cp "$file" "$out_file"
    else
         cp "${TEMP_DIR}/header.vcf" "$out_file"
         grep -v "^#" "$file" >> "$out_file"
    fi
    
    echo "$out_file" >> "$TMP_LIST1b"
    processed_files_count1b=$((processed_files_count1b + 1))
done
    # Attempt merging using available tools
merge_success1b=false
if command -v vcf-concat &> /dev/null && command -v vcf-sort &> /dev/null; then
    vcf-concat $(cat "$TMP_LIST1b") | vcf-sort -c > "$OUTPUT_FILE1b" && \
    merge_success1b=true
elif command -v bcftools &> /dev/null; then
    TEMP_BCF_OUT="${TEMP_DIR}/merged.bcf"
    bcftools concat --allow-overlaps -f "$TMP_LIST1b" -o "$TEMP_BCF_OUT" -O b && \
    bcftools sort "$TEMP_BCF_OUT" -o "$OUTPUT_FILE1b" -O v && \
    merge_success1b=true
else
    cp "${TEMP_DIR}/header.vcf" "$OUTPUT_FILE1b"
    grep -hv "^#" $(cat "$TMP_LIST1b") | sort -k1,1V -k2,2n >> "$OUTPUT_FILE1b" && \
    merge_success1b=true
fi

if ! ${merge_success1b}; then
    rm -rf "$TEMP_DIR"
    exit 1
fi

# --- Add SVTYPE annotation ---
ANNOTATED_FILE1b="${OUTPUT_FILE1b}.annotated"
awk 'BEGIN {OFS="\t"}
/^##/ {print; next}
/^#CHROM/ { 
    if ($0 !~ /##INFO=<ID=SVTYPE/) {
        print "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Variant type (SNP or INDEL)\">"
    }
    print; next
}
{ 
    ref=$4; alt=$5; info=$8; svtype="SNP";
    if (length(ref) != 1) {
        svtype="INDEL";
    } else {
        n = split(alt, alleles, ",");
        for (i=1; i<=n; i++) {
            if (length(alleles[i]) != 1) {
                svtype="INDEL";
                break;
            }
        }
    }
    if (info == "." || info == "") {
        $8 = "SVTYPE=" svtype
    } else if (info !~ /SVTYPE=/) {
        $8 = info ";SVTYPE=" svtype
    } else {
        $8 = info
    }
    $1=$1;
    print $0
}' "$OUTPUT_FILE1b" > "$ANNOTATED_FILE1b"

if [ $? -eq 0 ] && [ -s "$ANNOTATED_FILE1b" ]; then
    mv "$ANNOTATED_FILE1b" "$OUTPUT_FILE1b"
else
    rm -f "$ANNOTATED_FILE1b"
fi
# Compress with bgzip and index with tabix
compress_index_success1b=false
if command -v bgzip &> /dev/null; then
    bgzip -f "$OUTPUT_FILE1b" && COMPRESSED_FILE1b="${OUTPUT_FILE1b}.gz"
    if command -v tabix &> /dev/null && [ -f "$COMPRESSED_FILE1b" ]; then
        tabix -p vcf "$COMPRESSED_FILE1b" && compress_index_success1b=true
    else
        echo "bgzip succeeded but tabix failed or not found. Index not created." >> ${TIMING_LOG}
        compress_index_success1b=true
    fi
elif command -v gzip &> /dev/null; then
     gzip -f "$OUTPUT_FILE1b" && COMPRESSED_FILE1b="${OUTPUT_FILE1b}.gz"
     compress_index_success1b=true
else
    COMPRESSED_FILE1b="$OUTPUT_FILE1b"
    compress_index_success1b=true
fi

if ! ${compress_index_success1b}; then
    echo "Error during compression/indexing." >> ${TIMING_LOG}
fi
rm -rf "$TEMP_DIR"



# --- Step 2: Get Read Info ---
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting Step 2: Get Read Info"
STEP2_DIR="${OUTPUT_DIR}/2_get_read_info"
parallel -j ${THREAD_NUM} --joblog ${STEP2_DIR}/parallel_get_read_info.log < ${JOBLIST_STEP2}

current_time=$(date +%s)
duration=$((current_time - last_step_time))
echo "Step 2 (Get Read Info) completed in: $(date -u -d @${duration} +'%H hours %M minutes %S seconds')" >> ${TIMING_LOG}
last_step_time=${current_time}
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished Step 2."
echo "----------------------------------------" | tee -a ${TIMING_LOG}

# --- Dynamic Extension Rounds ---
# GET_FASTQ_THREAD_NUM=16
# GET_FASTQ_BATCH_SIZE=10000
# GET_FASTQ_BUFFER_SIZE=10000
# GET_FASTQ_QUEUE_SIZE=1000

for i in $(seq 1 ${NUM_EXTENSION_ROUNDS}); do
    step_num=$((i + 2))
    
    # Get FASTQ for all
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting Step ${step_num}.0: Get Fastq for all (Extension Round ${i})"
    GET_FASTQ_DIR="${GET_FASTQ_DIRS[$i]}"
    
    # Determine pickle directory and suffix
    if [ $i -eq 1 ]; then
        pickle_dir="${STEP2_DIR}"
        pickle_suffix="_read_ids.pkl"
    else
        pickle_dir="${EXTEND_DIRS[$((i-1))]}"
        pickle_suffix="_${i}_read_ids.pkl"
    fi
    echo "spliting fastq for all:"
    echo "python ${SRC_DIR}/get_fastq_for_all.py \
    --input_fastq ${FASTQ_FILES[$i]} \
    --output_dir ${GET_FASTQ_DIR} \
    --chunk_list ${CHUNK_LIST} \
    --pickle_dir ${pickle_dir} \
    --pickle_suffix ${pickle_suffix} \
    --threads 1" 
    python ${SRC_DIR}/get_fastq_for_all.py \
    --input_fastq ${FASTQ_FILES[$i]} \
    --output_dir ${GET_FASTQ_DIR} \
    --chunk_list ${CHUNK_LIST} \
    --pickle_dir ${pickle_dir} \
    --pickle_suffix ${pickle_suffix} \
    --threads 1

    # echo "[$(date '+%Y-%m-%d %H:%M:%S')] Indexing FASTQ files in ${GET_FASTQ_DIR}..." | tee -a ${TIMING_LOG}
    # find ${GET_FASTQ_DIR} -name "*.fastq" | parallel -j ${THREAD_NUM} --joblog ${GET_FASTQ_DIR}/parallel_bgzip.log "bgzip {}"
    # find ${GET_FASTQ_DIR} -name "*.fastq.gz" | parallel -j ${THREAD_NUM} --joblog ${GET_FASTQ_DIR}/parallel_fqidx.log "samtools fqidx {}"
    # echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished indexing FASTQ files." | tee -a ${TIMING_LOG}

    current_time=$(date +%s)
    duration=$((current_time - last_step_time))
    echo "Step ${step_num}.0 (Get Fastq Round ${i}) completed in: $(date -u -d @${duration} +'%H hours %M minutes %S seconds')" >> ${TIMING_LOG}
    last_step_time=${current_time}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished Step ${step_num}.0."
    echo "----------------------------------------" | tee -a ${TIMING_LOG}

    # Extend Reads
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting Step ${step_num}: Extend Reads (Round ${i})"
    EXTEND_DIR="${EXTEND_DIRS[$i]}"
    parallel -j ${THREAD_NUM} --joblog ${EXTEND_DIR}/parallel_extend_${i}.log < ${JOBLIST_EXTEND_FILES[$i]}

    current_time=$(date +%s)
    duration=$((current_time - last_step_time))
    echo "Step ${step_num} (Extend Reads Round ${i}) completed in: $(date -u -d @${duration} +'%H hours %M minutes %S seconds')" >> ${TIMING_LOG}
    last_step_time=${current_time}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished Step ${step_num}."
    echo "----------------------------------------" | tee -a ${TIMING_LOG}
done

# --- Step 5: Get Variant Info ---
variant_step_num=$((NUM_EXTENSION_ROUNDS + 3))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting Step ${variant_step_num}: Get Variant Info"
STEP5_DIR="${OUTPUT_DIR}/${variant_step_num}_variant_info_dp${VI_MIN_DEPTH}_freq${VI_MIN_FREQ}"

# Part 5a: Run variant info script per chunk
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running variant info script per chunk..."
parallel -j ${THREAD_NUM} --joblog ${STEP5_DIR}/parallel_variant_info.log < ${JOBLIST_STEP5A}

# Part 5b: Merge variant info per chromosome
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Merging variant info per chromosome..."
bash ${JOBLIST_STEP5B}

current_time=$(date +%s)
duration=$((current_time - last_step_time))
echo "Step ${variant_step_num} (Get Variant Info) completed in: $(date -u -d @${duration} +'%H hours %M minutes %S seconds')" >> ${TIMING_LOG}
last_step_time=${current_time}
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished Step ${variant_step_num}."
echo "----------------------------------------" | tee -a ${TIMING_LOG}

# --- Step 6: Variant Calling ---
calling_step_num=$((NUM_EXTENSION_ROUNDS + 4))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting Step ${calling_step_num}: Variant Calling"
STEP6_DIR="${OUTPUT_DIR}/${calling_step_num}_variant_calling_dp${VI_MIN_DEPTH}_freq${VI_MIN_FREQ}_new_caller"
mkdir -p ${STEP6_DIR}

parallel -j ${THREAD_NUM} --joblog ${STEP6_DIR}/parallel_variant_calling.log "${CALLER} --chr ${CHR_PREFIX}{1} \
--qual ${IV_QUAL} \
--fin_file ${STEP5_DIR}/${CHR_PREFIX}{1}_variant_info.txt \
--bed_file ${STEP6_DIR}/${CHR_PREFIX}{1}_candidatepos.txt \
--fout_vcf ${STEP6_DIR}/${CHR_PREFIX}{1}.vcf" ::: ${CHRS[@]}

current_time=$(date +%s)
duration=$((current_time - last_step_time))
echo "Step ${calling_step_num} (Variant Calling) completed in: $(date -u -d @${duration} +'%H hours %M minutes %S seconds')" >> ${TIMING_LOG}
last_step_time=${current_time}
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished Step ${calling_step_num}."
echo "----------------------------------------" | tee -a ${TIMING_LOG}

cat ${STEP5_DIR}/logs/tmpcandidates.bed ${STEP6_DIR}/*_candidatepos.txt | sort -k1.4n -k2,2n | uniq> ${OUTPUT_DIR}/candidates.bed

# --- VCF Merge and Sort ---
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting VCF Merge and Sort"
MERGE_OUTPUT_DIR="${STEP6_DIR}" 
OUTPUT_FILE="${OUTPUT_DIR}/output.vcf"
COMPRESSED_FILE="${OUTPUT_FILE}.gz"
TEMP_DIR=$(mktemp -d -p ${OUTPUT_DIR})

VCF_FILES=()
valid_chrs_for_merge=()
for chr_num in ${CHRS[@]}; do
    chr_name="${CHR_PREFIX}${chr_num}"
    vcf_file="${MERGE_OUTPUT_DIR}/${chr_name}.vcf"
    if [ -f "$vcf_file" ] && [ -s "$vcf_file" ]; then
        VCF_FILES+=("$vcf_file")
        valid_chrs_for_merge+=("$chr_name")
    else
         echo "Warning: VCF file ${vcf_file} not found or empty. Skipping from merge." >> ${TIMING_LOG}
    fi
done

if [ ${#VCF_FILES[@]} -eq 0 ]; then
    echo "Error: No valid VCF files found to merge. Check Step ${calling_step_num} output." | tee -a ${TIMING_LOG}
    rm -rf "$TEMP_DIR"
    exit 1
fi

echo "Merging ${#VCF_FILES[@]} VCF files into ${OUTPUT_FILE}" >> ${TIMING_LOG}
echo "Files to merge: ${VCF_FILES[@]}" >> ${TIMING_LOG}

# Create standard VCF header
cat > "${TEMP_DIR}/header.vcf" << EOF
##fileformat=VCFv4.2
##fileDate=$(date +%Y%m%d)
##source=run_realtime_v2_identity_flexible.sh_MergeScript
##reference=${REF_FASTA}
EOF

# Add contig lines only for chromosomes being merged
for chr_name in ${valid_chrs_for_merge[@]}; do
    chr_length=$(echo "${IDXSTATS_OUTPUT}" | awk -v chr="${chr_name}" '$1 == chr {print $2; exit}')
    if [[ -n "$chr_length" && "$chr_length" -gt 0 ]]; then
         echo "##contig=<ID=${chr_name},length=${chr_length}>" >> "${TEMP_DIR}/header.vcf"
    else
         echo "##contig=<ID=${chr_name}>" >> "${TEMP_DIR}/header.vcf"
    fi
done

cat >> "${TEMP_DIR}/header.vcf" << EOF
##INFO=<ID=PASS,Number=0,Type=Flag,Description="Passed variant filter">
##INFO=<ID=LowQual,Number=0,Type=Flag,Description="Low quality variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type (SNP or INDEL)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
EOF

# Process VCF files for merging
echo "Preparing VCF files for merging..." >> ${TIMING_LOG}
TMP_LIST="${TEMP_DIR}/vcf_list.txt"
> ${TMP_LIST}
processed_files_count=0
for file in "${VCF_FILES[@]}"; do
    base_name=$(basename "$file")
    out_file="${TEMP_DIR}/${base_name}.processed"
    if grep -q "^##fileformat=VCF" "$file"; then
         cp "$file" "$out_file"
    else
         echo "Adding header to ${base_name}" >> ${TIMING_LOG}
         cp "${TEMP_DIR}/header.vcf" "$out_file"
         grep -v "^#" "$file" >> "$out_file"
    fi
    
    echo "$out_file" >> "$TMP_LIST"
    # bgzip -f "$out_file"              # 压缩
    # out_file_bgz="${out_file}.gz"     # 压缩后的文件名
    # echo "$out_file_bgz" >> "$TMP_LIST"
    processed_files_count=$((processed_files_count + 1))
done

echo "Prepared ${processed_files_count} VCF files for merging." >> ${TIMING_LOG}

# Attempt merging using available tools
merge_success=false
if command -v vcf-concat &> /dev/null && command -v vcf-sort &> /dev/null; then
    echo "Using vcf-concat and vcf-sort..." >> ${TIMING_LOG}
    vcf-concat $(cat "$TMP_LIST") | vcf-sort -c > "$OUTPUT_FILE" && \
    merge_success=true
    if ! ${merge_success}; then echo "vcf-concat/sort failed." >> ${TIMING_LOG}; fi
elif command -v bcftools &> /dev/null; then
    echo "Using bcftools for merging and sorting..." >> ${TIMING_LOG}
    TEMP_BCF_OUT="${TEMP_DIR}/merged.bcf"
    bcftools concat --allow-overlaps -f "$TMP_LIST" -o "$TEMP_BCF_OUT" -O b && \
    bcftools sort "$TEMP_BCF_OUT" -o "$OUTPUT_FILE" -O v && \
    merge_success=true
    if ! ${merge_success}; then echo "bcftools merge failed." >> ${TIMING_LOG}; fi
else
    echo "Using basic grep and sort (less robust)..." >> ${TIMING_LOG}
    cp "${TEMP_DIR}/header.vcf" "$OUTPUT_FILE"
    grep -hv "^#" $(cat "$TMP_LIST") | sort -k1,1V -k2,2n >> "$OUTPUT_FILE" && \
    merge_success=true
    if ! ${merge_success}; then echo "grep/sort failed." >> ${TIMING_LOG}; fi
fi

if ! ${merge_success}; then
    echo "Error: VCF merging failed. Check logs in ${STEP6_DIR} and ${TEMP_DIR}." | tee -a ${TIMING_LOG}
    rm -rf "$TEMP_DIR"
    exit 1
fi

echo "VCF Merging completed." >> ${TIMING_LOG}

# --- Add SVTYPE annotation ---
echo "Annotating variants with SVTYPE (SNP/INDEL)..." >> ${TIMING_LOG}
ANNOTATED_FILE="${OUTPUT_FILE}.annotated"
awk 'BEGIN {OFS="\t"}
/^##/ {print; next}
/^#CHROM/ { 
    if ($0 !~ /##INFO=<ID=SVTYPE/) {
        print "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Variant type (SNP or INDEL)\">"
    }
    print; next
}
{ 
    ref=$4; alt=$5; info=$8; svtype="SNP";
    if (length(ref) != 1) {
        svtype="INDEL";
    } else {
        n = split(alt, alleles, ",");
        for (i=1; i<=n; i++) {
            if (length(alleles[i]) != 1) {
                svtype="INDEL";
                break;
            }
        }
    }
    if (info == "." || info == "") {
        $8 = "SVTYPE=" svtype
    } else if (info !~ /SVTYPE=/) {
        $8 = info ";SVTYPE=" svtype
    } else {
        $8 = info
    }
    $1=$1;
    print $0
}' "$OUTPUT_FILE" > "$ANNOTATED_FILE"

if [ $? -eq 0 ] && [ -s "$ANNOTATED_FILE" ]; then
    mv "$ANNOTATED_FILE" "$OUTPUT_FILE"
    echo "SVTYPE Annotation successful." >> ${TIMING_LOG}
else
    echo "Error during SVTYPE annotation. Original merged file kept." | tee -a ${TIMING_LOG}
    rm -f "$ANNOTATED_FILE"
fi

# Clean up temporary files and joblists
rm -rf "$TEMP_DIR"
rm -rf "$JOBLIST_DIR"

# Compress with bgzip and index with tabix
compress_index_success=false
if command -v bgzip &> /dev/null; then
    echo "Compressing with bgzip..." >> ${TIMING_LOG}
    bgzip -f "$OUTPUT_FILE" && COMPRESSED_FILE="${OUTPUT_FILE}.gz"
    if command -v tabix &> /dev/null && [ -f "$COMPRESSED_FILE" ]; then
        echo "Creating tabix index..." >> ${TIMING_LOG}
        tabix -p vcf "$COMPRESSED_FILE" && compress_index_success=true
    else
        echo "bgzip succeeded but tabix failed or not found. Index not created." >> ${TIMING_LOG}
        compress_index_success=true
    fi
elif command -v gzip &> /dev/null; then
     echo "bgzip not found. Compressing with gzip..." >> ${TIMING_LOG}
     gzip -f "$OUTPUT_FILE" && COMPRESSED_FILE="${OUTPUT_FILE}.gz"
     echo "Index cannot be created without tabix." >> ${TIMING_LOG}
     compress_index_success=true
else
    echo "Neither bgzip nor gzip found. Output file left uncompressed." >> ${TIMING_LOG}
    COMPRESSED_FILE="$OUTPUT_FILE"
    compress_index_success=true
fi

if ! ${compress_index_success}; then
    echo "Error during compression/indexing." >> ${TIMING_LOG}
fi

# --- Total Time Calculation ---
end_time_total=$(date +%s)
duration_total=$((end_time_total - start_time_total))
echo "----------------------------------------" >> ${TIMING_LOG}
echo "Pipeline finished." | tee -a ${TIMING_LOG}
echo "Total execution time: $(date -u -d @${duration_total} +'%H hours %M minutes %S seconds')" >> ${TIMING_LOG}
echo "Final output file: ${COMPRESSED_FILE}" >> ${TIMING_LOG}
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pipeline finished."

echo "Process completed."
echo "Final VCF file: ${COMPRESSED_FILE}" 





echo "${COMPRESSED_FILE}"
/home/user/cuimiao/software/rtg-tools/dist/rtg-tools-3.12.1-32d4c2d/rtg vcfeval --baseline=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz --bed-regions=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c ${COMPRESSED_FILE} -o ${OUTPUT_DIR}/rtg150_50_50_sp -t /home/user/cuimiao/HG/hg002_benchmark/sdf/ --squash-ploidy
/home/user/cuimiao/software/rtg-tools/dist/rtg-tools-3.12.1-32d4c2d/rtg vcfeval --baseline=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz --bed-regions=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed  -c ${COMPRESSED_FILE} -o ${OUTPUT_DIR}/rtg150_50_50_ -t /home/user/cuimiao/HG/hg002_benchmark/sdf/

bcftools view -v snps ${COMPRESSED_FILE} -Oz -o ${OUTPUT_DIR}/150_50_50snp.vcf.gz
tabix ${OUTPUT_DIR}/150_50_50snp.vcf.gz

bcftools view -v indels ${COMPRESSED_FILE} -Oz -o ${OUTPUT_DIR}/150_50_50indel.vcf.gz
tabix ${OUTPUT_DIR}/150_50_50indel.vcf.gz

echo "snp:"
/home/user/cuimiao/software/rtg-tools/dist/rtg-tools-3.12.1-32d4c2d/rtg vcfeval --baseline=/home/user/cuimiao/HG/hg002_benchmark/hg002_spl_snp_m.vcf.gz --bed-regions=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c ${OUTPUT_DIR}/150_50_50snp.vcf.gz -o ${OUTPUT_DIR}/rtg150_50_50_spsnp -t /home/user/cuimiao/HG/hg002_benchmark/sdf/ --squash-ploidy
/home/user/cuimiao/software/rtg-tools/dist/rtg-tools-3.12.1-32d4c2d/rtg vcfeval --baseline=/home/user/cuimiao/HG/hg002_benchmark/hg002_spl_snp_m.vcf.gz --bed-regions=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c ${OUTPUT_DIR}/150_50_50snp.vcf.gz -o ${OUTPUT_DIR}/rtg150_50_50_snp -t /home/user/cuimiao/HG/hg002_benchmark/sdf/

echo "indel:"
/home/user/cuimiao/software/rtg-tools/dist/rtg-tools-3.12.1-32d4c2d/rtg vcfeval --baseline=/home/user/cuimiao/HG/hg002_benchmark/HG002_indel.vcf.gz --bed-regions=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c ${OUTPUT_DIR}/150_50_50indel.vcf.gz -o ${OUTPUT_DIR}/rtg150_50_50_spind -t /home/user/cuimiao/HG/hg002_benchmark/sdf/ --squash-ploidy
/home/user/cuimiao/software/rtg-tools/dist/rtg-tools-3.12.1-32d4c2d/rtg vcfeval --baseline=/home/user/cuimiao/HG/hg002_benchmark/HG002_indel.vcf.gz --bed-regions=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c ${OUTPUT_DIR}/150_50_50indel.vcf.gz -o ${OUTPUT_DIR}/rtg150_50_50_ind -t /home/user/cuimiao/HG/hg002_benchmark/sdf/


# echo "First_part:"
# /home/user/cuimiao/software/rtg-tools/dist/rtg-tools-3.12.1-32d4c2d/rtg vcfeval --baseline=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz --bed-regions=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c ${OUTPUT_DIR}/First_part.vcf.gz -o ${OUTPUT_DIR}/rtgFirst_part_sp -t /home/user/cuimiao/HG/hg002_benchmark/sdf/ --squash-ploidy
# /home/user/cuimiao/software/rtg-tools/dist/rtg-tools-3.12.1-32d4c2d/rtg vcfeval --baseline=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz --bed-regions=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed  -c ${OUTPUT_DIR}/output.vcf.gz -o ${OUTPUT_DIR}/rtgFirst_part_ -t /home/user/cuimiao/HG/hg002_benchmark/sdf/

# bcftools view -v snps ${OUTPUT_DIR}/First_part.vcf.gz -Oz -o ${OUTPUT_DIR}/First_partsnp.vcf.gz
# tabix ${OUTPUT_DIR}/First_partsnp.vcf.gz

# bcftools view -v indels ${OUTPUT_DIR}/First_part.vcf.gz -Oz -o ${OUTPUT_DIR}/First_partindel.vcf.gz
# tabix ${OUTPUT_DIR}/First_partindel.vcf.gz

# echo "snp:"
# /home/user/cuimiao/software/rtg-tools/dist/rtg-tools-3.12.1-32d4c2d/rtg vcfeval --baseline=/home/user/cuimiao/HG/hg002_benchmark/hg002_spl_snp_m.vcf.gz --bed-regions=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c ${OUTPUT_DIR}/First_partsnp.vcf.gz -o ${OUTPUT_DIR}/rtgFirst_part_spsnp -t /home/user/cuimiao/HG/hg002_benchmark/sdf/ --squash-ploidy
# /home/user/cuimiao/software/rtg-tools/dist/rtg-tools-3.12.1-32d4c2d/rtg vcfeval --baseline=/home/user/cuimiao/HG/hg002_benchmark/hg002_spl_snp_m.vcf.gz --bed-regions=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c ${OUTPUT_DIR}/First_partsnp.vcf.gz -o ${OUTPUT_DIR}/rtgFirst_part_snp -t /home/user/cuimiao/HG/hg002_benchmark/sdf/

# echo "indel:"
# /home/user/cuimiao/software/rtg-tools/dist/rtg-tools-3.12.1-32d4c2d/rtg vcfeval --baseline=/home/user/cuimiao/HG/hg002_benchmark/HG002_indel.vcf.gz --bed-regions=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c ${OUTPUT_DIR}/First_partindel.vcf.gz -o ${OUTPUT_DIR}/rtgFirst_part_spind -t /home/user/cuimiao/HG/hg002_benchmark/sdf/ --squash-ploidy
# /home/user/cuimiao/software/rtg-tools/dist/rtg-tools-3.12.1-32d4c2d/rtg vcfeval --baseline=/home/user/cuimiao/HG/hg002_benchmark/HG002_indel.vcf.gz --bed-regions=/home/user/cuimiao/HG/hg002_benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -c ${OUTPUT_DIR}/First_partindel.vcf.gz -o ${OUTPUT_DIR}/rtgFirst_part_ind -t /home/user/cuimiao/HG/hg002_benchmark/sdf/

