# miniSNV

## Introduction
RVC is a real-time variant calling framework tailored for cycle-based NGS workflows. RVC incrementally processes partially sequenced reads and continuously updates variant evidence using a scanline-based alignment algorithm and a lightweight binomial scoring model.

---

## Installation
### 1.Environment Configuration
```
conda create -n RVC python=3.8 -y

conda activate RVC

conda install samtools bgzip parallel vcftools seqkit -y

pip install biopython

pip install pybind11

```
### 2.download RVC
```
git clone https://github.com/CuiMiao-HIT/RVC.git

cd RVC/Release
make clean && make -j 12

cd RVC/ksw2/python_bindings
pip install . (If error, please try 'python setup.py clean --all' then  'pip install .')

```

---	

## Usage
```
conda activate RVC

bash run_RVC.sh
```
### Options
**Required parameters:**  
```
	SRC_DIR               # Work-directory for distributed job.(eg : */RVC)
	BAM_FILE              # BAM file input. The input file must be samtools indexed.
	REF_FASTA             # FASTA reference file input. The input file must be samtools indexed.
	FASTQ_FILES=          # Each batch of fastq files input.(eg : ("The-first-cycles-batch-data.fastq"
                                                                    "The-second-cycles-batch-data.fastq"
                                                                    ...
                                                                    "The-last-cycles-batch-data.fastq"))
	READ_LENGTHS          # Each batch of fastq length.(eg : (130 50 50 20))
	OUTPUT_DIR            # Output-directory for distributed job.(default : */RVC/RUN_RESULT)
```
**Other parameters:**  
```
	CHRS                  # List of chromosomes to process.(default : chr1-chr22)
	THREAD_NUM            # Number of threads to use.(default : 22)
	CHUNK_SIZE            # Reference length to detect candidate in one loop.(default : 20000000)
	BATCH_SIZE            # The number of reads loaded from a FASTQ file into memory at one time.(default : 2000000)
	IDENTITY_THRESHOLD    # The consistency threshold used to filter the alignments.(default : 0.99)
	MIN_BQ                # Minimum base quality for screening candidate sites.(default : 0)
	MIN_MQ                # Minimum alignment quality for screening candidate sites.(default : 0)
	VI_MIN_DEPTH          # Minimum depth threshold for variant callers.(default : 8)
	VI_MIN_FREQ           # Minimum allele frequency threshold for variant caller.(default : 0.15)
	IV_QUAL               # Minimum variant calling quality for filtering.(default : 1)
```


## Contact
Please post on [Github Issue](https://github.com/CuiMiao-HIT/RVC/issues) or contact cuimiao@stu.hit.edu.cn.
