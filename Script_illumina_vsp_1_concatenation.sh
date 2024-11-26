
#!/bin/bash

### 1. Concatenate fastq files

# define directories
WORKING_DIR=$(pwd)
RAW_READS_DIR="$WORKING_DIR/raw_reads"
READS_DIR="$WORKING_DIR/reads_dir"

# concatenate fastq files
for sample in $(ls *_R1_001.fastq.gz | awk -F'_S' '{print $1}' | sort | uniq); do

  cat ${sample}_S*_L*_R1_001.fastq.gz > ${sample}_R1_combined.fastq.gz
  
  cat ${sample}_S*_L*_R2_001.fastq.gz > ${sample}_R2_combined.fastq.gz

  mv ${sample}_S*_L*_R1_001.fastq.gz $RAW_READS_DIR
  mv ${sample}_S*_L*_R2_001.fastq.gz $RAW_READS_DIR
done

# move files to their respective directories
mv *_combined.fastq.gz "$READS_DIR"

