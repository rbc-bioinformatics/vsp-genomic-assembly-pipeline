#!/bin/bash

### 2. Quality Control
# Trim adapters and filter low-quality bases using "fastp" for high-quality data processing.

# define number of threads to use
threads=4

# Define directories
WORKING_DIR=$(pwd)
READS_DIR="$WORKING_DIR/reads_dir"
CLEAN_READS="$WORKING_DIR/cleanReads"
FASTP_OUT_DIR="$WORKING_DIR/fastp_out"

./illumina_vsp_fastp.sh -i "$READS_DIR" -o "$FASTP_OUT_DIR" -c "$threads"

: << 'EOF'
- ```-i``` specifies the input reads directory.
- ```-o``` specifies the output directory for cleaned reads.
- ```-c``` specifies the number of threads to use.
EOF

# After running quality control, move and compress clean reads to the designated directory
mv "$FASTP_OUT_DIR"/*/*.fastq "$CLEAN_READS"
gzip "$CLEAN_READS"/*.fastq
chmod +rwx "$CLEAN_READS"/*.fastq.gz