
#!/bin/bash

### 3. Host Genome Removal with ```Bowtie2```

# define variables and directories needed
threads=4
WORKING_DIR=$(pwd)
CLEAN_READS="$WORKING_DIR/cleanReads"
HUMAN_REF_DIR="$WORKING_DIR/reference_genomes/human"

# Check if the human reference directory exists
if [[ ! -d $HUMAN_REF_DIR ]]; then
  mkdir -p "$HUMAN_REF_DIR"  # Create the directory if it doesn't exist
fi

# Define the path to the reference genome file
reference_genome="$HUMAN_REF_DIR/hg38.fa"

# define the index prefix variable
index_prefix="$HUMAN_REF_DIR/index"

# Check if the reference genome file exists
if [[ ! -f $reference_genome ]]; then
  echo "Reference genome file does not exist. Downloading the human reference genome..."
  
  # Download the human reference genome
  wget -P "$HUMAN_REF_DIR" http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

  # Unzip the genome
  gunzip "$HUMAN_REF_DIR/hg38.fa.gz"

  # Build Bowtie2 index for host genome
  bowtie2-build "$reference_genome" "$index_prefix" --threads "$threads"
fi

# Create output directory for non-host reads
NON_HOST_READS_DIR="$WORKING_DIR/nonHost"
mkdir -p "$NON_HOST_READS_DIR"

# create a directory to store all unmapped reads
UNMAPPED_READS_DIR="$NON_HOST_READS_DIR/unmapped_reads"
mkdir -p "$UNMAPPED_READS_DIR"

# Map reads to host genome and extract unmapped reads
for forward_read in "$CLEAN_READS"/*.fastp_1.fastq.gz; do
  sample=$(basename "$forward_read" .fastp_1.fastq.gz)
  reverse_read="${forward_read%.fastp_1.fastq.gz}.fastp_2.fastq.gz"

  # Create output directory for each sample
  SAMPLE_DIR="$WORKING_DIR/reference_genomes/human/${sample}"
  mkdir -p "$SAMPLE_DIR"

  # Output files
  sam_output="$SAMPLE_DIR/${sample}.sam"
  unmapped_output="$UNMAPPED_READS_DIR/${sample}_reads_unmapped.fastq"

  # Run Bowtie2
  bowtie2 -1 "$forward_read" -2 "$reverse_read" -S "$sam_output" --un-conc "$unmapped_output" --threads "$threads" -x "$index_prefix"

  echo "Mapping completed for ${sample}"

  # Print mapping statistics
  flagstat_output="$SAMPLE_DIR/${sample}.flagstat"
  samtools flagstat -@ "$threads" "$sam_output" > "${flagstat_output}"
done
