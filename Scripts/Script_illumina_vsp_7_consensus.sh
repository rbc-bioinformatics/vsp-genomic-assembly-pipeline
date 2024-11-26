
#!/bin/bash

### 7. Generate consensus using "iVar"

: << 'EOF'
At this step we now create the consensus fasta file which we will use for 
downstream analysis, here we call consensus for position with support from 
at least on read at that position, this is because we use metagenomic 
sequencing aproach the depth might be very low. if tilling approach was 
used for sequencing depth for calling consensus can be set to 5 or 10 
reads per position.

Generate consensus FASTA:
Optionally, you can set different parameters to define minimum thresholds 
for the consensus
EOF

# define directories
WORKING_DIR=$(pwd)
NON_HOST_READS_DIR="$WORKING_DIR/nonHost"

for sorted_mapped_bam_file in "$NON_HOST_READS_DIR"/*.sorted.mapped.bam; do
  sample=$(basename "$sorted_mapped_bam_file" .sorted.mapped.bam)

  # define variables
  consensus_file="${sorted_mapped_bam_file%.sorted.mapped.bam}.fa"

  samtools mpileup -A -Q 0 "${sorted_mapped_bam_file}" | awk '$4 >= 50' | ivar consensus -p "${consensus_file}" -q 20 -t 0.7 -m 50

  echo "${sample} consensus completed."
done
