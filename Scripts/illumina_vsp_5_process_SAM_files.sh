
#!/bin/bash

### 5. Processing of the mapped sequences using "Samtools"

### Sort SAM file

# define directories
WORKING_DIR=$(pwd)
NON_HOST_READS_DIR="$WORKING_DIR/nonHost"

for sam_file in "$NON_HOST_READS_DIR"/*.sam; do
  sample=$(basename "$sam_file" .sam)
  sorted_bam_file="${sam_file%.sam}.sorted.bam"
  sorted_mapped_bam_file="${sam_file%.sam}.sorted.mapped.bam"

  samtools sort "${sam_file}" > "${sorted_bam_file}"
  
  ### Discard un-mapped reads
  samtools view -F 0x04 -b "${sorted_bam_file}" > "${sorted_mapped_bam_file}"

  ### Index BAM file
  samtools index "${sorted_mapped_bam_file}"

  echo "- ${sample} processed."
done
