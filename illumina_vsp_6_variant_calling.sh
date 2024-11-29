
#!/bin/bash

### 6. Variant calling using "iVar"

: << 'EOF'
Variant calling is the process of identifying and cataloging the 
differences between the virus of interest sequencing reads and a 
reference genome. Variant calling enables us to know the amount of 
changes occurred on the genome of interest, we get the SNPs. MNPs, 
idels e.t.c. This step is crucial for knowing/detection of new variants 
of the virus in question.
EOF

# define directories
WORKING_DIR=$(pwd)
NON_HOST_READS_DIR="$WORKING_DIR/nonHost"
MVD_REF_DIR="$WORKING_DIR/reference_genomes/mvd"

# Define the reference genome:
mvd_reference="$MVD_REF_DIR/NC_001608.3.fna"

for sorted_mapped_bam_file in "$NON_HOST_READS_DIR"/*.sorted.mapped.bam; do
  sample=$(basename "$sorted_mapped_bam_file" .sorted.mapped.bam)

  # define the prefix
  prefix="${sorted_mapped_bam_file%.sorted.mapped.bam}_variants"

  samtools mpileup --reference "${mvd_reference}" "${sorted_mapped_bam_file}" | ivar variants -r "${mvd_reference}" -p "${prefix}"

  echo "${sample} variant calling completed."
done
