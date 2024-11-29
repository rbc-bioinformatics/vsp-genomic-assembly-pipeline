
#!/bin/bash

### 4. Mapping de-hosted reads to the reference sequence using "Minimap2"

: << 'EOF'
In this step the de-hosted reads are mapped to the reference 
sequence. Based on the sequence run or target virus the reference 
sequence can be downloaded from the NCBI-Virus database. The information 
on which reference (accession number) to be used for the selected virus 
can be searched in literature normally has a refseq tag in the NCBI search results.

### Mapping using minimap2 
Here is the syntax for mapping non-human reads to the reference sequence using minimap2:
EOF

# define directories
WORKING_DIR=$(pwd)
NON_HOST_READS_DIR="$WORKING_DIR/nonHost"
UNMAPPED_READS_DIR="$NON_HOST_READS_DIR/unmapped_reads"
MVD_REF_DIR="$WORKING_DIR/reference_genomes/mvd"

# Define the reference genome:
mvd_reference="$MVD_REF_DIR/NC_001608.3.fna"

for unmapped_forward_read in "$UNMAPPED_READS_DIR"/*_reads_unmapped.1.fastq; do
  sample=$(basename "$unmapped_forward_read" _reads_unmapped.1.fastq)
  unmapped_reverse_read="${unmapped_forward_read%_reads_unmapped.1.fastq}_reads_unmapped.2.fastq"

  # Define aligned SAM file output:
  sam_file="${NON_HOST_READS_DIR}/${sample}.sam"

  # Run Minimap2 command:
  minimap2 -ax sr "${mvd_reference}" "${unmapped_forward_read}" "${unmapped_reverse_read}" > "${sam_file}"
done
