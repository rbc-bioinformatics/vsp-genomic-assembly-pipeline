# Illumina VSP Pipeline

The **Illumina Viral Surveillance Panel (VSP) Pipeline** enables reference-based genome assembly for pathogens sequenced using the **Viral Surveillance Panel v2 Kit**. This kit leverages a hybrid-capture method, allowing for the sequencing of over 200 viral pathogens, including SARS-CoV-2, influenza, arboviruses, and hepatitis viruses. Its accuracy in identifying mutations makes it highly suitable for outbreak analysis and viral evolution studies.

---

## Table of Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Workflow](#workflow)
- [Installation](#installation)
- [Usage](#usage)
- [Notes](#notes)
- [License](#license)

---

## Overview

This pipeline provides an automated bioinformatics workflow for analyzing Illumina sequencing data generated by the Viral Surveillance Panel v2 Kit. It is specifically designed to:

- Perform quality control on raw sequencing reads.
- Remove host genome sequences.
- Assemble viral genomes.
- Identify variants and generate consensus sequences.

The pipeline is optimized for hybrid-capture sequencing, which requires less read depth compared to shotgun metagenomics and offers higher mutation detection accuracy than amplicon sequencing.

---

## Prerequisites

### Tools Required

- **fastp**: For trimming and quality control of raw sequencing reads.
- **Bowtie2**: For aligning reads to a host genome and filtering out host sequences.
- **Samtools**: For handling SAM/BAM files produced by `Bowtie2`.
- **minimap2**: For aligning the unmapped reads to the reference.
- **LoFred**: For variant calling.
- **ivar**: For variant calling and consensus sequence generation.
        
### Input Data

- Raw Illumina sequencing reads in FASTQ format.

---

## System Requirements

To run the pipeline, the following system specifications are recommended:

- **Operating System**: Linux or Windows Subsystem for Linux (WSL)
- **Processor**: Core i7 or higher
- **Memory**: 16 GB RAM or more
- **Storage**: 500 GB or more
- **Tools**: Ensure the required tools (e.g., fastp, Bowtie2, Samtools, minimap2, iVar or LoFreq) are installed and accessible.

---

## Workflow

The pipeline follows a stepwise process:

1.**FASTQ Files Preprocessing**
   - Concatenates FASTQ files to prepare them for quality control.
2. **Quality Control**:
   - Removes adapters and low-quality bases using **fastp** to ensure high-quality sequencing reads.
3. **Host Genome Removal**:
   - Filters out host sequences (e.g., human genome) using **Bowtie2**, retaining only viral reads for further analysis.
4. **Genome Assembly**:
   - Maps reads to a reference genome using **minimap2**.
5. **Processing Aligned Reads**:
   - Sorts and indexes mapped reads.
6. **Variant Calling**:
   - Identifies mutations using **iVar** or **LoFreq**, enabling insights into viral evolution.
7. **Consensus Sequence Generation**:
   - Produces a consensus sequence, critical for downstream analysis like phylogenetic studies or outbreak tracking.

---

## Installation

1. Clone this repository:
    ```bash
    git clone https://github.com/rbc-bioinformatics/vsp-genomic-assembly-pipeline.git
    cd vsp-genomic-assembly-pipeline
    ```

2. Ensure all scripts have executable permissions:
    ```bash
    chmod +x *.sh
    ```

---

## Usage

Place your raw sequencing data (FASTQ files) in the same directory as the pipeline scripts and reference genomes directory. To run the entire pipeline, execute the following command:

```bash
./illumina_vsp_fastp.sh
```
Or:
```bash
bash illumina_vsp_fastp.sh
```

# Notes
Target Reference Genomes: 

Depending on your study, download and specify the appropriate reference genome from databases like NCBI-Virus. Use accession numbers with the RefSeq tag for the most accurate sequences.

Depth Considerations:

The hybrid-capture approach is well-suited for low-depth sequencing typical of metagenomics.
For tiling-based sequencing, increase the minimum depth for consensus calling to 5 or 10 reads per position.

Variant Detection:

Use iVar for standard variant detection or LoFreq for high-sensitivity calls.
Variant data helps track mutations, single nucleotide polymorphisms (SNPs), and other genome alterations.

# License
This pipeline is open-source and distributed under the MIT License. Contributions and feedback are welcome!
