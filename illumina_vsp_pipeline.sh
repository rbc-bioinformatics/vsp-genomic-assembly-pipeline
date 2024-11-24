#!/bin/bash

: << 'EOF'
Note: Before running the pipeline you need to have your fastq files and all scripts in the same directory
EOF

# Define color variables
YELLOW='\033[1;33m'
GREEN='\033[1;32m'
RESET='\033[0m'

echo "===================================================="
START_DATE_TIME=$(date +'%Y-%m-%d %H:%M:%S')
START_TIME=$(date +%s)
echo "$START_DATE_TIME : Analysis Starts"
echo "===================================================="

### PREPROCESSING

: << 'EOF'
The pipeline begins by creating necessary directories to store intermediate and output files.
EOF

# Define directories
WORKING_DIR=$(pwd)
SCRIPTS_DIR="$WORKING_DIR/scripts"
CLEAN_READS="$WORKING_DIR/cleanReads"
FASTP_OUT_DIR="$WORKING_DIR/fastp_out"
READS_DIR="$WORKING_DIR/reads_dir"
RAW_READS_DIR="$WORKING_DIR/raw_reads"

# Create directories
mkdir -p "$SCRIPTS_DIR" "$CLEAN_READS" "$FASTP_OUT_DIR" "$READS_DIR" "$RAW_READS_DIR"

### STEP 1
echo ""
START_DATE_TIME_1=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "${YELLOW} *** $START_DATE_TIME_1 : (1/7) Concatenation Starts *** ${RESET}"
echo ""
./illumina_vsp_1_concatenation.sh
echo ""
END_DATE_TIME_1=$(date +'%Y-%m-%d %H:%M:%S')
echo " *** $END_DATE_TIME_1 : Concatenation Ends *** "
echo ""

### STEP 2
echo ""
START_DATE_TIME_2=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "${YELLOW} *** $START_DATE_TIME_2 : (2/7) QC Starts *** ${RESET}"
echo ""
./illumina_vsp_2_qc.sh
echo ""
END_DATE_TIME_2=$(date +'%Y-%m-%d %H:%M:%S')
echo " *** $END_DATE_TIME_2 : QC Ends *** "
echo ""

### STEP 3
echo ""
START_DATE_TIME_3=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "${YELLOW} *** $START_DATE_TIME_3 : (3/7) De-hosting Starts *** ${RESET}"
echo ""
./illumina_vsp_3_de-host.sh
echo ""
END_DATE_TIME_3=$(date +'%Y-%m-%d %H:%M:%S')
echo " *** $END_DATE_TIME_3 : De-hosting Ends *** "
echo ""

### STEP 4
echo ""
START_DATE_TIME_4=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "${YELLOW} *** $START_DATE_TIME_4 : (4/7) Alignment Starts *** ${RESET}"
echo ""
./illumina_vsp_4_alignment.sh
echo ""
END_DATE_TIME_4=$(date +'%Y-%m-%d %H:%M:%S')
echo " *** $END_DATE_TIME_4 : Alignment Ends *** "
echo ""

### STEP 5
echo ""
START_DATE_TIME_5=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "${YELLOW} *** $START_DATE_TIME_5 : (5/7) SAM Files Processing Starts *** ${RESET}"
echo ""
./illumina_vsp_5_process_SAM_files.sh
echo ""
END_DATE_TIME_5=$(date +'%Y-%m-%d %H:%M:%S')
echo " *** $END_DATE_TIME_5 : SAM Files Processing Ends *** "
echo ""

### STEP 6
echo ""
START_DATE_TIME_6=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "${YELLOW} *** $START_DATE_TIME_6 : (6/7) Variant Calling Starts *** ${RESET}"
echo ""
./illumina_vsp_6_variant_calling.sh
echo ""
END_DATE_TIME_6=$(date +'%Y-%m-%d %H:%M:%S')
echo " *** $END_DATE_TIME_6 : Variant Calling Ends ***"
echo ""

### STEP 7
echo ""
START_DATE_TIME_7=$(date +'%Y-%m-%d %H:%M:%S')
echo -e "${YELLOW} *** $START_DATE_TIME_7 : (7/7) Consensus Generation Starts *** ${RESET}"
echo ""
./illumina_vsp_7_consensus.sh
echo ""
END_DATE_TIME_7=$(date +'%Y-%m-%d %H:%M:%S')
echo " *** $END_DATE_TIME_7 : Consensus Generation Ends ***"
echo ""

echo "===================================================="
END_DATE_TIME=$(date +'%Y-%m-%d %H:%M:%S')
END_TIME=$(date +%s)
echo "$END_DATE_TIME : Analysis Ends"

# Calculate and display the duration of the analysis
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(( (DURATION % 3600) / 60 ))
SECONDS=$((DURATION % 60))
echo "===================================================="

echo -e "${GREEN} Analysis Total Duration: ${HOURS}h ${MINUTES}m ${SECONDS}s ${RESET}"

# Move scripts to the scripts directory
mv *sh "$SCRIPTS_DIR"