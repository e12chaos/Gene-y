#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Create output directory
mkdir -p ../output/alignment

# Define paths
REFERENCE="/Users/nikaelaaitken/gene-y/data/reference/hg38.fa"
TRIMMED_DIR="../output/trim"
ALIGNMENT_DIR="../output/alignment"

# Index reference if not already indexed
if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index "${REFERENCE}"
fi

# Loop through each R1 paired-end file
for R1 in ${TRIMMED_DIR}/*_R1_001_paired.fastq.gz; do
    # Derive corresponding R2 filename
    R2="${R1/_R1_/_R2_}"
    
    # Extract sample name
    SAMPLE=$(basename "${R1}" _R1_001_paired.fastq.gz)
    
    # Define output SAM file
    OUTPUT="${ALIGNMENT_DIR}/${SAMPLE}_aligned.sam"
    
    echo "Aligning sample: ${SAMPLE}"
    
    # Run bwa mem
    bwa mem -t 12 "${REFERENCE}" "${R1}" "${R2}" > "${OUTPUT}"
    
    echo "Alignment for ${SAMPLE} completed. Output: ${OUTPUT}"
done

echo "All alignments completed successfully."
