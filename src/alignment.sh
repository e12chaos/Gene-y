#!/bin/bash

# Create output directory
mkdir -p ../output/alignment

# Index reference if not already indexed
if [ ! -f ../data/hg38.fa.bwt ]; then
    bwa index ../data/hg38.fa
fi

# Run bwa mem for BM002B_A_S1_L001
bwa mem -t 12 ../data/hg38.fa ../output/trim/BM002B_A_S1_L001_R1_001_paired.fastq.gz ../output/trim/BM002B_A_S1_L001_R2_001_paired.fastq.gz > 