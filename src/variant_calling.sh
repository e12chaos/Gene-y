#!/bin/bash

# Usage: bash variant_calling.sh
# Exit immediately if a command exits with a non-zero status

set -e

# Load necessary modules
module load samtools
module load gatk
module load vcftools

# Input Arguments
SAMPLE_ID=$1
REFERENCE_GENOME=$2
KNOWN_SITES_VCF=$3

# Define directories
ALIGNMENT_DIR="../output/alignment"
VC_DIR="../output/VC"
VC_SORTING_DIR="$VC_DIR/VC_sorting"
VC_DEDUP_RECAL_DIR="$VC_DIR/VC_dedup_recal"
VC_BQSR_CALL_DIR="$VC_DIR/VC_bqsr_call"
VC_SORTING_OUTPUT_DIR="$VC_SORTING_DIR"
VC_DEDUP_RECAL_OUTPUT_DIR="$VC_DEDUP_RECAL_DIR"
VC_BQSR_CALL_OUTPUT_DIR="$VC_BQSR_CALL_DIR/vcf"

# Create necessary directories if they don't exist
mkdir -p $VC_SORTING_DIR
mkdir -p $VC_DEDUP_RECAL_DIR/dedup_recal_metrics
mkdir -p $VC_BQSR_CALL_OUTPUT_DIR
mkdir -p tmp
export TMPDIR=/scratch/aitken.n/data/bioinfo_fastq/tmp

# Step 1: Create Sequence Dictionary and Index Reference Genome
gatk CreateSequenceDictionary -R $REFERENCE_GENOME
samtools faidx $REFERENCE_GENOME

# Step 2: Sort SAM to BAM
INPUT_SAM="$ALIGNMENT_DIR/${SAMPLE_ID}.sam"
OUTPUT_SORTED_BAM="$VC_SORTING_DIR/${SAMPLE_ID}_sorted.bam"
gatk SortSam -I $INPUT_SAM -O $OUTPUT_SORTED_BAM -SO coordinate

# Step 3: Mark Duplicates
OUTPUT_DEDUP_BAM="$VC_DEDUP_RECAL_DIR/${SAMPLE_ID}_dedup.bam"
OUTPUT_METRICS="$VC_DEDUP_RECAL_DIR/dedup_recal_metrics/${SAMPLE_ID}_metrics.txt"
gatk --java-options "-Xmx60G" MarkDuplicates \
    -I $OUTPUT_SORTED_BAM \
    -O $OUTPUT_DEDUP_BAM \
    -M $OUTPUT_METRICS \
    --TMP_DIR $TMPDIR

# Step 4: Base Recalibration
RECAL_TABLE="$VC_DEDUP_RECAL_DIR/${SAMPLE_ID}_recal.table"
gatk BaseRecalibrator \
    -I $OUTPUT_DEDUP_BAM \
    -R $REFERENCE_GENOME \
    --known-sites $KNOWN_SITES_VCF \
    -O $RECAL_TABLE

# Step 5: Apply BQSR
OUTPUT_BQSR_BAM="$VC_BQSR_CALL_DIR/${SAMPLE_ID}_bqsr.bam"
gatk ApplyBQSR \
    -R $REFERENCE_GENOME \
    -I $OUTPUT_DEDUP_BAM \
    --bqsr-recal-file $RECAL_TABLE \
    -O $OUTPUT_BQSR_BAM

# Step 6: HaplotypeCaller
OUTPUT_VARIANTS_VCF="$VC_BQSR_CALL_OUTPUT_DIR/${SAMPLE_ID}_variants.vcf"
gatk HaplotypeCaller \
    -R $REFERENCE_GENOME \
    -I $OUTPUT_BQSR_BAM \
    -O $OUTPUT_VARIANTS_VCF

# Step 7: Filter VCF with VCFtools
OUTPUT_FILTERED_VCF="$VC_BQSR_CALL_OUTPUT_DIR/${SAMPLE_ID}_variants_filtered.vcf"
vcftools --vcf $OUTPUT_VARIANTS_VCF --minQ 30 --recode --out "${OUTPUT_FILTERED_VCF%.vcf}_filtered.vcf"

echo "Variant calling completed for sample: $SAMPLE_ID"
