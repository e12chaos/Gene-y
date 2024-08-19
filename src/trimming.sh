#!/bin/bash

# Create output directory if it doesn't exist
mkdir -p output/trim

# Adapter file path is passed as the first argument
ADAPTERS_PATH="$1"
shift  # Remove the first argument so that the remaining arguments are the input files

# Loop through input files and perform trimming
for file in "$@"
do
    base=$(basename "$file" "_R1_001.fastq.gz")
    R1="$file"
    R2="${file/R1/R2}"
    
    # Output files
    paired_R1="output/trim/${base}_R1_001_paired.fastq.gz"
    unpaired_R1="output/trim/${base}_R1_001_unpaired.fastq.gz"
    paired_R2="output/trim/${base}_R2_001_paired.fastq.gz"
    unpaired_R2="output/trim/${base}_R2_001_unpaired.fastq.gz"
    
    # Check if input files are not empty
    if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then
        echo "Error: One or both input files are empty or do not exist: $R1, $R2"
        continue
    fi
    
    # Run Trimmomatic with -phred33
    echo "Running Trimmomatic with -phred33..."
    trimmomatic PE -phred33 "$R1" "$R2" \
        "$paired_R1" "$unpaired_R1" \
        "$paired_R2" "$unpaired_R2" \
        ILLUMINACLIP:"$ADAPTERS_PATH":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    # Check if the trimming was successful
    if [ $? -ne 0 ]; then
        echo "Trimmomatic failed with -phred33, trying with -phred64..."
        
        # Retry with -phred64
        trimmomatic PE -phred64 "$R1" "$R2" \
            "$paired_R1" "$unpaired_R1" \
            "$paired_R2" "$unpaired_R2" \
            ILLUMINACLIP:"$ADAPTERS_PATH":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

        if [ $? -ne 0 ]; then
            echo "Error: Trimmomatic failed with both -phred33 and -phred64 for $file"
        else
            echo "Trimmomatic succeeded with -phred64 for $file"
        fi
    else
        echo "Trimmomatic succeeded with -phred33 for $file"
    fi
done
