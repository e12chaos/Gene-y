#!/bin/bash

# Create the output directory if it doesn't exist
mkdir -p output/QC_Reports

# Run FastQC on the provided input files
for file in "$@"
do
    fastqc "$file" -o output/QC_Reports
done
