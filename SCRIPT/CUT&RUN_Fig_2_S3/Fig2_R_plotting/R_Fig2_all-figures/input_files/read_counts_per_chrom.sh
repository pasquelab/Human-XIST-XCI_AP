#!/bin/bash

# Check if the 2 required arguments are given:
if [ $# -ne 2 ]; then
    echo "Usage: $0 <list_bam> <output_dir>"
    exit 1
fi

list_bam=$1
output_dir=$2

# Check if the output directory exists, if not, create it
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir" || { echo "Failed to create output directory: $output_dir"; exit 1; }
fi


output_file="${output_dir}/read_counts_per_chrom.tsv"


# Clear the output file if it already exists
> "$output_file"


for sample in $list_bam; do
    #sample_name=$(basename "$sample" "_downsampled.idxstats.txt")
    sample_name=$(basename "$sample" | cut -d'_' -f1,2)
    cat "$sample" | awk -v sample="$sample_name" '$1 !~ /_/{print sample "\t" $1 "\t" $3}' | grep -v 'chrM' | grep -v '*' >> "$output_file"
done
