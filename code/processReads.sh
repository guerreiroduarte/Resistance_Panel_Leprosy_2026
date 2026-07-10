#!/usr/bin/env bash

# Define optional varibles
MIN_LENGTH=100
MAX_LENGTH=1000
MIN_MEAN_Q=8

while getopts "i:o:l:L:q:" opt; do
    case "$opt" in
        i) READS_DIR="$OPTARG"  ;; # Directory with the input zipped reads in fastq format
        o) OUTPUT_DIR="$OPTARG" ;; # Directory in which the processed reads will be saved
        l) MIN_LENGTH="$OPTARG" ;; # Minimum read length required
        L) MAX_LENGTH="$OPTARG" ;; # Maximum read length required
        q) MIN_MEAN_Q="$OPTARG" ;; # Minimum mean Phre-scaled read score
        *) echo "Invalid option"; exit 1 ;;
    esac
done

mkdir -p "$OUTPUT_DIR"

for read in "$READS_DIR"/*.fastq.gz; do
    # Name of the samples (it is the same as the file name)
    base_name=$(basename "$read" .fastq.gz)
    trimmed_read="$OUTPUT_DIR/trimmed_${base_name}.fastq"
    filtered_read="$OUTPUT_DIR/filtered_${base_name}.fastq"
    
    porechop -i "$read" -o "$trimmed_read" --check_reads 1000

    filtlong -l "$MIN_LENGTH" -L "$MAX_LENGTH" -q "$MIN_MEAN_Q" | gzip > "$filtered_read"

    rm -rf "$trimmed_read"
done