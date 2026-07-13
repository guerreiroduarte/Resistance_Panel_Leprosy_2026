#!/usr/bin/env bash

MODEL_LINK="https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/r1041_e82_400bps_sup_v430_bacteria_finetuned/"
MODEL_NAME="r1041_e82_400bps_sup_v430_bacteria_finetuned"
SEED="20030407"
MAX_DEPTH="5000"

while getopts "i:o:r:t:d:s:" opt; do
    case "$opt" in
        i) ALIGNED_DIR="$OPTARG" ;; # Directory with the input zipped reads in fastq format
        o) OUTPUT_DIR="$OPTARG"  ;; # Directory in which the processed reads will be saved
        r) REF="$OPTARG"         ;; # Reference targets in fasta format
        t) TARGETS="$OPTARG"     ;; # Bed file of targets to call variants on
        d) MAX_DEPTH="$OPTARG"   ;;
        s) SEED="$OPTARG"        ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/Stats"

VAR_STATS="$OUTPUT_DIR/Stats/var_stats.txt"

rm -rf "$VAR_STATS"

echo -e "SAMPLE\tAMPLICON\tPOS\tREF\tALT\tQUAL\tDP\tDP4\tIs_Ti\tIs_Ts" > "$VAR_STATS"

for read in "$ALIGNED_DIR"/*.bam; do
    sample_name=$(basename "$read" .bam)
    var_file="$OUTPUT_DIR/${sample_name}_vars.vcf"
    
    if [ ! -f "$var_file" ]; then
        echo "Calling variants for $sample_name..."
        bcftools mpileup -Ov -B -q 15 -Q 10 -R "$TARGETS" -f "$REF" -d "$MAX_DEPTH" --seed "$SEED" "$read" | \
            bcftools call --ploidy 2 -V indels -mv -P 0.1 -o "$var_file"
    else
        echo "Variants for $sample_name were already called."
    fi

    paste \
        <(bcftools query -f "${sample_name}\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%DP\t%DP4\n" "$var_file") \
            <(bcftools stats "$var_file" | grep "^QUAL" | awk '{print $4"\t"$5}') >> "$VAR_STATS"
done