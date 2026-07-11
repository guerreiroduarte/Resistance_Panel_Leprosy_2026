#!/usr/bin/env bash

WGET_LINK="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/855/GCF_000195855.1_ASM19585v1/GCF_000195855.1_ASM19585v1_genomic.fna.gz"

while getopts "i:o:w:a:r:p:" opt; do
    case "$opt" in
        i) READS_DIR="$OPTARG"    ;; # Directory with the input zipped reads in fastq format
        o) OUTPUT_DIR="$OPTARG"   ;; # Directory in which the processed reads will be saved
        w) WGET_LINK="$OPTARG"    ;; # Link to download reference genome fasta file
        a) AMPLICON_BED="$OPTARG" ;; # Amplicon bed used to retain sequences from reference genome
        r) REGIONS_BED="$OPTARG"  ;; # Amplicon bed to compute respective coverage
        p) PRIMER_BED="$OPTARG"   ;; # Primers to clip them from aligned reads
        *) echo "Invalid option"; exit 1 ;;
    esac
done

ALIGNMENT_STATS="$OUTPUT_DIR/Stats"

mkdir -p "$OUTPUT_DIR" "$ALIGNMENT_STATS"

ref_genome="$OUTPUT_DIR/Mleprae_TN.fasta"
target_in_fasta="$OUTPUT_DIR/Targets.fasta"
faidx_format="$OUTPUT_DIR/tmp_regions.txt"

awk '{print $1":"$2"-"$3}' "$AMPLICON_BED" > "$faidx_format"

if [ ! -f "$ref_genome" ]; then
    echo "Downloading M. leprae TN reference genome from NCBI..."
    wget -qO- "$WGET_LINK" | gunzip > "$ref_genome"
else
    echo "M. leprae TN genome already existis in $OUTPUT_DIR"
fi

samtools faidx "$ref_genome" -r "$faidx_format" -o "$target_in_fasta"
samtools faidx "$target_in_fasta"

for read in "$READS_DIR"/filtered_*.fastq.gz; do
    file_name=$(basename "$read" .fastq.gz)
    sample_name=${file_name#filtered_}

    aligned_read="$OUTPUT_DIR"/${sample_name}_aln.sam
    amplicon_sorted="$OUTPUT_DIR"/${sample_name}.bam

    if [ ! -f "$aligned_read" ]; then

        echo "Aligninig reads to amplicons"
        minimap2 -ax map-ont -u b "$target_in_fasta" "$read" > "$aligned_read"

        echo "Filtering alignments and clipping primers"
        samtools view -b -F 4 -q 60 "$aligned_read" | \
        samtools ampliconclip -b "$PRIMER_BED" --hard-clip --both-ends --strand --clipped - -o - | \
        samtools sort -o "$amplicon_sorted" -

        samtools index "$amplicon_sorted"
    
    else
        echo "Reads were already aligned to reference targets..."
    fi

done

mkdir -p "$OUTPUT_DIR/Stats"

amplicon_stats="$OUTPUT_DIR/Stats/coverage.txt"

rm -rf "$amplicon_stats"

echo -e "METRIC\tSAMPLE\tRPOB\tGYRA\tGYRB\tFOLP1\tFOLP2\t23S_RNA_I\t23S_RNA_II" > "$amplicon_stats"

samtools ampliconstats "$PRIMER_BED" "$OUTPUT_DIR"/*.bam | \
    grep -E "^(FRPERC|FDEPTH|FVDEPTH|FREADS|FPCOV)" >> "$amplicon_stats"