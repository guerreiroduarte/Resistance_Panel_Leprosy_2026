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

rm -rf "$faidx_format"

aln_stats="$OUTPUT_DIR/Stats/alignment.txt"

rm -rf "$aln_stats"

echo -e "SAMPLE\tTOTAL\tALIGNED\tBOTH\tFORWARD\tREVERSE\tUNALIGNED" > "$aln_stats"

for read in "$READS_DIR"/filtered_*.fastq.gz; do
    file_name=$(basename "$read" .fastq.gz)
    sample_name=${file_name#filtered_}

    aligned_read="$OUTPUT_DIR"/${sample_name}_aln.sam
    amplicon_sorted="$OUTPUT_DIR"/${sample_name}.bam

    if [ ! -f "$aligned_read" ]; then

        echo "Aligning reads to amplicons"
        minimap2 -ax map-ont -u b "$target_in_fasta" "$read" > "$aligned_read"

    else
        echo "Reads were aligned"
    fi

    # Just to compute alignment statistics
    samtools view -b "$aligned_read" | samtools ampliconclip -b "$PRIMER_BED" --both-ends --strand - -o "$OUTPUT_DIR/tmp.bam" 2>&1 | \
        awk -v sample="$sample_name" 'BEGIN{OFS="\t"}
        /TOTAL READS:/     { total = $NF }
        /FORWARD CLIPPED:/ { total_f = $NF }
        /REVERSE CLIPPED:/ { total_r = $NF }
        /BOTH CLIPPED:/    { total_both = $NF }
        /EXCLUDED:/        { not_aligned = $NF }
        
        END {
            only_f = total_f - total_both
            only_r = total_r - total_both

            aligned = total - not_aligned
            
            print sample, total, aligned, total_both, only_f, only_r, not_aligned
        }' >> "$aln_stats"

    rm -rf "$OUTPUT_DIR/tmp.bam"

    if [ ! -f "$amplicon_sorted" ]; then

        echo "Filtering alignments and clipping primers"
        samtools view -b -F 4 -q 60 "$aligned_read" | \
        samtools ampliconclip -b "$PRIMER_BED" --both-ends --strand --clipped - -o - | \
        samtools sort -o "$amplicon_sorted" -

        samtools index "$amplicon_sorted"
    
    else
        echo "Primers were already clipped"
    fi

done

mkdir -p "$OUTPUT_DIR/Stats"

amplicon_stats="$OUTPUT_DIR/Stats/coverage.txt"

rm -rf "$amplicon_stats"

echo -e "METRIC\tSAMPLE\tRPOB\tGYRA\tGYRB\tFOLP1\tFOLP2\t23S_RNA_I\t23S_RNA_II" > "$amplicon_stats"

samtools ampliconstats "$PRIMER_BED" "$OUTPUT_DIR"/*.bam | \
    grep -E "^(FRPERC|FDEPTH|FVDEPTH|FREADS|FPCOV)" >> "$amplicon_stats"

REGIOES=($(awk 'OFS=":" {print $1, $2 "-" $NF}' "$REGIONS_BED"))

for bam in "$OUTPUT_DIR"/*.bam; do
    sample_name=$(basename "$bam" .bam)
    linha_metric="READSIZE\t${sample_name}"
    
    for regiao in "${REGIOES[@]}"; do
        bases_inside=$(samtools stats "$bam" "$regiao" | \
                       grep "average length:" | \
                       awk '{print $NF}')
    
        linha_metric="${linha_metric}\t${bases_inside}"
    done
    
    echo -e "$linha_metric" >> "$amplicon_stats"
done