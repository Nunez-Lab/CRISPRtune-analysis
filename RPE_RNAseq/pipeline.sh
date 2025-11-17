#!/usr/bin/env bash
################################################################################
## Configuration

# The path to the kallisto index to use for quantification
KALLISTO_INDEX="/mnt/Shared1/genomes/Ensembl_Transcriptome/Homo_sapiens.GRCh38.cdna.all.idx"

##
##############################################################################
## Main script

# Set input/output directory variables

INPUT_DIR="/media/peter/Archive01/Peter_Colias/2025-08-27_CRISPRtune_RPE_RNA-seq"
OUTPUT_DIR="/media/peter/Archive01/Peter_Colias/temporary/RNA-seq"
####################################################################

#Run FastQC to perform sequencing quality control checks

mkdir -p $OUTPUT_DIR/fastqc

fastqc \
    -t 48 \
    -o $OUTPUT_DIR/fastqc \
    $INPUT_DIR/*.fastq.gz

# Aggregate FastQC reports into a single page using MultiQC

multiqc \
    --filename $OUTPUT_DIR/fastqc/multiqc.html \
    $OUTPUT_DIR/fastqc

# Trim adapters (sequences obtained from NEB #E7760 manual)

mkdir -p $OUTPUT_DIR/trimmed
# 
# # Loop through the forward reads
for forward in $INPUT_DIR/*_R1_001.fastq.gz; do
    # Get the sample name
    name=$(basename ${forward} _R1_001.fastq.gz)
    echo "Running cutadapt on ${name}..."
    # Run cutadapt (trimming software) note: core=0 means use as many as available
    cutadapt \
        --cores=0 \
        -m 1 \
        --poly-a \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o "$OUTPUT_DIR/trimmed/${name}_R1.fastq.gz" \
        -p "$OUTPUT_DIR/trimmed/${name}_R2.fastq.gz" \
        "$INPUT_DIR/${name}_R1_001.fastq.gz" \
        "$INPUT_DIR/${name}_R2_001.fastq.gz"
done

##CHECK THE QUALITY OF POST TRIMMED FILES

mkdir -p $OUTPUT_DIR/fastqc-trimmed

fastqc \
    -t 48 \
    -o $OUTPUT_DIR/fastqc-trimmed \
    $OUTPUT_DIR/trimmed/*.fastq.gz

multiqc \
    --filename $OUTPUT_DIR/fastqc-trimmed/multiqc.html \
    $OUTPUT_DIR/fastqc-trimmed


## Quantify with kallisto

mkdir -p $OUTPUT_DIR/kallisto

for forward in $OUTPUT_DIR/trimmed/*_R1.fastq.gz; do
    name=$(basename ${forward} _R1.fastq.gz)
    mkdir -p $OUTPUT_DIR/kallisto/$name
    kallisto quant \
        -t 16 \
        -i $KALLISTO_INDEX \
        -o $OUTPUT_DIR/kallisto/$name \
        $OUTPUT_DIR/trimmed/${name}_R1.fastq.gz \
        $OUTPUT_DIR/trimmed/${name}_R2.fastq.gz
done
