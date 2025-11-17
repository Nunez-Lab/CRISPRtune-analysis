#!/usr/bin/env bash

################################################################################
## Configuration

# The name of the directory in raw_data where the sequencing files are stored
DIR="RNA_seq"

# The path to the kallisto index to use for quantification
KALLISTO_INDEX="/media/nunez/SSD-Data/genomes/Ensembl_Transcriptome/Homo_sapiens.GRCh38.cdna.all.idx"

##
##############################################################################
## Main script

# Set input/output directory variables

INPUT_DIR="raw_data/${DIR}"
OUTPUT_DIR="output/${DIR}"
####################################################################


# #Run FastQC to perform sequencing quality control checks

# mkdir -p $OUTPUT_DIR/fastqc

# fastqc \
#     -t 48 \
#     -o $OUTPUT_DIR/fastqc \
#     $INPUT_DIR/*.fastq.gz

# # Aggregate FastQC reports into a single page using MultiQC

# multiqc \
#     --filename $OUTPUT_DIR/fastqc/multiqc.html \
#     $OUTPUT_DIR/fastqc

# Trim adapters (sequences obtained from NEB #E7760 manual)

# mkdir -p $OUTPUT_DIR/trimmed
# # 
# # # Loop through the forward reads
# for forward in $INPUT_DIR/*_R1_001.fastq.gz; do
#     # Get the sample name
#     name=$(basename ${forward} _R1_001.fastq.gz)
#     echo "Running cutadapt on ${name}..."
#     # Run cutadapt (trimming software) note: core=0 means use as many as available
#     cutadapt \
#         --cores=0 \
#         -m 1 \
#         --poly-a \
#         -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
#         -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#         -o "$OUTPUT_DIR/trimmed/${name}_R1.fastq.gz" \
#         -p "$OUTPUT_DIR/trimmed/${name}_R2.fastq.gz" \
#         "$INPUT_DIR/${name}_R1_001.fastq.gz" \
#         "$INPUT_DIR/${name}_R2_001.fastq.gz"
# done


# # Run FastQC again to make sure trimming worked well

# mkdir -p $OUTPUT_DIR/fastqc-trimmed

# fastqc \
#     -t 48 \
#     -o $OUTPUT_DIR/fastqc-trimmed \
#     $OUTPUT_DIR/trimmed/*.fastq.gz

# multiqc \
#     --filename $OUTPUT_DIR/fastqc-trimmed/multiqc.html \
#     $OUTPUT_DIR/fastqc-trimmed

# Quantify with kallisto

# for forward in $OUTPUT_DIR/trimmed/*_R1.fastq.gz; do
#     name=$(basename ${forward} _R1.fastq.gz)
#     mkdir -p $OUTPUT_DIR/quant/kallisto_quant/v3_ensembl_kallisto_quant/$name
#     kallisto quant \
#         -t 4 \
#         -i $KALLISTO_INDEX \
#         -o $OUTPUT_DIR/quant/kallisto_quant/v3_ensembl_kallisto_quant/$name \
#         $OUTPUT_DIR/trimmed/${name}_R1.fastq.gz \
#         $OUTPUT_DIR/trimmed/${name}_R2.fastq.gz
# done
