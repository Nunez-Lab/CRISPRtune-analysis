# !/usr/bin/env bash
#adapted from Sinha et al., Molecular Cell 2024

# uv init

BOWTIE_INDEX="/mnt/Shared1/genomes/Bowtie_GRCh38_Index/GRCh38_noalt_as/GRCh38_noalt_as"
ENCODE_BLACKLIST="/mnt/Shared1/genomes/ENCODE_Exclusion_List_Regions/ENCFF356LFX.bed"
E_COLI_INDEX="/mnt/Shared1/genomes/Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Sequence/Bowtie2Index/genome"

CUT_INPUT_DIR="/mnt/Archive01/Peter_Colias/2025-08-27_CRISPRtune_RPE_ATAC_CUTRUN/Raw/CUTRUN"
ATAC_INPUT_DIR="/mnt/Archive01/Peter_Colias/2025-08-27_CRISPRtune_RPE_ATAC_CUTRUN/Raw/ATAC"

CUT_OUTPUT_DIR="/mnt/Archive01/Peter_Colias/2025-08-27_CRISPRtune_RPE_ATAC_CUTRUN/Processed/CUTRUN"
ATAC_OUTPUT_DIR="/mnt/Archive01/Peter_Colias/2025-08-27_CRISPRtune_RPE_ATAC_CUTRUN/Processed/ATAC"

CUT_OUTPUT_DIR_E_COLI="${CUT_OUTPUT_DIR}/E_COLI_NORM"

########################################
##FASTQC and Trimming of CUTRUN

# mkdir -p $CUT_OUTPUT_DIR/fastqc

# fastqc \
#     -t 48 \
#     -o $CUT_OUTPUT_DIR/fastqc \
#     $CUT_INPUT_DIR/*.fastq.gz

# Aggregate FastQC reports into a single page using MultiQC

# multiqc \
#     --filename $CUT_OUTPUT_DIR/fastqc/multiqc.html \
#     $CUT_OUTPUT_DIR/fastqc

# trim adaptors CUT&RUN

# mkdir -p $CUT_OUTPUT_DIR/trimmed

# # # Loop through the forward reads
# for forward in $CUT_INPUT_DIR/*_R1_001.fastq.gz; do
#     # Get the sample name
#     name=$(basename ${forward} _R1_001.fastq.gz)
#     echo "Running cutadapt on ${name}..."
#     # Run cutadapt (trimming software)
#     cutadapt \
#         --cores=32 \
#         -m 1 \
#         --poly-a \
#         -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
#         -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#         -o "$CUT_OUTPUT_DIR/trimmed/${name}_R1.fastq.gz" \
#         -p "$CUT_OUTPUT_DIR/trimmed/${name}_R2.fastq.gz" \
#         "$CUT_INPUT_DIR/${name}_R1_001.fastq.gz" \
#         "$CUT_INPUT_DIR/${name}_R2_001.fastq.gz"
# done

# # Run FastQC again on CUT&RUN to make sure trimming worked well

# mkdir -p $CUT_OUTPUT_DIR/fastqc-trimmed

# fastqc \
#     -t 48 \
#     -o $CUT_OUTPUT_DIR/fastqc-trimmed \
#     $CUT_OUTPUT_DIR/trimmed/*.fastq.gz

# multiqc \
#     --filename $CUT_OUTPUT_DIR/fastqc-trimmed/multiqc.html \
#     $CUT_OUTPUT_DIR/fastqc-trimmed

# ########################################
# ##Alignment
# mkdir -p $CUT_OUTPUT_DIR/aligned
# mkdir -p $CUT_OUTPUT_DIR/aligned_log

# for forward in $CUT_OUTPUT_DIR/trimmed/*_R1.fastq.gz; do
#     # Get the sample name
#     name=$(basename ${forward} _R1.fastq.gz)
#     echo "Running bowtie2 on ${name}..."
#     R1=$CUT_OUTPUT_DIR/trimmed/${name}_R1.fastq.gz
#     R2=$CUT_OUTPUT_DIR/trimmed/${name}_R2.fastq.gz
#     bowtie2 \
#         --threads 30 \
#         --end-to-end --very-sensitive \
#         --no-discordant \
#         --phred33 \
#         -I 10 \
#         -X 700 \
#         --no-mixed \
#         -x $BOWTIE_INDEX \
#         -1 "$R1" \
#         -2 "$R2" \
#         -S $CUT_OUTPUT_DIR/aligned/${name}_bowtie2.sam &> $CUT_OUTPUT_DIR/aligned_log/${name}_bowtie2.txt   
# done


# echo "Running quickcheck on alignment"
# err_count=0
# for forward in $CUT_OUTPUT_DIR/aligned/*_bowtie2.sam; do
#     name=$(basename ${forward} _bowtie2.sam)
#     if samtools \
#         quickcheck \
#         -v \
#         "$CUT_OUTPUT_DIR/aligned/${name}_bowtie2.sam"; then
#             :
#     else
#         ((err_count++))
#     fi
# done

# ########################################
# ##filtering and normalization

# # %% multimapping filter

# mkdir -p $CUT_OUTPUT_DIR_E_COLI/MAPQ_filter

# echo "MAPQ filtering started"

# for file in $CUT_OUTPUT_DIR/aligned/*_bowtie2.sam; do
#     name=$(basename ${file} _bowtie2.sam)
#     samtools \
#         view \
#         -h \
#         -q 30 \
#         -b \
#         -@ 32 \
#         -o $CUT_OUTPUT_DIR_E_COLI/MAPQ_filter/${name}_MAPQ_filter.bam \
#         $CUT_OUTPUT_DIR/aligned/${name}_bowtie2.sam
# done


# ## blacklist site removal
# mkdir -p $CUT_OUTPUT_DIR_E_COLI/DAC_exclude

# echo "ENCODE Blacklist Removal ..."
# for file in $CUT_OUTPUT_DIR_E_COLI/MAPQ_filter/*_MAPQ_filter.bam; do
#     name=$(basename ${file} _MAPQ_filter.bam)
#     echo "filtering ${name}..."
#     bedtools \
#         intersect \
#         -abam $file \
#         -b ${ENCODE_BLACKLIST} \
#         -v \
#         > "$CUT_OUTPUT_DIR_E_COLI/DAC_exclude/${name}_DAC_exclude.bam"
# done

# ##%% Ecoli norm deduplication

# mkdir -p $CUT_OUTPUT_DIR_E_COLI/SAM_name_sort

# for file in $CUT_OUTPUT_DIR_E_COLI/DAC_exclude/*_DAC_exclude.bam; do
#     name=$(basename ${file} _DAC_exclude.bam)
#     echo "Running samtools sort by name on ${name}_DAC_exclude.bam"
#     samtools \
#         sort \
#         -n \
#         -@ 64 \
#         -o "$CUT_OUTPUT_DIR_E_COLI/SAM_name_sort/${name}_sorted_names.bam" \
#         -O BAM \
#         "$CUT_OUTPUT_DIR_E_COLI/DAC_exclude/${name}_DAC_exclude.bam"
# done

# mkdir -p $CUT_OUTPUT_DIR_E_COLI/fixMate

# for file in $CUT_OUTPUT_DIR_E_COLI/SAM_name_sort/*_sorted_names.bam; do
#     name=$(basename ${file} _sorted_names.bam)
#     echo "Running samtools fixmate on ${name}_sorted_names.bam"
#     samtools \
#         fixmate \
#         -m \
#         -@ 64 \
#         "$CUT_OUTPUT_DIR_E_COLI/SAM_name_sort/${name}_sorted_names.bam" \
#         "$CUT_OUTPUT_DIR_E_COLI/fixMate/${name}_fixmate.bam"
# done

# mkdir -p $CUT_OUTPUT_DIR_E_COLI/coordinate_sort

# for file in $CUT_OUTPUT_DIR_E_COLI/fixMate/*_fixmate.bam; do
#     name=$(basename ${file} _fixmate.bam)
#     echo "Running samtools sort on ${name}_fixmate.bam"
#     samtools \
#         sort \
#         -@ 64 \
#         -o "$CUT_OUTPUT_DIR_E_COLI/coordinate_sort/${name}_coordinate_sort.bam" \
#         "$CUT_OUTPUT_DIR_E_COLI/fixMate/${name}_fixmate.bam"
# done

# mkdir -p $CUT_OUTPUT_DIR_E_COLI/deduplicated

# for file in $CUT_OUTPUT_DIR_E_COLI/coordinate_sort/*_coordinate_sort.bam; do
#     name=$(basename ${file} _coordinate_sort.bam)
#     echo "Running samtools markdup on ${name}_coordinate_sort.bam"
#     samtools \
#         markdup \
#         -@ 24 \
#         -r \
#         -S \
#         "$CUT_OUTPUT_DIR_E_COLI/coordinate_sort/${name}_coordinate_sort.bam" \
#         "$CUT_OUTPUT_DIR_E_COLI/deduplicated/${name}_deduplicated.bam"
# done

# for file in $CUT_OUTPUT_DIR_E_COLI/deduplicated/*_deduplicated.bam; do
#     name=$(basename ${file} _deduplicated.bam)
#     echo "running samtools index on ${name}_deduplicated.bam ..."
#     samtools \
#         index \
#         -b \
#         -@ 24 \
#         "$CUT_OUTPUT_DIR_E_COLI/deduplicated/${name}_deduplicated.bam"
# done

# mkdir -p $CUT_OUTPUT_DIR_E_COLI/post_alignmnet_QC

# for file in $CUT_OUTPUT_DIR_E_COLI/deduplicated/*_deduplicated.bam; do
#     name=$(basename ${file} _deduplicated.bam)
#     echo "Running QC on ${name}_deduplicated.bam..."
#     samtools \
#         flagstat \
#         -@ 24 \
#         $CUT_OUTPUT_DIR_E_COLI/deduplicated/${name}_deduplicated.bam > $CUT_OUTPUT_DIR_E_COLI/post_alignmnet_QC/${name}_post_alignment_QC.txt
# done

# # align to E. Coli genome and count reads
# mkdir -p $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/aligned
# mkdir -p $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/aligned_log

# for forward in $CUT_OUTPUT_DIR/trimmed/*_R1.fastq.gz; do
#     # Get the sample name
#     name=$(basename ${forward} _R1.fastq.gz)
#     echo "Running bowtie2 on ${name}..."
#     R1=$CUT_OUTPUT_DIR/trimmed/${name}_R1.fastq.gz
#     R2=$CUT_OUTPUT_DIR/trimmed/${name}_R2.fastq.gz
#     bowtie2 \
#         --threads 30 \
#         --end-to-end --very-sensitive \
#         --no-discordant \
#         --phred33 \
#         -I 10 \
#         -X 700 \
#         --no-mixed \
#         -x $E_COLI_INDEX \
#         -1 "$R1" \
#         -2 "$R2" \
#         -S $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/aligned/${name}_E_Coli_bowtie2.sam &> $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/aligned_log/${name}_E_Coli_bowtie2.txt   
# done

# # %%Filter multimapping
# mkdir -p $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/MAPQ_filter

# for file in $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/aligned/*_E_Coli_bowtie2.sam; do
#     name=$(basename ${file} _E_Coli_bowtie2.sam)
#     samtools \
#         view \
#         -h \
#         -q 30 \
#         -b \
#         -@ 32 \
#         -o $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/MAPQ_filter/${name}_E_Coli_MAPQ_filter.bam \
#         $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/aligned/${name}_E_Coli_bowtie2.sam
# done

# # %%deduplicate
# mkdir -p $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/SAM_name_sort

# for file in $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/MAPQ_filter/*_E_Coli_MAPQ_filter.bam; do
#     name=$(basename ${file} _E_Coli_MAPQ_filter.bam)
#     echo "Running samtools sort by name on ${name}_E_Coli_MAPQ_filter.bam"
#     samtools \
#         sort \
#         -n \
#         -@ 64 \
#         -o "$CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/SAM_name_sort/${name}_E_Coli_sorted_names.bam" \
#         -O BAM \
#         $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/MAPQ_filter/${name}_E_Coli_MAPQ_filter.bam
# done

# mkdir -p $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/fixMate

# for file in $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/SAM_name_sort/*_E_Coli_sorted_names.bam; do
#     name=$(basename ${file} _E_Coli_sorted_names.bam)
#     echo "Running samtools fixmate on ${name}_E_Coli_sorted_names.bam"
#     samtools \
#         fixmate \
#         -m \
#         -@ 64 \
#         "$CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/SAM_name_sort/${name}_E_Coli_sorted_names.bam" \
#         "$CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/fixMate/${name}_E_Coli_fixmate.bam"
# done

# mkdir -p $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/coordinate_sort

# for file in $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/fixMate/*_E_Coli_fixmate.bam; do
#     name=$(basename ${file} _E_Coli_fixmate.bam)
#     echo "Running samtools sort on ${name}_E_Coli_fixmate.bam"
#     samtools \
#         sort \
#         -@ 64 \
#         -o "$CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/coordinate_sort/${name}_E_Coli_coordinate_sort.bam" \
#         "$CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/fixMate/${name}_E_Coli_fixmate.bam"
# done

# mkdir -p $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/deduplicated

# for file in $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/coordinate_sort/*_E_Coli_coordinate_sort.bam; do
#     name=$(basename ${file} _E_Coli_coordinate_sort.bam)
#     echo "Running samtools markdup on ${name}_E_Coli_coordinate_sort.bam"
#     samtools \
#         markdup \
#         -@ 24 \
#         -r \
#         -S \
#         "$CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/coordinate_sort/${name}_E_Coli_coordinate_sort.bam" \
#         "$CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/deduplicated/${name}_E_Coli_deduplicated.bam"
# done

# # %%index
# for file in $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/deduplicated/*_E_Coli_deduplicated.bam; do
#     name=$(basename ${file} _E_Coli_deduplicated.bam)
#     echo "running samtools index on ${name}_E_Coli_deduplicated.bam ..."
#     samtools \
#         index \
#         -b \
#         -@ 24 \
#         "$CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/deduplicated/${name}_E_Coli_deduplicated.bam"
# done

# mkdir -p $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/post_alignmnet_QC

# for file in $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/deduplicated/*_E_Coli_deduplicated.bam; do
#     name=$(basename ${file} _E_Coli_deduplicated.bam)
#     echo "Running QC on ${name}_E_Coli_deduplicated.bam..."
#     samtools \
#         flagstat \
#         -@ 24 \
#         $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/deduplicated/${name}_E_Coli_deduplicated.bam > $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/post_alignmnet_QC/${name}_E_Coli_post_alignment_QC.txt
# done

# ########################################
# ##Generate bigwigs
# mkdir -p $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/dedup_bigwig_RPGC

# for file in $CUT_OUTPUT_DIR_E_COLI/deduplicated/*_deduplicated.bam; do
#     name=$(basename ${file} _deduplicated.bam)
#     #echo "$name"
#     total_reads=$CUT_OUTPUT_DIR_E_COLI/deduplicated/${name}_deduplicated.bam
#     ecoli_reads=$CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/deduplicated/${name}_E_Coli_deduplicated.bam
#     total_count=$(samtools \
#         view \
#         -c \
#         -F 4 \
#         -@ 48 \
#         $total_reads)
#     ecoli_count=$(samtools \
#         view \
#         -c \
#         -F 4 \
#         -@ 48 \
#         $ecoli_reads)
#     fraction_ecoli=$(awk "BEGIN {print $ecoli_count/$total_count}")
#     percent_ecoli=$(awk "BEGIN {print $fraction_ecoli*100}")
#     norm_factor=$(awk "BEGIN {print 1/$percent_ecoli}")
#     echo "${name}, ${norm_factor}"

#     bamCoverage \
#         --binSize 10 \
#         --extendReads \
#         --centerReads \
#         -p max \
#         --normalizeUsing RPGC \
#         --effectiveGenomeSize 2913022398 \
#         --scaleFactor $norm_factor \
#         -b $CUT_OUTPUT_DIR_E_COLI/deduplicated/${name}_deduplicated.bam \
#         -o $CUT_OUTPUT_DIR_E_COLI/E_COLI_Aligned/dedup_bigwig_RPGC/${name}_RPGC_deduplicated.bw
# done

# %%peak calling

# mkdir -p $CUT_OUTPUT_DIR_E_COLI/MACS3_broadpeaks
# mkdir -p $CUT_OUTPUT_DIR_E_COLI/MACS3_narrowpeaks

# uv venv
# source .venv/bin/activate
# uv add macs3

# igg1="$CUT_OUTPUT_DIR_E_COLI/deduplicated/PC_CR_IgG_WT_R1_S379_L006_deduplicated.bam"
# igg2="$CUT_OUTPUT_DIR_E_COLI/deduplicated/PC_CR_IgG_WT_R2_S380_L006_deduplicated.bam"

# subsample_list=("PC_CR_4" "PC_CR_9" "PC_CR_27")

# for editor in "${subsample_list[@]}"; do
#     echo "Running ${editor} samples ..."

#     for file in $CUT_OUTPUT_DIR_E_COLI/deduplicated/${editor}*_deduplicated.bam; do
#         name=$(basename ${file} _deduplicated.bam)
#         echo "running MACS3 broad on $name..."
#         macs3 \
#             callpeak \
#             -t "$CUT_OUTPUT_DIR_E_COLI/deduplicated/${name}_deduplicated.bam" \
#             -c  "${igg1}" "${igg2}" \
#             -f BAMPE \
#             -g hs \
#             -n "${name}" \
#             --keep-dup all \
#             --broad \
#             --broad-cutoff 0.1 \
#             --outdir $CUT_OUTPUT_DIR_E_COLI/MACS3_broadpeaks &
#     done
    
#     wait
#     echo "Finished ${editor} samples ..."

# done


# echo "started narrowpeak analysis"

# for editor in "${subsample_list[@]}"; do
#     echo "Running ${editor} samples ..."

#     for file in $CUT_OUTPUT_DIR_E_COLI/deduplicated/${editor}*_deduplicated.bam; do
#         name=$(basename ${file} _deduplicated.bam)
#         echo "running MACS3 narrow on ${name}..."
#         macs3 \
#             callpeak \
#             -t "$CUT_OUTPUT_DIR_E_COLI/deduplicated/${name}_deduplicated.bam" \
#             -c  "${igg1}" "${igg2}" \
#             -f BAMPE \
#             -g hs \
#             -n "${name}" \
#             --keep-dup all \
#             --qvalue 0.01 \
#             --call-summits \
#             --outdir $CUT_OUTPUT_DIR_E_COLI/MACS3_narrowpeaks &
#         done
    
#     wait
#     echo "Finished ${editor} samples ..."

# done


# deactivate
# echo "deactivated"

# mkdir -p $CUT_OUTPUT_DIR_E_COLI/DiffBind_broad
# Rscript "/mnt/Archive01/Peter_Colias/2025-08-27_CRISPRtune_RPE_ATAC_CUTRUN/scripts/DiffBind.R"
