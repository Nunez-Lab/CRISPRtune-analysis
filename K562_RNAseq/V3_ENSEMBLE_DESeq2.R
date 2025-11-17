library(DESeq2)
library(readr)

COUNT_FILE <- "/Users/p.coli/Desktop/tximport_quant/0A_KALLISTO QUANT V1-3/V3_ensembl_kallisto_quant/03_filtered/filtered_count.csv"
OUTPUT_DIR <- "/Users/p.coli/Desktop/tximport_quant/0A_KALLISTO QUANT V1-3/V3_ensembl_kallisto_quant/05_DESeq2/"

condition <- c("CD55","CD55","NT","NT")

dCas9_sample_name <- c("PCd551", "PCd552", "PCdN1", "PCdN2")
CRISPRi_sample_name <- c("PCi551", "PCi552", "PCiN1", "PCiN2")
CRISPRtune_sample_name <- c("PCM551", "PCM552", "PCMN1", "PCMN2")

dCas9_meta <- data.frame(row.names = dCas9_sample_name, condition = factor(condition))
CRISPRi_meta <- data.frame(row.names = CRISPRi_sample_name, condition = factor(condition))
CRISPRtune_meta <- data.frame(row.names = CRISPRtune_sample_name, condition = factor(condition))

counts<-read_csv(COUNT_FILE)
counts<- counts[,-1]
colnames(counts)[1] <- "GENEID"

GENEID <- counts$GENEID
count_matrix <- as.matrix(counts[,-1])

rownames(count_matrix) <- GENEID

dCas9_matrix <- count_matrix[,1:4]
CRISPRi_matrix <- count_matrix[,5:8]
CRISPRtune_matrix <- count_matrix[,9:12]

dCas9_dds = DESeqDataSetFromMatrix(
  countData = dCas9_matrix,
  colData = dCas9_meta,
  design = ~ condition
)

CRISPRi_dds = DESeqDataSetFromMatrix(
  countData = CRISPRi_matrix,
  colData = CRISPRi_meta,
  design = ~ condition
)

CRISPRtune_dds = DESeqDataSetFromMatrix(
  countData = CRISPRtune_matrix,
  colData = CRISPRtune_meta,
  design = ~ condition
)

dCas9_dds <- DESeq(dCas9_dds)
CRISPRi_dds <-DESeq(CRISPRi_dds)
CRISPRtune_dds <-DESeq(CRISPRtune_dds)


dCas9_res <-results(dCas9_dds)
CRISPRi_res <-results(CRISPRi_dds)
CRISPRtune_res <-results(CRISPRtune_dds)


write.csv(dCas9_res, file.path(OUTPUT_DIR, "dCas9_DeSEQ.csv"))
write.csv(CRISPRi_res, file.path(OUTPUT_DIR, "CRISPRi_DESeq2.csv"))
write.csv(CRISPRtune_res, file.path(OUTPUT_DIR, "CRISPRtune_DESeq2.csv"))
