library(tximportData)
library(tximport)
library(readr)
library(biomaRt)
library(rhdf5)
library(dplyr)


################################
GENE_METADATA_FILE = "/Users/p.coli/Desktop/tximport_quant/0A_KALLISTO QUANT V1-3/gene-metadata.tsv"

TPM_FILE ="/Users/p.coli/Desktop/tximport_quant/0A_KALLISTO QUANT V1-3/V3_ensembl_kallisto_quant/02_tximport_output/abundance.csv"
COUNT_FILE ="/Users/p.coli/Desktop/tximport_quant/0A_KALLISTO QUANT V1-3/V3_ensembl_kallisto_quant/02_tximport_output/counts.csv"


OUTPUT_DIR = "/Users/p.coli/Desktop/tximport_quant/0A_KALLISTO QUANT V1-3/V3_ensembl_kallisto_quant/03_filtered"

gene_metadata = read.table(
  GENE_METADATA_FILE,
  header=TRUE,
  sep="\t",
)

tpm_df = read.table(
  TPM_FILE,
  header=TRUE,
  sep=",",
)

count_df = read.table(
  COUNT_FILE,
  header=TRUE,
  sep=",",
)

###############################
biotype_count = gene_metadata %>%
  count(gene_biotype)


protein_coding = gene_metadata[
  gene_metadata$gene_biotype == "protein_coding",
  "ensembl_gene_id_version"
]
###############################

filtered_tpm = tpm_df[tpm_df[[1]] %in% protein_coding,]
filtered_count = count_df[count_df[[1]] %in% protein_coding,]
filtered_count[, -1] <- as.data.frame(lapply(filtered_count[, -1], function(x) as.integer(round(x))))


write.csv(filtered_tpm, file.path(OUTPUT_DIR, "filtered_tpm.csv"))
write.csv(filtered_count, file.path(OUTPUT_DIR, "filtered_count.csv"))

