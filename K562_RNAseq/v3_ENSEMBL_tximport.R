library(tximportData)
library(tximport)
library(readr)
library(biomaRt)
library(rhdf5)


################################
#where kallisto output is
dir <- "/Users/p.coli/Desktop/tximport_quant/0A_KALLISTO QUANT V1-3/V3_ensembl_kallisto_quant"
OUTPUT_DIR <- "/Users/p.coli/Desktop/tximport_quant/0A_KALLISTO QUANT V1-3/V3_ensembl_kallisto_quant/02_tximport_output"
#TX2GENE_FILE<- "/Users/p.coli/Desktop/tximport_quant/tx2gene_refseq.csv"

#kallisto output files
files <- file.path(dir, "01_V3_ENSEMBLE_Kallisto_QUANT", c("PCd551_S16","PCd552_S29","PCdN1_S19","PCdN2_S32","PCi551_S14", "PCi552_S27","PCiN1_S17","PCiN2_S30","PCM551_S15","PCM552_S28","PCMN1_S18","PCMN2_S31"),
                   "abundance.h5")

name_vector <- c("PCd551","PCd552","PCdN1","PCdN2","PCi551", "PCi552","PCiN1","PCiN2","PCM551","PCM552","PCMN1","PCMN2")

names(files) <- name_vector

################################
# %% Create tx2gene using biomaRt

mart = useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  version = 114
)


tx2gene = getBM(
  attributes =
    c("ensembl_transcript_id_version",
      "ensembl_gene_id_version"),
  mart = mart
)

#tx2gene = read_csv(TX2GENE_FILE)


# %% Aggregate abundances and save count file
txi = tximport(
  files,
  type = "kallisto",
  tx2gene = tx2gene,
  #ignoreTxVersion =TRUE
)

lenscaled_txi = tximport(
  files,
  type = "kallisto",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion =TRUE
  
)

write.csv(txi$counts, file.path(OUTPUT_DIR, "counts.csv"))
write.csv(txi$abundance, file.path(OUTPUT_DIR, "abundance.csv"))
write.csv(txi$countsFromAbundance, file.path(OUTPUT_DIR, "countsFromAbundance.csv"))

write.csv(lenscaled_txi$counts, file.path(OUTPUT_DIR, "lenscaled_counts.csv"))
write.csv(lenscaled_txi$abundance, file.path(OUTPUT_DIR, "lenscaled_abundance.csv"))
write.csv(lenscaled_txi$countsFromAbundance, file.path(OUTPUT_DIR, "lenscaled_countsFromAbundance.csv"))

#tx2gene <- read_csv("/media/nunez/HDD-Data/Da/CRISPRtune/tx2gene.csv")
#txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "scaledTPM", ignoreTxVersion = FALSE)

#setwd('/Users/p.coli/Desktop/tximport_quant/remake plots')

#for (name in names(txi)) {write.table(
#   txi[[name]], 
#  file = paste0(name, ".tsv"), 
# sep = "\t", 
#row.names = TRUE, 
#col.names = NA,
#quote = FALSE
#)
#}

