library(dplyr)
library(tidyverse)
library(DiffBind)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
SAMPLESHEET_FILE = args[1]
COMPARISONS= args[2]
OUTPUT_PATH = args[3]
  
samples=read.csv(SAMPLESHEET_FILE, header=TRUE,
)

comparisons=read.csv(COMPARISONS, header=TRUE,)
                 

# antibodies = unique(samples$Factor) #for downstream iteration

H3K27ac_samples = subset(samples, Factor == "H3K27ac")
H3K4me3_samples = subset(samples, Factor == "H3K4me3")
H3K9me3_samples = subset(samples, Factor == "H3K9me3")
# 
# analysis_list = list()
# report_list= list()


for (i in 1:nrow(comparisons)) {
  row = unname(unlist(comparisons[i,]))
  print(row[2])
  subset_df = subset(H3K27ac_samples, Condition == row[2] |Condition == row[3])
  CUTRUN.dba = dba(sampleSheet = subset_df)
  CUTRUN.counts = dba.count(CUTRUN.dba, minOverlap=1, summits=250) #prevents dropping windows that only appear in one sample
  CUTRUN.norm = dba.normalize(CUTRUN.counts)
  CUTRUN.contrast = dba.contrast(CUTRUN.norm, categories=DBA_CONDITION, minMembers=2 )
  CUTRUN.analyze = dba.analyze(CUTRUN.contrast, method=DBA_ALL_METHODS)
  # dba.show(CUTRUN.analyze, bContrasts=TRUE)
  CUTRUN.DB = dba.report(CUTRUN.analyze)
  CUTRUN.report.edgeR = dba.report(CUTRUN.analyze,
                             method = DBA_EDGER,
                             contrast = 1,
                             th = 1,
                             bUsePval = FALSE
  )
  CUTRUN.report.deseq2 = dba.report(CUTRUN.analyze,
                             method = DBA_DESEQ2, 
                             contrast = 1,
                             th = 1,
                             bUsePval = FALSE
  )
  CUTRUN.report.edgeR.df = as.data.frame(CUTRUN.report.edgeR)
  CUTRUN.report.deseq2.df = as.data.frame(CUTRUN.report.deseq2)
  
  png(file.path(OUTPUT_PATH, paste0("H3K27ac_",row[4], ".png")), width=2000, height=2000, res=300)
  dba.plotVolcano(CUTRUN.analyze,method = DBA_EDGER)
  dev.off()
  
  write.csv(CUTRUN.report.edgeR.df, file.path(OUTPUT_PATH, paste0("edgeR_H3K27ac_",row[4], ".csv")),)
  write.csv(CUTRUN.report.deseq2.df, file.path(OUTPUT_PATH, paste0("DESeq2_H3K27ac_",row[4], ".csv")),)
  
}



for (i in 1:nrow(comparisons)) {
  row = unname(unlist(comparisons[i,]))
  print(row[2])
  subset_df = subset(H3K4me3_samples, Condition == row[2] |Condition == row[3])
  CUTRUN.dba = dba(sampleSheet = subset_df)
  CUTRUN.counts = dba.count(CUTRUN.dba, minOverlap=1, summits=250) #prevents dropping windows that only appear in one sample
  CUTRUN.norm = dba.normalize(CUTRUN.counts)
  CUTRUN.contrast = dba.contrast(CUTRUN.norm, categories=DBA_CONDITION, minMembers=2 )
  CUTRUN.analyze = dba.analyze(CUTRUN.contrast, method=DBA_ALL_METHODS)
  # dba.show(CUTRUN.analyze, bContrasts=TRUE)
  CUTRUN.DB = dba.report(CUTRUN.analyze)
  CUTRUN.report.edgeR = dba.report(CUTRUN.analyze,
                                   method = DBA_EDGER,
                                   contrast = 1,
                                   th = 1,
                                   bUsePval = FALSE
  )
  CUTRUN.report.deseq2 = dba.report(CUTRUN.analyze,
                                    method = DBA_DESEQ2, 
                                    contrast = 1,
                                    th = 1,
                                    bUsePval = FALSE
  )
  CUTRUN.report.edgeR.df = as.data.frame(CUTRUN.report.edgeR)
  CUTRUN.report.deseq2.df = as.data.frame(CUTRUN.report.deseq2)
  png(file.path(OUTPUT_PATH, paste0("H3K4me3_",row[4], ".png")), width=2000, height=2000, res=300)
  dba.plotVolcano(CUTRUN.analyze,method = DBA_EDGER)
  dev.off()
  
  write.csv(CUTRUN.report.edgeR.df, file.path(OUTPUT_PATH, paste0("edgeR_H3K4me3_",row[4], ".csv")),)
  write.csv(CUTRUN.report.deseq2.df, file.path(OUTPUT_PATH, paste0("DESeq2_H3K4me3_",row[4], ".csv")),)
  
}



for (i in 1:nrow(comparisons)) {
  row = unname(unlist(comparisons[i,]))
  print(row[2])
  subset_df = subset(H3K9me3_samples, Condition == row[2] |Condition == row[3])
  CUTRUN.dba = dba(sampleSheet = subset_df)
  CUTRUN.counts = dba.count(CUTRUN.dba, minOverlap=1, summits=250) #prevents dropping windows that only appear in one sample
  CUTRUN.norm = dba.normalize(CUTRUN.counts)
  CUTRUN.contrast = dba.contrast(CUTRUN.norm, categories=DBA_CONDITION, minMembers=2 )
  CUTRUN.analyze = dba.analyze(CUTRUN.contrast, method=DBA_ALL_METHODS)
  # dba.show(CUTRUN.analyze, bContrasts=TRUE)
  CUTRUN.DB = dba.report(CUTRUN.analyze)
  CUTRUN.report.edgeR = dba.report(CUTRUN.analyze,
                                   method = DBA_EDGER,
                                   contrast = 1,
                                   th = 1,
                                   bUsePval = FALSE
  )
  CUTRUN.report.deseq2 = dba.report(CUTRUN.analyze,
                                    method = DBA_DESEQ2, 
                                    contrast = 1,
                                    th = 1,
                                    bUsePval = FALSE
  )
  CUTRUN.report.edgeR.df = as.data.frame(CUTRUN.report.edgeR)
  CUTRUN.report.deseq2.df = as.data.frame(CUTRUN.report.deseq2)
  png(file.path(OUTPUT_PATH, paste0("H3K9me3_", row[4], ".png")), width=2000, height=2000, res=300)
  dba.plotVolcano(CUTRUN.analyze,method = DBA_EDGER)
  dev.off()
  
  write.csv(CUTRUN.report.edgeR, file.path(OUTPUT_PATH, paste0("edgeR_H3K9me3_",row[4], ".csv")),)
  write.csv(CUTRUN.report.deseq2, file.path(OUTPUT_PATH, paste0("DESeq2_H3K9me3_",row[4], ".csv")),)
  
}





