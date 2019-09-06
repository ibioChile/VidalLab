source("https://bioconductor.org/biocLite.R")
biocLite("isomiRs", lib="/home/plantomics/pcamejo/Rlib")
library(isomiRs,lib.loc="/home/plantomics/pcamejo/Rlib")

fn_list <- c("/miraligner_fixed/1-150_S1_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/2-150_S2_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/3-150_S3_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/7-140_S4_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/8-150_S5_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/9-140_S6_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/13-140_S7_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/14-150_S8_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/15-150_S9_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/19-150_S10_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/20-150_S11_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/21-150_S12_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/25-150_S13_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/26-150_S14_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/27-140_S15_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/31-150_S16_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/32-150_S17_L000_R1_000.filtered.fixed.mirna",
             "/miraligner_fixed/33-150_S18_L000_R1_000.filtered.fixed.mirna")

de <- data.frame(row.names=c("f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13","f14","f15","f16","f17","f18"), condition = c("d5", "d5", "d5","d9","d9","d9","d13","d13","d13","d17","d17","d17","d21","d21","d21","d25","d25","d25"))
data_isomir <- IsomirDataSeqFromFiles(fn_list, coldata=de, header = TRUE)

ids <- counts(data_isomir)

write.table(ids,"miRNA_counts.tsv",quote = FALSE,sep="\t")
