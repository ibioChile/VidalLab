# microRNA Differential Expression


This pipeline explains how to find differentially expressed microRNAs among different conditions, including the processing of small RNA libraries. In this case, we will analyze an experiment where samples from *Arabidopsis thaliana* were collected at 6 different days (3 replicates each day). 

1. Remove adapters from microRNA libraries using cutadapt. Reads with length <18 bp and >28 bp are discarded. 

```for file in fastq/*.fastq; do base=${file##*/}; cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 18 -M 28 --discard-untrimmed -o adapter/${base%.*}.filtered.fastq $file; done```

2. Dowload miRBase [*hairpin.fa*](http://www.mirbase.org/ftp.shtml) and [*miRNA.str*](http://www.mirbase.org/ftp.shtml) files. The last one has to be downloaded from the FTP site. Create a folder ```DB/``` and store files there.

3. We use the [miraligner](https://code.google.com/p/seqbuster/wiki/miraligner) tool to map reads to microRNAs from Arabidopsis. The following script allows up to one mismatch during mapping.

```for file in adapter/*.fastq; do base=${file##*/}; java -jar miraligner.jar -sub 1 -trim 3 -add 3 -s ath -i $file -db DB/ -o miraligner_out/${base%.*}; done```

4. Format *.mirna* files to use file as input for next step. This script fixes 'Freq' column, changing 0 to 1 values.

```for file in /miraligner_out/*.mirna; do base=${file##*/}; awk -F$'\t' 'BEGIN {OFS = FS} {$3="1";print}' $file > /miraligner_fixed/${base%.*}.fixed.mirna; sed -i '' -e 's/name\t1/name\tfreq/g' /miraligner_fixed/${base%.*}.fixed.mirna; done```

5. A counts table can be generated using the [*isoCounts.R*](https://github.com/ibioChile/VidalLab/blob/master/Scripts/isoCounts.R) script, which takes the formatted *.mirna* files as input. This script uses the [isoCounts](http://lpantano.github.io/isomiRs/reference/isoCounts.html) function from [isomiRs](https://bioconductor.org/packages/release/bioc/html/isomiRs.html) package.

```Rscript isoCounts.R```

I recommend to run this script on a server, since it requires a long time of processing. The output of this script can be used for any DE analysis.

**We will explain 2 different methods for temporal DE analysis. The next steps will be run on RStudio** 

6. Import [counts](https://github.com/ibioChile/VidalLab/blob/master/Data/Arabidopsis_microRNA/miRNA_counts.tsv) file to R. 

```counts_isomirs <- read.table("miRNA_counts.tsv", header=TRUE, sep="\t")```

7. Filter counts: Only miRNA with a minimum of 5 counts in at least 25% of samples (5 in this case) are considered for DE analysis.

```
x <- counts_isomirs >= 5
counts_isomirs_add <- counts_isomirs[rowSums(x == TRUE) >= 5,]
```

## Approach 1: One-way analysis of variance

- Raw read counts were first median-normalized to adjust for the effect of library sizes and read count distributions *(Anders & Huber 2010)*. Normalized counts are converted to the log2 scale, using log2 (x + 1) for conversion.

```
library(DESeq2)
de <- data.frame(row.names=c("f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13","f14","f15","f16","f17","f18"), condition = c("d5", "d5", "d5","d9","d9","d9","d13","d13","d13","d17","d17","d17","d21","d21","d21","d25","d25","d25"))
dds = DESeqDataSetFromMatrix(counts_isomirs_add , de , design = ~ condition)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_log <- log2(normalized_counts + 1)
```


- To check that this type of normalization actually changed the variance of your data, we can plot the relative log distribution (rle), before and after normalization.

```
counts_isomirs_add_log <- log2(counts_isomirs_add  + 1)

mn_bef  <- apply(counts_isomirs_add_log, 1, median)
rle <- data.frame(sweep(counts_isomirs_add_log, MARGIN=1, STATS=mn_bef , FUN='-'))
boxplot(rle, xlab="Samples",ylab="RLE Before Normalization")

mn_aft  <- apply(normalized_counts_log, 1, median)
rle <- data.frame(sweep(normalized_counts_log, MARGIN=1, STATS=mn_aft, FUN='-'))
boxplot(rle, xlab="Samples",ylab="RLE After Normalization")
```

This is an example of these plots.

![Deseq2_normalization](https://user-images.githubusercontent.com/53570955/64458226-c2af1e00-d0c2-11e9-9a2f-bf194f16132e.jpg)

- Use the function *anova* to carry out an analysis of variance for each miRNA data.

```
data_in <- t(normalized_counts_log)
data_in <- cbind(time = de,data_in)

p_val_list <- character(0)
for (gene in rownames(normalized_counts_log)){
  data_to_aov <- data_in[,c("condition",gene)]
  colnames(data_to_aov) <- c("condition","miRNA")
  fmla <- as.formula("miRNA ~ condition")
  aov_fmla <- aov(fmla, data = data_to_aov)
  p_val  <- summary(aov_fmla)[[1]][["Pr(>F)"]][1]
  p_val_list <- as.numeric(append(p_val_list,p_val))
}
```

- The Benjamini-Hochberg procedure (Benjamini & Hochberg, 1995) is used to control the false discovery rate (FDR) based on the p-values obtained from the ANOVA analysis. 

```
pval_adjust <- p.adjust(p_val_list, method = "BH")
data_out_aov <- data.frame(cbind(normalized_counts_log,'p-val'=p_val_list, 'adjusted p-val' = pval_adjust))
```

- microRNAs having p-values with an FDR threshold < 0.05 are designated as differentially expressed.

```
DE_miRNA <- data_out_aov[data_out_aov$adjusted.p.val < 0.05,]
```

## Approach 2: RNentropy

RNentropy (*Zambelli et al. 2018*): identification of genes showing a significant variation of expression across all conditions studied. The samples, corresponding to different conditions sequenced in any number of replicates, are compared by taking into account the global gene expression level and at the same time the impact of biological variation across replicates. 

- Normalized counts are analyzed using the RN_calc tool of the RNentropy R package. Selection of DE miRNA is done with RN_select using the  Benjamini-Hochberg procedure (Benjamini & Hochberg, 1995)  to control the false discovery rate (FDR). Thresholds for global and local p-values are set to 0.01 (default). 

```
cond <- rbind(matrix(data=c(1,0,0,0,0,0),nrow=3,ncol=6,byrow = TRUE),
              matrix(data=c(0,1,0,0,0,0),nrow=3,ncol=6,byrow = TRUE),
              matrix(data=c(0,0,1,0,0,0),nrow=3,ncol=6,byrow = TRUE),
              matrix(data=c(0,0,0,1,0,0),nrow=3,ncol=6,byrow = TRUE),
              matrix(data=c(0,0,0,0,1,0),nrow=3,ncol=6,byrow = TRUE),
              matrix(data=c(0,0,0,0,0,1),nrow=3,ncol=6,byrow = TRUE))
cond <- data.frame(cond, row.names=c("f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13","f14","f15","f16","f17","f18"))
names(cond) <- c("d5","d9","d13","d17","d21","d25")

RNent_results <- RN_calc(as.data.frame(normalized_counts), as.matrix(cond))
RN_pmi(RNent_results)
D3_results_De <- RN_select(RNent_results, method = "BH")
DE_results <- D3_results_De$selected
```

- A post-processing filter was applied to determine DE miRNA. Only samples with an estimated expression value (non NA) in at least 5 of the 6 days analyzed were considered as DE over the time course. 

```
DE_results_filtered <- DE_results[rowSds(as.matrix(DE_results[3:8]),na.rm=TRUE) != 0,]
DE_results_filtered2 <- DE_results_filtered[rowSums(is.na(DE_results_filtered[3:8])) < 2,]

#Extract log data from DE miRNA
DE_miRNA <- row.names(DE_results_filtered2)
normalized_counts_log <- log2(normalized_counts + 1)
DE_miRNA_data <- normalized_counts_log[rownames(normalized_counts_log) %in% DE_miRNA,]
```

- Also, microRNA with a mean delta log2 among the maximum and minimum expression over time < 1 were discarded.

```
DE_miRNA_data_cond <- data.frame(t(DE_miRNA_data[,1:18]))
DE_miRNA_data_cond <- cbind('condition' = de$condition, DE_miRNA_data_cond)
mean_exp <- data.frame(DE_miRNA_data_cond %>% group_by(condition) %>% summarise_each(mean))
mean_exp <- mean_exp[,-1]
min_DE_miRNA_data <- apply(mean_exp, 2, FUN=min)
max_DE_miRNA_data <- apply(mean_exp, 2, FUN=max)
deltalog <- as.numeric(max_DE_miRNA_data) - as.numeric(min_DE_miRNA_data)

DE_miRNA_data2 <- DE_miRNA_data[deltalog > 1,]
DE_results_filtered3 <- DE_results_filtered2[row.names(DE_results_filtered2) %in% row.names(DE_miRNA_data2),]
```

*DE_miRNA_data2* : RNentropy data of under and over expressed microRNAs.
*DE_results_filtered3* : Expression of DE microRNAs.
