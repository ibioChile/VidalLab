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

6. Import counts file to R. 

```counts_isomirs <- read.table("miRNA_counts.tsv", header=TRUE, sep="\t")```

7. Counts from 5’ and 3’ strands of each miRNA should be added.

```
miRNA <- row.names(counts_isomirs)
miRNA <- gsub("-3p","",miRNA)
miRNA <- gsub("-5p","",miRNA)

counts_isomirs_add <- cbind("miRNA" = miRNA, counts_isomirs)
counts_isomirs_add <- as.data.frame(counts_isomirs_add %>% group_by(miRNA) %>% summarise_each(sum))
rownames(counts_isomirs_add) <- counts_isomirs_add[,1]
counts_isomirs_add <- counts_isomirs_add[,-1]
```

7. Filter counts: Only miRNA with a minimum of 5 counts in at least 25% of samples (5 samples in this case) are considered for the DE analysis.

```
x <- counts_isomirs_add >= 5
counts_isomirs_add <- counts_isomirs_add[rowSums(x == TRUE) >= 5,]
```

## Approach 1: One-way analysis of variance

- Raw read counts were first median-normalized to adjust for the effect of library sizes and read count distributions (Anders & Huber 2010). 
- Normalized counts were converted to the log2 scale, using log2 (x + 1) for the conversion.
- An analysis of variance (ANOVA) was carried out for normalized reads, the variables modeled were sampling day and miRNA abundance.
- The Benjamini-Hochberg procedure (Benjamini & Hochberg, 1995)  was used to control the false discovery rate (FDR) based on the p-values obtained from the ANOVA analysis. Genes having p-values with an FDR threshold < 0.05 were designated as differentially expressed.




