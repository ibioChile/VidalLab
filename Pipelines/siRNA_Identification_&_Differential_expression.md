# siRNA Identification & Differential expression

This pipeline explains how to identify clusters of siRNA and identify differentially expressed clusters among different conditions. In this case, we will analyze an experiment where samples from *Arabidopsis thaliana* were collected at 6 different days (3 replicates each day). 

1. Create environment and install programs required for this pipeline.

```conda create -n sirna```

```conda activate sirna```

```conda install -c bioconda shortstack infernal bedtools seqtk bbmap```

1. Remove adapters from microRNA libraries using cutadapt. Reads with length <18 bp and >28 bp are discarded. 

```for file in fastq/*.fastq; do base=${file##*/}; cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 18 -M 28 --discard-untrimmed -o adapter/${base%.*}.filtered.fastq $file; done```

2. - Download Rfam.cm and Rfam.clanin databases from the [Rfam FTP site](ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT). Create a folder ```DB/``` and store files there. Compress Rfam database.

```cmpress Rfam.cm```

3. We use the ```cmscan``` tool from [infernal](http://eddylab.org/infernal/) to search reads belonging to ncRNA, tRNA, mitRNA, miRNA, etc (Fastq files are first converted to fasta files using ```seqtk```).

```mkdir Rfam```

```for file in adapter/*.fastq; do base=${file##*/}; seqtk seq -A $file > ${file%.*}.fasta; cmscan --cpu 28 --tblout Rfam/${base%.*}.rfam -E 0.001 --noali --rfam --fmt 2 --clanin /home/pcamejo/databases/RFAM/Rfam.clanin /home/pcamejo/databases/RFAM/Rfam.cm ${file%.*}.fasta > Rfam/${base%.*}.cmscan; done```

4. Remove reads matching Rfam families, using the command ```filterbyname.sh``` from ```bbmap```.

```mkdir Rfam_filtered```

```for file in Rfam/*.rfam; do base=${file##*/}; awk '{print $4}’ $file >  ${base%.*}.tofilt.list; filterbyname.sh in=adapter/${base%.*}.fastq names=${base%.*}.tofilt.list out=Rfam_filtered/${base%.*}.rfam.fastq; done```

5. Run [Shortstack](https://github.com/MikeAxtell/ShortStack) with each libray to annotate small RNA-producing genes

```for file in ../Rfam_filtered/*.fastq; do base=${file##*/}; ShortStack --bowtie_cores 20 --readfile $file --genomefile GCF_000001735.4_TAIR10.1_genomic.fasta --outdir shortstack/${base%.*}_ssout --mismatches 0 --nohp; done```

6. We use ```multiIntersectBed``` from ```Bedtools``` to find shared clusters in the different libraries. We will create a bash file ```multiIntersectBed.sh``` containing the script to run ```multiIntersectBed``` with all libraries.

```cd shortstack/```

```touch multiIntersectBed.sh```

```echo -n "multiIntersectBed -header -i " > multiIntersectBed.sh```

```for folder in *_ssout/; do cd $folder; gff2bed  < ShortStack_D.gff3 > ${folder%_*_*_*_*}_ShortStack_D.bed; cd ..; echo -n "$folder${folder%_*_*_*_*}_ShortStack_D.bed " >> multiIntersectBed.sh; done```

```echo "> intersect_shortstack_DC.out" >> multiIntersectBed.sh```

```bash multiIntersectBed.sh```

7. The output ```intersect_shortstack_DC.out``` was filtered so that only regions that were present in at least 2 (out of 3) libraries for at least one day were used to serve as the final reference small RNA locus boundaries. Export the filtered regions as a bed file ```intersect_shortstack_DC_2out3lib.bed``` containing a list of filtered clusters genome, start and end information, separated by a tab (This step was carried out in excel).

8. Merge locus with 2nt of difference:

``` mergeBed -d 2 -i intersect_shortstack_DC_2out3lib.bed > output_merged_intervals_file.bed```

9. ShortStack-count mode under default settings was then used to find relative small RNA abundances on this reference list of each library.

```awk '{print $1":"$2"-"$3}' output_merged_intervals_file.bed > output_merged_intervals_file_loci.bed```

```cd Rfam_filtered/```

```ShortStack --bowtie_cores 10 --readfile  1-150_S1_L000_R1_000.filtered.rfam.fastq 20-150_S11_L000_R1_000.filtered.rfam.fastq 27-140_S15_L000_R1_000.filtered.rfam.fastq 7-140_S4_L000_R1_000.filtered.rfam.fastq 13-140_S7_L000_R1_000.filter ...(list all files)  --genomefile GCF_000001735.4_TAIR10.1_genomic.fasta --outdir shortstack_filt --locifile output_merged_intervals_file_loci.bed --sort_mem 1G --nohp```

10. Extract counts of small RNA clusters with DicerCall (DC) annotation.

```cd shortstack_filt```

```awk '$12 != "N" {print $1}' Results.txt > Results_DC_locus.txt```

```while read line; do grep $line Counts.txt >> Counts_DC.txt; done < Results_DC_locus.txt```


## RNentropy

- Raw read counts are first median-normalized to adjust for the effect of library sizes and read count distributions *(Anders & Huber 2010)*. Normalized counts are converted to the log2 scale, using log2 (x + 1) for conversion.

```
library(DESeq2)

Ara_siRNA <- read.table("Counts_DC.txt", header=TRUE, row.names = 1, sep="\t")

Ara_siRNA <- Ara_siRNA[,-1:-2]

de <- data.frame(row.names=colnames(Ara_siRNA), condition = c("d5", "d5", "d5","d9","d9","d9","d13","d13","d13","d17","d17","d17","d21","d21","d21","d25","d25","d25"))

dds = DESeqDataSetFromMatrix(Ara_siRNA , de , design = ~ condition)

dds <- estimateSizeFactors(dds)

sizeFactors(dds)

Ara_siRNA_norm <- counts(dds, normalized=TRUE)
```

- RNentropy (*Zambelli et al. 2018*): identification of genes showing a significant variation of expression across all conditions studied. The samples, corresponding to different conditions sequenced in any number of replicates, are compared by taking into account the global gene expression level and at the same time the impact of biological variation across replicates. 

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

*DE_miRNA_data2* : RNentropy data of under and over expre