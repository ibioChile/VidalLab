# microRNA Differential Expression

This pipeline explains how to find differentially expressed microRNAs among different conditions, including the processing of small RNA libraries. In this case, we will analyze an experiment where samples from *Arabidopsis thaliana* were collected at 6 different days (3 replicates each day). We will explain 2 different methods for temporal DE analysis. 

1. Remove adapters from miRNA libraries using cutadapt. Reads with length <18 bp and >28 bp are discarded. 

```for file in fastq/*.fastq; do base=${file##*/}; cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 18 -M 28 --discard-untrimmed -o adapter/${base%.*}.filtered.fastq $file; done```

2. Dowload miRBase [*hairpin.fa*](ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.zip) and [*miRNA.str*](ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.zip) files. Create a folder *DB/* and store files there.

2. We use the [miraligner](https://code.google.com/p/seqbuster/wiki/miraligner) tool to map reads to microRNAs from Arabidopsis. The following script allows up to one mismatch during mapping.

```for file in adapter/*.fastq; do base=${file##*/}; java -jar miraligner.jar -sub 1 -trim 3 -add 3 -s ath -i $file -db DB/ -o miraligner_out/${base%.*}; done```

3. Format *.mirna* files to use file as input for next step. This script fixes the 'Freq' columns, which is changed from 0 to 1.

```for file in /miraligner_out/*.mirna; do base=${file##*/}; awk -F$'\t' 'BEGIN {OFS = FS} {$3="1";print}' $file > /miraligner_fixed/${base%.*}.fixed.mirna; sed -i '' -e 's/name\t1/name\tfreq/g' /miraligner_fixed/${base%.*}.fixed.mirna; done```

3. The output *.mirna* files are then imported into R. A counts table can be generated using the [*isoCounts.R*](https://github.com/ibioChile/VidalLab/blob/master/Scripts/isoCounts.R) script. This R script uses the [isoCounts](http://lpantano.github.io/isomiRs/reference/isoCounts.html) function of the [isomiRs](https://bioconductor.org/packages/release/bioc/html/isomiRs.html) package.

```Rscript isomirs.R```


