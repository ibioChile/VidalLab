# microRNA Differential Expression

This pipeline explains how to determine differentially expressed microRNAs, starting from the libraries analysis. In this case, we will analyze an experiment were samples from *Arabidopsis thaliana* were collected during 6 different days (3 replicates each day). We will explain 2 different methods for temporal DE analysis. 

1. Remove adapters from miRNA libraries using cutadapt. Reads with length <18 bp and >28 bp are discarded. 

```for file in fastq/*.fastq; do base=${file##*/}; cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 18 -M 28 --discard-untrimmed -o adapter/${base%.*}.filtered.fastq $file; done```

2. Dowload miRBase [hairpin.fa](ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.zip) and [miRNA.str](ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.zip) files. Store on *DB* folder.

2. Reads are processed with the [miraligner](https://code.google.com/p/seqbuster/wiki/miraligner) tool. This tool map reads to the miRNA database of Arabidopsis, this cinfiguration allows up to one mismatch.

```for file in adapter/*.fastq; do base=${file##*/}; java -jar miraligner.jar -sub 1 -trim 3 -add 3 -s ath -i $file -db DB/ -o miraligner_out/${base%.*}; done```

3. The output .mirna files are then imported into R. A counts table can be generated using the [isoCounts](http://lpantano.github.io/isomiRs/reference/isoCounts.html) function of the [isomiRs](https://bioconductor.org/packages/release/bioc/html/isomiRs.html) package.


