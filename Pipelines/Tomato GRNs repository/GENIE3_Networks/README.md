# Tomato GENIE3 Pipeline for organ-level Gene Regulatory Networks

This repository provides a modular and reproducible pipeline to construct gene regulatory networks (GRNs) from bulk RNA-seq count data using the [GENIE3](https://bioconductor.org/packages/GENIE3) algorithm. 
It is optimized for *Solanum lycopersicum* and designed for execution on HPC systems with SLURM support.

---

## Overview

| Step | Script               | Description                                                                 |
|------|----------------------|-----------------------------------------------------------------------------|
| 1    | `01-run_GENIE3.R`    | Runs the GENIE3 algorithm on an expression matrix and list of TFs.          |
| 2    | `02-filter_network.R`| Filters the raw GENIE3 output to retain top-scoring TF-target interactions. |

---

## Requirements

- **OS**: Ubuntu 22.04+
- **RAM**: ≥ 60 GB recommended for large datasets
- **Disk space**: ≥ 2 TB for large-scale RNA-seq input
- **R version**: ≥ 4.1.0

### Required R packages

Install required R packages manually or via script:

```r
install.packages("data.table")
install.packages("dplyr")
install.packages("igraph")
install.packages("influential")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GENIE3")
---

## Input Format

Place your expression matrix and TF list inside the data/ directory.

data/
├── tissue1_counts.txt       # Raw counts table: rows = genes, columns = SRA library IDs
└── CuratedlistofTFs.txt     # Plain-text list of transcription factor gene IDs (one per line)

- The first column in the counts table must contain gene IDs
- No normalization is required prior to running GENIE3

---

## Output Structure

- results/raw_pGRN_tissue1.txt → Complete GENIE3 interaction list
- results/filtered-pGRN_tissue1.txt → Top 2% high-confidence interactions


---

##  Notes
- You can construct condition-specific pGRNs by inputting tissue-specific or treatment-specific count matrices.
- Filtering thresholds (e.g., top 1% or 2%) can be adjusted in the script.
- Suitable for any species with RNA-seq coverage and TF annotation

---

##  Master Script

You can execute the pipeline sequentially with the following script:

bash
#!/bin/bash

# Master pGRN Pipeline

echo "STEP 1: Running GENIE3 algorithm..."
Rscript scripts/01-run_GENIE3.R

echo "STEP 2: Filtering GENIE3 output to generate final pGRN..."
Rscript scripts/02-filter_network.R

---

## Citation

Organ-level Gene Regulatory Network models enable the identification of central transcription factors in Solanum lycopersicum (2025-04-01)
doi: https://doi.org/10.1101/2025.03.26.645553
