# ChIP-seq Processing Pipeline (from FASTQ to Peaks)

This repository provides a modular and reproducible pipeline to process ChIP-seq data from raw `.fq.gz` files to final peak calling using MACS2. It is designed to run efficiently on SLURM-based HPC systems.

---
## Overview

| Step | Script         | Description                                                               |
|------|----------------|---------------------------------------------------------------------------|
| 1    | `1-map.py`     | Maps reads to reference genome using Bowtie2.                            |
| 2    | `2-sortsam.py` | Sorts and converts SAM to BAM using Picard and Samtools.                 |
| 3    | `3-markdup.py`  | Removes PCR duplicates from BAM files using Picard MarkDuplicates.       |
| 4    | `4-qual.py`     | Filters BAM by mapping quality (MAPQ ≥ 10) using Samtools.               |
| 5    | `5-macs.sh`     | Calls peaks using MACS2 from filtered BAM files.                         |

---
## Requirements
- Python 3+
- Bowtie2
- Samtools
- Picard
- MACS2
- SLURM-compatible scheduler
---

---
Output Structure
1-SAM/         → Bowtie2 SAM files
2-SORT/        → Sorted and converted BAM files
3-Dedup/       → Deduplicated BAMs and metrics
4-Filtered/    → BAMs filtered by MAPQ ≥ 10
5-macs/        → MACS2 peak output (.narrowPeak)
---

---
Notes
Modify paths to Picard and Java modules based on your cluster.
The pipeline assumes paired reads are named with _R1 and _R2 or follow a consistent prefix convention.
Consider adding job dependencies for sequential control if needed
---
---
#!/bin/bash

# Master ChIP-seq Processing Pipeline
echo "STEP 1: Mapping reads with Bowtie2..."
python3 1-map2.py
echo "STEP 2: Sorting SAM files and converting to BAM..."
python3 2-sortsam2.py
echo "STEP 3: Removing PCR duplicates with Picard..."
python3 3-markdup.py
echo "STEP 4: Quality filtering (MAPQ ≥ 10)..."
python3 4-qual.py
echo "STEP 5: Peak calling using MACS2..."
bash 5-macs.sh
echo "Pipeline launched. Monitor job queue with 'squeue' or 'sacct'."
---
