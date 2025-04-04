# ATAC-seq Processing Pipeline (from FASTQ to Peaks)

This repository provides a modular and reproducible pipeline to process ATAC-seq data from raw `.fq.gz` files to final peak calling using MACS2. It is designed to run efficiently on SLURM-based HPC systems.

---

## Overview

| Step | Script           | Description                                                               |
|------|------------------|---------------------------------------------------------------------------|
| 1    | `1-mapping.py`       | Maps reads to reference genome using Bowtie2.                            |
| 2    | `2-sortSam.py`   | Sorts and converts SAM to BAM using Samtools.                 |
| 3    | `3-filterBAM.py` + `3.1-SendFiltering.sh`    | Marks duplicates and filters BAM files using Sambamba       |
| 4    | `4-Peakcalling.sh` + `4.1-SendPcalling.sh`     |  Calls peaks using MACS2 from filtered BAM files.               |

---

##  Requirements


- Python 3+
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/)
- [Samtools](http://www.htslib.org/)
- [Picard](https://broadinstitute.github.io/picard/)
- [Sambamba](https://github.com/biod/sambamba)
- [MACS2](https://github.com/macs3-project/MACS)
- SLURM-compatible job scheduler

---

## Input Format

All FASTQ files should be placed inside the `Files/` directory and named like:

```
Sample1_R1.fq.gz
Sample1_R2.fq.gz
Sample2.fq.gz    # single-end
```

---

## Output Structure

- `1-SAM/` → Bowtie2 SAM files
- `2-SORT/` → Sorted and converted BAM files
- `3-Filtered/` → Deduplicated BAMs filtered
- `4-Peaks/` → MACS2 peak output (.narrowPeak)

---

##  Notes
- Before starting, make sure `Sol4.fa` (the reference genome) is in the working directory
- Consider adding SLURM job dependencies (`--dependency=afterok:jobID`) for stricter execution order.
- Remember to make each .sh script executable (chmod +x script.sh).
- Modify the paths to Picard and Java modules based on your cluster environment.
- Adjust genome size in MACS2 for different species.


---

##  Master Script

You can execute the pipeline sequentially with the following script:

```bash
#!/bin/bash

# Master ATAC-seq Processing Pipeline

echo "STEP 1: Mapping reads with Bowtie2..."
python3 1-mapping.py

echo "STEP 2: Sorting SAM files and converting to BAM..."
python3 2-sortSam.py

echo "STEP 3: Removing  duplicates and filtering
bash 3-filterBAM.py 
sbatch 3.1-SendFiltering.sh

## Author

This pipeline is a modification of by Reynoso et al., (2019), , maintained by the Vidal Lab.
part of:
Organ-level Gene Regulatory Network models enable the identification of central transcription factors in Solanum lycopersicum (2025-04-01)
doi: https://doi.org/10.1101/2025.03.26.645553
