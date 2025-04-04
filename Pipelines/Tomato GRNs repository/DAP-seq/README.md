# Tomato DAP-seq Processing Pipeline (from FASTQ to Peaks)

This repository provides a modular and reproducible pipeline to process DAP-seq reads from raw `.fq.gz` FASTQ files to peak calling using MACS2. It is optimized for tomato (Solanum lycopersicum) and designed for use on SLURM-based HPC systems.

---

## Overview

| Step | Script           | Description                                                               |
|------|------------------|---------------------------------------------------------------------------|
| 1    | `1-map.py`       | Maps reads to reference genome using Bowtie2.                            |
| 2    | `2-filtering-sorting-deduplicating.py`   | Sorts and converts SAM to BAM, filter and depublicates using Samtools.                 |
| 3    | `3-MACSpeaks.sh` + `3.1-SendPcalling.sh`     |  Calls peaks using MACS2 from filtered BAM files.               |


##  Requirements


- Python 3+
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/)
- [Samtools](http://www.htslib.org/)
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
- `2-BAM/` → Deduplicated, sorted BAMs
- `3-Peaks/` → MACS2 peak files


---

##  Notes
- Before starting, make sure `Sol4.fa` (the reference genome) is in the working directory
- Consider adding SLURM job dependencies (`--dependency=afterok:jobID`) for stricter execution order.
- Remember to make each .sh script executable (chmod +x script.sh).
- Assumes `input.input.dedup.bam` is the control sample for MACS2.
- Adjust genome size (`-g 7.6e8`) in MACS2 for different species.

---

##  Master Script

You can execute the pipeline sequentially with the following script:

```bash
#!/bin/bash

# Master DNASE-seq Processing Pipeline

echo "STEP 1: Mapping reads with Bowtie2..."
python3 1-map.py

echo "STEP 2: Sorting SAM files and converting to BAM..."
python3 2-filtering-sorting-deduplicating.py

echo "STEP 3: Removing  duplicates and filtering
bash 3-MACSpeaks.sh
sbatch 3.1-SendPcalling.sh

## Author

This pipeline is a modification of by Hutin et al., (2023), maintained by the Vidal Lab.
part of:
Organ-level Gene Regulatory Network models enable the identification of central transcription factors in Solanum lycopersicum (2025-04-01)
doi: https://doi.org/10.1101/2025.03.26.645553
