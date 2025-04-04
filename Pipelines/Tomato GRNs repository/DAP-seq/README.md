# Tomato DAP-seq Processing Pipeline

This pipeline processes DAP-seq reads from raw `.fq.gz` FASTQ files to peak calling using MACS2. It is optimized for tomato (Solanum lycopersicum) and designed for use on SLURM-based HPC systems.

---

## 📌 Requirements

- `bowtie2`
- `samtools` (>=1.9)
- `macs2`
- Python 3.x
- SLURM job scheduler

---

## 🔁 Pipeline Overview

### 1. **Mapping Reads**
**Script**: `1-map.py`  
- Builds a Bowtie2 index from `Sol4.fa`
- Detects FASTQ files in `./Files/`
- Submits jobs via SLURM to align reads to the genome
- Output: SAM files → `1-SAM/`

```bash
python3 1-map.py
```

---

### 2. **Filtering, Sorting & Deduplication**
**Script**: `2-filtering-sorting-deduplicating2.py`  
- Filters reads (MAPQ > 30, XM:i tags)
- Converts to BAM
- Sorts, fixes mates, removes duplicates
- Output: Sorted, indexed BAM files → `2-BAM/`

```bash
python3 2-filtering-sorting-deduplicating2.py
```

---

### 3. **Peak Calling with MACS2**
**Scripts**:
- `3-MACSpeaks.sh`: Generates individual peak-calling scripts
- `3.1-SendPeakcalling.sh`: Submits them via SLURM array

Edit `path` and `peaks_dir` in `3-MACSpeaks.sh` accordingly.

```bash
bash 3-MACSpeaks.sh
sbatch 3.1-SendPeakcalling.sh
```

---

## 📂 Input Format

All FASTQ files should be placed inside the `Files/` directory and named like:

```
Sample1_R1.fq.gz
Sample1_R2.fq.gz
Sample2.fq.gz    # single-end
```

---

## 📤 Output

- Raw alignments: `1-SAM/`
- Deduplicated, sorted BAMs: `2-BAM/`
- MACS2 peak files: `3-Peaks/`

---

## ⚠️ Notes

- Assumes `input.input.dedup.bam` is the control sample for MACS2.
- Adjust genome size (`-g 7.6e8`) in MACS2 for different species.
- Ensure `Sol4.fa` is available for Bowtie2 index.

---

## 👨‍🔬 Author

This pipeline is maintained by the Vidal Lab.
