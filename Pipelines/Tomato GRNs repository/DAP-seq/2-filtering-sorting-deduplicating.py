from glob import glob
from os import system as s
import os

# Pipeline to filter, processes, and deduplicates alignments from SAM to sorted/indexed BAMs.
# Ensure output directories exist
s("mkdir -p 2-BAM")

if __name__ == "__main__":
    files = glob("./1-SAM/*.sam")
    for file in files:
        name2 = os.path.basename(file).replace(".sam", "")
        job_script = f"a_{name2}.sh"
        log_file = f"2-BAM/{name2}.out"
        err_file = f"2-BAM/{name2}.err"

        with open(job_script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH -c 10\n")
            f.write("#SBATCH --mem=30G\n")
            f.write(f"#SBATCH -J {name2}\n")
            f.write(f"#SBATCH -o {log_file}\n")
            f.write(f"#SBATCH -e {err_file}\n\n")

            # Step 1: Filter SAM
            f.write(f"samtools view -h {file} | awk '$1 ~ /^@/ || ($0 ~ /XM:i:[012][^0-9]/ && $5 > 30)' | grep -v 'XS:i:' > 1-SAM/{name2}.filtered.sam\n\n")

            # Step 2: Convert to BAM
            f.write(f"samtools view -bS 1-SAM/{name2}.filtered.sam > 2-BAM/{name2}.filtered.bam\n")

            # Step 3: Name sort
            f.write(f"samtools sort -n 2-BAM/{name2}.filtered.bam -o 2-BAM/{name2}.filtered.sorted.bam\n")

            # Step 4: Fixmate
            f.write(f"samtools fixmate -m 2-BAM/{name2}.filtered.sorted.bam 2-BAM/{name2}.fixmate.bam\n")

            # Step 5: Coordinate sort
            f.write(f"samtools sort 2-BAM/{name2}.fixmate.bam -o 2-BAM/{name2}.fixmate.sorted.bam\n")

            # Step 6: Mark duplicates
            f.write(f"samtools markdup -r 2-BAM/{name2}.fixmate.sorted.bam 2-BAM/{name2}.dedup.bam\n")

            # Step 7: Final sort and index
            f.write(f"samtools sort 2-BAM/{name2}.dedup.bam -o 2-BAM/{name2}.sorted.bam\n")
            f.write(f"samtools index 2-BAM/{name2}.sorted.bam\n")

            # Clean up intermediate BAMs 
            f.write(f"rm 2-BAM/{name2}.filtered.bam 2-BAM/{name2}.fixmate.bam 2-BAM/{name2}.dedup.bam\n")

        # Submit the SLURM script
        s(f"sbatch {job_script}")
