from glob import glob
import os

input_dir = "./2-SORT"
output_dir = "./3-Filtered"
os.makedirs(output_dir, exist_ok=True)

bam_files = sorted(glob(f"{input_dir}/*.bam"))

for idx, bam_file in enumerate(bam_files, start=1):
    sample = os.path.splitext(os.path.basename(bam_file))[0]
    markdup_bam = f"{input_dir}/{sample}-markdup.bam"
    filtered_bam = f"{output_dir}/{sample}-markdup-flt.bam"
    log_file = f"{input_dir}/{sample}_sambamba_log.log"

    with open(f"script-filt{idx}.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"# Mark duplicates and filter for {sample}\n")
        f.write("module load sambamba samtools\n")
        f.write(f"sambamba-1.0.1-linux-amd64-static markdup -t 2 {bam_file} {markdup_bam} > {log_file}\n")
        f.write(f"samtools view -@ 2 -bF 1804 -q 20 {markdup_bam} -o {filtered_bam}\n")
        f.write(f"rm {markdup_bam}\n")
    os.chmod(f"script-filt{idx}.sh", 0o755)
