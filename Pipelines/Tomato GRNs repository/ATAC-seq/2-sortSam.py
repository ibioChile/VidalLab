from glob import glob
from os import system as s
import os

s("mkdir -p 2-SORT")

if __name__ == "__main__":
    files = glob("./1-SAM/*.sam")
    
    for file in files:
        base_name = os.path.basename(file)
        sample_id = os.path.splitext(base_name)[0]
        script_name = f"sort_{sample_id}.sh"
        
        sorted_sam = f"2-SORT/{sample_id}.sorted.sam"
        bam_file = f"2-SORT/{sample_id}.bam"
        
        with open(script_name, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH -c 15\n")
            f.write("#SBATCH --mem=30G\n")
            f.write(f"#SBATCH -J {sample_id}\n")
            f.write("module load picard samtools\n")
            f.write(f"java -jar /media/data/opt/picard/picard.jar SortSam --TMP_DIR temp -I {file} -O {sorted_sam} --SO coordinate\n")
            f.write(f"samtools view -@ 15 -S -b {sorted_sam} -o {bam_file}\n")
            f.write(f"rm {sorted_sam}\n")
        
        print(f"Submitting {script_name}...")
        s(f"sbatch {script_name}")
