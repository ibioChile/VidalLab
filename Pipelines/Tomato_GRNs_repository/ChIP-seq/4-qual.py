from glob import glob
from os import system as s, makedirs
import os

makedirs("4-Filtered", exist_ok=True)

if __name__ == "__main__":
    for dedup in glob("3-Dedup/*_dedup.bam"):
        name = os.path.basename(dedup)
        with open("a.sh", "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH -J {name}
module load samtools

samtools view -bq 10 {dedup} > 4-Filtered/qFiltered_{name}
""")
        s("sbatch a.sh")
