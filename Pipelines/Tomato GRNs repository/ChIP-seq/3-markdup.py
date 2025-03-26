from glob import glob
from os import system as s, makedirs
import os

makedirs("3-Dedup", exist_ok=True)

if __name__ == "__main__":
    for bam in glob("2-SORT/*.bam"):
        name = os.path.basename(bam)
        with open("a.sh", "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH -J {name}
module load picard

java -jar /media/data/opt/picard/picard.jar MarkDuplicates -I {bam} -O 3-Dedup/{name.replace(".bam", "_dedup.bam")} -M 3-Dedup/{name.replace(".bam", "_metrics.txt")}
""")
        s("sbatch a.sh")
