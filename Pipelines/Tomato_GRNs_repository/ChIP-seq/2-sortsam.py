from glob import glob
from os import system as s, makedirs
import os

makedirs("2-SORT", exist_ok=True)

if __name__ == "__main__":
    for sam in glob("1-SAM/*.sam"):
        name = os.path.basename(sam)
        with open("a.sh", "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH -c 15
#SBATCH --mem=30G
#SBATCH -J {name}
module load picard
module load samtools

java -jar /media/data/opt/picard/picard.jar SortSam --TMP_DIR temp -I {sam} -O 2-SORT/{name.replace(".sam", ".sorted.sam")} --SO coordinate
samtools view -@ 15 -S -b 2-SORT/{name.replace(".sam", ".sorted.sam")} -o 2-SORT/{name.replace(".sam", ".bam")}
rm 2-SORT/{name.replace(".sam", ".sorted.sam")}
""")
        s("sbatch a.sh")
