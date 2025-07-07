#!/bin/bash

path=./3-Filtered
outdir=./4-Peaks
mkdir -p "$outdir"

m=1

for bam in "$path"/*.bam
do
    base=$(basename "$bam" .bam)

    echo "#!/bin/bash
#SBATCH --mem=30G
#SBATCH -c 4
#SBATCH -J peak-$base
#SBATCH -o $outdir/${base}_macs2.out

module load macs2

macs2 callpeak -t $bam -f BAMPE -n $base --outdir $outdir --nomodel  -q 0.05" > script-peaks$m.sh

    chmod +x script-peaks$m.sh
    m=$((m + 1))
done
