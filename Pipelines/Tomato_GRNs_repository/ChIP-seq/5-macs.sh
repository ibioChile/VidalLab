#!/bin/bash
mkdir -p 5-macs

for i in $(ls 4-Filtered/ | grep -v INPUT); do
    id=$(echo $i | perl -pe "s/.+(ID\d+).+/\1/")
    input=$(ls 4-Filtered/ | grep $id | grep INPUT)
    libname=$(echo $i | sed 's/.bam//; s/^.*-//')

    script="macs_$libname.sh"
    echo "#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH -J $libname
module load macs2

macs2 callpeak -t 4-Filtered/$i -c 4-Filtered/$input -f BAM -g 7.6e8 -p 0.01 -n peaks_${libname} --outdir 5-macs
" > $script

    chmod +x $script
    sbatch $script
done
