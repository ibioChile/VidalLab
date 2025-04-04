#!/bin/bash

####################
# Loop through all BAM files and generate MACS2 scripts
# (Exclude the input/control BAM file)
####################

path="/home2/jfernandez/DAP-seq/2-BAM"
peaks_dir="/home2/jfernandez/DAP-seq/3-Peaks"
control_bam="/home2/jfernandez/DAP-seq/2-BAM/input.input.dedup.bam"

mkdir -p "$peaks_dir"

m=1

for i in "$path"/*.bam; do
    base_name=$(basename "$i" .bam)
    # Skip the input control itself
    if [[ "$base_name" == "input.input" ]]; then
        continue
    fi

    script_name="script-peaks$m.sh"

    echo "#!/bin/bash

module load macs2
macs2 callpeak -t $i \
    -c $control_bam \
    -f BAMPE -g 7.6e8 -p 0.0001 --call-summits -B \
    --outdir $peaks_dir -n $base_name" > "$script_name"

    chmod +x "$script_name"
    m=$((m + 1))
done
