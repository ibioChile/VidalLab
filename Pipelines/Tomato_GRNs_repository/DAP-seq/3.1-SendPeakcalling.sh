#!/bin/bash

#SBATCH --mem=30G
#SBATCH -c 4
#SBATCH -J qual-$SLURM_ARRAY_TASK_ID
#SBATCH --array=1-10%10
#SBATCH --nice=214748364
#SBATCH --output=zpeak_%a.out
#SBATCH --error=zpeak_%a.err

cd /foldername/

./script-peaks$SLURM_ARRAY_TASK_ID.sh