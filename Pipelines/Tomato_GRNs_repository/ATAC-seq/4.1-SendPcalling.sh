#!/bin/bash

#SBATCH --mem=30G
#SBATCH -c 4
#SBATCH -J qual-$SLURM_ARRAY_TASK_ID
#SBATCH --array=1-N%X     # N = number of files, X = jobs in parallel
#SBATCH --nice=214748364
#SBATCH --output=zpeak_%a.out
#SBATCH --error=zpeak_%a.err

cd /Folder/

./script-peaks$SLURM_ARRAY_TASK_ID.sh
