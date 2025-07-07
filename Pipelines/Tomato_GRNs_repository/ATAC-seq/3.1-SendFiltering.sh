#!/bin/bash

#SBATCH --mem=15G
#SBATCH -c 2
#SBATCH -J qual-$SLURM_ARRAY_TASK_ID
#SBATCH --array=1-N%X     # N = number of files, X = jobs in parallel
#SBATCH --nice=214748364
#SBATCH --output=zdup_%a.out
#SBATCH --error=zdup_%a.err

cd /Folder/

./script-filt$SLURM_ARRAY_TASK_ID.sh
