#!/bin/bash
#SBATCH --job-name=fixedX
#SBATCH --mail-type=END
#SBATCH --mail-user=ssh.sentian@gmail.com
#SBATCH --output=/home/st1864/RAICc/output/%A_%a.out # master job id %A and array-task job id %a
#SBATCH --array=1-81
#SBATCH --time=02:00:00

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB


module purge
module load r/intel/3.6.0

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
cd /home/st1864/RAICc/code

Rscript run.R 'FALSE' 'TRUE' $SLURM_ARRAY_TASK_ID