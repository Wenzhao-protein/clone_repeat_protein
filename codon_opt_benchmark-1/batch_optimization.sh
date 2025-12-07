#!/bin/bash
#SBATCH -J step5_codon_opt
#SBATCH --cpus-per-task=8
#SBATCH -t 23:55:00
#SBATCH -c 8
#SBATCH --mem=24g
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wenzhao.dai@chem.uzh.ch
#SBATCH --error=/scratch/wdai/conda/Wenzhao-protein/cyclic_assembly_stacking/step5_mpnn/log/%A_%a.err
#SBATCH --output=/scratch/wdai/conda/Wenzhao-protein/cyclic_assembly_stacking/step5_mpnn/log/%A_%a.out

module load anaconda3
source activate codon_opt

# Determine if the script is running under SLURM and set TASK_ID
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Running outside of SLURM. It is a test run. Setting TASK_ID manually."
    TASK_ID=1
else
    echo "Running under SLURM.  It is a production run. Using SLURM_ARRAY_TASK_ID."
    TASK_ID=$(($SLURM_ARRAY_TASK_ID))
fi

export HYDRA_FULL_ERROR=1


# Get the task command from the tasks file
task=$(sed -n "${TASK_ID}p" tasks)

# Execute the task command
eval $task