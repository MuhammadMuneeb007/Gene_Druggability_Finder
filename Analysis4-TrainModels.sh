#!/bin/bash
#SBATCH --job-name=DrugML
#SBATCH --nodes=1
#SBATCH --partition=general
#SBATCH --time=24:00:00
#SBATCH --output=logs/ml_train.%A_%a.out
#SBATCH --error=logs/ml_train.%A_%a.err
#SBATCH --array=1-440
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8



#python Analysis4-TrainModels.py "$SLURM_ARRAY_TASK_ID" --no-permutation

python Analysis4-TrainModels.py "$SLURM_ARRAY_TASK_ID"   --no-permutation




 