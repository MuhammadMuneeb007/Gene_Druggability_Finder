#!/bin/bash
#SBATCH --job-name=Feature10FPocket
#SBATCH --nodes=1
#SBATCH --partition=general
#SBATCH --time=24:00:00
#SBATCH --output=feature10.%A_%a.out
#SBATCH --error=feature10.%A_%a.err
#SBATCH --array=1-200
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1

 
echo "============================================================"
echo "FEATURE 10 FPOCKET / POCKET GEOMETRY ARRAY JOB"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Working directory: $(pwd)"
echo "============================================================"
 

# Optional: check fpocket.
echo "Python: $(which python)"
echo "fpocket: $(which fpocket || true)"

# Main run.
# Default mode uses geometry proxies only.
# Add --run-fpocket if you want actual fpocket pocket scores.
#python Feature10-FPocket.py "${SLURM_ARRAY_TASK_ID}" \
#    --chunk-size 100 \
#    --hgnc databases/HGNC/hgnc_complete_set.txt \
#    --feature4-dir feature4_structure \
#    --outdir feature10_pocket_geometry

# If you want to run real fpocket, use this instead:
python Feature10-FPocket.py "${SLURM_ARRAY_TASK_ID}" \
     --chunk-size 100 \
     --hgnc databases/HGNC/hgnc_complete_set.txt \
     --feature4-dir feature4_structure \
     --outdir feature10_pocket_geometry \
     --run-fpocket

echo "============================================================"
echo "Finished: $(date)"
echo "============================================================"