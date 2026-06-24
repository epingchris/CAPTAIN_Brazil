#!/bin/bash
#SBATCH -A naiss2026-4-88
#SBATCH -N 1
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -t 1:00:00
#SBATCH -J sdm_pred_eval_batch
#SBATCH -a 1-5
#SBATCH --output=/cfs/klemming/projects/snic/brazilaf_sdm/logs/sdm_pred_eval_%a.log
#SBATCH --error=/cfs/klemming/projects/snic/brazilaf_sdm/logs/sdm_pred_eval_%a.err

# Load required modules
module load PDC/24.11
module load R/4.4.2-cpeGNU-24.11
module load gdal/3.10.0-cpeGNU-24.11
module load nano/7.2

# Export cluster environment variables
export SCRATCH_PATH="/cfs/klemming/scratch/e/epingrau"
export PROJECT_PATH="/cfs/klemming/projects/snic/brazilaf_sdm"
export HOME_PATH="/cfs/klemming/home/e/epingrau"

# Create necessary directories
mkdir -p "$SCRATCH_PATH/pred_endangered_60k"
mkdir -p "$PROJECT_PATH/logs"
mkdir -p "$PROJECT_PATH/CAPTAIN_Brazil_out"

# Log job info
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task: $SLURM_ARRAY_TASK_ID"
echo "Hostname: $(hostname)"
echo "Available CPUs: $SLURM_CPUS_PER_TASK"
echo "Allocated Memory: ${SLURM_MEM_PER_NODE}MB"
echo "Date: $(date)"
echo "Working dir: $(pwd)"
echo "=========================================="

# Change to project directory
cd "$HOME_PATH/CAPTAIN_Brazil" || exit 1

# Run the R fitting script
ls renv/   # will show if renv folder is visible
Rscript sdm_3_pred_eval.r $SLURM_ARRAY_TASK_ID
EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Job completed with exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "=========================================="

exit $EXIT_CODE
