#!/bin/bash
## ============================================================================
## submit_recovery.sh — SLURM array running the full RDM recovery sweep
##
## Aim:     Run the recovery pipeline as a 30-task array. Each task (seed =
##          array index) simulates one Study-2 dataset (50 participants, real
##          stimuli, all group-level parameters drawn from reasonable ranges;
##          design A ~ 1, 13 params), refits with the study2_rdm.R design +
##          priors, and saves group- and participant-level true-vs-recovered
##          estimates. Aggregate afterwards with recovery_aggregate.R.
## Inputs:  recovery_simulate.R, recovery_fit.R (and ../data/final_df_study2.rds)
## Outputs: recovery/{sim,fit,sum}_seedNN.* per task; logs/recov_%A_%a.{out,err}
## Usage:   sbatch submit_recovery.sh
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (RDM / EMC2 robustness
## analysis). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

#SBATCH --job-name=rdm_recov
#SBATCH --array=1-30
#SBATCH --qos=1day
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/recov_%A_%a.out
#SBATCH --error=logs/recov_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.com

set -euo pipefail

## CLUSTER/USER-SPECIFIC: adjust this path to your own R user library.
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs recovery

idx=$SLURM_ARRAY_TASK_ID
TAG=$(printf "seed%02d" "$idx")

SIM=recovery/sim_${TAG}.rds
FIT=recovery/fit_${TAG}.RData
SUM=recovery/sum_${TAG}.rds

echo "===================================================="
echo "Array idx : $idx"
echo "Tag       : $TAG"
echo "Started   : $(date)"
echo "===================================================="

if [ ! -f "$SIM" ]; then
  echo "[1/2] Simulating ${SIM} ..."
  Rscript --no-save recovery_simulate.R --seed "$idx" --out "$SIM"
fi

if [ ! -f "$SUM" ]; then
  echo "[2/2] Fitting ${FIT} ..."
  # Long chains (burnin 3000, samples 6000) for reliable convergence.
  Rscript --no-save recovery_fit.R \
      --sim "$SIM" --out "$FIT" --summary "$SUM" \
      --burnin 3000 --samples 6000
fi

echo "===================================================="
echo "Finished  : $(date)"
echo "===================================================="
