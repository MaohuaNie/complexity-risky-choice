#!/bin/bash
## ============================================================================
## scripts/submit_recovery.sh — SLURM wrapper for one recovery run
##
## Aim:     Submit a single parameter-recovery run (one model, one study) to a
##          SLURM cluster; wraps run_recovery.R with sensible resource requests.
## Inputs:  positional args <model_key> <n_datasets> <n_subjects> [study];
##          reads data/final_df_<study>.rds.
## Outputs: results/<model>_<study>_ds<N>_L<S>/ (run outputs); SLURM logs under
##          logs/.
## Usage:   sbatch scripts/submit_recovery.sh cpt_ccss_n_r_a_s 50 30 study2
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM parameter
## recovery). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================
#SBATCH --job-name=recovery
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --qos=6hours
#SBATCH --output=logs/recovery_%j.out
#SBATCH --error=logs/recovery_%j.err

set -euo pipefail

## Run from the recovery folder (05_ddm_recovery/) so relative paths resolve.
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"

MODEL=${1:-cpt_ccss_n_r_a_s}
N_DATASETS=${2:-10}
N_SUBJECTS=${3:-10}
STUDY=${4:-study2}    # study1 / study2 / study3  → data/final_df_<STUDY>.rds

DATA="data/final_df_${STUDY}.rds"
if [[ ! -f "$DATA" ]]; then
  echo "ERROR: data file not found: $DATA" >&2
  exit 1
fi

mkdir -p logs results stan_cache

## Cluster module lines: uncomment and edit to match your site's module names
## (or preload R + CmdStan in the submitting shell). The versions below are
## examples only.
# module load R/4.3.1-foss-2022a
# module load CmdStan/2.34.1

OUT_DIR="results/${MODEL}_${STUDY}_ds${N_DATASETS}_L${N_SUBJECTS}"

## Each chain uses 1 core → 4 chains × 2 parallel datasets = 8 cores.
Rscript run_recovery.R \
  --model "$MODEL" \
  --data  "$DATA" \
  --n_datasets "$N_DATASETS" \
  --n_subjects "$N_SUBJECTS" \
  --chains 4 \
  --parallel_chains 4 \
  --parallel_datasets 2 \
  --warmup 2000 \
  --sampling 2000 \
  --adapt_delta 0.95 \
  --seed 2026 \
  --out "$OUT_DIR"

echo "Finished. Results in: $OUT_DIR"
