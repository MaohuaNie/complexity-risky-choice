#!/bin/bash
## ============================================================================
## scripts/submit_recovery_long.sh — SLURM wrapper for a long recovery run
##
## Aim:     Same as submit_recovery.sh but with a larger time/QoS budget for the
##          slow CS-family recoveries (~17-29 h wall-clock at ds=50, L=30).
## Inputs:  positional args <model_key> <n_datasets> <n_subjects> [study];
##          reads data/final_df_<study>.rds.
## Outputs: results/<model>_<study>_ds<N>_L<S>/ (run outputs); SLURM logs under
##          logs/.
## Usage:   sbatch scripts/submit_recovery_long.sh cpt_cs_sp_dr 50 30 study2
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM parameter
## recovery). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================
#SBATCH --job-name=recovery_long
#SBATCH --time=2-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --qos=1week
#SBATCH --output=logs/recovery_%j.out
#SBATCH --error=logs/recovery_%j.err

set -euo pipefail

## Run from the recovery folder (05_ddm_recovery/) so relative paths resolve.
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"

MODEL=${1:-cpt_cs_sp_dr}
N_DATASETS=${2:-50}
N_SUBJECTS=${3:-30}
STUDY=${4:-study2}    # study1 / study2 / study3  → data/final_df_<STUDY>.rds

DATA="data/final_df_${STUDY}.rds"
if [[ ! -f "$DATA" ]]; then
  echo "ERROR: data file not found: $DATA" >&2
  exit 1
fi

mkdir -p logs results stan_cache

## Cluster modules: this variant expects R + CmdStan to be inherited from the
## submitting shell's already-loaded modules rather than reloaded here (some
## Lmod sites block same-name autoswap). Uncomment/edit to match your site if
## you prefer explicit loads. The versions below are examples only.
# module load R/4.3.1-foss-2022a
# module load CmdStan/2.34.1

OUT_DIR="results/${MODEL}_${STUDY}_ds${N_DATASETS}_L${N_SUBJECTS}"

## Each chain uses 1 core → 4 chains × 8 parallel datasets = 32 cores.
Rscript run_recovery.R \
  --model "$MODEL" \
  --data  "$DATA" \
  --n_datasets "$N_DATASETS" \
  --n_subjects "$N_SUBJECTS" \
  --chains 4 \
  --parallel_chains 4 \
  --parallel_datasets 8 \
  --warmup 2000 \
  --sampling 2000 \
  --adapt_delta 0.95 \
  --seed 2026 \
  --out "$OUT_DIR"

echo "Finished. Results in: $OUT_DIR"
