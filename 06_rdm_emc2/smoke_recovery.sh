#!/bin/bash
## ============================================================================
## smoke_recovery.sh — tiny end-to-end check of the RDM recovery pipeline
##
## Aim:     Sanity-check that simulate -> fit -> summary runs end-to-end
##          (including participant-level alpha extraction) with a tiny config
##          (n_subj = 10, chains = 2, burnin = 200, samples = 200) before
##          committing the full 30-task recovery array. NOT a real recovery;
##          the numbers are not meaningful.
## Inputs:  recovery_simulate.R, recovery_fit.R (and ../data/final_df_study2.rds)
## Outputs: recovery/smoke_*.{rds,RData}; logs/smoke_%j.{out,err}
## Usage:   sbatch smoke_recovery.sh
##          (watch: tail -f logs/smoke_*.out)
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (RDM / EMC2 robustness
## analysis). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

#SBATCH --job-name=rdm_smoke
#SBATCH --qos=30min
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/smoke_%j.out
#SBATCH --error=logs/smoke_%j.err

set -euo pipefail

## CLUSTER/USER-SPECIFIC: adjust this path to your own R user library.
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p recovery logs

echo "===================================================="
echo "Smoke test for recovery pipeline"
echo "Started   : $(date)"
echo "===================================================="

## --- step 1: simulate (tiny) -----------------------------------
echo
echo "[1/3] Simulating tiny dataset (n_subj = 10) ..."
Rscript --no-save recovery_simulate.R \
    --seed 999 \
    --n_subj 10 \
    --out recovery/smoke_sim.rds


## --- step 2: fit (short chains) --------------------------------
echo
echo "[2/3] Fitting (2 chains, 200 burnin, 200 samples) ..."
Rscript --no-save recovery_fit.R \
    --sim     recovery/smoke_sim.rds \
    --out     recovery/smoke_fit.RData \
    --summary recovery/smoke_sum.rds \
    --burnin  200 \
    --samples 200 \
    --chains  2


## --- step 3: print summary -------------------------------------
echo
echo "[3/3] Summary"
Rscript --no-save -e '
s <- readRDS("recovery/smoke_sum.rds")
cat("\n=== Smoke-test recovery (single tiny replicate) ===\n")
cat("Group-level true vs recovered:\n")
print(s$group, digits = 3)
cat(sprintf("\nParticipant-level rows extracted: %d (across %d parameters)\n",
            nrow(s$subject), length(unique(s$subject$parameter))))
cat("Head of participant-level table:\n")
print(utils::head(s$subject), digits = 3)
cat("\n(Tiny smoke run; recovery values are NOT meaningful with n=10 and 200 samples.\n")
cat(" What matters: the pipeline ran AND the participant-level alpha extraction worked.)\n")
'

echo
echo "===================================================="
echo "Finished  : $(date)"
echo "===================================================="
