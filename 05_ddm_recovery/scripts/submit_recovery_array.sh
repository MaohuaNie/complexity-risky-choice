#!/bin/bash
## ============================================================================
## scripts/submit_recovery_array.sh — SLURM array over (study, model)
##
## Aim:     Run parameter recovery for several models in one array submission;
##          each array task maps to one PRIORS entry (optionally filtered to the
##          CPT or MV subset) and runs run_recovery.R for that model.
## Inputs:  positional args STUDY [FAMILY_CLASS] [N_DATASETS] [N_SUBJECTS];
##          reads R/priors.R (model registry) and data/final_df_<STUDY>.rds.
## Outputs: one results/<model>_<study>_ds<N>_L<S>/ per task; SLURM logs under
##          logs/.
## Usage:   sbatch --array=0-5 scripts/submit_recovery_array.sh study2 all 50 30
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM parameter
## recovery). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================
##
## One array task per (study, model). Each task runs one parameter-recovery
## simulation for one PRIORS entry.
##
## Registry (PRIORS keys in R/priors.R) = 6 models total:
##   rows 1..3 : cpt_ccss_n_r_a_s, cpt_cs_sp_dr, cpt_cs_sp_dr_skew   (CPT class)
##   rows 4..6 : mv_ccss_n_r_a_s,  mv_cs_sp_dr,  mv_cs_sp_dr_skew    (MV  class)
##
## Runtime (at ds=50, L=30):
##   *_ccss_n_r_a_s : a few hours
##   *_cs_sp_dr[_skew] : ~17-29 h wall-clock
## The 2-day default below therefore covers the slowest CS-family recoveries;
## fast CCSS tasks finish well within that budget.
##
## Already-finished datasets are skipped inside run_recovery.R
## (file.exists(fit_dataset_{d}.rds) / sim_dataset_{d}.rds checks), so
## re-submitting a partially-done run is safe.
##
## USAGE:
##
##   sbatch scripts/submit_recovery_array.sh STUDY [FAMILY_CLASS] [N_DATASETS] [N_SUBJECTS]
##      STUDY         : study1 | study2 | study3                   (required)
##      FAMILY_CLASS  : cpt | mv | all                             (default: all)
##      N_DATASETS    : number of simulated datasets per model     (default: 50)
##      N_SUBJECTS    : number of participants per dataset         (default: 30)
##
##   When FAMILY_CLASS is cpt or mv, the script restricts the registry to
##   that subset and the array indices 0..(n-1) map to the filtered rows:
##      cpt : 3 rows per study
##      mv  : 3 rows per study
##      all : 6 rows per study
##
##   Examples (production settings, per class):
##
##     # CPT recoveries:
##     sbatch --array=0-2 --time=2-00:00:00 --qos=1week \
##         scripts/submit_recovery_array.sh study2 cpt 50 30
##
##     # MV recoveries:
##     sbatch --array=0-2 --time=2-00:00:00 --qos=1week \
##         scripts/submit_recovery_array.sh study2 mv 50 30
##
##     # Everything at once:
##     sbatch --array=0-5 --time=2-00:00:00 --qos=1week \
##         scripts/submit_recovery_array.sh study2 all 50 30
##
##   Replace --qos=1week with the real sciCORE qos name that covers the
##   requested --time (check with `sinfo -o "%P %q"`). Add `%N` at the end
##   of --array (e.g. `--array=0-5%3`) to cap concurrent tasks.
##
## SBATCH defaults below are placeholders; override --time and --qos on the
## command line as needed.
## ---------------------------------------------------------------------------
#SBATCH --job-name=recovery
#SBATCH --time=2-00:00:00                 # 2 d default; override per class
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --qos=1week                       # PLACEHOLDER — real sciCORE qos name
#SBATCH --output=logs/recovery_%A_%a.out
#SBATCH --error=logs/recovery_%A_%a.err

set -euo pipefail

## Run from the recovery folder (05_ddm_recovery/) so relative paths resolve.
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"

STUDY=${1:-study2}
FAMILY_CLASS=${2:-all}
N_DATASETS=${3:-50}
N_SUBJECTS=${4:-30}
export STUDY FAMILY_CLASS                  # Rscript subprocess reads via Sys.getenv()

case "$FAMILY_CLASS" in
  cpt|mv|all) ;;
  *) echo "FAMILY_CLASS must be one of: cpt | mv | all (got '$FAMILY_CLASS')" >&2; exit 2 ;;
esac

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

## Map SLURM_ARRAY_TASK_ID → model. The registry is the PRIORS list from
## R/priors.R, optionally restricted to the CPT or MV subset. Models are
## sorted to give a stable index ordering (cpt_* first, then mv_*).
TASK_FILE="/tmp/recovery_task_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"

Rscript -e '
  source("R/priors.R")
  models <- sort(names(PRIORS))
  fc <- Sys.getenv("FAMILY_CLASS", "all")
  if (fc == "cpt") models <- grep("^cpt_", models, value = TRUE)
  if (fc == "mv")  models <- grep("^mv_",  models, value = TRUE)

  idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) + 1L
  if (idx < 1 || idx > length(models))
    stop(sprintf("Array index %d out of range [1..%d] (class=%s, n_tasks=%d)",
                 idx, length(models), fc, length(models)))
  cat(models[idx], "\n")
' > "$TASK_FILE"

read -r TASK_MODEL < "$TASK_FILE"

## Rename the running task so SLURM's completion email says e.g.
##   "Name=recovery_study2_cpt_cs_sp_dr Ended, State=COMPLETED, ExitCode 0"
## instead of the abstract "Name=recovery".
scontrol update JobId="$SLURM_JOB_ID" \
  JobName="recovery_${STUDY}_${TASK_MODEL}" 2>/dev/null || true

OUT_DIR="results/${TASK_MODEL}_${STUDY}_ds${N_DATASETS}_L${N_SUBJECTS}"

echo "Task $SLURM_ARRAY_TASK_ID  →  $STUDY / $TASK_MODEL  (ds=$N_DATASETS, L=$N_SUBJECTS)"

## Each chain uses 1 core → 4 chains × 8 parallel datasets = 32 cores.
Rscript run_recovery.R \
  --model "$TASK_MODEL" \
  --data  "$DATA" \
  --n_datasets "$N_DATASETS" \
  --n_subjects "$N_SUBJECTS" \
  --chains 4 \
  --parallel_chains 4 \
  --parallel_datasets 8 \
  --warmup 2000 \
  --sampling 2000 \
  --adapt_delta 0.95 \
  --seed $((2026 + SLURM_ARRAY_TASK_ID)) \
  --out "$OUT_DIR"

echo "Finished. Results in: $OUT_DIR"
