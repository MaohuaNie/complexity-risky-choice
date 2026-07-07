#!/bin/bash
# ============================================================================
# submit_fit_array.sh — SLURM array driver for the DDM fitting stage
#
# Aim:     Submit one array task per (study, family, model). Each task maps its
#          SLURM_ARRAY_TASK_ID to a registry row (optionally restricted to the
#          CPT or MV family class), then calls run_fit.R to fit that one model
#          with 4 cores (one per chain). Already-fit models no-op almost
#          instantly, so re-submitting a partially-done study is safe.
# Inputs:  Positional args STUDY [FAMILY_CLASS]; the model registry (R/) and
#          Stan files (stan/); study RDS files (data/final_df_<study>.rds).
# Outputs: run_fit.R artefacts under results/<study>/<family>/<model>/; SLURM
#          logs under logs/; a transient task-mapping file under /tmp.
# Usage:   sbatch --array=0-39 --time=3-12:00:00 --qos=<your_qos> \
#              scripts/submit_fit_array.sh study1 cpt
#          (see the USAGE block below for the full per-class recipe)
#
# CLUSTER-SPECIFIC — you MUST adapt before running on your own cluster:
#   * Submit from the project root so run_fit.R's relative "R/", "stan/",
#     "data/", "results/" paths resolve (SLURM starts jobs in the submit dir).
#     If you submit from elsewhere, add: cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
#   * Set --qos to a real QoS on your scheduler (the value below is a placeholder).
#   * Uncomment / edit the `module load R/...` line to match your cluster's
#     module naming.
# ----------------------------------------------------------------------------
# Part of the complexity-under-risk replication package (DDM fitting stage).
# Pipeline order and dependencies are documented in ../../README.md.
# ============================================================================
## ---------------------------------------------------------------------------
## SLURM array job: one array task per (study, family, model). Each task
## fits one hierarchical model with 4 cores (one per chain).
##
## Registry (after skew removal) = 60 models per study:
##   rows  1..16 : cpt_ccss       (CPT class)
##   rows 17..32 : cpt_ccss_7o    (CPT class)
##   rows 33..48 : mv_ccss        (MV  class)
##   rows 49..52 : cpt_cs         (CPT class)
##   rows 53..56 : cpt_cs_7o      (CPT class)
##   rows 57..60 : mv_cs          (MV  class)
##
## CPT fits take ~3 days wall-clock; MV fits take ~1–1.5 days. Because of
## this, submit CPT and MV as SEPARATE sbatch jobs with different
## --time/--qos overrides — see the USAGE block below.
##
## Already-fit models are skipped by fit_one() (file.exists(fit.rds) check),
## so re-submitting a partially-done study is safe — finished tasks no-op
## almost instantly and only pending ones re-sample.
##
## USAGE:
##
##   ./submit_fit_array.sh STUDY [FAMILY_CLASS]
##      STUDY         : study1 | study2 | study3  (required)
##      FAMILY_CLASS  : cpt | mv | all            (default: all)
##
##   When FAMILY_CLASS is cpt or mv, the script restricts the registry to
##   that subset and the array indices 0..(n-1) map to the filtered rows:
##      cpt : 40 rows per study (cpt_ccss, cpt_ccss_7o, cpt_cs, cpt_cs_7o)
##      mv  : 20 rows per study (mv_ccss, mv_cs)
##      all : 60 rows per study
##
##   Submit each study twice — once per class — with the right time budget:
##
##     # CPT (~3 days):
##     sbatch --array=0-39 --time=3-12:00:00 --qos=1week \
##         scripts/submit_fit_array.sh study1 cpt
##     sbatch --array=0-39 --time=3-12:00:00 --qos=1week \
##         scripts/submit_fit_array.sh study2 cpt
##     sbatch --array=0-39 --time=3-12:00:00 --qos=1week \
##         scripts/submit_fit_array.sh study3 cpt
##
##     # MV (~1–1.5 days):
##     sbatch --array=0-19 --time=1-12:00:00 --qos=1week \
##         scripts/submit_fit_array.sh study1 mv
##     sbatch --array=0-19 --time=1-12:00:00 --qos=1week \
##         scripts/submit_fit_array.sh study2 mv
##     sbatch --array=0-19 --time=1-12:00:00 --qos=1week \
##         scripts/submit_fit_array.sh study3 mv
##
##   Replace --qos=1week with the real sciCORE qos name that covers the
##   requested --time (check with `sinfo -o "%P %q"`). Add `%N` at the end
##   of --array (e.g. `--array=0-39%20`) to cap concurrent tasks.
##
## SBATCH defaults below are placeholders for the "all-class" case; the
## recommended workflow is to override --time and --qos on the command
## line per class, as shown above.
## ---------------------------------------------------------------------------
#SBATCH --job-name=fit
#SBATCH --time=3-12:00:00                 # 3.5 d default; override per class
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --qos=1week                       # PLACEHOLDER — real sciCORE qos name
#SBATCH --output=logs/fit_%A_%a.out
#SBATCH --error=logs/fit_%A_%a.err

set -euo pipefail

STUDIES=${1:-study2}
FAMILY_CLASS=${2:-all}
export STUDIES FAMILY_CLASS               # Rscript subprocess reads via Sys.getenv()

case "$FAMILY_CLASS" in
  cpt|mv|all) ;;
  *) echo "FAMILY_CLASS must be one of: cpt | mv | all (got '$FAMILY_CLASS')" >&2; exit 2 ;;
esac

mkdir -p logs

## Load R on sciCORE (edit for your cluster's module naming):
# module load R/4.4.2-foss-2024a

## Map SLURM_ARRAY_TASK_ID → (study, family, model). The registry is
## optionally restricted to the CPT or MV subset; the filtered grid is
## study-major so indices 0..(N-1) cover all (n_filtered × n_studies) tasks.
TASK_FILE="/tmp/fit_task_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"

Rscript -e '
  suppressPackageStartupMessages(library(dplyr))
  source("R/model_registry.R")
  verify_registry("stan")

  reg <- MODEL_REGISTRY
  fc  <- Sys.getenv("FAMILY_CLASS", "all")
  if (fc == "cpt") reg <- filter(reg, grepl("^cpt_", family))
  if (fc == "mv")  reg <- filter(reg, grepl("^mv_",  family))

  studies <- strsplit(Sys.getenv("STUDIES", "study2"), ",")[[1]]
  grid <- expand.grid(i = seq_len(nrow(reg)),
                      study = studies,
                      stringsAsFactors = FALSE)           # study-major
  grid$family <- reg$family[grid$i]
  grid$model  <- reg$model[grid$i]

  ## Drop (study, family) combinations whose data structure does not match:
  ##   - study1 uses 7-outcome lotteries    → only *_7o CPT families apply
  ##   - study2, study3 use 2-outcome only → only non-*_7o CPT families apply
  ## MV families (mv_ccss, mv_cs) use moment-summary covariates and apply
  ## to all studies regardless of outcome count.
  studies_7o <- "study1"
  is_cpt <- grepl("^cpt_", grid$family)
  is_7o  <- grepl("_7o$",  grid$family)
  drop <- is_cpt & (
      (is_7o  & !(grid$study %in% studies_7o)) |
      (!is_7o &  (grid$study %in% studies_7o))
  )
  grid <- grid[!drop, , drop = FALSE]

  idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) + 1L
  if (idx < 1 || idx > nrow(grid))
    stop(sprintf("Array index %d out of range [1..%d] (class=%s, n_tasks=%d)",
                 idx, nrow(grid), fc, nrow(grid)))
  r <- grid[idx, ]
  cat(r$study, r$family, r$model, sep = "\t"); cat("\n")
' > "$TASK_FILE"

read -r TASK_STUDY TASK_FAMILY TASK_MODEL < "$TASK_FILE"

## Rename the running task so SLURM's completion email says e.g.
##   "Name=fit_study2_mv_cs_sp Ended, State=COMPLETED, ExitCode 0"
## instead of the abstract "Name=fit". Requires scontrol (no extra deps).
scontrol update JobId="$SLURM_JOB_ID" \
  JobName="fit_${TASK_STUDY}_${TASK_FAMILY}_${TASK_MODEL}" 2>/dev/null || true

echo "Task $SLURM_ARRAY_TASK_ID  →  $TASK_STUDY / $TASK_FAMILY / $TASK_MODEL"

Rscript run_fit.R \
  --study  "$TASK_STUDY" \
  --family "$TASK_FAMILY" \
  --model  "$TASK_MODEL" \
  --chains 4 --parallel_chains 4 \
  --warmup 3000 --sampling 3000 \
  --adapt_delta 0.95 --max_treedepth 12 \
  --seed $((2026 + SLURM_ARRAY_TASK_ID))
