#!/bin/bash
## ============================================================================
## submit_all_ppc.sh — submit PPC generation for every study/family/model
##
## Aim:     Fan out generate_ppc.R (via submit_one_ppc.sh) across all
##          (study, family, model) tuples needed for the manuscript figures.
##          Per (study, family): CS = baseline, sp, dr, sp_dr (4 models);
##          CCSS = baseline, n, r, a, s, <best> (6 models), where <best> is the
##          per-(study, family) parsimonious LOO pick (rank-1 unless rank-2 is
##          statistically indistinguishable by the 2*SE rule -> n_a_s). CPT
##          family is the _7o variant in study1 and the regular one in
##          studies 2-3. Total 6 study-family combos x 10 models = 60 jobs.
## Inputs:  scripts/submit_one_ppc.sh; fitted models under results/.
## Outputs: one SLURM job per tuple, each writing a
##          posterior_predictives.csv; SLURM logs under logs/.
## Usage:   bash 04_ppc/scripts/submit_all_ppc.sh
## ----------------------------------------------------------------------------
## Cluster-specific: submission uses SLURM (sbatch); adapt for your site.
## Part of the complexity-under-risk replication package (posterior predictive
## checks). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================

set -e
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"

mkdir -p logs

WRAPPER="scripts/submit_one_ppc.sh"

# CCSS family models: baseline + 4 single-shift + parsimonious best
ccss_models() {
  local STUDY="$1"
  local FAM="$2"
  # Parsimony-selected best (from the supplementary table)
  local BEST="n_a_s"
  if [ "$STUDY" = "study1" ] && [ "$FAM" = "mv_ccss" ]; then
    BEST="n_r_a_s"
  fi
  echo "baseline n r a s $BEST"
}

# CS family models: baseline + 2 single-mechanism + full
cs_models() {
  echo "baseline sp dr sp_dr"
}

# Per-study CPT family is _7o in study1, regular in studies 2 and 3
cpt_ccss_for_study() {
  if [ "$1" = "study1" ]; then echo "cpt_ccss_7o"; else echo "cpt_ccss"; fi
}
cpt_cs_for_study() {
  if [ "$1" = "study1" ]; then echo "cpt_cs_7o"; else echo "cpt_cs"; fi
}

N_SUBMIT=0

for STUDY in study1 study2 study3; do
  CPT_CCSS_FAM=$(cpt_ccss_for_study "$STUDY")
  CPT_CS_FAM=$(cpt_cs_for_study "$STUDY")

  # ---- CCSS condition: CPT + MV families ----
  for FAM in "$CPT_CCSS_FAM" "mv_ccss"; do
    for MODEL in $(ccss_models "$STUDY" "$FAM"); do
      echo "Submitting CCSS:  $STUDY / $FAM / $MODEL"
      sbatch "$WRAPPER" "$STUDY" "$FAM" "$MODEL"
      N_SUBMIT=$((N_SUBMIT + 1))
    done
  done

  # ---- CS condition: CPT + MV families ----
  for FAM in "$CPT_CS_FAM" "mv_cs"; do
    for MODEL in $(cs_models); do
      echo "Submitting CS:    $STUDY / $FAM / $MODEL"
      sbatch "$WRAPPER" "$STUDY" "$FAM" "$MODEL"
      N_SUBMIT=$((N_SUBMIT + 1))
    done
  done
done

echo
echo "Submitted $N_SUBMIT PPC jobs."
echo "Check with: squeue -u \"\$USER\" -h -o '%.12i %.6T %.10M %j' | grep ppc"
