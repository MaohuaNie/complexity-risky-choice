#!/bin/bash
## ============================================================================
## make_rtdist_best.sh — RT-distribution PPC plots for the 12 best models
##
## Aim:     Produce the 12 best-model RT-distribution figures (3 studies x 2
##          conditions x {MVS, CPT}). For each combo, if the model's
##          posterior_predictives.csv is missing it is built with
##          generate_ppc.R first, then make_rtdist_one.R renders the single
##          ppc_rt_distribution_<model> plot (no clobbering of other figures).
## Inputs:  the (study|family|best-model|condition) combos listed below;
##          results/<study>/<family>/<model>/posterior_predictives.csv (built
##          on demand).
## Outputs: results/<study>/<family>/ppc_rt_distribution_<model>.pdf/.png;
##          SLURM logs under logs/.
## Usage:   sbatch 04_ppc/scripts/make_rtdist_best.sh
## ----------------------------------------------------------------------------
## Cluster-specific: the `module load R/...` line and the --mail-user/--qos
## SBATCH directives target the scicore SLURM cluster; adapt them to your site.
## Part of the complexity-under-risk replication package (posterior predictive
## checks). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================
#SBATCH --job-name=rtdist
#SBATCH --output=logs/rtdist_%j.out
#SBATCH --error=logs/rtdist_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --qos=1day
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.com

set -e
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
module purge && module load R/4.4.2-foss-2024a
mkdir -p logs

## study | family | best model | condition
combos=(
  "study1|mv_ccss|n_r_a_s|ccss"
  "study1|mv_cs|sp_dr|cs"
  "study1|cpt_ccss_7o|n_a_s|ccss"
  "study1|cpt_cs_7o|sp_dr|cs"
  "study2|mv_ccss|n_a_s|ccss"
  "study2|mv_cs|sp_dr|cs"
  "study2|cpt_ccss|n_a_s|ccss"
  "study2|cpt_cs|sp_dr|cs"
  "study3|mv_ccss|n_a_s|ccss"
  "study3|mv_cs|sp_dr|cs"
  "study3|cpt_ccss|n_a_s|ccss"
  "study3|cpt_cs|sp_dr|cs"
)

for c in "${combos[@]}"; do
  IFS='|' read -r S FAM M COND <<< "$c"
  ppc="results/$S/$FAM/$M/posterior_predictives.csv"
  echo "==================================================================="
  echo "=== $S / $FAM / $M ($COND)   $(date +%H:%M:%S) ==="
  if [ ! -f "$ppc" ]; then
    echo "  posterior_predictives.csv missing -> generate_ppc.R"
    Rscript post/generate_ppc.R --study "$S" --family "$FAM" --model "$M" \
      --tail_per_chain 800 --thin 2 --cores "$SLURM_CPUS_PER_TASK"
  else
    echo "  posterior_predictives.csv present -> skip generate"
  fi
  Rscript post/make_rtdist_one.R --study "$S" --family "$FAM" --model "$M" --cond "$COND"
done

echo "=== all 12 done $(date +%H:%M:%S) ==="
