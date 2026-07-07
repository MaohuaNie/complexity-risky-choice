#!/bin/bash
## ============================================================================
## submit_one_ppc.sh — SLURM wrapper: generate PPC for one study/family/model
##
## Aim:     Submit one generate_ppc.R run as a SLURM job (tail 800 / thin 2).
## Inputs:  positional args STUDY FAMILY MODEL; the fitted model under
##          results/<study>/<family>/<model>/.
## Outputs: results/<study>/<family>/<model>/posterior_predictives.csv;
##          SLURM logs under logs/.
## Usage:   sbatch 04_ppc/scripts/submit_one_ppc.sh study1 mv_ccss n_a_s
## ----------------------------------------------------------------------------
## Cluster-specific: the `module load R/...` line and the --mail-user/--qos
## SBATCH directives target the scicore SLURM cluster; adapt them to your site.
## Part of the complexity-under-risk replication package (posterior predictive
## checks). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================
#SBATCH --job-name=ppc
#SBATCH --output=logs/ppc_%x_%j.out
#SBATCH --error=logs/ppc_%x_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --qos=1day
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.com

## Usage: sbatch submit_one_ppc.sh STUDY FAMILY MODEL
## Example: sbatch submit_one_ppc.sh study1 mv_ccss n_a_s

set -e
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
module purge && module load R/4.4.2-foss-2024a

STUDY="$1"; FAM="$2"; MODEL="$3"
echo "=== PPC $STUDY/$FAM/$MODEL at $(date +%H:%M:%S) ==="

Rscript post/generate_ppc.R \
  --study "$STUDY" --family "$FAM" --model "$MODEL" \
  --tail_per_chain 800 --thin 2 \
  --cores "$SLURM_CPUS_PER_TASK"

echo "=== Done at $(date +%H:%M:%S) ==="
