#!/bin/bash
## ============================================================================
## run_summarize_all_ppc.sh — summarize every posterior_predictives.csv
##
## Aim:     Driver job that finds every existing posterior_predictives.csv and
##          submits one SLURM job per file to run summarize_ppc.R, writing a
##          small ppc_summary.rds next to each CSV (skips files already
##          summarized). Jobs run in parallel.
## Inputs:  results/**/posterior_predictives.csv; summarize_ppc.R.
## Outputs: results/<study>/<family>/<model>/ppc_summary.rds per input CSV;
##          SLURM logs under logs/.
## Usage:   sbatch 04_ppc/scripts/run_summarize_all_ppc.sh
## ----------------------------------------------------------------------------
## Cluster-specific: the `module load R/...` line and SBATCH directives target
## the scicore SLURM cluster; the inner per-file wrapper is written to /tmp at
## submit time. Adapt for your site.
## Part of the complexity-under-risk replication package (posterior predictive
## checks). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================
#SBATCH --job-name=ppc_summarize_driver
#SBATCH --output=logs/ppc_summarize_driver_%j.out
#SBATCH --error=logs/ppc_summarize_driver_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --qos=30min

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"

# Inner wrapper, written on the fly so each submission gets the right env.
# It re-enters the submit dir via SLURM_SUBMIT_DIR (falls back to the current
# working directory, which the driver already set to the repo root above).
cat > /tmp/run_one_summarize.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=ppc_sum
#SBATCH --output=logs/ppc_sum_%x_%j.out
#SBATCH --error=logs/ppc_sum_%x_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --qos=6hours
set -e
cd "${SLURM_SUBMIT_DIR:-$PWD}"
module purge && module load R/4.4.2-foss-2024a
Rscript post/summarize_ppc.R --study "$1" --family "$2" --model "$3"
EOF
chmod +x /tmp/run_one_summarize.sh

N=0
for ppc_csv in $(find results -name "posterior_predictives.csv"); do
  # Parse path: results/<study>/<family>/<model>/posterior_predictives.csv
  dir=$(dirname "$ppc_csv")
  MODEL=$(basename "$dir")
  FAM=$(basename "$(dirname "$dir")")
  STUDY=$(basename "$(dirname "$(dirname "$dir")")")
  out="$dir/ppc_summary.rds"
  if [ -f "$out" ]; then
    echo "Skip (exists): $STUDY/$FAM/$MODEL"
    continue
  fi
  echo "Submit: $STUDY/$FAM/$MODEL"
  sbatch /tmp/run_one_summarize.sh "$STUDY" "$FAM" "$MODEL"
  N=$((N+1))
done

echo
echo "Submitted $N summarize jobs."
echo "Check: squeue -u \"\$USER\" -h -o '%.12i %.6T %.10M %j' | grep ppc_sum"
