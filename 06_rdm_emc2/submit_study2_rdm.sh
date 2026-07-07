#!/bin/bash
## ============================================================================
## submit_study2_rdm.sh — SLURM job that runs the real Study-2 RDM fit
##
## Aim:     Submit study2_rdm.R to a compute node (16 CPUs / 32G / <=24h) so the
##          expensive hierarchical RDM fit runs off the login node.
## Inputs:  study2_rdm.R (and its data ../data/final_df_study2.rds)
## Outputs: emc_RDM_CS_study2.RData, Table1_RDM_CS_study2.html; logs/<name>_<id>.{out,err}
## Usage:   sbatch submit_study2_rdm.sh
##          (status: squeue -u $USER ; cancel: scancel <JOBID>)
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (RDM / EMC2 robustness
## analysis). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

#SBATCH --job-name=rdm_study2
#SBATCH --qos=1day
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@example.com

set -euo pipefail

## --- environment -----------------------------------------------
## R is already on PATH in this sciCORE setup (vscode profile), so no
## module loading is needed. Uncomment if running from a plain login node.
#module purge
#module load R/4.4.2-gfbf-2024a   # must match the R under which EMC2 was installed

## Point R at the exact user library where EMC2 was installed.
## (Literal %p/%v placeholders are unreliable across R startup modes.)
## CLUSTER/USER-SPECIFIC: adjust this path to your own R user library.
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
export OMP_NUM_THREADS=1          # EMC2 handles its own parallelism
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

## --- diagnostics -----------------------------------------------
echo "===================================================="
echo "Job ID        : $SLURM_JOB_ID"
echo "Job name      : $SLURM_JOB_NAME"
echo "Node          : $SLURMD_NODENAME"
echo "CPUs per task : $SLURM_CPUS_PER_TASK"
echo "Partition     : $SLURM_JOB_PARTITION"
echo "Submit dir    : $SLURM_SUBMIT_DIR"
echo "Started       : $(date)"
echo "===================================================="

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs

## --- diagnostic: where will R look for packages? --------------
Rscript --no-save -e '
cat("R_LIBS_USER env:", Sys.getenv("R_LIBS_USER"), "\n")
cat(".libPaths():\n"); print(.libPaths())
cat("MASS:", "MASS" %in% rownames(installed.packages()),
    "  EMC2:", "EMC2" %in% rownames(installed.packages()), "\n")
'

## --- run -------------------------------------------------------
## Use --no-save (not --vanilla) so Rprofile.site from the R module still
## adds the module's site library to .libPaths() — that is where MASS and
## other "Recommended" packages live.
Rscript --no-save study2_rdm.R

echo "===================================================="
echo "Finished      : $(date)"
echo "===================================================="
