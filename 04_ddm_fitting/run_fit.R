#!/usr/bin/env Rscript
## ============================================================================
## run_fit.R — CLI entry point to fit one hierarchical model on one study
##
## Aim:     Parse command-line options, source the pipeline modules, verify the
##          registry, resolve default data/output paths, and call fit_one() for
##          a single (study, family, model). The atomic unit invoked by both the
##          SLURM array script and run_all_fits.R.
## Inputs:  CLI flags (--study/--family/--model plus sampler settings); the
##          study RDS (default data/final_df_<study>.rds); Stan files under
##          stan/. Run from the project root so the source("R/...") paths and
##          "stan"/"data"/"results" relative paths resolve.
## Outputs: Everything fit_one() writes under
##          results/<study>/<family>/<model>/ (draws, diagnostics, etc.).
## Usage:   Rscript run_fit.R --study study2 --family cpt_ccss --model n_r_a_s \
##              --chains 4 --parallel_chains 4 --warmup 2000 --sampling 2000 \
##              --adapt_delta 0.95 --max_treedepth 12 --seed 2026
##          Or source()'d from R with FIT_ARGS set, e.g.:
##            FIT_ARGS <- c("--study","study2","--family","cpt_ccss",
##                          "--model","n_r_a_s"); source("run_fit.R")
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM fitting stage).
## Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

source("R/model_registry.R")
source("R/preprocess.R")
source("R/inits.R")
source("R/fit.R")

verify_registry("stan")

opt_list <- list(
  make_option("--study",   type = "character"),
  make_option("--family",  type = "character"),
  make_option("--model",   type = "character"),
  make_option("--data",    type = "character", default = NULL,
              help = "Path to RDS; defaults to data/final_df_<study>.rds"),
  make_option("--out",     type = "character", default = NULL,
              help = "Output dir; defaults to results/<study>/<family>/<model>/"),
  make_option("--chains",           type = "integer", default = 4),
  make_option("--parallel_chains",  type = "integer", default = 4),
  make_option("--warmup",           type = "integer", default = 2000),
  make_option("--sampling",         type = "integer", default = 2000),
  make_option("--adapt_delta",      type = "double",  default = 0.95),
  make_option("--max_treedepth",    type = "integer", default = 12),
  make_option("--seed",             type = "integer", default = 2026),
  make_option("--refresh",          type = "integer", default = 200),
  make_option("--save_warmup",      action = "store_true", default = FALSE),
  make_option("--overwrite",        action = "store_true", default = FALSE)
)

.override_args <- if (exists("FIT_ARGS", envir = .GlobalEnv))
  get("FIT_ARGS", envir = .GlobalEnv) else NULL

opt <- parse_args(
  OptionParser(option_list = opt_list),
  args = if (is.null(.override_args)) commandArgs(trailingOnly = TRUE) else .override_args
)

stopifnot(!is.null(opt$study), !is.null(opt$family), !is.null(opt$model))

if (is.null(opt$data)) opt$data <- sprintf("data/final_df_%s.rds", opt$study)
if (is.null(opt$out))  opt$out  <- file.path("results", opt$study, opt$family, opt$model)

fit_one(
  study = opt$study, family = opt$family, model = opt$model,
  data_rds = opt$data, out_dir = opt$out,
  chains = opt$chains, parallel_chains = opt$parallel_chains,
  iter_warmup = opt$warmup, iter_sampling = opt$sampling,
  adapt_delta = opt$adapt_delta, max_treedepth = opt$max_treedepth,
  seed = opt$seed, refresh = opt$refresh,
  save_warmup = opt$save_warmup, overwrite = opt$overwrite
)
