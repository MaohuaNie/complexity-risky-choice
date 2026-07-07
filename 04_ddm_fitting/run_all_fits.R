#!/usr/bin/env Rscript
## ============================================================================
## run_all_fits.R — fit every registered model on every requested study
##
## Aim:     Serial driver that loops over the selected studies x registry rows
##          and calls fit_one() on each, catching per-model errors so one
##          failure doesn't abort the batch. Resumable: a (study, family, model)
##          already fit under results/ is skipped unless --overwrite is set.
##          The single-machine alternative to the SLURM array script.
## Inputs:  CLI flags (--studies, optional --families/--models subsets, sampler
##          settings); study RDS files (data/final_df_<study>.rds); Stan files
##          under stan/. Run from the project root.
## Outputs: fit_one() artefacts under results/<study>/<family>/<model>/, plus an
##          appended run log at logs/run_all_fits.log.
## Usage:   Rscript run_all_fits.R --studies study1,study2,study3
##          Rscript run_all_fits.R --studies study2 --families cpt_ccss,mv_ccss
##          Rscript run_all_fits.R --studies study2 --models baseline,n_r_a_s
##          Or source()'d from R with RUN_ALL_ARGS set.
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM fitting stage).
## Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

source("R/model_registry.R")
source("R/preprocess.R")
source("R/inits.R")
source("R/fit.R")

verify_registry("stan")

opt_list <- list(
  make_option("--studies",  type = "character", default = "study2",
              help = "Comma-separated list: study1, study2, study3"),
  make_option("--families", type = "character", default = NULL,
              help = "Comma-separated subset of families. Default: all four."),
  make_option("--models",   type = "character", default = NULL,
              help = "Comma-separated subset of model suffixes. Default: all."),
  make_option("--chains",           type = "integer", default = 4),
  make_option("--parallel_chains",  type = "integer", default = 4),
  make_option("--warmup",           type = "integer", default = 2000),
  make_option("--sampling",         type = "integer", default = 2000),
  make_option("--adapt_delta",      type = "double",  default = 0.95),
  make_option("--max_treedepth",    type = "integer", default = 12),
  make_option("--seed",             type = "integer", default = 2026),
  make_option("--overwrite",        action = "store_true", default = FALSE),
  make_option("--stop_on_error",    action = "store_true", default = FALSE)
)

.override_args <- if (exists("RUN_ALL_ARGS", envir = .GlobalEnv))
  get("RUN_ALL_ARGS", envir = .GlobalEnv) else NULL
opt <- parse_args(
  OptionParser(option_list = opt_list),
  args = if (is.null(.override_args)) commandArgs(trailingOnly = TRUE) else .override_args
)

studies  <- trimws(strsplit(opt$studies, ",")[[1]])
families <- if (is.null(opt$families)) {
  unique(MODEL_REGISTRY$family)
} else {
  trimws(strsplit(opt$families, ",")[[1]])
}
models <- if (is.null(opt$models)) NULL else trimws(strsplit(opt$models, ",")[[1]])

reg <- MODEL_REGISTRY |> filter(family %in% families)
if (!is.null(models)) reg <- reg |> filter(model %in% models)

if (nrow(reg) == 0) stop("No models selected.")

dir.create("logs", showWarnings = FALSE)
LOG_PATH <- "logs/run_all_fits.log"

.append_log <- function(..., sep = " ") {
  msg <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., sep = sep)
  cat(msg, "\n", file = LOG_PATH, append = TRUE)
  cat(msg, "\n")
}

.append_log(sprintf("START  studies=%s  families=%s  models=%s  total_tasks=%d",
                    paste(studies, collapse=","),
                    paste(families, collapse=","),
                    if (is.null(models)) "ALL" else paste(models, collapse=","),
                    length(studies) * nrow(reg)))

total <- length(studies) * nrow(reg)
done <- 0
skipped <- 0
failed <- 0

for (study in studies) {
  data_rds <- sprintf("data/final_df_%s.rds", study)
  if (!file.exists(data_rds)) {
    .append_log("SKIP  missing data:", data_rds)
    skipped <- skipped + nrow(reg)
    next
  }
  for (i in seq_len(nrow(reg))) {
    family <- reg$family[i]; model <- reg$model[i]
    out_dir <- file.path("results", study, family, model)
    fit_path <- file.path(out_dir, "fit.rds")

    if (file.exists(fit_path) && !opt$overwrite) {
      skipped <- skipped + 1
      .append_log(sprintf("SKIP  %-8s %-8s %-15s  already fit (%s)",
                          study, family, model, fit_path))
      next
    }

    result <- tryCatch({
      fit_one(
        study = study, family = family, model = model,
        data_rds = data_rds, out_dir = out_dir,
        chains = opt$chains, parallel_chains = opt$parallel_chains,
        iter_warmup = opt$warmup, iter_sampling = opt$sampling,
        adapt_delta = opt$adapt_delta, max_treedepth = opt$max_treedepth,
        seed = opt$seed, refresh = 0, overwrite = opt$overwrite
      )
      "OK"
    }, error = function(e) conditionMessage(e))

    if (identical(result, "OK")) {
      done <- done + 1
      .append_log(sprintf("OK    %-8s %-8s %-15s", study, family, model))
    } else {
      failed <- failed + 1
      .append_log(sprintf("FAIL  %-8s %-8s %-15s  err=%s",
                          study, family, model, result))
      if (isTRUE(opt$stop_on_error)) stop(result)
    }
  }
}

.append_log(sprintf("END    done=%d  skipped=%d  failed=%d  total=%d",
                    done, skipped, failed, total))
