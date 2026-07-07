#!/usr/bin/env Rscript
## ============================================================================
## run_recovery.R — main entry point for one parameter-recovery run
##
## Aim:     Simulate n_datasets synthetic experiments from known raw-scale
##          parameters (drawn over plausible ranges), refit each with the same
##          Stan model used in the fitting stage (folder 03), and compare
##          recovered vs true parameters. This validates that the hierarchical
##          DDM models can recover their parameters from data of the size and
##          trial structure of the real study.
## Inputs:  --data <preprocessed Study RDS> (trial template); Stan sources in
##          stan/; helper R/ modules; CLI options (see below).
## Outputs: Under <out>/: sim_dataset_{d}.rds, fit_dataset_{d}.rds,
##          recovery_long.csv, recovery_stats.csv, recovery.pdf/.png,
##          recovery_caption.txt, run_meta.rds.
## Usage:
##   Rscript run_recovery.R \
##       --model   cpt_ccss_n_r_a_s | mv_ccss_n_r_a_s | cpt_cs_sp_dr_skew | mv_cs_sp_dr_skew \
##       --data    data/final_df_study2.rds \
##       --n_datasets 5 --n_subjects 10 \
##       --chains  4 --warmup 1000 --sampling 1000 --seed 2026 \
##       --parallel_datasets 1 \
##       --out     results/cpt_ccss_n_r_a_s_study2
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM parameter
## recovery). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(parallel)
  library(cmdstanr)
  library(posterior)
  library(ggplot2)
})

## Working directory must be this recovery folder (05_ddm_recovery/): set it
## with `cd` before invoking Rscript. All source() and file paths below are
## relative to this folder.

source("R/priors.R")
source("R/transforms.R")
source("R/trials.R")
source("R/simulate.R")
source("R/fit.R")
source("R/summarize.R")

## ---- CLI ------------------------------------------------------------------
opt_list <- list(
  make_option("--model",    type = "character"),
  make_option("--data",     type = "character"),
  make_option("--family",   type = "character", default = NULL,
              help = "Override auto-detected family (cpt_ccss/mv_ccss/cpt_cs/mv_cs)"),
  make_option("--n_datasets", type = "integer", default = 5),
  make_option("--n_subjects", type = "integer", default = 10),
  make_option("--chains",   type = "integer", default = 4),
  make_option("--parallel_chains", type = "integer", default = 4),
  make_option("--warmup",   type = "integer", default = 2000),
  make_option("--sampling", type = "integer", default = 2000),
  make_option("--adapt_delta", type = "double", default = 0.95),
  make_option("--max_treedepth", type = "integer", default = 12),
  make_option("--seed",     type = "integer", default = 2026),
  make_option("--parallel_datasets", type = "integer", default = 1),
  make_option("--template_subject", type = "character", default = NULL),
  make_option("--out",      type = "character", default = "results/recovery_run"),
  make_option("--overwrite", action = "store_true", default = FALSE,
              help = "Re-run datasets that already have sim+fit RDS files cached.")
)
## When invoked via Rscript, args come from commandArgs().
## When source()d from an R session, set RECOVERY_ARGS (character vector)
## in the global env before sourcing — e.g.
##     RECOVERY_ARGS <- c("--model", "cpt_ccss_n_r_a_s",
##                        "--data",  "data/final_df_study2.rds",
##                        "--n_datasets", "2", "--n_subjects", "10",
##                        "--out",   "results/smoke")
##     source("run_recovery.R")
.override_args <- if (exists("RECOVERY_ARGS", envir = .GlobalEnv))
  get("RECOVERY_ARGS", envir = .GlobalEnv) else NULL
opt <- parse_args(OptionParser(option_list = opt_list),
                  args = if (is.null(.override_args))
                           commandArgs(trailingOnly = TRUE)
                         else .override_args)

stopifnot(!is.null(opt$model), !is.null(opt$data))
if (!opt$model %in% names(PRIORS)) {
  stop("Unknown --model: ", opt$model, "\nKnown: ",
       paste(names(PRIORS), collapse = ", "))
}
if (is.null(opt$family)) {
  opt$family <- sub("_(n_r_a_s|sp_dr_skew|sp_dr)$", "", opt$model)
}

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

cat("Recovery run\n")
cat("  Model:              ", opt$model, "\n")
cat("  Family:             ", opt$family, "\n")
cat("  Data:               ", opt$data, "\n")
cat("  N datasets:         ", opt$n_datasets, "\n")
cat("  N subjects/dataset: ", opt$n_subjects, "\n")
cat("  Output dir:         ", opt$out, "\n\n")

## ---- load trials template once --------------------------------------------
trials_tpl <- load_trials_template(opt$data, opt$family,
                                   template_subject = opt$template_subject)
cat("Trial template:", nrow(trials_tpl), "trials from subject ",
    unique(trials_tpl$subject), "\n")

## ---- compile once, reuse for all datasets ---------------------------------
mod <- compile_model(opt$model, cache_dir = "stan_cache")

## ---- set up parallelism over datasets -------------------------------------
## Each dataset already runs `parallel_chains` chains. Only set
## parallel_datasets > 1 if you have (parallel_datasets * parallel_chains) cores.
##
## We use parallel::mclapply (fork-based) rather than future_lapply because
## forked workers inherit the entire parent R environment automatically —
## so helper functions like softplus, pweight, and EFFECTIVE[[]] work
## without any manual globals-export. mclapply works on Linux/macOS (sciCORE
## is Linux). On Windows it falls back to sequential with a warning.
use_fork <- opt$parallel_datasets > 1 && .Platform$OS.type != "windows"

param_names <- PRIORS[[opt$model]]$params

## ---- one-dataset workflow --------------------------------------------------
run_one <- function(d) {
  seed_d   <- opt$seed + d * 997L
  sim_path <- file.path(opt$out, sprintf("sim_dataset_%d.rds", d))
  fit_path <- file.path(opt$out, sprintf("fit_dataset_%d.rds", d))

  ## Resumability: if both cached files exist and --overwrite not set,
  ## reload them and rebuild the merged table without re-simulating or re-fitting.
  if (!opt$overwrite && file.exists(sim_path) && file.exists(fit_path)) {
    cat(sprintf("[dataset %d] SKIP - using cached sim + fit\n", d))
    sim    <- readRDS(sim_path)
    cached <- readRDS(fit_path)
    post_tbl <- posterior_table(cached$draws, param_names)
    merged <- merge_truth(post_tbl,
                          true_params = sim$true_params,
                          mu_true     = sim$mu_true,
                          sigma_true  = sim$sigma_true) %>%
      mutate(dataset = d)
    return(merged)
  }

  t0 <- Sys.time()
  cat(sprintf("[dataset %d] simulating (seed=%d)...\n", d, seed_d))
  sim <- simulate_dataset(model_key = opt$model, family = opt$family,
                          trials = trials_tpl, L = opt$n_subjects,
                          seed = seed_d)
  saveRDS(sim, sim_path, compress = "xz")

  stan_data <- build_stan_data(sim, family = opt$family, model_key = opt$model)
  cat(sprintf("[dataset %d] fitting (N=%d, L=%d)...\n",
              d, stan_data$N, stan_data$L))
  fit_result <- fit_recovery(
    mod = mod, stan_data = stan_data, model_key = opt$model,
    chains = opt$chains, parallel_chains = opt$parallel_chains,
    iter_warmup = opt$warmup, iter_sampling = opt$sampling,
    adapt_delta = opt$adapt_delta, max_treedepth = opt$max_treedepth,
    seed = seed_d, refresh = 0, show_messages = FALSE
  )
  saveRDS(list(draws = fit_result$draws,
               summary = fit_result$summary,
               diagnostics = fit_result$diagnostics),
          fit_path, compress = "xz")

  post_tbl <- posterior_table(fit_result$draws, param_names)
  merged <- merge_truth(post_tbl,
                        true_params = sim$true_params,
                        mu_true     = sim$mu_true,
                        sigma_true  = sim$sigma_true) %>%
    mutate(dataset = d)

  t1 <- Sys.time()
  cat(sprintf("[dataset %d] done in %.1fs (divergences=%d)\n",
              d, as.numeric(difftime(t1, t0, units = "secs")),
              sum(fit_result$diagnostics$num_divergent)))
  merged
}

## ---- run all datasets ------------------------------------------------------
t_start <- Sys.time()
if (use_fork) {
  per_dataset <- parallel::mclapply(
    seq_len(opt$n_datasets), run_one,
    mc.cores = opt$parallel_datasets,
    mc.preschedule = FALSE,     # dynamic scheduling for heterogeneous runtimes
    mc.set.seed = TRUE
  )
  ## mclapply returns error objects rather than throwing; surface the first one
  errs <- sapply(per_dataset, inherits, what = "try-error")
  if (any(errs)) stop("Dataset ", which(errs)[1], " failed: ",
                      attr(per_dataset[[which(errs)[1]]], "condition")$message)
} else {
  per_dataset <- lapply(seq_len(opt$n_datasets), run_one)
}
t_end <- Sys.time()
cat(sprintf("\nAll %d datasets done in %.1fs.\n",
            opt$n_datasets,
            as.numeric(difftime(t_end, t_start, units = "secs"))))

merged_all <- bind_rows(per_dataset)
readr::write_csv(merged_all, file.path(opt$out, "recovery_long.csv"))

## per-dataset stats
stats <- merged_all %>%
  group_by(dataset, level, param) %>%
  summarise(n = n(),
            r = if (n() > 2) suppressWarnings(cor(true_value, est_mean)) else NA_real_,
            mae  = mean(abs(true_value - est_mean)),
            rmse = sqrt(mean((true_value - est_mean)^2)),
            .groups = "drop")
readr::write_csv(stats, file.path(opt$out, "recovery_stats.csv"))

## aggregated recovery plot (participant level only, all datasets pooled)
## APA 7: no title inside the figure — title and note go in the caption.
p <- plot_participant_recovery(merged_all)

## write APA caption to accompany the figure
writeLines(
  recovery_caption(opt$model,
                   n_datasets = opt$n_datasets,
                   n_subjects = opt$n_subjects,
                   figure_number = 1),
  con = file.path(opt$out, "recovery_caption.txt")
)
## dynamic figure size: 3 columns wide, one row per 3 params
n_params <- length(PRIORS[[opt$model]]$params)
plot_ncol <- min(3, n_params)
plot_nrow <- ceiling(n_params / plot_ncol)
fig_w <- 3.5 * plot_ncol + 1
fig_h <- 3.5 * plot_nrow + 0.6
pdf_dev <- tryCatch({
  grDevices::cairo_pdf(tempfile(fileext = ".pdf"), width = 1, height = 1)
  grDevices::dev.off()
  cairo_pdf
}, error = function(e) "pdf", warning = function(e) "pdf")
ggsave(file.path(opt$out, "recovery.pdf"), p,
       width = fig_w, height = fig_h, device = pdf_dev)
ggsave(file.path(opt$out, "recovery.png"), p,
       width = fig_w, height = fig_h, dpi = 300)

## meta
saveRDS(list(
  options   = opt,
  timestamp = Sys.time(),
  elapsed   = as.numeric(difftime(t_end, t_start, units = "secs"))
), file.path(opt$out, "run_meta.rds"))

## compact summary to stdout
cat("\n=== Recovery summary (pooled across datasets) ===\n")
merged_all %>%
  group_by(level, param) %>%
  summarise(r = suppressWarnings(cor(true_value, est_mean)),
            mae = mean(abs(true_value - est_mean)),
            .groups = "drop") %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\nFiles written to:", opt$out, "\n")
