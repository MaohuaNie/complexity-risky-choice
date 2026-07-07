#!/usr/bin/env Rscript
## ============================================================================
## scripts/aggregate_partial.R — rebuild recovery summaries from partial runs
##
## Aim:     Reconstruct the pooled recovery artifacts from whatever
##          sim_dataset_*.rds + fit_dataset_*.rds pairs exist in a results dir,
##          for use when run_recovery.R was killed (wall-time, scancel, OOM)
##          before its post-loop aggregation ran.
## Inputs:  <results_dir> containing matched sim/fit RDS pairs; R/priors.R and
##          R/summarize.R.
## Outputs: in <results_dir>/: recovery_long.csv, recovery_stats.csv,
##          recovery.pdf/.png, recovery_caption.txt.
## Usage:   Rscript scripts/aggregate_partial.R results/cpt_ccss_n_r_a_s_study2_ds50_L30
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM parameter
## recovery). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/aggregate_partial.R <results_dir>")
out_dir <- args[1]
if (!dir.exists(out_dir)) stop("Not a directory: ", out_dir)

source("R/priors.R")
source("R/summarize.R")

basename_dir <- basename(normalizePath(out_dir, mustWork = TRUE))
model_key <- names(PRIORS)[sapply(names(PRIORS),
                                  function(k) startsWith(basename_dir, k))]
if (length(model_key) != 1)
  stop("Could not infer a unique model_key from '", basename_dir, "'")

cat(sprintf("Model: %s    Dir: %s\n", model_key, out_dir))

param_names <- PRIORS[[model_key]]$params

## ---- find finished datasets -----------------------------------------------
sim_files <- list.files(out_dir, pattern = "^sim_dataset_\\d+\\.rds$",
                         full.names = TRUE)
sim_idx   <- as.integer(sub(".*sim_dataset_(\\d+)\\.rds$", "\\1", sim_files))

complete <- Filter(function(d) {
  file.exists(file.path(out_dir, sprintf("sim_dataset_%d.rds", d))) &&
  file.exists(file.path(out_dir, sprintf("fit_dataset_%d.rds", d)))
}, sort(sim_idx))

if (length(complete) == 0)
  stop("No matched sim/fit pairs found in ", out_dir)

cat(sprintf("Found %d finished datasets: %s\n",
            length(complete),
            paste(range(complete), collapse = "..")))

## ---- assemble merged truth+estimate per dataset ---------------------------
per_dataset <- lapply(complete, function(d) {
  sim    <- readRDS(file.path(out_dir, sprintf("sim_dataset_%d.rds", d)))
  cached <- readRDS(file.path(out_dir, sprintf("fit_dataset_%d.rds", d)))
  post_tbl <- posterior_table(cached$draws, param_names)
  merge_truth(post_tbl,
              true_params = sim$true_params,
              mu_true     = sim$mu_true,
              sigma_true  = sim$sigma_true) %>%
    mutate(dataset = d)
})

merged_all <- bind_rows(per_dataset)

## ---- write artifacts ------------------------------------------------------
write_csv(merged_all, file.path(out_dir, "recovery_long.csv"))

stats <- merged_all %>%
  group_by(dataset, level, param) %>%
  summarise(
    n    = n(),
    r    = if (n() > 2) suppressWarnings(cor(true_value, est_mean)) else NA_real_,
    mae  = mean(abs(true_value - est_mean)),
    rmse = sqrt(mean((true_value - est_mean)^2)),
    .groups = "drop"
  )
write_csv(stats, file.path(out_dir, "recovery_stats.csv"))

p <- plot_participant_recovery(merged_all)

writeLines(
  recovery_caption(model_key,
                   n_datasets    = length(complete),
                   n_subjects    = length(unique(merged_all$id[merged_all$level == "participant"])),
                   figure_number = 1),
  con = file.path(out_dir, "recovery_caption.txt")
)

n_params <- length(param_names)
plot_ncol <- min(3, n_params)
plot_nrow <- ceiling(n_params / plot_ncol)
fig_w <- 3.5 * plot_ncol + 1
fig_h <- 3.5 * plot_nrow + 0.6

pdf_dev <- tryCatch({
  grDevices::cairo_pdf(tempfile(fileext = ".pdf"), width = 1, height = 1)
  grDevices::dev.off()
  cairo_pdf
}, error = function(e) "pdf", warning = function(e) "pdf")

ggsave(file.path(out_dir, "recovery.pdf"), p,
       width = fig_w, height = fig_h, device = pdf_dev)
ggsave(file.path(out_dir, "recovery.png"), p,
       width = fig_w, height = fig_h, dpi = 300)

cat("\nWrote:\n")
cat("  ", file.path(out_dir, "recovery_long.csv"),  "\n")
cat("  ", file.path(out_dir, "recovery_stats.csv"), "\n")
cat("  ", file.path(out_dir, "recovery.pdf"),       "\n")
cat("  ", file.path(out_dir, "recovery.png"),       "\n")

## ---- pooled summary -------------------------------------------------------
cat("\n=== Pooled recovery (across", length(complete), "datasets) ===\n")
merged_all %>%
  group_by(level, param) %>%
  summarise(
    n   = n(),
    r   = suppressWarnings(cor(true_value, est_mean, use = "complete.obs")),
    mae = mean(abs(true_value - est_mean), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(level, param) %>%
  as.data.frame() %>%
  print(digits = 3)
