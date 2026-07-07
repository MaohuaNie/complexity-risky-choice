#!/usr/bin/env Rscript
## ============================================================================
## scripts/check_one_dataset.R — inspect recovery for a single dataset
##
## Aim:     Quick per-dataset diagnostic: load one cached sim + fit pair, print
##          sampler diagnostics and an ordered recovery table (r/mae/rmse per
##          parameter), and optionally save a true-vs-estimated scatter plot.
## Inputs:  <results_dir> <dataset_index> [--plot]; R/priors.R, R/summarize.R.
## Outputs: table printed to stdout; with --plot, check_dataset_<index>.pdf in
##          <results_dir>.
## Usage:   Rscript scripts/check_one_dataset.R results/cpt_cs_sp_dr_study2_ds50_L30 20 --plot
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM parameter
## recovery). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript scripts/check_one_dataset.R <results_dir> <dataset_index> [--plot]")
}
out_dir    <- args[1]
ds_idx     <- as.integer(args[2])
want_plot  <- "--plot" %in% args

## Infer model_key from the results dir name (first token, e.g.
## "cpt_cs_sp_dr_study2_ds50_L30" → model starts with "cpt_cs_sp_dr").
## We match against the PRIORS registry so misnamed dirs fail loudly.
source("R/priors.R")
source("R/summarize.R")

basename_dir <- basename(normalizePath(out_dir, mustWork = TRUE))
model_key <- names(PRIORS)[sapply(names(PRIORS), function(k) startsWith(basename_dir, k))]
if (length(model_key) != 1) {
  stop(sprintf("Could not infer a unique model_key from '%s'. Matches: %s",
               basename_dir, paste(model_key, collapse = ", ")))
}
cat(sprintf("Model: %s    Dataset: %d    Dir: %s\n", model_key, ds_idx, out_dir))

sim_path <- file.path(out_dir, sprintf("sim_dataset_%d.rds", ds_idx))
fit_path <- file.path(out_dir, sprintf("fit_dataset_%d.rds", ds_idx))
if (!file.exists(sim_path)) stop("Missing: ", sim_path)
if (!file.exists(fit_path)) stop("Missing: ", fit_path, " — dataset not yet finished")

sim    <- readRDS(sim_path)
cached <- readRDS(fit_path)

param_names <- PRIORS[[model_key]]$params

## ---- diagnostics -----------------------------------------------------------
diag <- cached$diagnostics
n_div <- sum(diag$num_divergent %||% NA_integer_)
max_td <- if (!is.null(diag$num_max_treedepth)) max(diag$num_max_treedepth) else NA_integer_
ebfmi  <- if (!is.null(diag$ebfmi)) round(min(diag$ebfmi), 3) else NA_real_

cat(sprintf("\nDiagnostics:\n  divergences (sum over chains): %d\n  max_treedepth hits (max)     : %s\n  min E-BFMI                   : %s\n\n",
            n_div, as.character(max_td), as.character(ebfmi)))

## ---- recovery table --------------------------------------------------------
post_tbl <- posterior_table(cached$draws, param_names)
merged <- merge_truth(
  post_tbl,
  true_params = sim$true_params,
  mu_true     = sim$mu_true,
  sigma_true  = sim$sigma_true
)

## ---- canonical parameter ordering -----------------------------------------
## Drift block first (the substantive parameters) → response block → offset.
##   1. theta          drift scaler
##   2. beta           utility curvature (CPT) / risk sensitivity (MV)
##   3. gamma / eta    probability weighting (CPT) / skew sensitivity (MV)
##   4. threshold (a)  response caution
##   5. ndt (tau)      non-decision time
##   6. sp             starting-point bias
##   7. zeta           additive drift offset
##   Deltas (context shifts) follow their base parameter.
canonical_order <- c(
  "theta_raw",       "delta_theta",
  "beta_raw",        "beta",         "delta_beta",
  "gamma_raw",       "delta_gamma",
  "eta",             "delta_eta",
  "threshold_raw",   "delta_threshold",
  "ndt_raw",
  "sp_raw",
  "zeta"
)
order_params <- function(x) {
  x <- as.character(x)
  idx <- match(x, canonical_order)
  idx[is.na(idx)] <- length(canonical_order) + seq_len(sum(is.na(idx)))
  x[order(idx)]
}

param_order_used <- order_params(unique(merged$param))

## ---- recovery table (ordered, greek-labelled) -----------------------------
stats <- recovery_stats(merged) %>%
  mutate(
    label = greek_label(param),
    param = factor(param, levels = param_order_used),
    level = factor(level, levels = c("participant", "mu", "sigma"))
  ) %>%
  arrange(level, param) %>%
  select(level, param, label, n, r, mae, rmse)

cat("Recovery (raw scale, ordered):\n")
print(stats, n = Inf)

## ---- optional scatter plot -------------------------------------------------
if (want_plot) {
  df <- merged %>%
    filter(level == "participant") %>%
    mutate(
      covers = (true_value >= q2.5) & (true_value <= q97.5),
      param  = factor(param, levels = param_order_used)
    )

  panel_letters <- letters[seq_along(param_order_used)]
  fac_labels <- setNames(
    vapply(seq_along(param_order_used), function(i) {
      sprintf('bold("(%s)")~%s',
              panel_letters[i],
              greek_label(param_order_used[i]))
    }, character(1)),
    param_order_used
  )

  badges <- df %>%
    group_by(param) %>%
    summarise(
      rho = suppressWarnings(cor(true_value, est_mean, use = "complete.obs")),
      pct = 100 * mean(covers, na.rm = TRUE),
      lo  = min(c(true_value, q2.5, est_mean), na.rm = TRUE),
      hi  = max(c(true_value, q97.5, est_mean), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      x        = lo + 0.03 * (hi - lo),
      y_rho    = hi - 0.05 * (hi - lo),
      y_cov    = hi - 0.15 * (hi - lo),
      label_rho = sprintf("\u03C1 = %.2f", rho),
      label_cov = sprintf("%.0f%% in 95%% CrI", pct)
    )

  p <- ggplot(df, aes(true_value, est_mean)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5,
                       colour = covers), alpha = 0.35, linewidth = 0.35) +
    geom_point(aes(colour = covers), size = 1.2) +
    geom_text(data = badges, aes(x = x, y = y_rho, label = label_rho),
              hjust = 0, size = 3, inherit.aes = FALSE) +
    geom_text(data = badges, aes(x = x, y = y_cov, label = label_cov),
              hjust = 0, size = 2.8, colour = "grey25", inherit.aes = FALSE) +
    facet_wrap(~ param, scales = "free",
               labeller = as_labeller(fac_labels, default = label_parsed),
               ncol = min(3, length(param_order_used))) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "#C1272D"),
                        guide = "none") +
    labs(x = "true value (raw)", y = "posterior mean (raw)",
         title = sprintf("%s — dataset %d (divs = %d)",
                         model_key, ds_idx, n_div)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "plain"))

  plot_path <- file.path(out_dir, sprintf("check_dataset_%d.pdf", ds_idx))
  ggsave(plot_path, p, width = 9, height = 6.5)
  cat(sprintf("\nSaved plot: %s\n", plot_path))
}
