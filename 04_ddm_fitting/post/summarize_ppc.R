#!/usr/bin/env Rscript
## ============================================================================
## summarize_ppc.R — reduce one posterior_predictives.csv to a small RDS
##
## Aim:     For each model fit the CSV is ~1.25 GB (draws x trials); the
##          plotting code only needs per-(sample, grouping-var) aggregates.
##          This computes those aggregates once (group/EVD-bin choice
##          proportions, RT differences, per-participant obs-vs-pred summaries,
##          and the CCSS EVD x complexity GLM interaction) and stores them in a
##          small list-RDS consumed by plot_ppc_local.R and
##          plot_ppc_individual.R.
## Inputs:  results/<study>/<family>/<model>/posterior_predictives.csv.
## Outputs: results/<study>/<family>/<model>/ppc_summary.rds.
## Usage:   Rscript 04_ppc/summarize_ppc.R --study study1 --family mv_cs \
##            --model sp_dr
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (posterior predictive
## checks). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tidyr); library(readr)
})

opt_list <- list(
  make_option("--study",  type = "character", default = NULL),
  make_option("--family", type = "character", default = NULL),
  make_option("--model",  type = "character", default = NULL),
  make_option("--results_dir", type = "character", default = "results")
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$study), !is.null(opt$family), !is.null(opt$model))

ppc_dir  <- file.path(opt$results_dir, opt$study, opt$family, opt$model)
ppc_path <- file.path(ppc_dir, "posterior_predictives.csv")
out_path <- file.path(ppc_dir, "ppc_summary.rds")
stopifnot(file.exists(ppc_path))

cat(sprintf("Reading: %s\n", ppc_path))
d <- read_csv(ppc_path, show_col_types = FALSE)
cat(sprintf("  %d rows, %d sample_ids, %d participants\n",
            nrow(d), length(unique(d$sample_id)), length(unique(d$participant))))

is_cs   <- grepl("_cs", opt$family) & !grepl("ccss", opt$family)
is_ccss <- grepl("ccss", opt$family)

EVD_BREAKS <- c(-Inf, -15, -5, 5, 15, Inf)
EVD_LABELS <- c("-21 to -19", "-11 to -9", "-1 to 1", "9 to 11", "19 to 21")

summary_list <- list(study = opt$study, family = opt$family, model = opt$model)

# ============================================================================
# CS family
# ============================================================================
if (is_cs && "chose_complex_obs" %in% names(d)) {
  cat("Aggregating CS summaries...\n")

  # Predicted chose-complex on each trial: derived from observed direction
  d <- d %>% mutate(
    oa_complex = as.integer(sign(chose_complex_obs * observed_choice)),
    chose_complex_pred = as.integer(oa_complex * pred_choice)
  )

  # ---- Observed (constant across samples) ----
  obs_first <- d %>% filter(sample_id == min(sample_id))

  obs_overall <- tibble(
    obs_prop_complex = mean(obs_first$chose_complex_obs == 1, na.rm = TRUE)
  )

  obs_per_subj <- obs_first %>%
    group_by(participant) %>%
    summarise(
      obs_prop_complex = mean(chose_complex_obs == 1, na.rm = TRUE),
      obs_med_rt_complex = median(observed_rt[chose_complex_obs ==  1], na.rm = TRUE),
      obs_med_rt_simple  = median(observed_rt[chose_complex_obs == -1], na.rm = TRUE),
      .groups = "drop"
    ) %>% mutate(obs_rt_diff = obs_med_rt_complex - obs_med_rt_simple)

  obs_rt_overall <- tibble(
    obs_rt_diff = mean(obs_per_subj$obs_rt_diff, na.rm = TRUE)
  )

  # ---- Predicted, per sample_id (overall) ----
  pred_overall <- d %>%
    group_by(sample_id) %>%
    summarise(prop_complex = mean(chose_complex_pred == 1, na.rm = TRUE),
              .groups = "drop")

  # ---- Predicted RT diff, per (sample_id, participant) then averaged ----
  pred_rt_per_subj <- d %>%
    group_by(sample_id, participant, chose_complex_pred) %>%
    summarise(med_rt = median(pred_rt, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = chose_complex_pred, values_from = med_rt,
                names_prefix = "rt_") %>%
    mutate(rt_diff = `rt_1` - `rt_-1`)

  pred_rt_overall <- pred_rt_per_subj %>%
    group_by(sample_id) %>%
    summarise(mean_rt_diff = mean(rt_diff, na.rm = TRUE), .groups = "drop")

  # ---- Predicted per (participant) — averaged over sample_id for scatter ----
  pred_per_subj <- d %>%
    group_by(sample_id, participant) %>%
    summarise(prop_complex = mean(chose_complex_pred == 1, na.rm = TRUE),
              .groups = "drop") %>%
    group_by(participant) %>%
    summarise(pred_prop_mean = mean(prop_complex),
              pred_prop_lo   = quantile(prop_complex, 0.025),
              pred_prop_hi   = quantile(prop_complex, 0.975),
              .groups = "drop")

  pred_rt_per_subj_summary <- pred_rt_per_subj %>%
    group_by(participant) %>%
    summarise(pred_rt_diff_mean = mean(rt_diff, na.rm = TRUE),
              pred_rt_diff_lo   = quantile(rt_diff, 0.025, na.rm = TRUE),
              pred_rt_diff_hi   = quantile(rt_diff, 0.975, na.rm = TRUE),
              .groups = "drop")

  summary_list$cs <- list(
    obs_overall    = obs_overall,
    obs_rt_overall = obs_rt_overall,
    obs_per_subj   = obs_per_subj,
    pred_overall   = pred_overall,
    pred_rt_overall= pred_rt_overall,
    pred_per_subj  = pred_per_subj,
    pred_rt_per_subj_summary = pred_rt_per_subj_summary
  )
}

# ============================================================================
# CCSS family
# ============================================================================
if (is_ccss) {
  cat("Aggregating CCSS summaries...\n")

  if ("sdd" %in% names(d)) {
    d <- d %>% mutate(
      is_risky_obs  = as.integer((sdd > 0 & observed_choice == 1) |
                                  (sdd < 0 & observed_choice == -1)),
      is_risky_pred = as.integer((sdd > 0 & pred_choice == 1) |
                                  (sdd < 0 & pred_choice == -1)),
      condition = ifelse(con == 1, "CC", "SS")
    )
  }
  if ("skew" %in% names(d)) {
    d <- d %>% mutate(
      is_rskew_obs  = as.integer((skew > 0 & observed_choice == 1) |
                                  (skew < 0 & observed_choice == -1)),
      is_rskew_pred = as.integer((skew > 0 & pred_choice == 1) |
                                  (skew < 0 & pred_choice == -1))
    )
  }
  if ("evd" %in% names(d)) {
    d <- d %>% mutate(
      is_ev_consistent_obs  = as.integer((evd > 0 & observed_choice == 1) |
                                          (evd < 0 & observed_choice == -1)),
      is_ev_consistent_pred = as.integer((evd > 0 & pred_choice == 1) |
                                          (evd < 0 & pred_choice == -1)),
      evd_bin = cut(evd, breaks = EVD_BREAKS, labels = EVD_LABELS)
    )
  }

  # ---- Observed ----
  obs_first <- d %>% filter(sample_id == min(sample_id))

  obs_by_cond <- obs_first %>%
    group_by(condition) %>%
    summarise(obs_prop_risky = mean(is_risky_obs, na.rm = TRUE),
              obs_prop_rskew = mean(is_rskew_obs, na.rm = TRUE),
              .groups = "drop")

  obs_by_evd <- obs_first %>%
    filter(!is.na(evd_bin)) %>%
    group_by(evd_bin, condition) %>%
    summarise(obs_prop_ev_cons = mean(is_ev_consistent_obs, na.rm = TRUE),
              obs_prop_risky   = mean(is_risky_obs, na.rm = TRUE),
              .groups = "drop")

  obs_per_subj <- obs_first %>%
    group_by(participant, con) %>%
    summarise(med_rt = median(observed_rt, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = con, values_from = med_rt, names_prefix = "rt_con") %>%
    mutate(obs_rt_diff = `rt_con1` - `rt_con-1`)

  obs_rt_overall <- tibble(obs_rt_diff = mean(obs_per_subj$obs_rt_diff, na.rm = TRUE))

  # ---- Predicted, per sample_id × condition ----
  pred_by_cond <- d %>%
    group_by(sample_id, condition) %>%
    summarise(prop_risky = mean(is_risky_pred, na.rm = TRUE),
              prop_rskew = mean(is_rskew_pred, na.rm = TRUE),
              .groups = "drop")

  # ---- Predicted, per sample_id × evd_bin × condition ----
  pred_by_evd <- d %>%
    filter(!is.na(evd_bin)) %>%
    group_by(sample_id, evd_bin, condition) %>%
    summarise(prop_ev_cons = mean(is_ev_consistent_pred, na.rm = TRUE),
              prop_risky   = mean(is_risky_pred, na.rm = TRUE),
              .groups = "drop")

  # ---- Predicted RT diff (CC - SS), per sample_id ----
  pred_rt_per_subj <- d %>%
    group_by(sample_id, participant, con) %>%
    summarise(med_rt = median(pred_rt, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = con, values_from = med_rt, names_prefix = "rt_con") %>%
    mutate(rt_diff = `rt_con1` - `rt_con-1`)

  pred_rt_overall <- pred_rt_per_subj %>%
    group_by(sample_id) %>%
    summarise(mean_rt_diff = mean(rt_diff, na.rm = TRUE), .groups = "drop")

  # ---- Per-participant predicted (for scatter) ----
  pred_per_subj_choice <- d %>%
    group_by(sample_id, participant, condition) %>%
    summarise(prop_risky = mean(is_risky_pred, na.rm = TRUE),
              prop_rskew = mean(is_rskew_pred, na.rm = TRUE),
              .groups = "drop") %>%
    group_by(participant, condition) %>%
    summarise(pred_prop_risky_mean = mean(prop_risky),
              pred_prop_risky_lo   = quantile(prop_risky, 0.025),
              pred_prop_risky_hi   = quantile(prop_risky, 0.975),
              pred_prop_rskew_mean = mean(prop_rskew),
              pred_prop_rskew_lo   = quantile(prop_rskew, 0.025),
              pred_prop_rskew_hi   = quantile(prop_rskew, 0.975),
              .groups = "drop")

  pred_rt_per_subj_summary <- pred_rt_per_subj %>%
    group_by(participant) %>%
    summarise(pred_rt_diff_mean = mean(rt_diff, na.rm = TRUE),
              pred_rt_diff_lo   = quantile(rt_diff, 0.025, na.rm = TRUE),
              pred_rt_diff_hi   = quantile(rt_diff, 0.975, na.rm = TRUE),
              .groups = "drop")

  # ---- Per-participant OBSERVED choice proportions, by condition ----
  # Feeds the individual-level PPC scatterplots (risky, right-skewed) and the
  # observed EV-consistency (proxy for SNR: proportion choosing the higher-EV
  # option). EV-consistency is model-free; the CC-SS difference indexes how much
  # complexity degrades choice consistency.
  obs_choice_per_subj <- obs_first %>%
    group_by(participant, condition) %>%
    summarise(obs_prop_ev_cons = mean(is_ev_consistent_obs, na.rm = TRUE),
              obs_prop_risky   = mean(is_risky_obs, na.rm = TRUE),
              obs_prop_rskew   = mean(is_rskew_obs, na.rm = TRUE),
              .groups = "drop")

  # ---- Per-participant OBSERVED EV-consistency CC-SS difference ----
  obs_consistency_per_subj <- obs_choice_per_subj %>%
    select(participant, condition, obs_prop_ev_cons) %>%
    pivot_wider(names_from = condition, values_from = obs_prop_ev_cons) %>%
    mutate(obs_cons_diff = CC - SS) %>%
    select(participant, obs_cons_diff)

  # ---- Per-participant PREDICTED EV-consistency CC-SS difference (with CrI) ----
  # Per posterior sample, compute each participant's consistency in CC and SS,
  # take the CC-SS difference, then summarise across samples.
  pred_consistency_per_subj <- d %>%
    group_by(sample_id, participant, condition) %>%
    summarise(prop_ev_cons = mean(is_ev_consistent_pred, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = prop_ev_cons) %>%
    mutate(cons_diff = CC - SS) %>%
    group_by(participant) %>%
    summarise(pred_cons_diff_mean = mean(cons_diff, na.rm = TRUE),
              pred_cons_diff_lo   = quantile(cons_diff, 0.025, na.rm = TRUE),
              pred_cons_diff_hi   = quantile(cons_diff, 0.975, na.rm = TRUE),
              .groups = "drop")

  # ---- Per-participant PREDICTED CC-SS differences for risky / right-skew ----
  # Difference of the per-condition proportions, computed per posterior sample
  # (so the CrI reflects the joint CC/SS posterior), then summarised.
  pred_choice_diff_per_subj <- d %>%
    group_by(sample_id, participant, condition) %>%
    summarise(prop_risky = mean(is_risky_pred, na.rm = TRUE),
              prop_rskew = mean(is_rskew_pred, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = condition,
                values_from = c(prop_risky, prop_rskew)) %>%
    mutate(risky_diff = prop_risky_CC - prop_risky_SS,
           rskew_diff = prop_rskew_CC - prop_rskew_SS) %>%
    group_by(participant) %>%
    summarise(pred_risky_diff_mean = mean(risky_diff, na.rm = TRUE),
              pred_risky_diff_lo   = quantile(risky_diff, 0.025, na.rm = TRUE),
              pred_risky_diff_hi   = quantile(risky_diff, 0.975, na.rm = TRUE),
              pred_rskew_diff_mean = mean(rskew_diff, na.rm = TRUE),
              pred_rskew_diff_lo   = quantile(rskew_diff, 0.025, na.rm = TRUE),
              pred_rskew_diff_hi   = quantile(rskew_diff, 0.975, na.rm = TRUE),
              .groups = "drop")

  # ---- Per-participant EVD x complexity interaction (GLM SNR proxy) ----------
  # Logistic GLM  choice(option 1) ~ evd * con.  The evd:con coefficient indexes
  # how complexity modulates EV-sensitivity: a magnitude-weighted consistency
  # measure. Fit on observed choices, and on each of up to 200 posterior samples
  # (for the CrI). This is the measure used for the individual-level PPC.
  suppressPackageStartupMessages(library(data.table))
  DT <- as.data.table(d)
  # Bayesian logistic (weakly-informative Cauchy prior, arm::bayesglm) instead of
  # plain glm: bounds the interaction coefficient under separation, which is
  # common in per-participant / per-sample fits on ~90 trials and otherwise
  # produces |coef| ~ 100 outliers that make the estimate and plot unusable.
  .fit_int <- function(sub, ycol) {
    y <- as.integer(sub[[ycol]] == 1)
    tryCatch(unname(coef(suppressWarnings(
      arm::bayesglm(y ~ evd * con, data = sub, family = binomial())))["evd:con"]),
      error = function(e) NA_real_)
  }
  obs_consistency_glm <- DT[sample_id == min(sample_id),
                            .(obs_int = .fit_int(.SD, "observed_choice")),
                            by = participant]
  .sids <- sort(unique(DT$sample_id))
  if (length(.sids) > 200) .sids <- .sids[round(seq(1, length(.sids), length.out = 200))]
  pred_consistency_glm <- DT[sample_id %in% .sids,
                             .(int = .fit_int(.SD, "pred_choice")),
                             by = .(sample_id, participant)
                            ][, .(pred_int_mean = median(int, na.rm = TRUE),  # posterior median: robust to per-sample GLM separation (mean is not)
                                  pred_int_lo   = quantile(int, 0.025, na.rm = TRUE),
                                  pred_int_hi   = quantile(int, 0.975, na.rm = TRUE)),
                              by = participant]

  summary_list$ccss <- list(
    obs_by_cond    = obs_by_cond,
    obs_by_evd     = obs_by_evd,
    obs_per_subj   = obs_per_subj,
    obs_rt_overall = obs_rt_overall,
    obs_choice_per_subj      = obs_choice_per_subj,
    obs_consistency_per_subj = obs_consistency_per_subj,
    obs_consistency_glm      = as.data.frame(obs_consistency_glm),
    pred_by_cond   = pred_by_cond,
    pred_by_evd    = pred_by_evd,
    pred_rt_overall= pred_rt_overall,
    pred_per_subj_choice = pred_per_subj_choice,
    pred_consistency_per_subj = pred_consistency_per_subj,
    pred_consistency_glm      = as.data.frame(pred_consistency_glm),
    pred_choice_diff_per_subj = pred_choice_diff_per_subj,
    pred_rt_per_subj_summary = pred_rt_per_subj_summary
  )
}

# ============================================================================
saveRDS(summary_list, out_path)
cat(sprintf("\nWrote: %s (%s)\n", out_path,
            format(file.info(out_path)$size, big.mark = ",")))
