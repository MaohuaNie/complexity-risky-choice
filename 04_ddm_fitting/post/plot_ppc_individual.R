#!/usr/bin/env Rscript
## ============================================================================
## plot_ppc_individual.R — per-participant (individual-level) PPC scatters
##
## Aim:     For each metric, plot observed (x) vs posterior-predictive mean (y)
##          per participant with faint 95% CrI whiskers, a 45-degree identity
##          line and a per-panel obs-vs-pred correlation, one square facet per
##          model. Styled for APA 7 (Helvetica, thin
##          axes, no in-figure title; the caption is supplied in LaTeX).
##          Metrics whose fields are absent from ppc_summary.rds are skipped
##          (the CCSS consistency/risky/skew per-subject fields are produced by
##          summarize_ppc.R).
## Inputs:  results/<study>/<family>/<model>/ppc_summary.rds for the requested
##          CS and CCSS families/models.
## Outputs: results/<study>/individual_ppcs_<mv|cpt>/<metric>.{pdf,png}, with
##          metrics: CS complex-choice proportion + RT difference; CCSS RT
##          difference + EV-consistency (EVD x complexity GLM interaction).
## Usage:   Rscript 04_ppc/plot_ppc_individual.R --study study2 \
##            --families mv_cs,mv_ccss \
##            --models_cs baseline,sp,dr,sp_dr \
##            --models_ccss baseline,n,r,a,s,n_a_s \
##            --best_cs sp_dr --best_ccss n_a_s
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (posterior predictive
## checks). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tidyr); library(ggplot2)
})

opt_list <- list(
  make_option("--study",       type = "character", default = NULL),
  make_option("--families",    type = "character", default = "mv_cs,mv_ccss"),
  make_option("--models_cs",   type = "character", default = "baseline,sp,dr,sp_dr"),
  make_option("--models_ccss", type = "character", default = "baseline,n,r,a,s,n_a_s"),
  make_option("--best_cs",     type = "character", default = "sp_dr"),
  make_option("--best_ccss",   type = "character", default = "n_a_s"),
  make_option("--results_dir", type = "character", default = "results"),
  make_option("--out_dir",     type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$study))

study_dir <- file.path(opt$results_dir, opt$study)

families    <- trimws(strsplit(opt$families, ",")[[1]])
cs_family   <- families[grepl("cs", families) & !grepl("ccss", families)][1]
ccss_family <- families[grepl("ccss", families)][1]
models_cs   <- trimws(strsplit(opt$models_cs, ",")[[1]])
models_ccss <- trimws(strsplit(opt$models_ccss, ",")[[1]])

## Family class (mv | cpt) -> output folder individual_ppcs_<class>.
cls <- if (any(grepl("^cpt_", families))) "cpt" else "mv"
out_dir <- if (!is.null(opt$out_dir)) opt$out_dir else
             file.path(study_dir, paste0("individual_ppcs_", cls))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Output: ", out_dir)

## --- display names / order (match group figure) ----------------------------
display_cs <- c(baseline = "Baseline", sp = "Starting Point Bias",
                dr = "Drift Rate Adjustment", sp_dr = "Full Model")
## CCSS single-shift labels differ by family class: the MVS r/s shifts are on
## risk (beta) and skewness (eta); the CPT r/s shifts are on utility curvature
## (beta) and probability weighting (gamma).
display_ccss <- if (cls == "cpt") {
  c(baseline = "Baseline", n = "Signal-to-noise Ratio",
    r = "Utility Curvature", a = "Threshold", s = "Probability Weighting")
} else {
  c(baseline = "Baseline", n = "Signal-to-noise Ratio",
    r = "Risk Preference", a = "Threshold", s = "Skewness Preference")
}
disp <- function(code, type, best) {
  if (type == "cs")   { if (!is.null(best) && code == best) return("Full Model"); return(unname(display_cs[code])) }
  if (type == "ccss") { if (!is.null(best) && code == best) return("Best");       return(unname(display_ccss[code])) }
  code
}
cs_levels   <- vapply(models_cs,   disp, character(1), type = "cs",   best = opt$best_cs,   USE.NAMES = FALSE)
ccss_levels <- vapply(models_ccss, disp, character(1), type = "ccss", best = opt$best_ccss, USE.NAMES = FALSE)

## --- load summaries ---------------------------------------------------------
load_summaries <- function(family, model_list) {
  out <- list()
  for (m in model_list) {
    p <- file.path(study_dir, family, m, "ppc_summary.rds")
    if (!file.exists(p)) { message("  missing: ", p); next }
    out[[m]] <- readRDS(p)
  }
  out
}
cs_sum   <- if (!is.na(cs_family))   load_summaries(cs_family,   models_cs)   else list()
ccss_sum <- if (!is.na(ccss_family)) load_summaries(ccss_family, models_ccss) else list()

num <- function(x) as.numeric(unname(x))

## --- theme (APA 7) -----------------------------------
theme_ppc <- theme_classic(base_size = 9) +
  theme(
    text             = element_text(family = "Helvetica", colour = "black"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text       = element_text(size = 8, face = "bold", margin = margin(2.5, 0, 2.5, 0)),
    axis.line        = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.4, colour = "black"),
    axis.title       = element_text(size = 9.5),
    axis.text        = element_text(size = 7.5, colour = "black"),
    panel.spacing    = unit(6, "pt"),
    plot.margin      = margin(4, 6, 2, 2),
    aspect.ratio     = 1,
    legend.position  = "bottom",
    legend.title     = element_text(size = 8.5),
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(9, "pt"),
    legend.margin    = margin(0, 0, 0, 0)
  )

## Combine per-participant obs + pred across models into one long df.
combine_indiv <- function(summaries, get_obs, get_pred, type, best, levs, by = "participant") {
  rows <- lapply(names(summaries), function(m) {
    s <- summaries[[m]]
    ob <- tryCatch(get_obs(s),  error = function(e) NULL)
    pr <- tryCatch(get_pred(s), error = function(e) NULL)
    if (is.null(ob) || is.null(pr)) return(NULL)
    inner_join(ob, pr, by = by) %>%
      filter(!is.na(obs), !is.na(pred_mean)) %>%
      mutate(model_display = factor(disp(m, type, best), levels = levs))
  })
  bind_rows(rows)
}

## CC-SS difference per participant for a per-condition choice metric.
## Observed diff from obs_choice_per_subj; predicted diff prefers the extended
## per-sample field (pred_choice_diff_per_subj, with CrI) and otherwise falls
## back to difference-of-the-means (point only, no whiskers) so it still renders
## before a re-summarize.
combine_ccss_diff <- function(summaries, obs_col, prefix, levs, best) {
  rows <- lapply(names(summaries), function(m) {
    s <- summaries[[m]]$ccss
    if (is.null(s) || is.null(s$obs_choice_per_subj)) return(NULL)
    obs <- s$obs_choice_per_subj %>%
      select(participant, condition, v = all_of(obs_col)) %>%
      pivot_wider(names_from = condition, values_from = v) %>%
      transmute(participant, obs = CC - SS)
    ext <- s$pred_choice_diff_per_subj
    mcol <- paste0("pred_", prefix, "_diff_mean")
    if (!is.null(ext) && mcol %in% names(ext)) {
      pr <- ext %>% transmute(participant,
              pred_mean = .data[[mcol]],
              pred_lo   = num(.data[[paste0("pred_", prefix, "_diff_lo")]]),
              pred_hi   = num(.data[[paste0("pred_", prefix, "_diff_hi")]]))
    } else {
      pr <- s$pred_per_subj_choice %>%
        select(participant, condition, v = all_of(paste0("pred_prop_", prefix, "_mean"))) %>%
        pivot_wider(names_from = condition, values_from = v) %>%
        transmute(participant, pred_mean = CC - SS, pred_lo = CC - SS, pred_hi = CC - SS)
    }
    inner_join(obs, pr, by = "participant") %>%
      filter(!is.na(obs), !is.na(pred_mean)) %>%
      mutate(model_display = factor(disp(m, "ccss", best), levels = levs))
  })
  bind_rows(rows)
}

## Faceted obs-vs-pred scatter. No baked-in title (APA caption is separate).
indiv_plot <- function(df, xlab, ylab, ncol, colour_by = NULL) {
  lims <- range(c(df$obs, df$pred_mean, df$pred_lo, df$pred_hi), na.rm = TRUE)

  ## observed-vs-predicted correlation (rho) per model panel, top-left
  rho_df <- df %>%
    group_by(model_display) %>%
    summarise(rho = suppressWarnings(cor(obs, pred_mean, use = "complete.obs")),
              .groups = "drop") %>%
    mutate(label = ifelse(is.finite(rho), sprintf("rho == %.2f", rho), "rho == '--'"),
           x = lims[1] + 0.03 * diff(lims),
           y = lims[2] - 0.01 * diff(lims))

  base_aes <- if (is.null(colour_by)) aes(obs, pred_mean) else
                aes(obs, pred_mean, colour = .data[[colour_by]])
  p <- ggplot(df, base_aes) +
    geom_abline(slope = 1, intercept = 0, linetype = "22",
                colour = "grey55", linewidth = 0.3) +
    geom_errorbar(aes(ymin = pred_lo, ymax = pred_hi), width = 0,
                  linewidth = 0.25, alpha = 0.35) +
    geom_point(size = 0.7, alpha = 0.7) +
    geom_text(data = rho_df, aes(x = x, y = y, label = label),
              parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 1,
              size = 2.6, colour = "black") +
    facet_wrap(~ model_display, ncol = ncol,
               labeller = label_wrap_gen(width = 15)) +
    scale_x_continuous(limits = lims, expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = lims, expand = expansion(mult = 0.02)) +
    labs(x = xlab, y = ylab) +
    theme_ppc
  if (!is.null(colour_by))
    p <- p + scale_colour_manual(values = c("CC" = "#298c8c", "SS" = "#f1a226"),
                                 name = "Trial type")
  p
}

save_plot <- function(p, name, ncol, nrow, has_legend = FALSE) {
  panel <- 3.7                                   # cm, square panel
  w <- 1.8 + ncol * panel                        # y-axis title/labels + panels
  h <- 1.6 + nrow * (panel + 1.0) + (if (has_legend) 0.9 else 0)
  for (dev in c("pdf", "png")) {
    ggsave(file.path(out_dir, paste0(name, ".", dev)), p,
           width = w, height = h, units = "cm", dpi = 600,
           device = dev, limitsize = FALSE)
  }
  message("  wrote ", name, "  (", ncol, " x ", nrow, ")")
}

## Facet layout: both conditions use two rows (CS = 2x2, CCSS = 3x2).
cs_ncol   <- ceiling(length(cs_levels)   / 2); cs_nrow   <- 2
ccss_ncol <- ceiling(length(ccss_levels) / 2); ccss_nrow <- 2

## ===========================================================================
## CS metrics
## ===========================================================================
if (length(cs_sum) > 0) {
  df <- combine_indiv(cs_sum,
    function(s) s$cs$obs_per_subj  %>% transmute(participant, obs = obs_prop_complex),
    function(s) s$cs$pred_per_subj %>% transmute(participant,
        pred_mean = pred_prop_mean, pred_lo = num(pred_prop_lo), pred_hi = num(pred_prop_hi)),
    "cs", opt$best_cs, cs_levels)
  if (nrow(df)) save_plot(indiv_plot(df, "Observed P(complex)", "Predicted P(complex)", cs_ncol),
                          "cs_complex_choice", cs_ncol, cs_nrow)

  df <- combine_indiv(cs_sum,
    function(s) s$cs$obs_per_subj %>% transmute(participant, obs = obs_rt_diff),
    function(s) s$cs$pred_rt_per_subj_summary %>% transmute(participant,
        pred_mean = pred_rt_diff_mean, pred_lo = num(pred_rt_diff_lo), pred_hi = num(pred_rt_diff_hi)),
    "cs", opt$best_cs, cs_levels)
  if (nrow(df)) save_plot(indiv_plot(df, "Observed RT difference (s)", "Predicted RT difference (s)", cs_ncol),
                          "cs_rt_diff", cs_ncol, cs_nrow)
}

## ===========================================================================
## CCSS metrics
## ===========================================================================
if (length(ccss_sum) > 0) {
  ## RT difference (CC - SS)  [available now]
  df <- combine_indiv(ccss_sum,
    function(s) s$ccss$obs_per_subj %>% transmute(participant, obs = obs_rt_diff),
    function(s) s$ccss$pred_rt_per_subj_summary %>% transmute(participant,
        pred_mean = pred_rt_diff_mean, pred_lo = num(pred_rt_diff_lo), pred_hi = num(pred_rt_diff_hi)),
    "ccss", opt$best_ccss, ccss_levels)
  if (nrow(df)) save_plot(indiv_plot(df, "Observed RT difference (s)", "Predicted RT difference (s)", ccss_ncol),
                          "ccss_rt_diff", ccss_ncol, ccss_nrow)

  ## EV-consistency: EVD x complexity interaction (GLM SNR proxy)
  ## [needs extended summary fields obs_consistency_glm / pred_consistency_glm]
  df <- combine_indiv(ccss_sum,
    function(s) s$ccss$obs_consistency_glm %>% transmute(participant, obs = obs_int),
    function(s) s$ccss$pred_consistency_glm %>% transmute(participant,
        pred_mean = pred_int_mean, pred_lo = num(pred_int_lo), pred_hi = num(pred_int_hi)),
    "ccss", opt$best_ccss, ccss_levels)
  if (!is.null(df) && nrow(df))
    save_plot(indiv_plot(df, "Observed EVD x complexity interaction",
                         "Predicted EVD x complexity interaction", ccss_ncol),
              "ccss_consistency", ccss_ncol, ccss_nrow)
  else message("  [skip] ccss_consistency - GLM summary field absent")

  ## risky and right-skewed CC-SS differences are intentionally NOT plotted at
  ## the individual level: risk preference is weakly identified / near-null, and
  ## the main text reports the SNR (consistency) and threshold (RT) mechanisms.
  ## They remain available in the group-level figure.
}

message("Done: ", out_dir)
