#!/usr/bin/env Rscript
## ============================================================================
## plot_ppc_local.R — build the 6-panel PPC figure from ppc_summary.rds files
##
## Aim:     Assemble the manuscript 6-panel combined PPC figure from the small
##          per-model summary RDS files (produced by summarize_ppc.R), avoiding
##          a re-read of the large posterior_predictives.csv files. Layout:
##          observed value as a background bar per x-group; each model's
##          prediction as a dodged point + 95% CrI error bar; per-panel dynamic
##          y-axis; panel (d) shows the P(chose risky) EVD psychometric curve.
## Inputs:  results/<study>/<family>/<model>/ppc_summary.rds for the requested
##          CS and CCSS families/models.
## Outputs: the combined figure at --out (.pdf) plus a sibling .png.
## Usage:   Rscript 04_ppc/plot_ppc_local.R --study study1 \
##            --families mv_cs,mv_ccss \
##            --models_cs   baseline,sp,dr,sp_dr \
##            --models_ccss baseline,n,r,a,s,n_r_a_s \
##            --best_cs sp_dr --best_ccss n_r_a_s \
##            --out results/study1/ppc_combined_mv_local.pdf
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (posterior predictive
## checks). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork)
})

opt_list <- list(
  make_option("--study",       type = "character", default = NULL),
  make_option("--families",    type = "character", default = NULL),
  make_option("--models_cs",   type = "character", default = NULL),
  make_option("--models_ccss", type = "character", default = NULL),
  make_option("--best_cs",     type = "character", default = NULL),
  make_option("--best_ccss",   type = "character", default = NULL),
  make_option("--results_dir", type = "character", default = "results"),
  make_option("--out",         type = "character", default = "ppc_combined_local.pdf")
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$study), !is.null(opt$families))

families <- trimws(strsplit(opt$families, ",")[[1]])
ccss_families <- families[grepl("ccss", families)]
cs_families   <- setdiff(families[grepl("_cs", families)], ccss_families)

models_cs   <- if (!is.null(opt$models_cs))   trimws(strsplit(opt$models_cs,   ",")[[1]]) else NULL
models_ccss <- if (!is.null(opt$models_ccss)) trimws(strsplit(opt$models_ccss, ",")[[1]]) else NULL

study_dir <- file.path(opt$results_dir, opt$study)

# ---------------------------------------------------------------------------
# Display names
# ---------------------------------------------------------------------------
display_cs <- c(baseline = "Baseline", sp = "Starting Point Bias",
                dr = "Drift Rate Adjustment", sp_dr = "Full Model")
display_ccss_fixed <- c(baseline = "Baseline", n = "Signal-to-noise Ratio",
                        r = "Risk Preference", a = "Threshold",
                        s = "Skewness Preference")
get_display <- function(code, family_type, best_code = NULL) {
  if (family_type == "cs") {
    if (!is.null(best_code) && code == best_code) return("Full Model")
    return(unname(display_cs[code]))
  }
  if (family_type == "ccss") {
    if (!is.null(best_code) && code == best_code) return("Best")
    return(unname(display_ccss_fixed[code]))
  }
  code
}

# Colours: model palette (no Observed; bar is plain grey)
cs_pal   <- c("Baseline" = "#D55E00", "Starting Point Bias" = "#0072B2",
              "Drift Rate Adjustment" = "#7A68A6", "Full Model" = "#009E73")
ccss_pal <- c("Baseline" = "#D55E00", "Signal-to-noise Ratio" = "#0072B2",
              "Risk Preference" = "#009E73", "Threshold" = "#7A68A6",
              "Skewness Preference" = "#56B4E9", "Best" = "#E69F00")

# CC/SS trial-type colours: two greys, chosen so they don't collide with any
# of the saturated model colours (which include blue, gold/orange, green, etc.).
trial_type_pal <- c("CC" = "#6E6E6E",  # dark grey
                    "SS" = "#C8C8C8")  # light grey

# Order models in the display order (left → right within bar)
cs_levels   <- sapply(models_cs,   get_display,
                       family_type = "cs",   best_code = opt$best_cs,   USE.NAMES = FALSE)
ccss_levels <- sapply(models_ccss, get_display,
                       family_type = "ccss", best_code = opt$best_ccss, USE.NAMES = FALSE)

# ---------------------------------------------------------------------------
# Load summary RDS files
# ---------------------------------------------------------------------------
load_summaries <- function(family, model_list) {
  out <- list()
  for (mdl in model_list) {
    p <- file.path(study_dir, family, mdl, "ppc_summary.rds")
    if (!file.exists(p)) { message("Missing: ", p); next }
    s <- readRDS(p)
    out[[mdl]] <- s
  }
  out
}

cs_sum   <- if (length(cs_families)   > 0 && length(models_cs)   > 0)
              load_summaries(cs_families[1],   models_cs)   else NULL
ccss_sum <- if (length(ccss_families) > 0 && length(models_ccss) > 0)
              load_summaries(ccss_families[1], models_ccss) else NULL

# ---------------------------------------------------------------------------
# Dynamic y-axis helper
# ---------------------------------------------------------------------------
dyn_ylim <- function(values, pad_frac = 0.10, min_pad = 0,
                      hard_min = NULL, hard_max = NULL,
                      include_zero = FALSE) {
  r <- range(values, na.rm = TRUE)
  if (include_zero) r <- range(c(r, 0), na.rm = TRUE)
  span <- diff(r)
  if (span < 1e-6) span <- abs(r[1]) * 0.1 + 0.05
  # Padding is the larger of `pad_frac * span` and `min_pad` (absolute floor)
  pad <- max(span * pad_frac, min_pad)
  lo <- r[1] - pad
  hi <- r[2] + pad
  if (!is.null(hard_min)) lo <- max(lo, hard_min)
  if (!is.null(hard_max)) hi <- min(hi, hard_max)
  c(lo, hi)
}

# ---------------------------------------------------------------------------
# Build per-model predicted summaries from a list of summaries
# ---------------------------------------------------------------------------
# pred_values_fn returns a numeric vector of per-sample_id predicted values
build_pred_simple <- function(summaries, pred_values_fn,
                               family_type, best_code, levs) {
  lapply(names(summaries), function(mdl) {
    v <- pred_values_fn(summaries[[mdl]])
    tibble(
      model_display = factor(get_display(mdl, family_type, best_code), levels = levs),
      mean_val = mean(v, na.rm = TRUE),
      lo       = quantile(v, 0.025, na.rm = TRUE),
      hi       = quantile(v, 0.975, na.rm = TRUE)
    )
  }) %>% bind_rows()
}

build_pred_grouped <- function(summaries, pred_df_fn, group_cols,
                                family_type, best_code, levs) {
  lapply(names(summaries), function(mdl) {
    pred_df_fn(summaries[[mdl]]) %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(mean_val = mean(prop, na.rm = TRUE),
                lo = quantile(prop, 0.025, na.rm = TRUE),
                hi = quantile(prop, 0.975, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(model_display = factor(get_display(mdl, family_type, best_code), levels = levs))
  }) %>% bind_rows()
}

# ---------------------------------------------------------------------------
# Theme
# ---------------------------------------------------------------------------
theme_panel <- function() {
  theme_classic(base_size = 11) +
    theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          panel.grid.minor   = element_blank(),
          plot.title    = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9.5, colour = "grey30"),
          legend.position = "right",
          legend.key.size = unit(0.9, "lines"),
          axis.title.x = element_text(margin = margin(t = 6)))
}

# ---------------------------------------------------------------------------
# Single-category panel: 1 observed bar + N model dots dodged on top
# ---------------------------------------------------------------------------
panel_single <- function(obs_val, pred_df, palette, title, xlab, ylab,
                          x_label = "Group",
                          legend_name = "Model",
                          hline_dashed = NULL, hline_solid = NULL,
                          dynamic_y = TRUE,
                          y_hard_min = NULL, y_hard_max = NULL,
                          y_min_pad = 0) {
  obs_df  <- tibble(x_var = x_label, y_val = obs_val)
  pred_df <- pred_df %>% mutate(x_var = x_label)

  # Difference / single-category panel: use a warm off-white fill that is
  # clearly distinct from the CC (dark grey) and SS (light grey) bars in
  # panels c/d/f, so the observed bar in (a), (b), (e) reads as a separate
  # "single overall value" rather than as a third grey shade.
  diff_bar_fill <- "#E8E0D2"  # warm light beige
  p <- ggplot() +
    geom_col(data = obs_df, aes(x = x_var, y = y_val),
             fill = diff_bar_fill, colour = "grey45",
             width = 0.65, alpha = 0.95) +
    geom_errorbar(data = pred_df,
                  aes(x = x_var, ymin = lo, ymax = hi, colour = model_display),
                  position = position_dodge(width = 0.55),
                  width = 0.28, linewidth = 1.0) +
    geom_point(data = pred_df,
               aes(x = x_var, y = mean_val, colour = model_display,
                   shape = model_display),
               position = position_dodge(width = 0.55), size = 2.8,
               stroke = 1.0) +
    scale_colour_manual(values = palette, name = legend_name, drop = FALSE) +
    scale_shape_manual(values = setNames(c(16, 17, 15, 18, 8, 4)[seq_along(palette)],
                                          names(palette)),
                       name = legend_name, drop = FALSE) +
    labs(title = title, x = xlab, y = ylab) +
    theme_panel()

  if (!is.null(hline_dashed))
    p <- p + geom_hline(yintercept = hline_dashed, linetype = "dashed",
                        colour = "grey50", linewidth = 0.4)
  if (!is.null(hline_solid))
    p <- p + geom_hline(yintercept = hline_solid, colour = "black", linewidth = 0.4)

  if (dynamic_y) {
    all_vals <- c(obs_df$y_val, pred_df$lo, pred_df$hi)
    yl <- dyn_ylim(all_vals, pad_frac = 0.15, min_pad = y_min_pad,
                    hard_min = y_hard_min, hard_max = y_hard_max,
                    include_zero = !is.null(hline_solid))
    p <- p + coord_cartesian(ylim = yl)
  }
  p
}

# ---------------------------------------------------------------------------
# Multi-category panel: observed bars per category + model dots dodged within
# ---------------------------------------------------------------------------
panel_multi <- function(obs_df, pred_df, palette, title, xlab, ylab,
                         legend_name = "Model",
                         paired_bars = FALSE,  # bars dodged by condition within each x?
                         hline_dashed = NULL,
                         dynamic_y = TRUE,
                         y_hard_min = NULL, y_hard_max = NULL,
                         y_min_pad = 0) {
  # obs_df and pred_df MUST have a `condition` column (CC/SS).
  # paired_bars=TRUE: condition dodged within each x_var (panel d)
  # paired_bars=FALSE: x_var IS the condition; one bar per x (panels c, f)
  stopifnot("condition" %in% names(obs_df), "condition" %in% names(pred_df))

  dodge_width <- if (paired_bars) 0.8 else 0.7
  bar_width   <- if (paired_bars) 0.75 else 0.65

  p <- ggplot() +
    geom_col(data = obs_df,
             aes(x = x_var, y = y_val, fill = condition),
             position = position_dodge(width = dodge_width, preserve = "single"),
             width = bar_width, alpha = 0.85, colour = "grey50") +
    geom_errorbar(data = pred_df,
                  aes(x = x_var, ymin = lo, ymax = hi,
                      colour = model_display,
                      group = interaction(model_display, condition)),
                  position = position_dodge(width = dodge_width),
                  # thinner error bars in paired-bar panels (panel d) where
                  # 12 dots per EVD bin would otherwise visually collide
                  width = if (paired_bars) 0.18 else 0.22,
                  linewidth = if (paired_bars) 0.6 else 1.0) +
    geom_point(data = pred_df,
               aes(x = x_var, y = mean_val,
                   colour = model_display, shape = model_display,
                   group = interaction(model_display, condition)),
               position = position_dodge(width = dodge_width),
               size = if (paired_bars) 1.8 else 2.6,
               stroke = if (paired_bars) 0.7 else 1.0) +
    scale_fill_manual(values = trial_type_pal, name = "Trial Type") +
    scale_colour_manual(values = palette, name = legend_name, drop = FALSE) +
    scale_shape_manual(values = setNames(c(16, 17, 15, 18, 8, 4)[seq_along(palette)],
                                          names(palette)),
                       name = legend_name, drop = FALSE) +
    labs(title = title, x = xlab, y = ylab) +
    theme_panel()

  if (!is.null(hline_dashed))
    p <- p + geom_hline(yintercept = hline_dashed, linetype = "dashed",
                        colour = "grey50", linewidth = 0.4)

  if (dynamic_y) {
    all_vals <- c(obs_df$y_val, pred_df$lo, pred_df$hi)
    yl <- dyn_ylim(all_vals, pad_frac = 0.12, min_pad = y_min_pad,
                    hard_min = y_hard_min, hard_max = y_hard_max)
    p <- p + coord_cartesian(ylim = yl)
  }
  p
}

# =============================================================================
# Panel (a) — Complex Choice Proportion (CS)
# =============================================================================
p_a <- NULL
if (!is.null(cs_sum)) {
  obs_val <- cs_sum[[1]]$cs$obs_overall$obs_prop_complex
  pred_df <- build_pred_simple(cs_sum,
                                function(s) s$cs$pred_overall$prop_complex,
                                family_type = "cs", best_code = opt$best_cs,
                                levs = cs_levels)
  p_a <- panel_single(obs_val, pred_df, cs_pal,
                      title = "(a) Complex Choice Proportion",
                      xlab = NULL, ylab = "P(chose complex)",
                      x_label = "Complex vs. Simple",
                      legend_name = "Model (CS)",
                      hline_dashed = 0.5, y_hard_min = 0, y_hard_max = 1,
                      y_min_pad = 0.04)
}

# =============================================================================
# Panel (b) — RT Difference (Complex - Simple), CS
# =============================================================================
p_b <- NULL
if (!is.null(cs_sum)) {
  obs_val <- cs_sum[[1]]$cs$obs_rt_overall$obs_rt_diff
  pred_df <- build_pred_simple(cs_sum,
                                function(s) s$cs$pred_rt_overall$mean_rt_diff,
                                family_type = "cs", best_code = opt$best_cs,
                                levs = cs_levels)
  p_b <- panel_single(obs_val, pred_df, cs_pal,
                      title = "(b) RT Difference (Complex - Simple)",
                      xlab = NULL, ylab = "RT difference (s)",
                      x_label = "Complex vs. Simple",
                      legend_name = "Model (CS)",
                      hline_solid = 0)
}

# =============================================================================
# Panel (c) — Risky Choice Proportion (CCSS, CC vs SS)
# =============================================================================
p_c <- NULL
if (!is.null(ccss_sum)) {
  obs_df <- ccss_sum[[1]]$ccss$obs_by_cond %>%
    transmute(x_var = condition, condition = condition, y_val = obs_prop_risky)
  pred_df <- build_pred_grouped(ccss_sum,
                                 pred_df_fn = function(s) s$ccss$pred_by_cond %>%
                                   transmute(x_var = condition, condition = condition,
                                             prop = prop_risky),
                                 group_cols = c("x_var", "condition"),
                                 family_type = "ccss", best_code = opt$best_ccss,
                                 levs = ccss_levels)
  p_c <- panel_multi(obs_df, pred_df, ccss_pal,
                     title = "(c) Risky Choice Proportion",
                     xlab = "Trial Type", ylab = "P(chose risky)",
                     legend_name = "Model (CCSS)",
                     paired_bars = FALSE,
                     hline_dashed = 0.5, y_hard_min = 0, y_hard_max = 1,
                     y_min_pad = 0.04)
}

# =============================================================================
# Panel (d) — Risky Choice by EVD bin (CCSS), paired CC/SS bars
# =============================================================================
p_d <- NULL
if (!is.null(ccss_sum)) {
  obs_df <- ccss_sum[[1]]$ccss$obs_by_evd %>%
    transmute(x_var = evd_bin, condition, y_val = obs_prop_risky)
  pred_df <- build_pred_grouped(ccss_sum,
                                 pred_df_fn = function(s) s$ccss$pred_by_evd %>%
                                   transmute(x_var = evd_bin, condition, prop = prop_risky),
                                 group_cols = c("x_var", "condition"),
                                 family_type = "ccss", best_code = opt$best_ccss,
                                 levs = ccss_levels)
  p_d <- panel_multi(obs_df, pred_df, ccss_pal,
                     title = "(d) Risky Choice across EVD Levels",
                     xlab = "EVD (EV_risky - EV_safe)",
                     ylab = "P(chose risky)",
                     legend_name = "Model (CCSS)",
                     paired_bars = TRUE,
                     hline_dashed = 0.5,
                     y_hard_min = 0, y_hard_max = 1) +
         theme(axis.text.x = element_text(size = 8))
}

# =============================================================================
# Panel (e) — RT Difference (CC - SS), CCSS
# =============================================================================
p_e <- NULL
if (!is.null(ccss_sum)) {
  obs_val <- ccss_sum[[1]]$ccss$obs_rt_overall$obs_rt_diff
  pred_df <- build_pred_simple(ccss_sum,
                                function(s) s$ccss$pred_rt_overall$mean_rt_diff,
                                family_type = "ccss", best_code = opt$best_ccss,
                                levs = ccss_levels)
  p_e <- panel_single(obs_val, pred_df, ccss_pal,
                      title = "(e) RT Difference (CC - SS)",
                      xlab = NULL, ylab = "RT difference (s)",
                      x_label = "CC vs. SS",
                      legend_name = "Model (CCSS)",
                      hline_solid = 0)
}

# =============================================================================
# Panel (f) — Skewness Preference (CCSS, CC vs SS)
# =============================================================================
p_f <- NULL
if (!is.null(ccss_sum)) {
  obs_df <- ccss_sum[[1]]$ccss$obs_by_cond %>%
    transmute(x_var = condition, condition = condition, y_val = obs_prop_rskew)
  pred_df <- build_pred_grouped(ccss_sum,
                                 pred_df_fn = function(s) s$ccss$pred_by_cond %>%
                                   transmute(x_var = condition, condition = condition,
                                             prop = prop_rskew),
                                 group_cols = c("x_var", "condition"),
                                 family_type = "ccss", best_code = opt$best_ccss,
                                 levs = ccss_levels)
  p_f <- panel_multi(obs_df, pred_df, ccss_pal,
                     title = "(f) Skewness Preference",
                     xlab = "Trial Type", ylab = "P(chose right-skewed)",
                     legend_name = "Model (CCSS)",
                     paired_bars = FALSE,
                     hline_dashed = 0.5, y_hard_min = 0, y_hard_max = 1,
                     y_min_pad = 0.04)
}

# =============================================================================
# Assemble
# =============================================================================
panels <- list(p_a, p_b, p_c, p_d, p_e, p_f)
available <- !sapply(panels, is.null)
n_avail <- sum(available)

cat(sprintf("Panels assembled: %d/6\n", n_avail))

if (n_avail >= 2) {
  # No figure-level title — caption is added in LaTeX.
  # Force exactly 3 legends:
  #   Panel a contributes "Model (CS)"
  #   Panel c contributes "Model (CCSS)" + "Trial Type"
  #   All other panels drop their guides so patchwork doesn't duplicate.
  drop_guides <- function(p) {
    if (is.null(p)) return(NULL)
    p + guides(colour = "none", shape = "none", fill = "none")
  }
  panels[[2]] <- drop_guides(panels[[2]])   # b
  panels[[4]] <- drop_guides(panels[[4]])   # d
  panels[[5]] <- drop_guides(panels[[5]])   # e
  panels[[6]] <- drop_guides(panels[[6]])   # f

  combined <- wrap_plots(panels[available], ncol = 2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.spacing.y = unit(0.3, "cm"))
  ggsave(opt$out, combined, width = 13, height = 12.5)
  ggsave(sub("\\.pdf$", ".png", opt$out), combined,
         width = 13, height = 12.5, dpi = 200)
  cat("Wrote:\n  ", opt$out, "\n  ",
      sub("\\.pdf$", ".png", opt$out), "\n", sep = "")
} else {
  message("Not enough panels (need >= 2)")
}
