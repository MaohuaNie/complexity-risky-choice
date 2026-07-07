#!/usr/bin/env Rscript
## ============================================================================
## plot_ppc.R — build PPC figures directly from posterior_predictives.csv
##
## Aim:     Read the full per-model PPC CSVs for one or more CS and/or CCSS
##          families and produce the APA-7 manuscript figures: the 6-panel
##          combined figure, RT-quantile plots, RT-distribution densities, and
##          individual-level obs-vs-pred scatters. (For large runs prefer the
##          summarize_ppc.R -> plot_ppc_local.R path, which avoids re-reading
##          the giant CSVs.)
## Inputs:  results/<study>/<family>/<model>/posterior_predictives.csv for each
##          requested family/model.
## Outputs: results/<study>/posterior_predictive_combined_<fam_tag>.pdf/.png and
##          per-family results/<study>/<family>/ppc_rt_quantiles.*,
##          ppc_rt_distribution_<model>.*, ppc_individual_choice.*,
##          ppc_individual_rt.*, ppc_individual_complex.* (CS),
##          ppc_individual_rt_diff.*, ppc_individual_consistency.* (CCSS),
##          ppc_complex_by_evd.* (CS).
## Usage:   Rscript 04_ppc/plot_ppc.R --study study2 \
##            --families mv_cs,mv_ccss \
##            --models_cs baseline,sp,dr,sp_dr \
##            --models_ccss baseline,n,r,a,s,n_r_a_s \
##            --best_cs sp_dr --best_ccss n_r_a_s
##
##          # single family, backward compatible:
##          Rscript 04_ppc/plot_ppc.R --study study2 --family mv_ccss \
##            --models baseline,n,r,a,s,n_r_a_s
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (posterior predictive
## checks). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

opt_list <- list(
  make_option("--study",       type = "character", default = NULL),
  ## New multi-family interface
  make_option("--families",    type = "character", default = NULL,
              help = "Comma-separated family names (e.g., mv_cs,mv_ccss)"),
  make_option("--models_cs",   type = "character", default = NULL,
              help = "Comma-separated CS model names"),
  make_option("--models_ccss", type = "character", default = NULL,
              help = "Comma-separated CCSS model names"),
  make_option("--best_cs",     type = "character", default = NULL,
              help = "LOO-CV best CS model name"),
  make_option("--best_ccss",   type = "character", default = NULL,
              help = "LOO-CV best CCSS model name"),
  ## Backward-compatible single-family interface
  make_option("--family",      type = "character", default = NULL),
  make_option("--models",      type = "character", default = NULL),
  make_option("--results_dir", type = "character", default = "results")
)

.override <- if (exists("PLOT_PPC_ARGS", envir = .GlobalEnv))
  get("PLOT_PPC_ARGS", envir = .GlobalEnv) else NULL
opt <- parse_args(OptionParser(option_list = opt_list),
                  args = if (is.null(.override)) commandArgs(trailingOnly = TRUE) else .override)

## ---------------------------------------------------------------------------
## Resolve families and models
## ---------------------------------------------------------------------------
if (!is.null(opt$families)) {
  families <- trimws(strsplit(opt$families, ",")[[1]])
} else if (!is.null(opt$family)) {
  families <- opt$family
} else {
  stop("Must specify --families or --family")
}

## Separate CS vs CCSS families
## CS = any family that contains "_cs" but isn't the CCSS variant
## (catches cpt_cs, mv_cs, and cpt_cs_7o). CCSS is the disambiguator.
ccss_families <- families[grepl("ccss", families)]
cs_families   <- setdiff(families[grepl("_cs", families)], ccss_families)

if (!is.null(opt$models_cs)) {
  models_cs <- trimws(strsplit(opt$models_cs, ",")[[1]])
} else if (!is.null(opt$models) && length(cs_families) > 0) {
  models_cs <- trimws(strsplit(opt$models, ",")[[1]])
} else {
  models_cs <- character(0)
}

if (!is.null(opt$models_ccss)) {
  models_ccss <- trimws(strsplit(opt$models_ccss, ",")[[1]])
} else if (!is.null(opt$models) && length(ccss_families) > 0) {
  models_ccss <- trimws(strsplit(opt$models, ",")[[1]])
} else {
  models_ccss <- character(0)
}

study_dir <- file.path(opt$results_dir, opt$study)

## ---------------------------------------------------------------------------
## APA 7 theme
## ---------------------------------------------------------------------------
theme_apa <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text              = element_text(colour = "black"),
      strip.background  = element_blank(),
      strip.text        = element_text(size = 12, face = "bold", hjust = 0),
      panel.grid        = element_blank(),
      axis.line         = element_line(colour = "black", linewidth = 0.4),
      axis.ticks        = element_line(colour = "black", linewidth = 0.4),
      axis.text         = element_text(colour = "black", size = 9),
      axis.title        = element_text(size = 11),
      axis.title.x      = element_text(margin = margin(t = 8)),
      axis.title.y      = element_text(margin = margin(r = 8)),
      legend.position   = "top",
      legend.title      = element_text(size = 11, face = "bold"),
      legend.text       = element_text(size = 10),
      plot.margin       = margin(8, 12, 8, 8)
    )
}

## ---------------------------------------------------------------------------
## Model display names and colours
## ---------------------------------------------------------------------------
CS_DISPLAY <- c(
  "baseline" = "Baseline",
  "sp"       = "Starting Point Bias",
  "dr"       = "Drift Rate Adjustment",
  "sp_dr"    = "Full Model",
  "skew"     = "Skewness"
)

CCSS_DISPLAY <- c(
  "baseline" = "Baseline",
  "n"        = "Signal-to-noise Ratio",
  "r"        = "Risk Preference",
  "a"        = "Threshold",
  "s"        = "Skewness Preference"
)

get_display_name <- function(model_code, family_type, best_code = NULL) {
  if (!is.null(best_code) && model_code == best_code) return("Best")
  tbl <- if (family_type == "cs") CS_DISPLAY else CCSS_DISPLAY
  nm <- tbl[model_code]
  if (is.na(nm)) return(model_code)
  nm
}

CS_COLOURS <- c(
  "Behavioral"              = "black",
  "Baseline"                = "#E41A1C",
  "Starting Point Bias"     = "#377EB8",
  "Drift Rate Adjustment"   = "#4DAF4A",
  "Full Model"              = "#984EA3",
  "Skewness"                = "#00CED1"
)

CCSS_COLOURS <- c(
  "Behavioral"              = "black",
  "Baseline"                = "#E41A1C",
  "Signal-to-noise Ratio"   = "#377EB8",
  "Risk Preference"         = "#4DAF4A",
  "Threshold"               = "#984EA3",
  "Skewness Preference"     = "#00CED1",
  "Best"                    = "#FF7F00"
)

## ---------------------------------------------------------------------------
## Helper: load posterior predictive CSVs for a family
## ---------------------------------------------------------------------------
load_ppc <- function(family, model_list) {
  fam_dir <- file.path(study_dir, family)
  all_ppc <- list()
  for (mdl in model_list) {
    ppc_path <- file.path(fam_dir, mdl, "posterior_predictives.csv")
    if (!file.exists(ppc_path)) {
      cat("  WARNING: not found for", family, "/", mdl, "- skipping\n")
      next
    }
    d <- read_csv(ppc_path, show_col_types = FALSE) %>% mutate(model = mdl)
    all_ppc[[mdl]] <- d
    cat(sprintf("  %s/%s: %d rows\n", family, mdl, nrow(d)))
  }
  if (length(all_ppc) == 0) return(NULL)
  bind_rows(all_ppc) %>% mutate(model = factor(model, levels = model_list))
}

## ---------------------------------------------------------------------------
## Helper: summarise a binary metric (group-level)
## ---------------------------------------------------------------------------
summarise_binary <- function(ppc, obs_var, pred_var, group_var = NULL) {
  obs_data <- ppc %>%
    filter(model == levels(model)[1]) %>%
    distinct(participant, trial_in_participant, .keep_all = TRUE)

  if (!is.null(group_var)) {
    obs_summary <- obs_data %>%
      group_by(!!sym(group_var)) %>%
      summarise(obs_prop = mean(!!sym(obs_var), na.rm = TRUE), .groups = "drop")

    pred_summary <- ppc %>%
      group_by(model, sample_id, !!sym(group_var)) %>%
      summarise(prop = mean(!!sym(pred_var), na.rm = TRUE), .groups = "drop") %>%
      group_by(model, !!sym(group_var)) %>%
      summarise(mean_prop = mean(prop),
                lower = quantile(prop, 0.025),
                upper = quantile(prop, 0.975), .groups = "drop")
  } else {
    obs_summary <- obs_data %>%
      summarise(obs_prop = mean(!!sym(obs_var), na.rm = TRUE))

    pred_summary <- ppc %>%
      group_by(model, sample_id) %>%
      summarise(prop = mean(!!sym(pred_var), na.rm = TRUE), .groups = "drop") %>%
      group_by(model) %>%
      summarise(mean_prop = mean(prop),
                lower = quantile(prop, 0.025),
                upper = quantile(prop, 0.975), .groups = "drop")
  }
  list(obs = obs_summary, pred = pred_summary)
}

## ---------------------------------------------------------------------------
## Fixed EVD bins
## ---------------------------------------------------------------------------
EVD_BREAKS <- c(-Inf, -15, -5, 5, 15, Inf)
EVD_LABELS <- c("-21 to -19", "-11 to -9", "-1 to 1", "9 to 11", "19 to 21")

## =========================================================================
## Load data
## =========================================================================
cat("Loading posterior predictives...\n")

ppc_cs   <- NULL
ppc_ccss <- NULL

if (length(cs_families) > 0 && length(models_cs) > 0) {
  ppc_cs <- load_ppc(cs_families[1], models_cs)
}
if (length(ccss_families) > 0 && length(models_ccss) > 0) {
  ppc_ccss <- load_ppc(ccss_families[1], models_ccss)
}

if (is.null(ppc_cs) && is.null(ppc_ccss))
  stop("No posterior predictive data found.")

## =========================================================================
## Compute derived columns
## =========================================================================

if (!is.null(ppc_cs)) {
  ## For CS: chose_complex for predicted trials
  ## chose_complex_obs is in the data; for predicted, map pred_choice through same logic
  ## In CS preprocessing: cho=1 means chose option A (key "f")
  ## chose_complex = 1 if (oa_complex==1 & cho==1) | (oa_complex==-1 & cho==-1)
  ## Since we stored chose_complex_obs directly, for predictions we use the same

  ## mapping: the DDM upper boundary is always option A, so pred_choice follows
  ## the same cho encoding. We need to know which side was complex per trial.
  ## chose_complex_obs == 1  when participant chose complex
  ## cho (observed_choice) == 1 when participant chose option A
  ## So oa_complex = sign(chose_complex_obs * observed_choice) tells us if A=complex
  ## Then pred_chose_complex = oa_complex * pred_choice (but capped to {-1, 1})
  if ("chose_complex_obs" %in% names(ppc_cs)) {
    ppc_cs <- ppc_cs %>%
      mutate(
        oa_complex = as.integer(sign(chose_complex_obs * observed_choice)),
        chose_complex_pred = as.integer(oa_complex * pred_choice)
      )
  }

  ## Display names
  best_cs_code <- opt$best_cs
  ppc_cs <- ppc_cs %>%
    mutate(model_display = sapply(as.character(model),
                                  get_display_name, family_type = "cs",
                                  best_code = best_cs_code))
  cs_model_order <- unique(ppc_cs$model_display)
  ppc_cs$model_display <- factor(ppc_cs$model_display, levels = cs_model_order)
}

if (!is.null(ppc_ccss)) {
  ## Derived binary indicators for CCSS
  if ("sdd" %in% names(ppc_ccss)) {
    ppc_ccss <- ppc_ccss %>%
      mutate(
        ## Risky choice: higher-variance option
        is_risky_obs  = as.integer((sdd > 0 & observed_choice == 1) |
                                     (sdd < 0 & observed_choice == -1)),
        is_risky_pred = as.integer((sdd > 0 & pred_choice == 1) |
                                     (sdd < 0 & pred_choice == -1)),
        ## Condition label
        condition = ifelse(con == 1, "CC", "SS")
      )
  }
  if ("skew" %in% names(ppc_ccss)) {
    ppc_ccss <- ppc_ccss %>%
      mutate(
        ## Right-skewed choice
        is_rskew_obs  = as.integer((skew > 0 & observed_choice == 1) |
                                     (skew < 0 & observed_choice == -1)),
        is_rskew_pred = as.integer((skew > 0 & pred_choice == 1) |
                                     (skew < 0 & pred_choice == -1))
      )
  }
  if ("evd" %in% names(ppc_ccss)) {
    ppc_ccss <- ppc_ccss %>%
      mutate(
        ## EV-consistent choice
        is_ev_consistent_obs  = as.integer((evd > 0 & observed_choice == 1) |
                                            (evd < 0 & observed_choice == -1)),
        is_ev_consistent_pred = as.integer((evd > 0 & pred_choice == 1) |
                                            (evd < 0 & pred_choice == -1)),
        evd_bin = cut(evd, breaks = EVD_BREAKS, labels = EVD_LABELS)
      )
  }

  ## Display names
  best_ccss_code <- opt$best_ccss
  ppc_ccss <- ppc_ccss %>%
    mutate(model_display = sapply(as.character(model),
                                   get_display_name, family_type = "ccss",
                                   best_code = best_ccss_code))
  ccss_model_order <- unique(ppc_ccss$model_display)
  ppc_ccss$model_display <- factor(ppc_ccss$model_display, levels = ccss_model_order)
}

## =========================================================================
## Panel (a): Complex Choice Proportion (CS)
## =========================================================================
p_a <- NULL
if (!is.null(ppc_cs) && "chose_complex_obs" %in% names(ppc_cs)) {
  cat("Panel (a): Complex choice proportion...\n")

  obs_complex <- ppc_cs %>%
    filter(model == levels(model)[1]) %>%
    distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
    summarise(obs_prop = mean(chose_complex_obs == 1, na.rm = TRUE)) %>%
    mutate(condition = "Complex vs. Simple")

  pred_complex <- ppc_cs %>%
    group_by(model_display, sample_id) %>%
    summarise(prop = mean(chose_complex_pred == 1, na.rm = TRUE), .groups = "drop") %>%
    group_by(model_display) %>%
    summarise(mean_prop = mean(prop),
              lower = quantile(prop, 0.025),
              upper = quantile(prop, 0.975), .groups = "drop") %>%
    mutate(condition = "Complex vs. Simple")

  p_a <- ggplot() +
    geom_col(data = obs_complex, aes(x = condition, y = obs_prop),
             fill = "#B2D8D8", width = 0.6, alpha = 0.7) +
    geom_point(data = pred_complex,
               aes(x = condition, y = mean_prop, colour = model_display),
               position = position_dodge(width = 0.5), size = 2.5) +
    geom_errorbar(data = pred_complex,
                  aes(x = condition, ymin = lower, ymax = upper,
                      colour = model_display),
                  position = position_dodge(width = 0.5), width = 0.3, linewidth = 0.6) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = CS_COLOURS, name = "Model Type") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "Complex Choice Proportion",
         x = NULL, y = "Complex Choice Proportion") +
    theme_apa()
}

## =========================================================================
## Panel (b): RT Difference (CS)
## =========================================================================
p_b <- NULL
if (!is.null(ppc_cs) && "chose_complex_obs" %in% names(ppc_cs)) {
  cat("Panel (b): RT difference (CS)...\n")

  obs_rt_cs <- ppc_cs %>%
    filter(model == levels(model)[1]) %>%
    distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
    group_by(participant, chose_complex_obs) %>%
    summarise(med_rt = median(observed_rt, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = chose_complex_obs, values_from = med_rt,
                names_prefix = "rt_") %>%
    mutate(rt_diff = `rt_1` - `rt_-1`) %>%
    summarise(obs_rt_diff = mean(rt_diff, na.rm = TRUE)) %>%
    mutate(condition = "Complex vs. Simple")

  pred_rt_cs <- ppc_cs %>%
    group_by(model_display, sample_id, participant, chose_complex_pred) %>%
    summarise(med_rt = median(pred_rt, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = chose_complex_pred, values_from = med_rt,
                names_prefix = "rt_") %>%
    mutate(rt_diff = `rt_1` - `rt_-1`) %>%
    group_by(model_display, sample_id) %>%
    summarise(mean_rt_diff = mean(rt_diff, na.rm = TRUE), .groups = "drop") %>%
    group_by(model_display) %>%
    summarise(mean_diff = mean(mean_rt_diff),
              lower = quantile(mean_rt_diff, 0.025),
              upper = quantile(mean_rt_diff, 0.975), .groups = "drop") %>%
    mutate(condition = "Complex vs. Simple")

  p_b <- ggplot() +
    geom_col(data = obs_rt_cs, aes(x = condition, y = obs_rt_diff),
             fill = "#B2D8D8", width = 0.6, alpha = 0.7) +
    geom_point(data = pred_rt_cs,
               aes(x = condition, y = mean_diff, colour = model_display),
               position = position_dodge(width = 0.5), size = 2.5) +
    geom_errorbar(data = pred_rt_cs,
                  aes(x = condition, ymin = lower, ymax = upper,
                      colour = model_display),
                  position = position_dodge(width = 0.5), width = 0.3, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "solid", colour = "black", linewidth = 0.3) +
    scale_colour_manual(values = CS_COLOURS, name = "Model Type") +
    labs(title = "RT Difference (Complex - Simple)",
         x = NULL, y = "RT Difference (s)") +
    theme_apa()
}

## =========================================================================
## Panel (c): Risky Choice Proportion (CCSS)
## =========================================================================
p_c <- NULL
if (!is.null(ppc_ccss) && "is_risky_obs" %in% names(ppc_ccss)) {
  cat("Panel (c): Risky choice proportion...\n")

  obs_risky <- ppc_ccss %>%
    filter(model == levels(model)[1]) %>%
    distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
    group_by(condition) %>%
    summarise(obs_prop = mean(is_risky_obs, na.rm = TRUE), .groups = "drop")

  pred_risky <- ppc_ccss %>%
    group_by(model_display, sample_id, condition) %>%
    summarise(prop = mean(is_risky_pred, na.rm = TRUE), .groups = "drop") %>%
    group_by(model_display, condition) %>%
    summarise(mean_prop = mean(prop),
              lower = quantile(prop, 0.025),
              upper = quantile(prop, 0.975), .groups = "drop")

  p_c <- ggplot() +
    geom_col(data = obs_risky, aes(x = condition, y = obs_prop, fill = condition),
             width = 0.6, alpha = 0.9) +
    geom_point(data = pred_risky,
               aes(x = condition, y = mean_prop, colour = model_display),
               position = position_dodge(width = 0.5), size = 2.5) +
    geom_errorbar(data = pred_risky,
                  aes(x = condition, ymin = lower, ymax = upper,
                      colour = model_display),
                  position = position_dodge(width = 0.5), width = 0.3, linewidth = 0.6) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    scale_fill_manual(values = c("CC" = "#5F5F5F", "SS" = "#BFBFBF"),
                      name = "Trial Type") +
    scale_colour_manual(values = CCSS_COLOURS, name = "Model Type") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "Risky Choice Proportion",
         x = "Trial Type", y = "Risky Choice Proportion") +
    theme_apa()
}

## =========================================================================
## Panel (d): EV Consistency (CCSS)
## =========================================================================
p_d <- NULL
if (!is.null(ppc_ccss) && "evd_bin" %in% names(ppc_ccss)) {
  cat("Panel (d): EV consistency...\n")

  ## Only include trials with non-NA evd_bin
  ppc_evd <- ppc_ccss %>% filter(!is.na(evd_bin))

  ## Paired CC/SS bars per EVD bin (matches manuscript reference figure).
  obs_ev <- ppc_evd %>%
    filter(model == levels(model)[1]) %>%
    distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
    group_by(evd_bin, condition) %>%
    summarise(obs_prop = mean(is_ev_consistent_obs, na.rm = TRUE),
              .groups = "drop")

  pred_ev <- ppc_evd %>%
    group_by(model_display, sample_id, evd_bin, condition) %>%
    summarise(prop = mean(is_ev_consistent_pred, na.rm = TRUE), .groups = "drop") %>%
    group_by(model_display, evd_bin, condition) %>%
    summarise(mean_prop = mean(prop),
              lower = quantile(prop, 0.025),
              upper = quantile(prop, 0.975), .groups = "drop")

  ## Position dodging: bars separate CC/SS per bin; dots for each
  ## model within each (bin, condition) cell.
  ## Bar dodge is controlled purely by the `condition` aesthetic.
  ## Dots nest further: outer dodge by condition, inner dodge by model.
  bar_dodge <- position_dodge(width = 0.8)
  dot_dodge <- position_jitterdodge(dodge.width = 0.8,
                                    jitter.width = 0.15, jitter.height = 0)

  p_d <- ggplot() +
    geom_col(data = obs_ev,
             aes(x = evd_bin, y = obs_prop, fill = condition),
             position = bar_dodge, width = 0.75, alpha = 0.7) +
    geom_point(data = pred_ev,
               aes(x = evd_bin, y = mean_prop,
                   colour = model_display, group = condition),
               position = bar_dodge, size = 2) +
    geom_errorbar(data = pred_ev,
                  aes(x = evd_bin, ymin = lower, ymax = upper,
                      colour = model_display, group = condition),
                  position = bar_dodge, width = 0.2, linewidth = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    scale_fill_manual(values = c("CC" = "#5F5F5F", "SS" = "#BFBFBF"),
                      name = "Trial Type") +
    scale_colour_manual(values = CCSS_COLOURS, name = "Model Type") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "EV Consistency across EVD Levels",
         x = "EVD (EV_risky - EV_safe)",
         y = "Higher-EV Choice Proportion") +
    theme_apa() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8))
}

## =========================================================================
## Panel (e): RT Difference CC vs SS (CCSS)
## =========================================================================
p_e <- NULL
if (!is.null(ppc_ccss)) {
  cat("Panel (e): RT difference (CCSS)...\n")

  obs_rt_ccss <- ppc_ccss %>%
    filter(model == levels(model)[1]) %>%
    distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
    group_by(participant, con) %>%
    summarise(med_rt = median(observed_rt, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = con, values_from = med_rt, names_prefix = "rt_con") %>%
    mutate(rt_diff = `rt_con1` - `rt_con-1`) %>%
    summarise(obs_rt_diff = mean(rt_diff, na.rm = TRUE)) %>%
    mutate(condition = "CC vs. SS")

  pred_rt_ccss <- ppc_ccss %>%
    group_by(model_display, sample_id, participant, con) %>%
    summarise(med_rt = median(pred_rt, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = con, values_from = med_rt, names_prefix = "rt_con") %>%
    mutate(rt_diff = `rt_con1` - `rt_con-1`) %>%
    group_by(model_display, sample_id) %>%
    summarise(mean_rt_diff = mean(rt_diff, na.rm = TRUE), .groups = "drop") %>%
    group_by(model_display) %>%
    summarise(mean_diff = mean(mean_rt_diff),
              lower = quantile(mean_rt_diff, 0.025),
              upper = quantile(mean_rt_diff, 0.975), .groups = "drop") %>%
    mutate(condition = "CC vs. SS")

  p_e <- ggplot() +
    geom_col(data = obs_rt_ccss, aes(x = condition, y = obs_rt_diff),
             fill = "#B2D8D8", width = 0.6, alpha = 0.7) +
    geom_point(data = pred_rt_ccss,
               aes(x = condition, y = mean_diff, colour = model_display),
               position = position_dodge(width = 0.5), size = 2.5) +
    geom_errorbar(data = pred_rt_ccss,
                  aes(x = condition, ymin = lower, ymax = upper,
                      colour = model_display),
                  position = position_dodge(width = 0.5), width = 0.3, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.3) +
    scale_colour_manual(values = CCSS_COLOURS, name = "Model Type") +
    labs(title = "RT Difference (CC - SS)",
         x = NULL, y = "RT Difference (s)") +
    theme_apa()
}

## =========================================================================
## Panel (f): Skewness Preference (CCSS)
## =========================================================================
p_f <- NULL
if (!is.null(ppc_ccss) && "is_rskew_obs" %in% names(ppc_ccss)) {
  cat("Panel (f): Skewness preference...\n")

  obs_skew <- ppc_ccss %>%
    filter(model == levels(model)[1]) %>%
    distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
    group_by(condition) %>%
    summarise(obs_prop = mean(is_rskew_obs, na.rm = TRUE), .groups = "drop")

  pred_skew <- ppc_ccss %>%
    group_by(model_display, sample_id, condition) %>%
    summarise(prop = mean(is_rskew_pred, na.rm = TRUE), .groups = "drop") %>%
    group_by(model_display, condition) %>%
    summarise(mean_prop = mean(prop),
              lower = quantile(prop, 0.025),
              upper = quantile(prop, 0.975), .groups = "drop")

  p_f <- ggplot() +
    geom_col(data = obs_skew, aes(x = condition, y = obs_prop, fill = condition),
             width = 0.6, alpha = 0.9) +
    geom_point(data = pred_skew,
               aes(x = condition, y = mean_prop, colour = model_display),
               position = position_dodge(width = 0.5), size = 2.5) +
    geom_errorbar(data = pred_skew,
                  aes(x = condition, ymin = lower, ymax = upper,
                      colour = model_display),
                  position = position_dodge(width = 0.5), width = 0.3, linewidth = 0.6) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    scale_fill_manual(values = c("CC" = "#5F5F5F", "SS" = "#BFBFBF"),
                      name = "Trial Type") +
    scale_colour_manual(values = CCSS_COLOURS, name = "Model Type") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "Skewness Preference",
         x = "Trial Type", y = "Right-Skewed Choice Proportion") +
    theme_apa()
}

## =========================================================================
## Assemble 6-panel combined figure
## =========================================================================
panels <- list(p_a, p_b, p_c, p_d, p_e, p_f)
available <- !sapply(panels, is.null)

if (sum(available) >= 2) {
  cat("Assembling combined figure...\n")

  ## Build the combined figure
  ## Top row: CS legend; Bottom 4 panels: CCSS legend
  ## Use patchwork with panel labels

## Build the final combined figure with legends
  ## CS panels (a, b) share a legend on top row
  ## CCSS panels (c-f) share a legend on middle row
  has_cs   <- !is.null(p_a)
  has_ccss <- !is.null(p_c)

  ## Legend strategy:
  ##   * Model Type colour legend appears ONCE on panel (a) (CS) and ONCE on
  ##     panel (c) (CCSS) — CS and CCSS can have different model sets.
  ##   * Trial Type fill legend (CC/SS) appears ONCE on panel (c); panels (d)
  ##     and (f) share that same legend via patchwork's `guides = "collect"`.
  ##   * Panels (b), (d), (e), (f) suppress their own legends.
  strip_legend <- function(p) p + theme(legend.position = "none")

  if (has_cs && has_ccss) {
    ## Row 1: CS (a, b) with CS Model Type legend on (a)
    top_row <- (p_a + theme(legend.position = "top",
                             legend.box = "horizontal") +
                  guides(colour = guide_legend(nrow = 1, title = "Model Type"))) |
               strip_legend(p_b)

    ## Row 2: CCSS (c, d) with CCSS Model Type + Trial Type legends on (c)
    mid_row <- (p_c + theme(legend.position = "top",
                             legend.box = "horizontal") +
                  guides(colour = guide_legend(nrow = 2, title = "Model Type",
                                               order = 1),
                         fill   = guide_legend(nrow = 1, title = "Trial Type",
                                               order = 2))) |
               strip_legend(p_d)

    bot_row <- strip_legend(p_e) | strip_legend(p_f)

    combined_final <- top_row / mid_row / bot_row
  } else {
    combined_final <- wrap_plots(panels[available], ncol = 2) &
      theme(legend.position = "top")
  }

  combined_final <- combined_final +
    plot_annotation(
      tag_levels = list(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")[available]),
      theme = theme(plot.margin = margin(4, 6, 4, 6))
    ) &
    theme(plot.tag = element_text(size = 12, face = "bold", hjust = 0),
          plot.tag.position = c(0.02, 0.98),
          plot.title = element_text(size = 11, face = "bold"))

  ## Tag filename with family kind so mv vs cpt (vs cpt_7o) figures don't
  ## overwrite each other when plot_ppc.R is called multiple times per study.
  is_mv    <- any(grepl("^mv_",    families))
  is_cpt_o <- any(grepl("_7o",     families))
  is_cpt   <- any(grepl("^cpt_",   families)) && !is_cpt_o
  fam_tag <- if (is_mv && !is_cpt && !is_cpt_o) "mv"
             else if (is_cpt && !is_mv)         "cpt"
             else if (is_cpt_o && !is_mv)      "cpt_7o"
             else paste(families, collapse = "_")

  out_final <- file.path(study_dir,
                          sprintf("posterior_predictive_combined_%s.pdf", fam_tag))
  ## Tighter figure: ~10 x 12 inches instead of 14 x 18
  ggsave(out_final, combined_final, width = 10.5, height = 12)
  ggsave(sub("\\.pdf$", ".png", out_final), combined_final,
         width = 10.5, height = 12, dpi = 300)
  cat("  Saved:", out_final, "\n")
}

## =========================================================================
## RT quantile plot (per family)
## =========================================================================
plot_rt_quantiles <- function(ppc, fam_dir, is_ccss) {
  cat("Plotting RT quantiles...\n")
  models <- levels(ppc$model)

  ## Drop trials where simulation produced NA pred_choice/pred_rt so the
  ## response factor doesn't get a phantom "NA" level in facets.
  ppc <- ppc %>% filter(!is.na(pred_choice), !is.na(pred_rt))

  obs_quant <- ppc %>%
    filter(model == models[1]) %>%
    distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
    filter(!is.na(observed_choice), !is.na(observed_rt))

  if (is_ccss) {
    ## In CCSS, 'Option A' is by convention the one with higher-variance
    ## outcomes after .reorder_ccss(). Label as A / B to match Stan naming.
    obs_quant <- obs_quant %>%
      mutate(condition = ifelse(con == 1, "CC", "SS"),
             response  = ifelse(observed_choice == 1, "Option A", "Option B"))
  } else {
    obs_quant <- obs_quant %>%
      mutate(condition = "CS",
             response  = ifelse(observed_choice == 1, "Chose A", "Chose B"))
  }

  obs_q <- obs_quant %>%
    group_by(condition, response) %>%
    summarise(
      across(observed_rt,
             list(q10 = ~quantile(., 0.1, na.rm = TRUE),
                  q30 = ~quantile(., 0.3, na.rm = TRUE),
                  q50 = ~quantile(., 0.5, na.rm = TRUE),
                  q70 = ~quantile(., 0.7, na.rm = TRUE),
                  q90 = ~quantile(., 0.9, na.rm = TRUE))),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = starts_with("observed_rt_q"),
                 names_to = "quantile", values_to = "rt") %>%
    mutate(quantile = as.numeric(sub("observed_rt_q", "", quantile)))

  pred_q_raw <- ppc
  if (is_ccss) {
    pred_q_raw <- pred_q_raw %>%
      mutate(condition = ifelse(con == 1, "CC", "SS"),
             response  = ifelse(pred_choice == 1, "Option A", "Option B"))
  } else {
    pred_q_raw <- pred_q_raw %>%
      mutate(condition = "CS",
             response  = ifelse(pred_choice == 1, "Chose A", "Chose B"))
  }

  pred_q <- pred_q_raw %>%
    group_by(model_display, sample_id, condition, response) %>%
    summarise(
      q10 = quantile(pred_rt, 0.1, na.rm = TRUE),
      q30 = quantile(pred_rt, 0.3, na.rm = TRUE),
      q50 = quantile(pred_rt, 0.5, na.rm = TRUE),
      q70 = quantile(pred_rt, 0.7, na.rm = TRUE),
      q90 = quantile(pred_rt, 0.9, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = starts_with("q"), names_to = "quantile",
                 values_to = "rt", names_prefix = "q") %>%
    mutate(quantile = as.numeric(quantile)) %>%
    group_by(model_display, condition, response, quantile) %>%
    summarise(
      mean_rt = mean(rt, na.rm = TRUE),
      lower   = quantile(rt, 0.025, na.rm = TRUE),
      upper   = quantile(rt, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  colour_pal <- if (is_ccss) CCSS_COLOURS else CS_COLOURS

  p_quant <- ggplot() +
    geom_ribbon(data = pred_q,
                aes(x = quantile, ymin = lower, ymax = upper,
                    fill = model_display, group = model_display),
                alpha = 0.15) +
    geom_line(data = pred_q,
              aes(x = quantile, y = mean_rt, colour = model_display,
                  group = model_display),
              linewidth = 0.7) +
    geom_point(data = obs_q,
               aes(x = quantile, y = rt), colour = "black", size = 2.2) +
    facet_grid(condition ~ response) +
    scale_x_continuous(breaks = c(10, 30, 50, 70, 90),
                       labels = c("10%", "30%", "50%", "70%", "90%")) +
    scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50),
                  labels = c("0.5", "1", "2", "5", "10", "20", "50")) +
    scale_colour_manual(values = colour_pal, name = "Model") +
    scale_fill_manual(values = colour_pal, name = "Model") +
    labs(title = "Posterior Predictive Check: RT Quantiles",
         x = "RT percentile",
         y = "response time (seconds, log scale)",
         caption = paste0("Black points: observed RT quantiles (10th-90th). ",
                          "Lines: posterior predictive mean; ",
                          "bands: 95% credible interval across draws.")) +
    theme_apa() +
    theme(plot.caption = element_text(size = 8, hjust = 0, colour = "grey40",
                                       margin = margin(t = 6)),
          plot.title   = element_text(size = 12, face = "bold"))

  n_cond <- length(unique(obs_q$condition))
  ggsave(file.path(fam_dir, "ppc_rt_quantiles.pdf"), p_quant,
         width = 8, height = 3 + 2.5 * n_cond)
  ggsave(file.path(fam_dir, "ppc_rt_quantiles.png"), p_quant,
         width = 8, height = 3 + 2.5 * n_cond, dpi = 300)
  cat("  Saved ppc_rt_quantiles\n")
}

## =========================================================================
## RT distribution density plot (per family, per model)
## =========================================================================
plot_rt_distribution <- function(ppc, fam_dir, is_ccss) {
  cat("Plotting RT distribution...\n")
  models <- levels(ppc$model)

  for (mdl in models) {
    ppc_m <- ppc %>% filter(model == mdl)
    mdl_display <- unique(ppc_m$model_display)

    obs_data <- ppc_m %>%
      distinct(participant, trial_in_participant, .keep_all = TRUE)

    if (is_ccss) {
      obs_data <- obs_data %>%
        mutate(Condition = ifelse(observed_choice == 1, "Choosing A", "Choosing B"))
      pred_data <- ppc_m %>%
        mutate(Condition = ifelse(pred_choice == 1, "Choosing A", "Choosing B"))
    } else {
      has_cc <- "chose_complex_obs" %in% names(obs_data)
      if (has_cc) {
        obs_data <- obs_data %>%
          mutate(Condition = ifelse(chose_complex_obs == 1,
                                    "Choosing Complex", "Choosing Simple"))
        pred_data <- ppc_m %>%
          mutate(Condition = ifelse(chose_complex_pred == 1,
                                    "Choosing Complex", "Choosing Simple"))
      } else {
        obs_data <- obs_data %>%
          mutate(Condition = ifelse(observed_choice == 1, "Choosing A", "Choosing B"))
        pred_data <- ppc_m %>%
          mutate(Condition = ifelse(pred_choice == 1, "Choosing A", "Choosing B"))
      }
    }

    ## Aggregate predicted densities across draws
    p_dist <- ggplot() +
      geom_histogram(data = obs_data,
                     aes(x = observed_rt, y = after_stat(density), fill = Condition),
                     alpha = 0.4, bins = 50, position = "identity") +
      geom_density(data = pred_data,
                   aes(x = pred_rt, colour = Condition),
                   linewidth = 0.8) +
      scale_fill_manual(values = c("#FF9999", "#66CCCC")) +
      scale_colour_manual(values = c("#CC3333", "#009999")) +
      coord_cartesian(xlim = c(0, 20)) +
      labs(title = paste("Empirical vs. Predicted (", mdl_display, ") RT Distributions"),
           x = "RT", y = "Density") +
      theme_apa()

    ggsave(file.path(fam_dir, paste0("ppc_rt_distribution_", mdl, ".pdf")),
           p_dist, width = 7, height = 5)
    ggsave(file.path(fam_dir, paste0("ppc_rt_distribution_", mdl, ".png")),
           p_dist, width = 7, height = 5, dpi = 300)
  }
  cat("  Saved ppc_rt_distribution\n")
}

## =========================================================================
## Individual-level PPC plots (per family)
## =========================================================================
plot_individual <- function(ppc, fam_dir, is_ccss) {
  cat("Plotting individual-level...\n")
  models <- levels(ppc$model)

  ## ---- Choice proportion scatter ----
  indiv_obs <- ppc %>%
    filter(model == models[1]) %>%
    distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
    group_by(participant) %>%
    summarise(obs_choice_prop = mean(observed_choice == 1, na.rm = TRUE),
              obs_median_rt   = median(observed_rt, na.rm = TRUE),
              .groups = "drop")

  indiv_pred <- ppc %>%
    group_by(model_display, participant) %>%
    summarise(pred_choice_prop = mean(pred_choice == 1, na.rm = TRUE),
              pred_median_rt   = median(pred_rt, na.rm = TRUE),
              .groups = "drop") %>%
    left_join(indiv_obs, by = "participant")

  ## Choice proportion scatter
  p_choice <- ggplot(indiv_pred,
                      aes(x = obs_choice_prop, y = pred_choice_prop)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = "blue", alpha = 0.5) +
    geom_vline(xintercept = 0.5, linetype = "dotted", colour = "blue", alpha = 0.5) +
    geom_point(alpha = 0.6, size = 1.5, colour = "black") +
    facet_wrap(~model_display, ncol = 2) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = "Observed vs. Predicted Choice Proportion",
         x = "Observed Choice Proportion",
         y = "Predicted Choice Proportion") +
    theme_apa() + theme(legend.position = "none")

  n_models <- length(unique(indiv_pred$model_display))
  plot_nrow <- ceiling(n_models / 2)
  ggsave(file.path(fam_dir, "ppc_individual_choice.pdf"), p_choice,
         width = 8, height = 4 * plot_nrow + 0.5)
  ggsave(file.path(fam_dir, "ppc_individual_choice.png"), p_choice,
         width = 8, height = 4 * plot_nrow + 0.5, dpi = 300)

  ## ---- RT scatter ----
  p_rt <- ggplot(indiv_pred,
                  aes(x = obs_median_rt, y = pred_median_rt)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
    geom_point(alpha = 0.6, size = 1.5, colour = "black") +
    facet_wrap(~model_display, ncol = 2) +
    labs(title = "Observed vs. Predicted Median RT",
         x = "Observed Median RT (s)",
         y = "Predicted Median RT (s)") +
    theme_apa() + theme(legend.position = "none")

  ggsave(file.path(fam_dir, "ppc_individual_rt.pdf"), p_rt,
         width = 8, height = 4 * plot_nrow + 0.5)
  ggsave(file.path(fam_dir, "ppc_individual_rt.png"), p_rt,
         width = 8, height = 4 * plot_nrow + 0.5, dpi = 300)

  ## ---- CS-specific: Complex choice per participant ----
  if (!is_ccss && "chose_complex_obs" %in% names(ppc)) {
    indiv_complex_obs <- ppc %>%
      filter(model == models[1]) %>%
      distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
      group_by(participant) %>%
      summarise(obs_complex_prop = mean(chose_complex_obs == 1, na.rm = TRUE),
                .groups = "drop")

    indiv_complex_pred <- ppc %>%
      group_by(model_display, participant) %>%
      summarise(pred_complex_prop = mean(chose_complex_pred == 1, na.rm = TRUE),
                .groups = "drop") %>%
      left_join(indiv_complex_obs, by = "participant")

    p_complex <- ggplot(indiv_complex_pred,
                         aes(x = obs_complex_prop, y = pred_complex_prop)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
      geom_hline(yintercept = 0.5, linetype = "dotted", colour = "blue", alpha = 0.5) +
      geom_point(alpha = 0.6, size = 1.5, colour = "black") +
      facet_wrap(~model_display, ncol = 2) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(title = "Observed vs. Predicted Complex Choice Proportion",
           x = "Observed Complex Choice Proportion",
           y = "Predicted Complex Choice Proportion") +
      theme_apa() + theme(legend.position = "none")

    ggsave(file.path(fam_dir, "ppc_individual_complex.pdf"), p_complex,
           width = 8, height = 4 * plot_nrow + 0.5)
    ggsave(file.path(fam_dir, "ppc_individual_complex.png"), p_complex,
           width = 8, height = 4 * plot_nrow + 0.5, dpi = 300)
  }

  ## ---- RT difference per participant ----
  if (is_ccss) {
    ## RT diff = median RT(CC) - median RT(SS) per participant
    indiv_rt_diff_obs <- ppc %>%
      filter(model == models[1]) %>%
      distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
      group_by(participant, con) %>%
      summarise(med_rt = median(observed_rt, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = con, values_from = med_rt, names_prefix = "rt_") %>%
      mutate(obs_rt_diff = `rt_1` - `rt_-1`)

    indiv_rt_diff_pred <- ppc %>%
      group_by(model_display, sample_id, participant, con) %>%
      summarise(med_rt = median(pred_rt, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = con, values_from = med_rt, names_prefix = "rt_") %>%
      mutate(pred_rt_diff = `rt_1` - `rt_-1`) %>%
      group_by(model_display, participant) %>%
      summarise(pred_rt_diff = mean(pred_rt_diff, na.rm = TRUE), .groups = "drop") %>%
      left_join(indiv_rt_diff_obs %>% select(participant, obs_rt_diff),
                by = "participant")

    p_rt_diff <- ggplot(indiv_rt_diff_pred,
                         aes(x = obs_rt_diff, y = pred_rt_diff)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "blue", alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "blue", alpha = 0.5) +
      geom_point(alpha = 0.6, size = 1.5, colour = "black") +
      facet_wrap(~model_display, ncol = 2) +
      labs(title = "Observed vs. Predicted RT Difference",
           x = "Observed RT Difference (CC - SS, s)",
           y = "Predicted RT Difference (s)") +
      theme_apa() + theme(legend.position = "none")

    ggsave(file.path(fam_dir, "ppc_individual_rt_diff.pdf"), p_rt_diff,
           width = 8, height = 4 * plot_nrow + 0.5)
    ggsave(file.path(fam_dir, "ppc_individual_rt_diff.png"), p_rt_diff,
           width = 8, height = 4 * plot_nrow + 0.5, dpi = 300)
  } else if ("chose_complex_obs" %in% names(ppc)) {
    ## CS: RT diff = median RT(chose complex) - median RT(chose simple) per participant
    indiv_rt_diff_obs <- ppc %>%
      filter(model == models[1]) %>%
      distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
      group_by(participant, chose_complex_obs) %>%
      summarise(med_rt = median(observed_rt, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = chose_complex_obs, values_from = med_rt,
                  names_prefix = "rt_") %>%
      mutate(obs_rt_diff = `rt_1` - `rt_-1`)

    indiv_rt_diff_pred <- ppc %>%
      group_by(model_display, sample_id, participant, chose_complex_pred) %>%
      summarise(med_rt = median(pred_rt, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = chose_complex_pred, values_from = med_rt,
                  names_prefix = "rt_") %>%
      mutate(pred_rt_diff = `rt_1` - `rt_-1`) %>%
      group_by(model_display, participant) %>%
      summarise(pred_rt_diff = mean(pred_rt_diff, na.rm = TRUE), .groups = "drop") %>%
      left_join(indiv_rt_diff_obs %>% select(participant, obs_rt_diff),
                by = "participant")

    p_rt_diff <- ggplot(indiv_rt_diff_pred,
                         aes(x = obs_rt_diff, y = pred_rt_diff)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "blue", alpha = 0.5) +
      geom_point(alpha = 0.6, size = 1.5, colour = "black") +
      facet_wrap(~model_display, ncol = 2) +
      labs(title = "Observed vs. Predicted RT Difference",
           x = "Observed RT Difference (Complex - Simple, s)",
           y = "Predicted RT Difference (s)") +
      theme_apa() + theme(legend.position = "none")

    ggsave(file.path(fam_dir, "ppc_individual_rt_diff.pdf"), p_rt_diff,
           width = 8, height = 4 * plot_nrow + 0.5)
    ggsave(file.path(fam_dir, "ppc_individual_rt_diff.png"), p_rt_diff,
           width = 8, height = 4 * plot_nrow + 0.5, dpi = 300)
  }

  ## ---- CCSS: EV consistency per participant ----
  if (is_ccss && "evd" %in% names(ppc)) {
    ## EV consistency = proportion of EV-consistent choices per participant
    indiv_ev_obs <- ppc %>%
      filter(model == models[1]) %>%
      distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
      mutate(is_ev_consistent = as.integer(
        (evd > 0 & observed_choice == 1) | (evd < 0 & observed_choice == -1))) %>%
      group_by(participant) %>%
      summarise(obs_ev_prop = mean(is_ev_consistent, na.rm = TRUE), .groups = "drop")

    indiv_ev_pred <- ppc %>%
      mutate(is_ev_consistent = as.integer(
        (evd > 0 & pred_choice == 1) | (evd < 0 & pred_choice == -1))) %>%
      group_by(model_display, participant) %>%
      summarise(pred_ev_prop = mean(is_ev_consistent, na.rm = TRUE),
                .groups = "drop") %>%
      left_join(indiv_ev_obs, by = "participant")

    p_ev <- ggplot(indiv_ev_pred,
                    aes(x = obs_ev_prop, y = pred_ev_prop)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
      geom_hline(yintercept = 0.5, linetype = "dotted", colour = "blue", alpha = 0.5) +
      geom_point(alpha = 0.6, size = 1.5, colour = "black") +
      facet_wrap(~model_display, ncol = 2) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(title = "Observed vs. Predicted EV Consistency",
           x = "Observed EV-Consistent Choice Proportion",
           y = "Predicted EV-Consistent Choice Proportion") +
      theme_apa() + theme(legend.position = "none")

    ggsave(file.path(fam_dir, "ppc_individual_consistency.pdf"), p_ev,
           width = 8, height = 4 * plot_nrow + 0.5)
    ggsave(file.path(fam_dir, "ppc_individual_consistency.png"), p_ev,
           width = 8, height = 4 * plot_nrow + 0.5, dpi = 300)
  }

  cat("  Saved individual-level plots\n")
}

## =========================================================================
## CS-only: complex choice proportion across 5 EVD bins (standalone plot)
## =========================================================================
plot_complex_by_evd <- function(ppc, fam_dir) {
  cat("Plotting complex-choice by EVD...\n")

  if (!"evd" %in% names(ppc)) {
    cat("  SKIP - `evd` not in posterior_predictives.csv; re-run generate_ppc.R\n")
    return(invisible(NULL))
  }
  ppc <- ppc %>%
    mutate(evd_bin = cut(evd, breaks = EVD_BREAKS, labels = EVD_LABELS)) %>%
    filter(!is.na(evd_bin))

  ## Observed complex choice proportion per EVD bin (one row per bin).
  obs_cbe <- ppc %>%
    filter(model == levels(model)[1]) %>%
    distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
    group_by(evd_bin) %>%
    summarise(obs_prop = mean(chose_complex_obs == 1, na.rm = TRUE),
              .groups = "drop")

  ## Predicted: per draw, compute proportion within each bin, then CrI.
  pred_cbe <- ppc %>%
    group_by(model_display, sample_id, evd_bin) %>%
    summarise(prop = mean(chose_complex_pred == 1, na.rm = TRUE),
              .groups = "drop") %>%
    group_by(model_display, evd_bin) %>%
    summarise(mean_prop = mean(prop),
              lower = quantile(prop, 0.025),
              upper = quantile(prop, 0.975),
              .groups = "drop")

  p <- ggplot() +
    geom_col(data = obs_cbe, aes(x = evd_bin, y = obs_prop),
             fill = "#B2D8D8", width = 0.6, alpha = 0.7) +
    geom_point(data = pred_cbe,
               aes(x = evd_bin, y = mean_prop, colour = model_display),
               position = position_dodge(width = 0.6), size = 2) +
    geom_errorbar(data = pred_cbe,
                  aes(x = evd_bin, ymin = lower, ymax = upper,
                      colour = model_display),
                  position = position_dodge(width = 0.6),
                  width = 0.2, linewidth = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = CS_COLOURS, name = "Model Type") +
    labs(title = "Posterior Predictive Check: Complex Choice across EVD Levels",
         subtitle = "Observed vs. predicted complex-choice proportion by EV_complex - EV_simple",
         x = "EVD Levels  (EV_complex - EV_simple)",
         y = "Complex Choice Proportion") +
    theme_apa() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8))

  ggsave(file.path(fam_dir, "ppc_complex_by_evd.pdf"), p,
         width = 8, height = 5)
  ggsave(file.path(fam_dir, "ppc_complex_by_evd.png"), p,
         width = 8, height = 5, dpi = 300)
  cat("  Saved ppc_complex_by_evd\n")
}

## =========================================================================
## Generate per-family plots
## =========================================================================
for (fam in families) {
  fam_dir <- file.path(study_dir, fam)
  is_ccss <- grepl("ccss", fam)

  if (is_ccss && !is.null(ppc_ccss)) {
    plot_rt_quantiles(ppc_ccss, fam_dir, is_ccss = TRUE)
    plot_rt_distribution(ppc_ccss, fam_dir, is_ccss = TRUE)
    plot_individual(ppc_ccss, fam_dir, is_ccss = TRUE)
  } else if (!is_ccss && !is.null(ppc_cs)) {
    plot_rt_quantiles(ppc_cs, fam_dir, is_ccss = FALSE)
    plot_rt_distribution(ppc_cs, fam_dir, is_ccss = FALSE)
    plot_individual(ppc_cs, fam_dir, is_ccss = FALSE)
    plot_complex_by_evd(ppc_cs, fam_dir)
  }
}

cat("\nAll PPC plots complete.\n")
