## ============================================================================
## R/summarize.R — merge truth with posterior, compute stats, and plot
##
## Aim:     Turn cmdstanr draws into a tidy true-vs-estimated table, join the
##          generating parameters, compute recovery statistics (r/mae/rmse and
##          credible-interval coverage), and build the publication-quality
##          participant-level recovery scatter grid plus its APA caption.
## Inputs:  posterior draws + the true params/hyperparams saved by simulate.R.
## Outputs: in-memory tibbles and a ggplot; run_recovery.R / the scripts write
##          recovery_long.csv, recovery_stats.csv, recovery.pdf/.png and the
##          caption text.
## Usage:   source("R/summarize.R"); posterior_table(); merge_truth();
##          plot_participant_recovery()  (called from run_recovery.R).
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM parameter
## recovery). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================
##
## Two comparisons are reported per recovery dataset:
##   (1) Group-level: posterior mean of mu and sigma vs. the mu_true / sigma_true
##       that generated the dataset.
##   (2) Participant-level: posterior mean of participant_params[k, l] vs. the
##       raw true parameter that generated participant l.
##
## All comparisons are done on the RAW (unconstrained) scale — the space Stan
## samples in — so true and estimated values live on the same scale.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(posterior)
})

## ---------------------------------------------------------------------------
## Pull posterior means and quantiles for mu, sigma, participant_params.
## Returns a tidy long data.frame with columns:
##   level     -- "mu" | "sigma" | "participant"
##   param     -- raw-scale parameter name
##   id        -- participant index (NA for mu/sigma)
##   est_mean  -- posterior mean
##   est_sd    -- posterior sd
##   q2.5, q97.5 -- 95% credible interval
## ---------------------------------------------------------------------------
posterior_table <- function(draws, param_names) {
  K <- length(param_names)
  summ <- summarise_draws(draws, mean, sd,
                          ~quantile2(.x, probs = c(0.025, 0.975)))
  ## posterior >=1.5 names these "q2.5"/"q97.5"; older versions used
  ## "2.5%"/"97.5%". Handle both.
  names(summ)[names(summ) %in% c("2.5%",  "q2.5")]  <- "q2.5"
  names(summ)[names(summ) %in% c("97.5%", "q97.5")] <- "q97.5"

  ## parse variable names like "mu[3]", "sigma[1]", "participant_params[2,5]"
  parse_one <- function(v) {
    if (grepl("^mu\\[", v)) {
      k <- as.integer(sub("mu\\[(\\d+)\\]", "\\1", v))
      list(level = "mu", param = param_names[k], id = NA_integer_)
    } else if (grepl("^sigma\\[", v)) {
      k <- as.integer(sub("sigma\\[(\\d+)\\]", "\\1", v))
      list(level = "sigma", param = param_names[k], id = NA_integer_)
    } else if (grepl("^participant_params\\[", v)) {
      m <- regmatches(v, regexec("participant_params\\[(\\d+),(\\d+)\\]", v))[[1]]
      k <- as.integer(m[2]); l <- as.integer(m[3])
      list(level = "participant", param = param_names[k], id = l)
    } else {
      list(level = NA, param = NA, id = NA_integer_)
    }
  }
  parsed <- do.call(rbind, lapply(summ$variable, function(v) as.data.frame(parse_one(v))))
  out <- cbind(parsed, summ[, c("mean", "sd", "q2.5", "q97.5")])
  names(out)[names(out) == "mean"] <- "est_mean"
  names(out)[names(out) == "sd"]   <- "est_sd"
  dplyr::as_tibble(out) %>% dplyr::filter(!is.na(level))
}

## ---------------------------------------------------------------------------
## Merge estimates with true values.
##
## Inputs:
##   post_tbl    -- output of posterior_table()
##   true_params -- data.frame [L, K] of raw true participant params
##   mu_true     -- named numeric (length K)
##   sigma_true  -- named numeric (length K)
##
## Returns a tibble ready to plot/summarise.
## ---------------------------------------------------------------------------
merge_truth <- function(post_tbl, true_params, mu_true, sigma_true) {
  ## participant-level truth in long form
  pl_truth <- true_params %>%
    tibble::rownames_to_column("id") %>%
    mutate(id = as.integer(id)) %>%
    tidyr::pivot_longer(-id, names_to = "param", values_to = "true_value") %>%
    mutate(level = "participant")

  group_truth <- dplyr::bind_rows(
    tibble::tibble(level = "mu",
                   param = names(mu_true),
                   id = NA_integer_,
                   true_value = as.numeric(mu_true)),
    tibble::tibble(level = "sigma",
                   param = names(sigma_true),
                   id = NA_integer_,
                   true_value = as.numeric(sigma_true))
  )

  truth <- dplyr::bind_rows(pl_truth, group_truth)
  dplyr::left_join(post_tbl, truth, by = c("level", "param", "id"))
}

## ---------------------------------------------------------------------------
## Recovery correlations per-parameter, per-level.
## Returns tibble with columns: level, param, r, mae, rmse, n.
## ---------------------------------------------------------------------------
recovery_stats <- function(merged) {
  merged %>%
    dplyr::group_by(level, param) %>%
    dplyr::summarise(
      n    = dplyr::n(),
      r    = if (dplyr::n() > 2) suppressWarnings(cor(true_value, est_mean)) else NA_real_,
      mae  = mean(abs(true_value - est_mean), na.rm = TRUE),
      rmse = sqrt(mean((true_value - est_mean)^2, na.rm = TRUE)),
      .groups = "drop"
    )
}

## ---------------------------------------------------------------------------
## Plotting helpers — publication-quality recovery scatter grid.
## Only participant-level parameters are shown (group-level mu/sigma require
## many datasets for r to be interpretable, and if individual recovery is
## good the hyperparameters inherit that information).
## ---------------------------------------------------------------------------

## Greek / mathematical labels for every raw parameter in every model.
## Values are plotmath expressions used with label_parsed / parse(text=).
PARAM_LABELS <- c(
  ## common
  beta_raw         = "beta",
  beta             = "beta",
  theta_raw        = "theta",
  threshold_raw    = "alpha",
  ndt_raw          = "tau",
  gamma_raw        = "gamma",
  eta              = "eta",
  zeta             = "zeta",
  sp_raw           = "z[sp]",
  ## deltas
  delta_beta       = "Delta*beta",
  delta_theta      = "Delta*theta",
  delta_threshold  = "Delta*alpha",
  delta_gamma      = "Delta*gamma",
  delta_eta        = "Delta*eta"
)

greek_label <- function(name) {
  out <- PARAM_LABELS[name]
  out[is.na(out)] <- name[is.na(out)]
  out
}

## ---------------------------------------------------------------------------
## Participant parameter-recovery plot
## Style: Baribault & Collins (2025, Fig. 9), formatted to APA 7 standards.
##
## APA 7 figure conventions applied:
##   - Sans-serif font throughout (family = "sans"); caps on family via theme.
##   - No figure title inside the graphic (APA puts the title in the caption).
##   - Lowercase panel letters (a, b, c, ...) in the top-left of each panel.
##   - Axis labels in sentence case ("true value", "estimated value").
##   - Tick marks pointing outward; axis lines visible; no panel borders
##     beyond the two axes.
##   - Colour used only to convey data (red = interval misses the true value).
##   - No gridlines; no drop shadows; no unnecessary decoration.
##   - Figure-quality text sizes (8 – 12 pt in-figure).
##
## Use `recovery_caption()` to generate the accompanying APA-style caption.
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## APA-style caption string for the recovery figure.
## Returns a single string you can paste into a manuscript; the figure number
## and model name should be customised at call time.
## ---------------------------------------------------------------------------
recovery_caption <- function(model_label, n_datasets, n_subjects,
                             figure_number = 1) {
  sprintf(
paste0(
"Figure %d\n\n",
"Parameter recovery for the %s model\n\n",
"Note. Each point represents one simulated participant from %d simulated ",
"data sets (%d participants per data set). The x-axis shows the true value ",
"used to generate the data; the y-axis shows the posterior mean recovered by ",
"the hierarchical Bayesian model. Vertical lines indicate 95%% Bayesian ",
"credible intervals. The dashed diagonal marks perfect recovery. Intervals ",
"that exclude the true value are shown in red with an \u00d7 marker; ",
"intervals that cover the true value are shown in black. Each panel's ",
"top-left badge reports the proportion of 95%% credible intervals that ",
"covered the true value (nominal coverage \u2248 95%%, following Rubin, 1984). ",
"Recovery is evaluated by the clustering of posterior means around the ",
"identity line and by the calibration of the 95%% credible intervals."
),
    figure_number, model_label, n_datasets, n_subjects
  )
}

plot_participant_recovery <- function(merged,
                                      ncol            = NULL,
                                      colour_cover    = "black",
                                      colour_miss     = "#C1272D",
                                      colour_identity = "grey55",
                                      pad_frac        = 0.04,
                                      base_family     = "sans") {

  df <- merged %>% dplyr::filter(level == "participant")

  ## determine coverage per observation
  df <- df %>%
    dplyr::mutate(
      covers = (true_value >= q2.5) & (true_value <= q97.5)
    )

  ## preserve parameter order; wrap greek labels for facet strips,
  ## prefixing each with an APA-style lowercase panel letter "(a)  β" etc.
  param_order <- unique(df$param)
  df <- df %>% dplyr::mutate(param = factor(param, levels = param_order))
  n_params <- length(param_order)
  panel_letters <- letters[seq_len(n_params)]

  ## plotmath expressions: e.g.  bold("(a)")~beta
  fac_labels <- setNames(
    vapply(seq_along(param_order), function(i) {
      sprintf('bold("(%s)")~%s',
              panel_letters[i],
              greek_label(as.character(param_order[i])))
    }, character(1)),
    param_order
  )

  ## common per-facet xy range so panels are square with equal limits
  ranges <- df %>%
    dplyr::group_by(param) %>%
    dplyr::summarise(
      lo = min(c(true_value, q2.5,  est_mean), na.rm = TRUE),
      hi = max(c(true_value, q97.5, est_mean), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      span = hi - lo,
      lo   = lo - pad_frac * span,
      hi   = hi + pad_frac * span
    )

  ## coverage-% + rho badge per facet
  badges <- df %>%
    dplyr::group_by(param) %>%
    dplyr::summarise(
      n      = dplyr::n(),
      pct    = 100 * mean(covers, na.rm = TRUE),
      rho    = suppressWarnings(cor(true_value, est_mean, use = "complete.obs")),
      .groups = "drop"
    ) %>%
    dplyr::left_join(ranges, by = "param") %>%
    dplyr::mutate(
      x      = lo + 0.02 * (hi - lo),
      y_cov  = hi - 0.02 * (hi - lo),
      y_rho  = hi - 0.08 * (hi - lo),
      label_cov = sprintf("%.0f%% in 95%% CrI", pct),
      label_rho = sprintf("\u03C1 = %.3f", rho)
    )

  ## facet-specific identity lines that match each panel's range exactly
  identity_df <- ranges %>%
    dplyr::mutate(x0 = lo, x1 = hi, y0 = lo, y1 = hi)

  if (is.null(ncol)) {
    ncol <- min(3, length(param_order))
  }

  ## "blank" ghost points at the four corners of each facet's square range:
  ## forces x and y scale-free panels to use matched limits so each facet
  ## looks square (identity line is truly 45°).
  corners <- do.call(rbind, lapply(seq_len(nrow(ranges)), function(i) {
    r <- ranges[i, ]
    data.frame(param   = r$param,
               xcorner = c(r$lo, r$lo, r$hi, r$hi),
               ycorner = c(r$lo, r$hi, r$lo, r$hi))
  }))
  corners$param <- factor(corners$param, levels = param_order)

  ## adaptive transparency: more points → more transparent
  ## tuned so 20 pts looks solid, 200 pts is semi-transparent, 1500 is light
  ## but still clearly visible as a dark diagonal band
  n_per_panel <- nrow(df) / n_params
  alpha_ci    <- pmin(0.70, pmax(0.08, 8  / n_per_panel^0.50))
  alpha_point <- pmin(0.85, pmax(0.15, 12 / n_per_panel^0.50))
  point_size  <- pmin(1.8,  pmax(0.6,  12 / n_per_panel^0.40))
  ci_width    <- pmin(0.35, pmax(0.12, 4  / n_per_panel^0.40))

  ## horizontal jitter to separate overlapping points; proportional
  ## to each facet's range — small enough to preserve the diagonal pattern
  jitter_amount <- ranges$span * 0.012
  jitter_lookup <- setNames(jitter_amount, ranges$param)
  df <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      jit = jitter_lookup[as.character(param)],
      x_jit = true_value + runif(1, -jit, jit)
    ) %>%
    dplyr::ungroup()

  p <- ggplot(df, aes(x = x_jit, y = est_mean)) +
    geom_blank(data = corners, aes(x = xcorner, y = ycorner),
               inherit.aes = FALSE) +
    geom_segment(
      data = identity_df,
      aes(x = x0, xend = x1, y = y0, yend = y1),
      colour = colour_identity, linetype = "dashed", linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    ## covered CrIs: transparent, background layer
    geom_linerange(
      data = function(d) d[d$covers, ],
      aes(ymin = q2.5, ymax = q97.5),
      colour = colour_cover, linewidth = ci_width, alpha = alpha_ci
    ) +
    geom_point(
      data = function(d) d[d$covers, ],
      aes(x = x_jit, y = est_mean),
      colour = colour_cover, shape = 16, fill = colour_cover,
      size = point_size, stroke = 0.3, alpha = alpha_point
    ) +
    ## missed CrIs: red on top — visible but not overpowering
    geom_linerange(
      data = function(d) d[!d$covers, ],
      aes(ymin = q2.5, ymax = q97.5),
      colour = colour_miss, linewidth = ci_width * 1.0, alpha = 0.35
    ) +
    geom_point(
      data = function(d) d[!d$covers, ],
      aes(x = x_jit, y = est_mean),
      colour = colour_miss, shape = 4, fill = "white",
      size = point_size * 1.2, stroke = 0.45, alpha = 0.50
    ) +
    geom_text(
      data = badges,
      aes(x = x, y = y_cov, label = label_cov),
      hjust = 0, vjust = 1, size = 4.5, colour = "grey25",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = badges,
      aes(x = x, y = y_rho, label = label_rho),
      hjust = 0, vjust = 1, size = 4.5, colour = "grey25",
      fontface = "italic",
      inherit.aes = FALSE
    ) +
    facet_wrap(~param, scales = "free", ncol = ncol,
               labeller = labeller(param = as_labeller(fac_labels, label_parsed))) +
    labs(x = "true value", y = "estimated value") +
    theme_classic(base_size = 11, base_family = base_family) +
    theme(
      text              = element_text(family = base_family, colour = "black"),
      strip.background  = element_blank(),
      strip.placement   = "outside",
      strip.text        = element_text(size = 12, hjust = 0,
                                       margin = margin(b = 4, t = 2)),
      panel.border      = element_blank(),
      panel.grid        = element_blank(),
      axis.line         = element_line(colour = "black", linewidth = 0.4),
      axis.ticks        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(-0.12, "cm"),    # APA: ticks point outward
      axis.text         = element_text(colour = "black", size = 9,
                                       margin = margin(t = 4, r = 4)),
      axis.text.x       = element_text(margin = margin(t = 6)),
      axis.text.y       = element_text(margin = margin(r = 6)),
      axis.title        = element_text(size = 11),
      axis.title.x      = element_text(margin = margin(t = 8)),
      axis.title.y      = element_text(margin = margin(r = 8)),
      plot.margin       = margin(8, 12, 8, 8)
    )

  ## apply per-facet xy limits after coord_fixed via expand_limits hack is not
  ## needed: coord_fixed with scales="free" already lets each facet use its
  ## own range. The padded lo/hi are enforced implicitly through the identity
  ## segment and the point/CI data.
  p
}
