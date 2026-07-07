#!/usr/bin/env Rscript
## ============================================================================
## generate_ppc.R — simulate posterior-predictive choices/RTs from a fitted DDM
##
## Aim:     For each subsampled posterior draw, extract per-participant raw
##          params, map them to effective DDM params via the same transform
##          functions used in fitting, and simulate choice + RT with
##          rtdists::rdiffusion. Produces the large per-model PPC table that all
##          downstream plot_ppc*/summarize_ppc scripts consume.
## Inputs:  results/<study>/<family>/<model>/draws.rds (posterior draws),
##          results/<study>/<family>/<model>/meta.rds (records the data file),
##          the trial data RDS (from meta.rds, --data_rds, or
##          data/final_df_<study>.rds), and R/ helpers (model_registry,
##          preprocess, transforms, param_names).
## Outputs: results/<study>/<family>/<model>/posterior_predictives.csv
## Usage:   Rscript 04_ppc/generate_ppc.R \
##            --study study2 --family mv_ccss --model baseline \
##            --tail_per_chain 800 --thin 2 --cores 8
##
##          # or from R:
##          PPC_ARGS <- c("--study", "study2", "--family", "mv_ccss",
##                        "--model", "baseline")
##          source("04_ppc/generate_ppc.R")
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (posterior predictive
## checks). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(posterior)
  library(rtdists)
  library(parallel); library(pbmcapply)
})

source("R/model_registry.R")
source("R/preprocess.R")
source("R/transforms.R")
source("R/param_names.R")

## ---------------------------------------------------------------------------
## Helper: build trial-level list from stan_data for effective_* functions
## ---------------------------------------------------------------------------
.build_trials_df <- function(family, sd) {
  if (family == "cpt_ccss") {
    list(
      con      = sd$con,
      o_risky1 = sd$o_risky[, 1], o_risky2 = sd$o_risky[, 2],
      p_risky1 = sd$p_risky[, 1], p_risky2 = sd$p_risky[, 2],
      o_safe1  = sd$o_safe[, 1],  o_safe2  = sd$o_safe[, 2],
      p_safe1  = sd$p_safe[, 1],  p_safe2  = sd$p_safe[, 2]
    )
  } else if (family == "cpt_ccss_7o") {
    ## Pass the full 9-col matrices so effective_cpt_ccss_7o_* can branch on con.
    list(
      con = sd$con,
      o_risky = sd$o_risky, o_safe = sd$o_safe,
      p_risky = sd$p_risky, p_safe = sd$p_safe
    )
  } else if (family == "mv_ccss") {
    list(con = sd$con, evd = sd$evd, sdd = sd$sdd, skew = sd$skew)
  } else if (family == "cpt_cs") {
    list(
      o_complex1 = sd$o_complex[, 1], o_complex2 = sd$o_complex[, 2],
      p_complex1 = sd$p_complex[, 1], p_complex2 = sd$p_complex[, 2],
      o_simple1  = sd$o_simple[, 1],  o_simple2  = sd$o_simple[, 2],
      p_simple1  = sd$p_simple[, 1],  p_simple2  = sd$p_simple[, 2]
    )
  } else if (family == "cpt_cs_7o") {
    ## 7-outcome complex + 2-outcome simple.
    list(
      o_complex = sd$o_complex, p_complex = sd$p_complex,
      o_simple  = sd$o_simple,  p_simple  = sd$p_simple
    )
  } else if (family == "mv_cs") {
    list(evd = sd$evd, sdd = sd$sdd, skew_c = sd$skew_c, skew_s = sd$skew_s)
  } else {
    stop("Unknown family: ", family)
  }
}

opt_list <- list(
  make_option("--study",   type = "character", default = NULL),
  make_option("--family",  type = "character", default = NULL),
  make_option("--model",   type = "character", default = NULL),
  ## Draw-selection strategy: take the last `tail_per_chain` iterations of
  ## each chain (most converged portion) then thin by `thin` to decorrelate.
  ## Default 800 tail, thin=2 → 400 draws per chain × n_chains (=1600 on the
  ## standard 4-chain/2000-iter fit). Raise tail_per_chain for more coverage
  ## or thin for tighter decorrelation.
  make_option("--tail_per_chain", type = "integer", default = 800,
              help = "Final iterations to keep from each chain (pre-thin)"),
  make_option("--thin", type = "integer", default = 2,
              help = "Keep every `thin`-th iteration within the tail"),
  make_option("--cores",   type = "integer", default = 4),
  make_option("--results_dir", type = "character", default = "results"),
  make_option("--data_rds", type = "character", default = NULL,
              help = "Override data file (defaults to meta.rds-recorded path, or data/final_df_<study>.rds)"),
  make_option("--max_rt",  type = "double", default = 30,
              help = "Cap simulated RTs at this value")
)

.override <- if (exists("PPC_ARGS", envir = .GlobalEnv))
  get("PPC_ARGS", envir = .GlobalEnv) else NULL
opt <- parse_args(OptionParser(option_list = opt_list),
                  args = if (is.null(.override)) commandArgs(trailingOnly = TRUE) else .override)

model_dir <- file.path(opt$results_dir, opt$study, opt$family, opt$model)
draws_path <- file.path(model_dir, "draws.rds")
if (!file.exists(draws_path)) stop("draws.rds not found: ", draws_path)

## ---------------------------------------------------------------------------
## 1. Load trial data (same preprocessing as fitting)
##
## Resolution order for the data file:
##   1. --data_rds CLI flag if provided
##   2. data_rds recorded in meta.rds (written by fit_one)
##   3. data/final_df_<study>.rds (legacy fallback)
##
## Using the meta.rds-recorded path ensures PPC uses the EXACT same data
## file that the fit was trained on — critical when the fit was produced
## from a subset (e.g. smoke test using only 10 participants).
## ---------------------------------------------------------------------------
data_rds <- opt$data_rds
if (is.null(data_rds)) {
  meta_path <- file.path(model_dir, "meta.rds")
  if (file.exists(meta_path)) {
    meta <- readRDS(meta_path)
    if (!is.null(meta$data_rds) && file.exists(meta$data_rds)) {
      data_rds <- meta$data_rds
      cat(sprintf("  Using data path from meta.rds: %s\n", data_rds))
    }
  }
  if (is.null(data_rds)) {
    data_rds <- sprintf("data/final_df_%s.rds", opt$study)
  }
}
if (!file.exists(data_rds)) stop("Data file not found: ", data_rds)
pp <- PREPROCESSOR[[opt$family]](data_rds)
stan_data <- pp$stan_data
id_map    <- pp$id_map
N <- stan_data$N
L <- stan_data$L

## Build a trial-level list matching what the effective_* functions need
trials <- .build_trials_df(opt$family, stan_data)

cat(sprintf("PPC generation: %s / %s / %s\n  N=%d trials, L=%d participants, %d draws\n",
            opt$study, opt$family, opt$model, N, L, opt$n_draws))

## ---------------------------------------------------------------------------
## 2. Load posterior draws and subsample
## ---------------------------------------------------------------------------
draws_raw <- readRDS(draws_path)

## ---------------------------------------------------------------------------
## Tail-and-thin subsampling (chain-aware):
##   1. Load as a [iteration, chain, variable] array to preserve chain info.
##   2. For each chain, keep the last `tail_per_chain` iterations — these
##      come from deeper in sampling where the chain is most settled.
##   3. Within that tail, thin by `thin` (stride=2 by default → every other
##      iteration) to decorrelate adjacent highly-correlated MCMC states.
##   4. Flatten the subset to a matrix for downstream indexing.
##
## Output has length(iter_idx) × n_chains total draws. This replaces the
## earlier flat random/evenly-spaced thinning — the tail step biases the
## selection toward the converged region of each chain rather than mixing
## early and late iterations uniformly, and because we apply the same
## (tail, thin) recipe to every chain, the resulting sample is balanced
## across chains and deterministic (no seed needed).
## ---------------------------------------------------------------------------
draws_arr     <- posterior::as_draws_array(draws_raw)
n_iter_chain  <- posterior::niterations(draws_arr)
n_chains      <- posterior::nchains(draws_arr)

tail_n <- min(opt$tail_per_chain, n_iter_chain)
if (tail_n < opt$tail_per_chain) {
  message(sprintf(
    "  --tail_per_chain=%d exceeds iterations/chain=%d; clamping to %d.",
    opt$tail_per_chain, n_iter_chain, tail_n))
}
if (opt$thin < 1) stop("--thin must be >= 1")

## Iteration indices within each chain: last `tail_n` iters, every `thin`-th.
iter_idx <- seq(n_iter_chain - tail_n + 1, n_iter_chain, by = opt$thin)
if (length(iter_idx) == 0)
  stop("No iterations selected — check --tail_per_chain and --thin.")

## Subset array then flatten to matrix. Row order: iteration-major, chain-minor
## per posterior::as_draws_matrix convention (rows = draw, cols = variable).
draws_arr <- draws_arr[iter_idx, , , drop = FALSE]
draws_mat <- posterior::as_draws_matrix(draws_arr)
n_total_draws <- posterior::ndraws(posterior::as_draws_array(draws_raw))  # full (pre-subset) count
draw_idx <- seq_len(nrow(draws_mat))  # row indices into the *subset* matrix

cat(sprintf(
  "  Draws: %d total (%d chains x %d iter) -> tail %d/chain, thin %d -> %d per chain x %d chains = %d draws.\n",
  n_total_draws, n_chains, n_iter_chain,
  tail_n, opt$thin, length(iter_idx), n_chains, nrow(draws_mat)))

## Get model's number of params
reg_row <- MODEL_REGISTRY %>% filter(family == opt$family, model == opt$model)
K <- reg_row$n_params

## Pre-extract participant_params for all selected draws
## participant_params[k,l] in Stan → variable name "participant_params[k,l]"
pp_names <- paste0("participant_params[",
                   rep(1:K, L), ",",
                   rep(1:L, each = K), "]")
pp_mat <- draws_mat[draw_idx, pp_names, drop = FALSE]

## ---------------------------------------------------------------------------
## 3. Get the effective-params function for this family
## ---------------------------------------------------------------------------
eff_fn <- EFFECTIVE[[opt$family]]
if (is.null(eff_fn)) stop("No effective function for family: ", opt$family)

## ---------------------------------------------------------------------------
## 4. Simulate for each draw × participant
## ---------------------------------------------------------------------------
simulate_one_draw <- function(draw_i) {
  row_idx <- draw_i  # index into pp_mat (1..n_draws)

  results <- vector("list", L)
  for (l in seq_len(L)) {
    ## Extract this participant's K raw params for this draw
    col_start <- (l - 1) * K + 1
    col_end   <- l * K
    raw_k <- as.numeric(pp_mat[row_idx, col_start:col_end])

    ## Pad to full param vector (missing deltas → 0), mapping by name so
    ## intermediate-shift models slot into the correct canonical positions.
    raw_full <- build_full_raw(opt$family, raw_k, model = opt$model)

    ## Get per-trial effective params
    ## Filter trials to this participant (handle both vectors and row-wise matrices)
    trial_mask <- (stan_data$participant == l)
    trials_l <- lapply(trials, function(x) {
      if (is.matrix(x) && nrow(x) == N) x[trial_mask, , drop = FALSE]
      else if (length(x) == N)            x[trial_mask]
      else                                 x
    })

    eff <- tryCatch(eff_fn(raw_full, trials_l), error = function(e) NULL)
    if (is.null(eff)) {
      results[[l]] <- NULL
      next
    }

    n_trials_l <- sum(trial_mask)

    ## Simulate choice + RT via rdiffusion
    drift     <- as.numeric(eff$drift)
    threshold <- as.numeric(eff$threshold)
    ndt       <- as.numeric(eff$ndt)
    rel_sp    <- as.numeric(eff$rel_sp)

    pred_rt     <- rep(NA_real_, n_trials_l)
    pred_choice <- rep(NA_integer_, n_trials_l)

    for (t in seq_len(n_trials_l)) {
      if (!is.finite(drift[t]) || !is.finite(threshold[t]) ||
          threshold[t] <= 0 || !is.finite(ndt[t]) || ndt[t] <= 0 ||
          !is.finite(rel_sp[t]) || rel_sp[t] <= 0 || rel_sp[t] >= 1) next

      sim <- tryCatch(
        rtdists::rdiffusion(1, a = threshold[t], v = drift[t],
                            t0 = ndt[t], z = rel_sp[t] * threshold[t]),
        error = function(e) NULL
      )
      if (!is.null(sim) && is.finite(sim$rt) && sim$rt <= opt$max_rt) {
        pred_rt[t]     <- sim$rt
        pred_choice[t] <- if (sim$response == "upper") 1L else -1L
      }
    }

    results[[l]] <- data.frame(
      sample_id   = draw_idx[draw_i],
      participant = l,
      trial_in_participant = seq_len(n_trials_l),
      pred_choice = pred_choice,
      pred_rt     = pred_rt
    )
  }
  do.call(rbind, results)
}

cat("  Simulating...")
t0 <- Sys.time()

if (opt$cores > 1 && .Platform$OS.type != "windows") {
  ppc_list <- pbmclapply(seq_along(draw_idx), simulate_one_draw,
                       mc.cores = opt$cores, mc.set.seed = TRUE)
} else {
  ppc_list <- lapply(seq_along(draw_idx), simulate_one_draw)
}

ppc_df <- do.call(rbind, ppc_list)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf(" done in %.0fs (%d rows)\n", elapsed, nrow(ppc_df)))

## ---------------------------------------------------------------------------
## 5. Merge observed data
## ---------------------------------------------------------------------------
## Build a trial index per participant matching the Stan data order
obs_df <- data.frame(
  participant = stan_data$participant,
  observed_choice = stan_data$cho,
  observed_rt     = stan_data$rt
)
obs_df <- obs_df %>%
  group_by(participant) %>%
  mutate(trial_in_participant = row_number()) %>%
  ungroup()

## Add condition info
if ("con" %in% names(stan_data)) {
  obs_df$con <- stan_data$con
}

## ---------------------------------------------------------------------------
## Add trial-level behavioural variables for PPC metrics
## ---------------------------------------------------------------------------
is_ccss <- grepl("ccss", opt$family)
## CS: any family containing "_cs" that isn't CCSS (catches cpt_cs, mv_cs, cpt_cs_7o)
is_cs   <- !is_ccss && grepl("_cs", opt$family)

if (is_ccss) {
  ## MV families store evd/sdd/skew directly; CPT families need computation
  if ("evd" %in% names(stan_data)) {
    obs_df$evd  <- stan_data$evd
    obs_df$sdd  <- stan_data$sdd
    obs_df$skew <- stan_data$skew
  } else {
    ## CPT CCSS: compute EV, SD, Skew from outcomes/probabilities
    o_r <- stan_data$o_risky; p_r <- stan_data$p_risky
    o_s <- stan_data$o_safe;  p_s <- stan_data$p_safe
    ev_a <- o_r[,1] * p_r[,1] + o_r[,2] * p_r[,2]
    ev_b <- o_s[,1] * p_s[,1] + o_s[,2] * p_s[,2]
    var_a <- p_r[,1] * (o_r[,1] - ev_a)^2 + p_r[,2] * (o_r[,2] - ev_a)^2
    var_b <- p_s[,1] * (o_s[,1] - ev_b)^2 + p_s[,2] * (o_s[,2] - ev_b)^2
    obs_df$evd  <- ev_a - ev_b
    obs_df$sdd  <- sqrt(var_a) - sqrt(var_b)
    ## Skewness: (1-2p)/sqrt(p(1-p)) * sign(o_hi - o_lo)
    .skew <- function(p_hi, o_hi, o_lo) {
      p <- pmin(pmax(p_hi, 1e-6), 1 - 1e-6)
      ((1 - 2*p) / sqrt(p * (1-p))) * sign(o_hi - o_lo)
    }
    obs_df$skew <- .skew(p_r[,1], o_r[,1], o_r[,2]) -
                   .skew(p_s[,1], o_s[,1], o_s[,2])
  }
}

if (is_cs) {
  ## chose_complex: 1 if participant chose the complex option, -1 otherwise
  ## accuracy_flipped == 1 means chose simple, so chose_complex = -accuracy_flipped mapped
  if ("accuracy_flipped" %in% names(stan_data)) {
    obs_df$chose_complex_obs <- ifelse(stan_data$accuracy_flipped == 1L, -1L, 1L)
  }
  ## Moment-difference columns (complex - simple direction) for behavioural PPC metrics.
  ## Only mv_cs / mv_cs_* families populate these in stan_data; CPT CS families can compute
  ## moments from the o_complex / p_complex matrices if needed.
  if ("evd" %in% names(stan_data))  obs_df$evd  <- stan_data$evd
  if ("sdd" %in% names(stan_data))  obs_df$sdd  <- stan_data$sdd
  if ("skew_c" %in% names(stan_data) && "skew_s" %in% names(stan_data))
    obs_df$skew <- stan_data$skew_c - stan_data$skew_s
  ## CPT CS: compute evd from o_complex / p_complex matrices
  if (is.null(obs_df$evd) && "o_complex" %in% names(stan_data)) {
    oc <- stan_data$o_complex; pc <- stan_data$p_complex
    os <- stan_data$o_simple;  ps <- stan_data$p_simple
    ## Works for both 2-col and 7-col complex (cpt_cs / cpt_cs_7o).
    ev_complex <- rowSums(oc * pc)
    ev_simple  <- rowSums(os * ps)
    obs_df$evd <- as.numeric(ev_complex - ev_simple)
  }
}

ppc_full <- ppc_df %>%
  left_join(obs_df, by = c("participant", "trial_in_participant")) %>%
  left_join(id_map, by = "participant")

## ---------------------------------------------------------------------------
## 6. Save
## ---------------------------------------------------------------------------
out_path <- file.path(model_dir, "posterior_predictives.csv")
write_csv(ppc_full, out_path)
cat(sprintf("  Saved: %s (%.1f MB)\n", out_path,
            file.info(out_path)$size / 1e6))


