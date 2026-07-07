## ============================================================================
## R/simulate.R — DDM data simulation + Stan-data builder for one dataset
##
## Aim:     Generate one synthetic recovery dataset: draw participant raw
##          parameters from the hierarchical prior, map them to per-trial DDM
##          parameters (transforms.R), simulate choices/RTs with rtdists, apply
##          rejection sampling to discard pathological draws, and package the
##          result into the list a Stan model expects.
## Inputs:  a model key/family, a trial template (trials.R), L (participants),
##          a seed; relies on EFFECTIVE (transforms.R) and
##          sample_hierarchical_params (priors.R).
## Outputs: in-memory list (simulated trials + true raw params + hyperparams);
##          run_recovery.R persists it as sim_dataset_{d}.rds.
## Usage:   source("R/simulate.R"); simulate_dataset(...) / build_stan_data(...)
##          (called from run_recovery.R).
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM parameter
## recovery). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================
##
## rtdists::rdiffusion is C-backed and called per trial (n = 1) with the
## participant's per-trial drift / threshold / ndt / rel_sp. A "dataset" = one
## synthetic experiment of L participants sharing the same trial template;
## participant-level raw params are drawn from the hierarchical prior
## (priors.R).

suppressPackageStartupMessages({
  library(rtdists)
  library(dplyr)
})

## priors.R and transforms.R must be source()d by the top-level script
## before simulate.R is used. We don't source them here to avoid double-
## loading and path-dependence.

## ---------------------------------------------------------------------------
## Simulate choice + rt for one participant given a vector of effective
## (drift, threshold, ndt, rel_sp) per trial.
##
## rtdists::rdiffusion takes scalar a, v, t0, z; vectorise over trials by
## calling it in a loop but with n=1 per call. For typical trial counts
## (~90) this is <5ms per participant.
##
## Returns a data.frame with rt, cho (±1) of length length(drift).
## ---------------------------------------------------------------------------
simulate_participant <- function(eff, max_rt = 30) {
  ## unname to avoid rdiffusion complaining about named numerics
  drift     <- as.numeric(unname(eff$drift))
  threshold <- as.numeric(unname(eff$threshold))
  ndt       <- as.numeric(unname(eff$ndt))
  rel_sp    <- as.numeric(unname(eff$rel_sp))

  n <- length(drift)
  rt <- rep(NA_real_, n)
  cho <- rep(NA_integer_, n)
  for (i in seq_len(n)) {
    ## Skip a trial whose params are non-finite rather than crashing.
    ## Result: the trial is flagged `bad` and the whole participant is
    ## rejected above (>10% bad trials).
    if (!is.finite(drift[i]) || !is.finite(threshold[i]) || threshold[i] <= 0 ||
        !is.finite(ndt[i])   || ndt[i] <= 0 ||
        !is.finite(rel_sp[i])|| rel_sp[i] <= 0 || rel_sp[i] >= 1) {
      next
    }
    out <- tryCatch(
      rtdists::rdiffusion(
        n  = 1,
        a  = threshold[i],
        v  = drift[i],
        t0 = ndt[i],
        z  = rel_sp[i] * threshold[i]
      ),
      error = function(e) NULL
    )
    if (is.null(out) || !is.finite(out$rt)) next
    rt[i]  <- out$rt
    cho[i] <- if (out$response == "upper") 1L else -1L
  }
  bad <- is.na(rt) | rt <= 0 | rt > max_rt | !is.finite(rt)
  list(rt = rt, cho = cho, bad = bad)
}

## ---------------------------------------------------------------------------
## Generate one recovery dataset for a model.
##
## Arguments:
##   model_key  — one of names(EFFECTIVE)
##   family     — "cpt_ccss" | "mv_ccss" | "cpt_cs" | "mv_cs"
##   trials     — trial template (output of load_trials_template)
##   L          — number of synthetic participants
##   seed       — RNG seed
##   max_reject — max attempts per participant before giving up
##
## Returns a list:
##   trials:          long-format data.frame (participant, trial, covariates,
##                    rt, cho, ...ready for stan_data builder)
##   true_params:     data.frame [L x K] of RAW params per participant
##   mu_true, sigma_true, raw_matrix:  raw-scale hyperparameters
## ---------------------------------------------------------------------------
simulate_dataset <- function(model_key, family, trials, L = 10,
                             seed = 1, max_reject = 20, verbose = TRUE) {
  set.seed(seed)
  sampler <- sample_hierarchical_params(model_key, L = L, rng = seed)
  eff_fn  <- EFFECTIVE[[model_key]]
  if (is.null(eff_fn)) stop("Unknown model_key: ", model_key)

  sim_list <- vector("list", L)
  raw_accepted <- matrix(NA_real_, nrow = L, ncol = ncol(sampler$raw))
  colnames(raw_accepted) <- colnames(sampler$raw)

  ## tally of why rejections happen (diagnostic)
  reject_reasons <- integer(7)
  names(reject_reasons) <- c("nonfinite", "threshold_le0", "ndt_le0",
                             "rel_sp_bad", "drift_too_large",
                             "too_many_bad_trials", "stuck_responder")

  for (l in seq_len(L)) {
    attempts <- 0
    accepted <- FALSE
    while (attempts < max_reject) {
      attempts <- attempts + 1
      if (attempts == 1) {
        raw_l <- sampler$raw[l, ]
      } else {
        extra <- sample_hierarchical_params(model_key, L = 1,
                                            rng = seed * 1000 + l * 100 + attempts)
        raw_l <- extra$raw[1, ]
      }
      eff <- eff_fn(raw_l, trials)
      all_eff <- c(as.numeric(eff$drift), as.numeric(eff$threshold),
                   as.numeric(eff$ndt),   as.numeric(eff$rel_sp))
      if (any(!is.finite(all_eff))) { reject_reasons["nonfinite"] <- reject_reasons["nonfinite"] + 1; next }
      if (any(as.numeric(eff$threshold) <= 0)) { reject_reasons["threshold_le0"] <- reject_reasons["threshold_le0"] + 1; next }
      if (any(as.numeric(eff$ndt)       <= 0)) { reject_reasons["ndt_le0"] <- reject_reasons["ndt_le0"] + 1; next }
      if (any(as.numeric(eff$rel_sp) <= 0 | as.numeric(eff$rel_sp) >= 1)) { reject_reasons["rel_sp_bad"] <- reject_reasons["rel_sp_bad"] + 1; next }
      if (max(abs(as.numeric(eff$drift))) > 50) { reject_reasons["drift_too_large"] <- reject_reasons["drift_too_large"] + 1; next }

      sim <- simulate_participant(eff)
      if (mean(sim$bad, na.rm = TRUE) > 0.1) { reject_reasons["too_many_bad_trials"] <- reject_reasons["too_many_bad_trials"] + 1; next }
      ch <- sim$cho[!is.na(sim$cho)]
      if (length(ch) < 0.8 * length(sim$cho)) { reject_reasons["too_many_bad_trials"] <- reject_reasons["too_many_bad_trials"] + 1; next }
      if (mean(ch == 1) < 0.02 || mean(ch == 1) > 0.98) { reject_reasons["stuck_responder"] <- reject_reasons["stuck_responder"] + 1; next }
      accepted <- TRUE
      raw_accepted[l, ] <- raw_l
      sim_list[[l]] <- tibble::tibble(
        participant = l,
        trial       = trials$trial_idx,
        rt          = sim$rt,
        cho         = sim$cho
      )
      break
    }
    if (!accepted) {
      if (verbose) {
        cat("  Rejection tally after ", attempts, " attempts (participant ", l, "):\n", sep = "")
        print(reject_reasons)
        ## diagnostic snapshot of last eff
        cat("  Last draw summary:\n")
        cat(sprintf("    drift:     min=%.3g max=%.3g finite=%d/%d\n",
                    suppressWarnings(min(eff$drift, na.rm = TRUE)),
                    suppressWarnings(max(eff$drift, na.rm = TRUE)),
                    sum(is.finite(eff$drift)), length(eff$drift)))
        cat(sprintf("    threshold: min=%.3g max=%.3g\n",
                    suppressWarnings(min(eff$threshold, na.rm = TRUE)),
                    suppressWarnings(max(eff$threshold, na.rm = TRUE))))
        cat(sprintf("    ndt:       min=%.3g max=%.3g\n",
                    suppressWarnings(min(eff$ndt, na.rm = TRUE)),
                    suppressWarnings(max(eff$ndt, na.rm = TRUE))))
        cat(sprintf("    rel_sp:    min=%.3g max=%.3g\n",
                    suppressWarnings(min(eff$rel_sp, na.rm = TRUE)),
                    suppressWarnings(max(eff$rel_sp, na.rm = TRUE))))
      }
      stop(sprintf("Participant %d: could not find acceptable params after %d attempts",
                   l, max_reject))
    }
  }

  if (verbose && sum(reject_reasons) > 0) {
    cat("  Rejection tally across all participants:\n")
    print(reject_reasons[reject_reasons > 0])
  }

  ## merge trial covariates into simulated long table
  sim_df <- dplyr::bind_rows(sim_list) %>%
    dplyr::left_join(trials, by = c("trial" = "trial_idx"),
                     suffix = c("", ".tpl"))

  list(
    trials       = sim_df,
    true_params  = as.data.frame(raw_accepted),
    mu_true      = sampler$mu_true,
    sigma_true   = sampler$sigma_true
  )
}

## ---------------------------------------------------------------------------
## Build a Stan data list from a simulated long dataset.
## ---------------------------------------------------------------------------
build_stan_data <- function(sim, family, model_key = NULL) {
  df <- sim$trials
  base <- list(
    N = nrow(df),
    L = max(df$participant),
    participant = as.integer(df$participant),
    cho = as.integer(df$cho),
    rt  = as.numeric(df$rt)
  )

  if (family == "cpt_ccss") {
    c(base, list(
      o_risky = cbind(df$o_risky1, df$o_risky2),
      o_safe  = cbind(df$o_safe1,  df$o_safe2),
      p_risky = cbind(df$p_risky1, df$p_risky2),
      p_safe  = cbind(df$p_safe1,  df$p_safe2),
      con = as.integer(df$con),
      starting_point = 0.5
    ))
  } else if (family == "mv_ccss") {
    c(base, list(
      evd = as.numeric(df$evd),
      sdd = as.numeric(df$sdd),
      skew = as.numeric(df$skew),
      con = as.integer(df$con),
      starting_point = 0.5
    ))
  } else if (family == "cpt_cs") {
    ## accuracy_flipped: 1 if choice was "lower" (cho == -1), else 0  (matches Stan)
    c(base, list(
      accuracy_flipped = as.integer(df$cho == -1L),
      o_complex = cbind(df$o_complex1, df$o_complex2),
      o_simple  = cbind(df$o_simple1,  df$o_simple2),
      p_complex = cbind(df$p_complex1, df$p_complex2),
      p_simple  = cbind(df$p_simple1,  df$p_simple2)
    ))
  } else if (family == "mv_cs") {
    ## mv_cs_sp_dr.stan takes a single `skew = skew_c - skew_s`;
    ## mv_cs_sp_dr_skew.stan takes `skew_c` and `skew_s` separately.
    mv_cs_extra <- if (!is.null(model_key) && model_key == "mv_cs_sp_dr") {
      list(skew = as.numeric(df$skew_c) - as.numeric(df$skew_s))
    } else {
      list(skew_c = as.numeric(df$skew_c),
           skew_s = as.numeric(df$skew_s))
    }
    c(base, list(
      accuracy_flipped = as.integer(df$cho == -1L),
      evd = as.numeric(df$evd),
      sdd = as.numeric(df$sdd)
    ), mv_cs_extra)
  } else {
    stop("Unknown family: ", family)
  }
}
