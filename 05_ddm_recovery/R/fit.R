## ============================================================================
## R/fit.R — compile Stan models and fit one recovery dataset (cmdstanr)
##
## Aim:     Compile (once, cached) the Stan model matching a model key, build
##          prior-centred inits, and refit a simulated dataset with the same
##          hierarchical DDM model used in the fitting stage (folder 03).
## Inputs:  a model key, a Stan-data list (from build_stan_data), sampling
##          controls; Stan sources in stan/; PRIORS from priors.R.
## Outputs: in-memory list (draws, summary, diagnostics); compiled binaries
##          cached in stan_cache/. run_recovery.R persists draws as
##          fit_dataset_{d}.rds.
## Usage:   source("R/fit.R"); compile_model(key); fit_recovery(...)
##          (called from run_recovery.R).
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM parameter
## recovery). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================
##
## cmdstanr is used instead of rstan for faster compilation (binaries cached on
## disk in stan_cache/), cleaner per-chain parallelism, and a draws() interface
## that returns the posterior directly. All Stan sources live in this folder's
## stan/ so the recovery stage is self-contained.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

## Map model_key → path to the Stan source file.
## All four Stan files live inside recovery/stan/ so the folder is
## self-contained (copy the whole `recovery/` directory to the cluster).
STAN_PATHS <- list(
  cpt_ccss_n_r_a_s  = "stan/cpt_ccss_n_r_a_s.stan",
  mv_ccss_n_r_a_s   = "stan/mv_ccss_n_r_a_s.stan",
  cpt_cs_sp_dr_skew = "stan/cpt_cs_sp_dr_skew.stan",
  mv_cs_sp_dr_skew  = "stan/mv_cs_sp_dr_skew.stan",
  cpt_cs_sp_dr      = "stan/cpt_cs_sp_dr.stan",
  mv_cs_sp_dr       = "stan/mv_cs_sp_dr.stan"
)

## ---------------------------------------------------------------------------
## Compile (or fetch from cache) the Stan model for a given model_key.
## cmdstanr auto-caches by hash; we also pin the output dir so files live
## inside recovery/stan_cache.
## ---------------------------------------------------------------------------
compile_model <- function(model_key, cache_dir = "stan_cache") {
  stan_file <- STAN_PATHS[[model_key]]
  if (is.null(stan_file) || !file.exists(stan_file)) {
    stop("Stan file for '", model_key, "' not found at: ", stan_file)
  }
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  cmdstan_model(
    stan_file       = stan_file,
    dir             = cache_dir,
    compile         = TRUE,
    force_recompile = FALSE,
    cpp_options     = list(stan_threads = FALSE)
  )
}

## ---------------------------------------------------------------------------
## Initialisation helpers. For each model we need inits for mu, sigma, L_corr,
## z matching the Stan parameter block. We use mildly perturbed prior-centred
## starts to keep warmup short.
## ---------------------------------------------------------------------------
make_inits <- function(model_key, L, chains = 4, jitter = 0.3,
                       min_rt = NULL) {
  P <- PRIORS[[model_key]]
  K <- length(P$params)

  ## Derive a raw-scale prior centre from the uniform-on-real schema:
  ## midpoint of (real_lo, real_hi) → invert to raw via the declared transform.
  real_mid  <- (P$real_lo + P$real_hi) / 2
  mu_centre <- vapply(seq_len(K),
                      function(k) .inv_transform(real_mid[k], P$transform[k]),
                      numeric(1))

  ## Find the index of ndt_raw (if present in this model) and pin its init
  ## to softplus(ndt_raw) < 0.5 * min_rt, so the wiener density's
  ## rt > ndt constraint is satisfied from the very first log-prob eval.
  ndt_idx <- which(P$params == "ndt_raw")

  lapply(seq_len(chains), function(ch) {
    mu_init <- mu_centre + rnorm(K, 0, jitter)
    sigma_init <- pmax(runif(K, P$sigma_lo, P$sigma_hi) + rnorm(K, 0, 0.05),
                       0.05)
    ## Cap ndt init so softplus(mu_ndt) < 0.5 * min_rt (plenty of headroom
    ## even with z-based per-subject offsets). softplus^-1(x) = log(exp(x)-1).
    if (length(ndt_idx) == 1 && !is.null(min_rt) && is.finite(min_rt)) {
      target_ndt <- max(0.02, 0.3 * min_rt)            # 30% of fastest RT
      cap <- log(exp(target_ndt) - 1)                  # inverse softplus
      mu_init[ndt_idx]    <- min(mu_init[ndt_idx], cap)
      sigma_init[ndt_idx] <- min(sigma_init[ndt_idx], 0.15)
    }
    list(
      mu     = mu_init,
      sigma  = sigma_init,
      L_corr = diag(K),
      z      = matrix(rnorm(K * L, 0, 1), nrow = K, ncol = L)
    )
  })
}

## ---------------------------------------------------------------------------
## Fit one recovery dataset.
##
## Returns a list:
##   draws      — posterior draws (rvar format via posterior::as_draws_rvars)
##   summary    — convergence summary (rhat, ess)
##   diagnostics — divergences/treedepth counts
## ---------------------------------------------------------------------------
fit_recovery <- function(mod, stan_data, model_key,
                         chains = 4, parallel_chains = 4,
                         iter_warmup = 1000, iter_sampling = 1000,
                         adapt_delta = 0.9, max_treedepth = 12,
                         seed = 1, refresh = 0, show_messages = FALSE) {

  inits <- make_inits(model_key, L = stan_data$L, chains = chains,
                      min_rt = suppressWarnings(min(stan_data$rt, na.rm = TRUE)))

  fit <- mod$sample(
    data            = stan_data,
    seed            = seed,
    chains          = chains,
    parallel_chains = parallel_chains,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth,
    init            = inits,
    refresh         = refresh,
    show_messages   = show_messages
  )

  ## drop heavy arrays we don't need for recovery (log_lik, per-trial drifts)
  keep_vars <- c("mu", "sigma", "participant_params")
  draws <- fit$draws(variables = keep_vars)

  list(
    fit         = fit,              # keep so user can inspect if needed
    draws       = draws,
    summary     = fit$summary(variables = keep_vars),
    diagnostics = fit$diagnostic_summary()
  )
}
