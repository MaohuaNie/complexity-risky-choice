## ============================================================================
## R/fit.R — compile and sample one hierarchical Stan model on real data
##
## Aim:     Core fitting engine. compile_model() caches a cmdstanr binary from
##          the registry; fit_one() preprocesses a study's data, generates
##          inits, runs HMC, and writes lightweight posterior + diagnostic
##          artefacts. Called by run_fit.R and run_all_fits.R.
## Inputs:  MODEL_REGISTRY (R/model_registry.R), PREPROCESSOR family functions
##          (R/preprocess.R), generate_inits() (R/inits.R); a study data RDS;
##          Stan files under stan/. Uses cmdstanr (binaries cached in
##          stan_cache/).
## Outputs: Into results/<study>/<family>/<model>/:
##            draws.rds         — posterior draws of mu, sigma,
##                                participant_params, Omega
##            diagnostics.csv   — Rhat / ESS / divergence summary per param
##            params_long.csv   — per-participant posterior mean + 95% CrI (tidy)
##            omega_summary.csv — correlation matrix (Omega) posterior summary
##                                (upper triangle only; diagonal = 1)
##            log_lik.rds       — per-trial log-likelihood draws (input for LOO)
##            meta.rds          — data/stan hashes, seeds, timing, Stan version
##            id_map.csv        — Stan participant index <-> subject ID mapping
## Usage:   Not run directly; source()'d by run_fit.R / run_all_fits.R.
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM fitting stage).
## Pipeline order and dependencies are documented in ../README.md.
## ============================================================================
##
## Note: the full CmdStanMCMC envelope (fit.rds) is intentionally NOT saved.
## It was several GB per model and its serialisation caused OOM kills during
## post-sampling. All downstream scripts consume draws.rds + log_lik.rds + the
## CSVs, so only the per-parameter draws are kept.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

.timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")

.hash_sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (!requireNamespace("digest", quietly = TRUE))
    return(as.character(file.info(path)$mtime))
  digest::digest(file = path, algo = "sha256")
}

## ---------------------------------------------------------------------------
## Compile (cached) a model from the registry.
## ---------------------------------------------------------------------------
compile_model <- function(family, model, cache_dir = "stan_cache",
                          stanc_options = list(), force_recompile = FALSE) {
  row <- MODEL_REGISTRY |> filter(family == !!family, model == !!model)
  if (nrow(row) != 1) stop("No registry row for ", family, "/", model)
  stan_file <- row$stan_file
  if (!file.exists(stan_file)) stop("Stan file missing: ", stan_file)
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  cmdstan_model(stan_file = stan_file, dir = cache_dir,
                compile = TRUE, force_recompile = force_recompile,
                stanc_options = stanc_options)
}

## ---------------------------------------------------------------------------
## Fit one (study, family, model) combination.
## ---------------------------------------------------------------------------
fit_one <- function(study, family, model, data_rds, out_dir,
                    chains = 4, parallel_chains = 4,
                    iter_warmup = 2000, iter_sampling = 2000,
                    adapt_delta = 0.95, max_treedepth = 12,
                    seed = 2026, refresh = 200, save_warmup = FALSE,
                    overwrite = FALSE) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ## Skip marker: draws.rds is the canonical "fit is done" artefact. A legacy
  ## fit.rds (from older runs) also counts as done so we don't re-fit a
  ## completed model.
  draws_path <- file.path(out_dir, "draws.rds")
  legacy_fit <- file.path(out_dir, "fit.rds")
  if ((file.exists(draws_path) || file.exists(legacy_fit)) && !overwrite) {
    message("[SKIP] fit already exists in ", out_dir,
            " — use overwrite=TRUE to re-fit")
    return(invisible(NULL))
  }

  ## registry lookup
  row <- MODEL_REGISTRY |> filter(family == !!family, model == !!model)
  if (nrow(row) != 1) stop("No registry row for ", family, "/", model)
  n_params <- row$n_params
  stan_file <- row$stan_file

  ## preprocess
  pp <- PREPROCESSOR[[family]]
  if (is.null(pp)) stop("No preprocessor for family ", family)
  pd <- pp(data_rds)
  stan_data <- pd$stan_data
  write_csv(pd$id_map, file.path(out_dir, "id_map.csv"))

  ## compile + init
  mod <- compile_model(family, model)
  min_rt <- suppressWarnings(min(stan_data$rt, na.rm = TRUE))
  inits <- generate_inits(family, n_params, L = stan_data$L,
                          chains = chains, min_rt = min_rt)

  ## sample
  message(sprintf("[%s] Fitting %s/%s/%s  N=%d  L=%d  (%d chains x %d iter)",
                  .timestamp(), study, family, model,
                  stan_data$N, stan_data$L,
                  chains, iter_warmup + iter_sampling))
  t0 <- Sys.time()
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
    save_warmup     = save_warmup,
    refresh         = refresh,
    show_messages   = TRUE
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ## Save lightweight posterior draws (mu, sigma, participant_params, Omega) —
  ## what every downstream script consumes. We deliberately do NOT call
  ## fit$save_object(): serialising the full CmdStanMCMC envelope is several GB
  ## per model and caused OOM kills during post-sampling. Keeping only the
  ## per-parameter draws cuts peak memory and per-model disk usage sharply.
  saveRDS(fit$draws(variables = c("mu", "sigma", "participant_params", "Omega")),
          draws_path)
  gc(verbose = FALSE)

  ## Correlation-matrix (Omega) posterior summary: mean + 95% CrI per entry.
  ## Upper triangle only; diagonal = 1 by construction.
  .write_omega_summary(fit, out_dir, study, family, model)
  gc(verbose = FALSE)

  ## diagnostics
  diag <- .diagnose(fit)
  write_csv(diag$summary, file.path(out_dir, "diagnostics.csv"))

  ## participant-level posterior summary
  params_long <- .summarize_participants(fit, n_params, pd$id_map, family, model)
  write_csv(params_long, file.path(out_dir, "params_long.csv"))
  gc(verbose = FALSE)

  ## save log_lik for later LOO (post/run_loo.R consumes this)
  log_lik_path <- file.path(out_dir, "log_lik.rds")
  saveRDS(fit$draws(variables = "log_lik", format = "draws_matrix"),
          log_lik_path)
  gc(verbose = FALSE)

  ## metadata
  meta <- list(
    study           = study,
    family          = family,
    model           = model,
    n_params        = n_params,
    data_rds        = normalizePath(data_rds, mustWork = TRUE),
    data_hash       = .hash_sha256(data_rds),
    stan_file       = normalizePath(stan_file, mustWork = TRUE),
    stan_hash       = .hash_sha256(stan_file),
    cmdstan_version = cmdstan_version(),
    n_chains        = chains,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth,
    seed            = seed,
    N               = stan_data$N,
    L               = stan_data$L,
    timestamp       = .timestamp(),
    elapsed_sec     = elapsed,
    divergences     = diag$n_divergent,
    max_rhat        = diag$max_rhat,
    min_ess         = diag$min_ess
  )
  saveRDS(meta, file.path(out_dir, "meta.rds"))

  message(sprintf("[%s] Done %s/%s/%s in %.0fs  div=%d  max_rhat=%.3f",
                  .timestamp(), study, family, model, elapsed,
                  diag$n_divergent, diag$max_rhat))

  invisible(meta)
}

## ---------------------------------------------------------------------------
## Correlation-matrix (Omega) posterior summary.
## Writes omega_summary.csv with one row per upper-triangle entry Omega[i,j]
## (i <= j), columns: study, family, model, i, j, mean, sd, q2.5, q97.5.
## The diagonal is always 1 (symmetric correlation matrix), included for
## completeness. Matches the format produced by post/extract_omega.R so
## downstream aggregators can consume both interchangeably.
## ---------------------------------------------------------------------------
.write_omega_summary <- function(fit, out_dir, study, family, model) {
  omega_draws <- tryCatch(
    fit$draws("Omega", format = "draws_matrix"),
    error = function(e) { message("  Omega extraction failed: ",
                                   conditionMessage(e)); NULL }
  )
  if (is.null(omega_draws)) return(invisible(NULL))

  summ <- posterior::summarise_draws(
    omega_draws, mean, sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975))
  )
  names(summ)[names(summ) %in% c("2.5%",  "q2.5")]  <- "q2.5"
  names(summ)[names(summ) %in% c("97.5%", "q97.5")] <- "q97.5"

  out <- summ |>
    tidyr::extract(variable, into = c("i", "j"),
                   regex = "Omega\\[(\\d+),(\\d+)\\]",
                   convert = TRUE, remove = FALSE) |>
    dplyr::filter(i <= j) |>
    dplyr::mutate(study = study, family = family, model = model) |>
    dplyr::select(study, family, model, i, j, mean, sd, q2.5, q97.5)

  write_csv(out, file.path(out_dir, "omega_summary.csv"))
  invisible(out)
}

## ---------------------------------------------------------------------------
## Diagnostics extractor — returns convergence summary and top-line numbers.
## ---------------------------------------------------------------------------
.diagnose <- function(fit) {
  summ <- fit$summary(variables = c("mu", "sigma", "participant_params"),
                      .num_args = list(sigfig = 4))
  n_divergent <- sum(fit$diagnostic_summary()$num_divergent)
  max_rhat <- suppressWarnings(max(summ$rhat, na.rm = TRUE))
  min_ess  <- suppressWarnings(min(c(summ$ess_bulk, summ$ess_tail), na.rm = TRUE))
  list(summary = summ,
       n_divergent = n_divergent,
       max_rhat = max_rhat,
       min_ess = min_ess)
}

## ---------------------------------------------------------------------------
## Map (family, model) → canonical parameter name at each raw-vector position.
## Replicates .active_full_names() in R/param_names.R without requiring a
## source() — keeps fit.R self-contained.
## ---------------------------------------------------------------------------
.param_names_for <- function(family, model) {
  is_mv    <- grepl("^mv_",  family)
  is_ccss  <- grepl("ccss",  family)

  ## Base names (positions 1-5) — depends on family class, not on the
  ## _7o suffix (7o variants share parameterisation with their base family).
  base <- if (is_mv && is_ccss) {
    c("beta", "theta_raw", "threshold_raw", "ndt_raw", "eta")
  } else if (is_mv && !is_ccss) {
    c("beta", "theta_raw", "threshold_raw", "ndt_raw", "eta")
  } else if (!is_mv && is_ccss) {
    c("beta_raw", "theta_raw", "threshold_raw", "ndt_raw", "gamma_raw")
  } else {
    c("beta_raw", "theta_raw", "threshold_raw", "ndt_raw", "gamma_raw")
  }

  if (model == "baseline") return(base)

  if (is_ccss) {
    letters_present <- strsplit(model, "_")[[1]]
    canonical <- c("r", "n", "a", "s")
    active    <- canonical[canonical %in% letters_present]
    dm <- if (is_mv)
      c(r = "delta_beta", n = "delta_theta", a = "delta_threshold", s = "delta_eta")
    else
      c(r = "delta_beta", n = "delta_theta", a = "delta_threshold", s = "delta_gamma")
    extras <- unname(dm[active])
  } else {
    words_present <- strsplit(model, "_")[[1]]
    canonical <- c("sp", "dr", "skew")
    active    <- canonical[canonical %in% words_present]
    em <- if (is_mv)
      c(sp = "sp_raw", dr = "zeta", skew = "delta_eta")
    else
      c(sp = "sp_raw", dr = "zeta", skew = "delta_gamma")
    extras <- unname(em[active])
  }
  c(base, extras)
}

## ---------------------------------------------------------------------------
## Per-participant posterior summary, tidy for downstream analysis.
## Output columns: family, model, subject, participant, param, param_name,
##                 param_idx, mean, sd, q2.5, q97.5
## ---------------------------------------------------------------------------
.summarize_participants <- function(fit, n_params, id_map, family, model) {
  ## Stan's `participant_params` is [n_params, L]. We extract, rename by
  ## position (1..n_params), attach canonical param names, and join subject IDs.
  pp <- fit$draws("participant_params", format = "draws_matrix")
  summ <- posterior::summarise_draws(
    pp, mean, sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975))
  )
  names(summ)[names(summ) %in% c("2.5%",  "q2.5")]  <- "q2.5"
  names(summ)[names(summ) %in% c("97.5%", "q97.5")] <- "q97.5"

  parse_one <- function(v) {
    m <- regmatches(v, regexec("participant_params\\[(\\d+),(\\d+)\\]", v))[[1]]
    if (length(m) < 3) return(data.frame(param_idx = NA, participant = NA))
    data.frame(param_idx = as.integer(m[2]), participant = as.integer(m[3]))
  }
  idx <- do.call(rbind, lapply(summ$variable, parse_one))

  ## Readable param names, e.g. "beta_raw", "delta_theta".
  param_names <- .param_names_for(family, model)

  out <- cbind(idx, summ[, c("mean", "sd", "q2.5", "q97.5")]) |>
    as_tibble() |>
    filter(!is.na(participant)) |>
    mutate(family = family, model = model,
           param      = sprintf("p%02d", param_idx),
           param_name = param_names[param_idx]) |>
    left_join(id_map, by = "participant") |>
    select(family, model, subject, participant,
           param, param_name, param_idx,
           mean, sd, q2.5, q97.5)
  out
}
