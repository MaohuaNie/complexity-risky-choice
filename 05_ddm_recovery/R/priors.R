## ============================================================================
## R/priors.R — generative hierarchical priors for the recovery simulation
##
## Aim:     Define, per model, the ranges from which true parameters are drawn
##          and the sampler that produces raw-scale (Stan-space) hyper-means,
##          group SDs, and per-participant parameters for one dataset.
## Inputs:  none (self-contained PRIORS registry + sampler).
## Outputs: PRIORS list and sample_hierarchical_params(); consumed by
##          simulate.R and fit.R (make_inits). No files written.
## Usage:   source("R/priors.R") from run_recovery.R (do not run standalone).
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM parameter
## recovery). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================
##
## Design:
##   * GROUP MEAN  mu_real ~ Uniform(real_lo, real_hi)  on the REAL scale.
##     This sweeps the full psychologically-meaningful range across datasets
##     rather than clustering every dataset near one hyper-mean point.
##   * mu_true (raw scale, what Stan samples) = inverse_transform(mu_real).
##   * PARTICIPANT SD  sigma_raw ~ Uniform(sigma_lo, sigma_hi) on the RAW scale.
##     Participants within a dataset are Normal(mu_true, sigma_raw) on raw.
##
## Transforms (must match transforms.R / each Stan model's TP block):
##   "softplus" : real = log(1 + exp(raw)),     real > 0
##   "phi"      : real = Phi(raw),              real in (0, 1)
##   "identity" : real = raw                    (unconstrained)
##
## Target REAL-scale ranges encoded below:
##   threshold (alpha)          1.0 – 5.0
##   ndt (tau)                  0.2 – 0.6 s
##   theta (drift scaler)       0.1 – 0.4
##   beta, gamma (CPT)          0.4 – 1.2   (utility / prob-weighting curv.)
##   beta (MV)                 -0.8 – 0.8   (unconstrained weight on SD)
##   eta  (MV)                 -3   – 3     (unconstrained weight on skew)
##   rel_sp                     0.3 – 0.7   (Phi(sp_raw))
##   zeta (additive offset)   -15   – 15
##   delta_* (context shifts)  -1   – 1     (unconstrained)
## ---------------------------------------------------------------------------

PRIORS <- list()

## ---- CPT CCSS n_r_a_s (9 params) ------------------------------------------
PRIORS[["cpt_ccss_n_r_a_s"]] <- list(
  params    = c("beta_raw", "theta_raw", "threshold_raw", "ndt_raw", "gamma_raw",
                "delta_beta", "delta_theta", "delta_threshold", "delta_gamma"),
  transform = c("softplus", "softplus",  "softplus",      "softplus","softplus",
                "identity", "identity",  "identity",      "identity"),
  real_lo   = c(0.4, 0.1, 1.0, 0.2, 0.4,  -1, -1, -1, -1),
  real_hi   = c(1.2, 0.4, 5.0, 0.6, 1.2,   1,  1,  1,  1),
  sigma_lo  = c(0.15, 0.15, 0.20, 0.15, 0.15, 0.20, 0.20, 0.20, 0.20),
  sigma_hi  = c(0.35, 0.35, 0.40, 0.30, 0.35, 0.40, 0.40, 0.40, 0.40)
)

## ---- MV CCSS n_r_a_s (9 params) -------------------------------------------
##   beta, eta are unconstrained (identity transform).
PRIORS[["mv_ccss_n_r_a_s"]] <- list(
  params    = c("beta", "theta_raw", "threshold_raw", "ndt_raw", "eta",
                "delta_beta", "delta_theta", "delta_threshold", "delta_eta"),
  transform = c("identity", "softplus", "softplus", "softplus", "identity",
                "identity", "identity", "identity", "identity"),
  real_lo   = c(-0.8, 0.1, 1.0, 0.2, -3,  -1, -1, -1, -1),
  real_hi   = c( 0.8, 0.4, 5.0, 0.6,  3,   1,  1,  1,  1),
  sigma_lo  = c(0.15, 0.15, 0.20, 0.15, 0.50, 0.20, 0.20, 0.20, 0.20),
  sigma_hi  = c(0.35, 0.35, 0.40, 0.30, 1.00, 0.40, 0.40, 0.40, 0.40)
)

## ---- CPT CS sp_dr_skew (8 params) -----------------------------------------
##   sp_raw -> Phi() -> rel_sp in (0, 1); zeta unconstrained.
PRIORS[["cpt_cs_sp_dr_skew"]] <- list(
  params    = c("beta_raw", "theta_raw", "threshold_raw", "ndt_raw", "gamma_raw",
                "sp_raw", "zeta", "delta_gamma"),
  transform = c("softplus", "softplus", "softplus", "softplus", "softplus",
                "phi",      "identity", "identity"),
  real_lo   = c(0.4, 0.1, 1.0, 0.2, 0.4, 0.3, -15, -1),
  real_hi   = c(1.2, 0.4, 5.0, 0.6, 1.2, 0.7,  15,  1),
  sigma_lo  = c(0.15, 0.15, 0.20, 0.15, 0.15, 0.15, 4.0, 0.20),
  sigma_hi  = c(0.35, 0.35, 0.40, 0.30, 0.35, 0.30, 6.0, 0.40)
)

## ---- MV CS sp_dr_skew (8 params) ------------------------------------------
PRIORS[["mv_cs_sp_dr_skew"]] <- list(
  params    = c("beta", "theta_raw", "threshold_raw", "ndt_raw", "eta",
                "sp_raw", "zeta", "delta_eta"),
  transform = c("identity", "softplus", "softplus", "softplus", "identity",
                "phi",      "identity", "identity"),
  real_lo   = c(-0.8, 0.1, 1.0, 0.2, -3, 0.3, -15, -1),
  real_hi   = c( 0.8, 0.4, 5.0, 0.6,  3, 0.7,  15,  1),
  sigma_lo  = c(0.15, 0.15, 0.20, 0.15, 0.50, 0.15, 4.0, 0.20),
  sigma_hi  = c(0.35, 0.35, 0.40, 0.30, 1.00, 0.30, 6.0, 0.40)
)

## ---- CPT CS sp_dr (7 params; no delta_gamma) ------------------------------
PRIORS[["cpt_cs_sp_dr"]] <- list(
  params    = c("beta_raw", "theta_raw", "threshold_raw", "ndt_raw", "gamma_raw",
                "sp_raw", "zeta"),
  transform = c("softplus", "softplus", "softplus", "softplus", "softplus",
                "phi",      "identity"),
  real_lo   = c(0.4, 0.1, 1.0, 0.2, 0.4, 0.3, -15),
  real_hi   = c(1.2, 0.4, 5.0, 0.6, 1.2, 0.7,  15),
  sigma_lo  = c(0.15, 0.15, 0.20, 0.15, 0.15, 0.15, 4.0),
  sigma_hi  = c(0.35, 0.35, 0.40, 0.30, 0.35, 0.30, 6.0)
)

## ---- MV CS sp_dr (7 params; no delta_eta) ---------------------------------
PRIORS[["mv_cs_sp_dr"]] <- list(
  params    = c("beta", "theta_raw", "threshold_raw", "ndt_raw", "eta",
                "sp_raw", "zeta"),
  transform = c("identity", "softplus", "softplus", "softplus", "identity",
                "phi",      "identity"),
  real_lo   = c(-0.8, 0.1, 1.0, 0.2, -3, 0.3, -15),
  real_hi   = c( 0.8, 0.4, 5.0, 0.6,  3, 0.7,  15),
  sigma_lo  = c(0.15, 0.15, 0.20, 0.15, 0.50, 0.15, 4.0),
  sigma_hi  = c(0.35, 0.35, 0.40, 0.30, 1.00, 0.30, 6.0)
)

## ---------------------------------------------------------------------------
## Inverse transforms: real-scale mu -> raw-scale mu.
## softplus_inv(y) = log(exp(y) - 1); numerically stable at large y.
## phi_inv(p)      = qnorm(p); p must be in (0, 1).
## ---------------------------------------------------------------------------
.softplus_inv <- function(y) {
  ## For large y, log(exp(y) - 1) ~ y; use log(expm1(y)) for small y.
  ifelse(y > 20, y, log(expm1(y)))
}
.inv_transform <- function(real_value, type) {
  switch(type,
    softplus = .softplus_inv(real_value),
    phi      = qnorm(real_value),
    identity = real_value,
    stop("Unknown transform: ", type)
  )
}

## ---------------------------------------------------------------------------
## Sample raw-scale participant parameters from the hierarchical prior.
##
## Pipeline:
##   1. mu_real[k] ~ Uniform(real_lo[k], real_hi[k])
##   2. mu_true[k] = inverse_transform(mu_real[k], transform[k])   # raw scale
##   3. sigma_true[k] ~ Uniform(sigma_lo[k], sigma_hi[k])           # raw scale
##   4. raw[l, k]  ~ Normal(mu_true[k], sigma_true[k])
##
## Returns a list with:
##   mu_true     [K]    — raw-scale group means (what Stan's `mu` samples)
##   sigma_true  [K]    — raw-scale group SDs
##   raw         [L, K] — per-participant raw parameters
## ---------------------------------------------------------------------------
sample_hierarchical_params <- function(model_key, L, rng = NULL) {
  if (!is.null(rng)) set.seed(rng)
  P <- PRIORS[[model_key]]
  if (is.null(P)) stop("Unknown model_key: ", model_key)
  K <- length(P$params)

  ## sanity: all per-param vectors must be length K
  stopifnot(length(P$transform) == K, length(P$real_lo) == K,
            length(P$real_hi)   == K, length(P$sigma_lo) == K,
            length(P$sigma_hi)  == K)

  mu_real <- runif(K, P$real_lo, P$real_hi)
  mu_true <- vapply(seq_len(K),
                    function(k) .inv_transform(mu_real[k], P$transform[k]),
                    numeric(1))
  mu_true <- setNames(mu_true, P$params)

  sigma_true <- setNames(runif(K, P$sigma_lo, P$sigma_hi), P$params)

  z <- matrix(rnorm(L * K), nrow = L, ncol = K)
  raw <- sweep(z %*% diag(sigma_true), 2, mu_true, `+`)
  colnames(raw) <- P$params

  list(mu_true    = mu_true,
       sigma_true = sigma_true,
       raw        = raw)
}
