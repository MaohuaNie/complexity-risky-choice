## ============================================================================
## R/inits.R — initial values for the hierarchical Stan models
##
## Aim:     generate_inits() builds per-chain init lists for any model in the
##          registry (one parametric generator handling every n_params). It
##          centres each group-mean at a psychologically plausible raw-scale
##          value, jitters per chain, and caps the non-decision-time init so
##          the first log-prob evaluation satisfies rt > ndt everywhere
##          (avoiding wiener_lpdf init failures).
## Inputs:  family, n_params, L (number of participants), chains, and the
##          data's minimum RT (min_rt) — supplied by fit.R.
## Outputs: A list of `chains` init lists (mu, sigma, L_corr, z), returned to
##          the caller; nothing is written to disk.
## Usage:   Not run directly; source()'d and called by R/fit.R.
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM fitting stage).
## Pipeline order and dependencies are documented in ../README.md.
## ============================================================================
##
## Param-index contract: position 4 in Stan's mu/sigma vectors is always ndt.
## This matches every Stan file in stan/; change NDT_IDX if that ever changes.

NDT_IDX <- 4L

## raw-scale centres for the first 5 (base) parameters per family.
BASE_INIT_CENTRE <- list(
  cpt_ccss    = c(0.5414, -1.0,  1.3, -0.8, 0.5414),  # beta, theta, threshold, ndt, gamma
  cpt_ccss_7o = c(0.5414, -1.0,  1.3, -0.8, 0.5414),
  cpt_cs      = c(0.5414, -1.0,  1.3, -0.8, 0.5414),
  cpt_cs_7o   = c(0.5414, -1.0,  1.3, -0.8, 0.5414),
  mv_ccss     = c(0.0,    -2.0,  1.3, -0.8, 0.0),     # beta, theta, threshold, ndt, eta
  mv_cs       = c(0.0,    -2.0,  1.3, -0.8, 0.0)
)

## Extra (post-5) parameters: their init centres default to 0. CS-family
## adds sp/dr/skew extras; CCSS adds deltas. 0 is a neutral, defensible
## starting point for all of them.

generate_inits <- function(family, n_params, L, chains = 4,
                           jitter_mu = 0.25, min_rt = NULL) {
  base <- BASE_INIT_CENTRE[[family]]
  if (is.null(base)) stop("Unknown family: ", family)

  mu_centre <- c(base, rep(0, max(0, n_params - 5L)))

  lapply(seq_len(chains), function(ch) {
    mu_init    <- mu_centre + rnorm(n_params, 0, jitter_mu)
    sigma_init <- pmax(runif(n_params, 0.15, 0.4), 0.05)

    ## Cap ndt init so softplus(mu_ndt) < 0.3 * min_rt.
    ## Inverse softplus: log(exp(x) - 1).
    if (!is.null(min_rt) && is.finite(min_rt) && n_params >= NDT_IDX) {
      target_ndt <- max(0.02, 0.3 * min_rt)
      cap        <- log(exp(target_ndt) - 1)
      mu_init[NDT_IDX]    <- min(mu_init[NDT_IDX], cap)
      sigma_init[NDT_IDX] <- min(sigma_init[NDT_IDX], 0.15)
    }

    list(
      mu     = mu_init,
      sigma  = sigma_init,
      L_corr = diag(n_params),
      z      = matrix(rnorm(n_params * L, 0, 1), n_params, L)
    )
  })
}
