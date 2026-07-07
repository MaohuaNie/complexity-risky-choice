## ============================================================================
## R/transforms.R — raw→effective parameter transforms (mirror the Stan models)
##
## Aim:     Reproduce, in R, each Stan model's transformed-parameters block so
##          that simulated data are generated under exactly the same drift /
##          threshold / ndt / starting-point mapping the model will later fit.
##          Provides one effective_*() per model plus the EFFECTIVE dispatch
##          table keyed by model.
## Inputs:  raw parameter vectors (from priors.R) + a trial template
##          (from trials.R).
## Outputs: per-trial lists (drift, threshold, ndt, rel_sp); no files written.
## Usage:   source("R/transforms.R") from run_recovery.R (do not run standalone).
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM parameter
## recovery). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================
## Each branch below must mirror the corresponding Stan model exactly:
## raw scale → constrained ("effective") scale used inside the drift formula.

softplus <- function(x) log1p(exp(x))
phi_norm <- function(x) pnorm(x)               # Stan's Phi()

## Prelec probability weighting, matching pweight() in the CPT Stan files.
pweight <- function(p, gamma) exp(-(-log(p))^gamma)

## ---------------------------------------------------------------------------
## Apply the per-trial transforms for one participant, one family.
##
## Input: raw parameters (named numeric vector) + trial-level covariates.
## Output: a list of per-trial effective parameters needed for DDM sim:
##           drift, threshold, ndt, rel_sp (starting point in [0,1]).
##
## Each branch mirrors the transformed-parameters block of the matching
## Stan file.
## ---------------------------------------------------------------------------

effective_cpt_ccss_n_r_a_s <- function(raw, trials) {
  con <- trials$con
  n   <- length(con)
  beta      <- softplus(raw["beta_raw"]      + raw["delta_beta"]      * con)
  theta     <- softplus(raw["theta_raw"]     + raw["delta_theta"]     * con)
  threshold <- softplus(raw["threshold_raw"] + raw["delta_threshold"] * con)
  ndt       <- rep(as.numeric(softplus(raw["ndt_raw"])), n)  # broadcast to n trials
  gamma     <- softplus(raw["gamma_raw"]     + raw["delta_gamma"]     * con)

  w_r1 <- pweight(trials$p_risky1, gamma)
  w_s1 <- pweight(trials$p_safe1,  gamma)
  u_risky <- w_r1 * trials$o_risky1^beta + (1 - w_r1) * trials$o_risky2^beta
  u_safe  <- w_s1 * trials$o_safe1^beta  + (1 - w_s1) * trials$o_safe2^beta
  drift <- theta * (u_risky^(1 / beta) - u_safe^(1 / beta))

  list(drift = as.numeric(drift),
       threshold = as.numeric(threshold),
       ndt = ndt,
       rel_sp = rep(0.5, n))
}

effective_mv_ccss_n_r_a_s <- function(raw, trials) {
  con <- trials$con
  n   <- length(con)
  beta      <- raw["beta"]      + raw["delta_beta"]      * con
  theta     <- softplus(raw["theta_raw"]     + raw["delta_theta"]     * con)
  threshold <- softplus(raw["threshold_raw"] + raw["delta_threshold"] * con)
  ndt       <- rep(as.numeric(softplus(raw["ndt_raw"])), n)
  eta       <- raw["eta"]       + raw["delta_eta"]       * con

  drift <- theta * (trials$evd + beta * trials$sdd + eta * trials$skew)

  list(drift = as.numeric(drift),
       threshold = as.numeric(threshold),
       ndt = ndt,
       rel_sp = rep(0.5, n))
}

effective_cpt_cs_sp_dr_skew <- function(raw, trials) {
  beta      <- softplus(raw["beta_raw"])
  theta     <- softplus(raw["theta_raw"])
  threshold <- softplus(raw["threshold_raw"])
  ndt       <- softplus(raw["ndt_raw"])
  gamma     <- softplus(raw["gamma_raw"])
  rel_sp    <- phi_norm(raw["sp_raw"])
  zeta      <- raw["zeta"]
  delta_g   <- raw["delta_gamma"]
  gamma_c   <- softplus(gamma + delta_g)
  gamma_s   <- softplus(gamma - delta_g)

  w_c1 <- pweight(trials$p_complex1, gamma_c)
  w_s1 <- pweight(trials$p_simple1,  gamma_s)
  u_c <- w_c1 * trials$o_complex1^beta + (1 - w_c1) * trials$o_complex2^beta
  u_s <- w_s1 * trials$o_simple1^beta  + (1 - w_s1) * trials$o_simple2^beta
  drift <- theta * (u_c^(1 / beta) - u_s^(1 / beta) + zeta)
  n <- length(drift)

  list(drift = as.numeric(drift),
       threshold = rep(as.numeric(threshold), n),
       ndt = rep(as.numeric(ndt), n),
       rel_sp = rep(as.numeric(rel_sp), n))
}

effective_mv_cs_sp_dr_skew <- function(raw, trials) {
  beta      <- raw["beta"]
  theta     <- softplus(raw["theta_raw"])
  threshold <- softplus(raw["threshold_raw"])
  ndt       <- softplus(raw["ndt_raw"])
  eta       <- raw["eta"]
  rel_sp    <- phi_norm(raw["sp_raw"])
  zeta      <- raw["zeta"]
  delta_eta <- raw["delta_eta"]

  drift <- theta * (trials$evd + beta * trials$sdd +
                    ((eta + delta_eta) * trials$skew_c -
                     (eta - delta_eta) * trials$skew_s) + zeta)
  n <- length(drift)

  list(drift = as.numeric(drift),
       threshold = rep(as.numeric(threshold), n),
       ndt = rep(as.numeric(ndt), n),
       rel_sp = rep(as.numeric(rel_sp), n))
}

effective_cpt_cs_sp_dr <- function(raw, trials) {
  beta      <- softplus(raw["beta_raw"])
  theta     <- softplus(raw["theta_raw"])
  threshold <- softplus(raw["threshold_raw"])
  ndt       <- softplus(raw["ndt_raw"])
  gamma     <- softplus(raw["gamma_raw"])
  rel_sp    <- phi_norm(raw["sp_raw"])
  zeta      <- raw["zeta"]

  w_c1 <- pweight(trials$p_complex1, gamma)
  w_s1 <- pweight(trials$p_simple1,  gamma)
  u_c <- w_c1 * trials$o_complex1^beta + (1 - w_c1) * trials$o_complex2^beta
  u_s <- w_s1 * trials$o_simple1^beta  + (1 - w_s1) * trials$o_simple2^beta
  drift <- theta * (u_c^(1 / beta) - u_s^(1 / beta) + zeta)
  n <- length(drift)

  list(drift = as.numeric(drift),
       threshold = rep(as.numeric(threshold), n),
       ndt = rep(as.numeric(ndt), n),
       rel_sp = rep(as.numeric(rel_sp), n))
}

effective_mv_cs_sp_dr <- function(raw, trials) {
  beta      <- raw["beta"]
  theta     <- softplus(raw["theta_raw"])
  threshold <- softplus(raw["threshold_raw"])
  ndt       <- softplus(raw["ndt_raw"])
  eta       <- raw["eta"]
  rel_sp    <- phi_norm(raw["sp_raw"])
  zeta      <- raw["zeta"]

  ## mv_cs_sp_dr.stan uses a single `skew = skew_c - skew_s` input
  skew <- trials$skew_c - trials$skew_s
  drift <- theta * (trials$evd + beta * trials$sdd + eta * skew + zeta)
  n <- length(drift)

  list(drift = as.numeric(drift),
       threshold = rep(as.numeric(threshold), n),
       ndt = rep(as.numeric(ndt), n),
       rel_sp = rep(as.numeric(rel_sp), n))
}

## Dispatch table: model_key → effective-params function
EFFECTIVE <- list(
  cpt_ccss_n_r_a_s  = effective_cpt_ccss_n_r_a_s,
  mv_ccss_n_r_a_s   = effective_mv_ccss_n_r_a_s,
  cpt_cs_sp_dr_skew = effective_cpt_cs_sp_dr_skew,
  mv_cs_sp_dr_skew  = effective_mv_cs_sp_dr_skew,
  cpt_cs_sp_dr      = effective_cpt_cs_sp_dr,
  mv_cs_sp_dr       = effective_mv_cs_sp_dr
)
