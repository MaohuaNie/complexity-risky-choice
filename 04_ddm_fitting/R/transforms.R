## ============================================================================
## R/transforms.R — raw-to-effective parameter transforms (R mirror of Stan)
##
## Aim:     R re-implementation of each Stan model's transformed-parameters
##          block. The effective_*() functions take a named raw-parameter
##          vector plus trial-level covariates and return per-trial effective
##          quantities (drift, threshold, ndt, rel_sp) for the DDM. Used by the
##          PPC-generation stage to simulate choices/RTs; these MUST stay
##          numerically identical to their matching stan/ files.
## Inputs:  raw (named numeric vector) and trials (per-trial covariates); no
##          files read.
## Outputs: In-memory per-trial parameter lists via the EFFECTIVE dispatch
##          table; nothing written to disk.
## Usage:   Not run directly; source()'d by the PPC-generation stage.
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM fitting stage).
## Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

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

## ---------------------------------------------------------------------------
## Study-1 (7-outcome) variants.
## Mirror the Stan likelihood: rank-dependent CPT for CC (con == +1) over
## outcome columns 3..9 (7-outcome) in CCSS, or over the 7-col o_complex matrix
## in CS. SS (con == -1) uses 2-outcome CPT over cols 1..2 in CCSS, or 2-col
## o_simple in CS. Outcomes are assumed sorted in descending order upstream.
## ---------------------------------------------------------------------------

## Helper: 7-outcome rank-dependent CPT utility for one trial row.
##   o, p: length-7 vectors, o sorted descending
##   gamma, beta: scalars
.u_rdu7 <- function(o, p, gamma, beta) {
  ## Cumulative-weighted decision weights:
  ##   w_k = pweight(sum p[1..k]) - pweight(sum p[1..k-1])   for k = 1..6
  ##   w_7 = 1 - sum(w_1..w_6)
  cp <- cumsum(p[1:6])
  W  <- pweight(cp, gamma)
  w1 <- W[1]
  w2 <- W[2] - W[1]
  w3 <- W[3] - W[2]
  w4 <- W[4] - W[3]
  w5 <- W[5] - W[4]
  w6 <- W[6] - W[5]
  w7 <- 1 - (w1 + w2 + w3 + w4 + w5 + w6)
  w1 * o[1]^beta + w2 * o[2]^beta + w3 * o[3]^beta + w4 * o[4]^beta +
    w5 * o[5]^beta + w6 * o[6]^beta + w7 * o[7]^beta
}

## Helper: 2-outcome CPT utility (vectorised).
.u_cpt2 <- function(o1, o2, p1, gamma, beta) {
  w1 <- pweight(p1, gamma)
  w1 * o1^beta + (1 - w1) * o2^beta
}

effective_cpt_ccss_7o_n_r_a_s <- function(raw, trials) {
  con <- trials$con
  n   <- length(con)
  beta      <- softplus(raw["beta_raw"]      + raw["delta_beta"]      * con)
  theta     <- softplus(raw["theta_raw"]     + raw["delta_theta"]     * con)
  threshold <- softplus(raw["threshold_raw"] + raw["delta_threshold"] * con)
  ndt       <- rep(as.numeric(softplus(raw["ndt_raw"])), n)
  gamma     <- softplus(raw["gamma_raw"]     + raw["delta_gamma"]     * con)

  u_risky <- numeric(n)
  u_safe  <- numeric(n)
  for (i in seq_len(n)) {
    if (con[i] == 1L) {
      ## 7-outcome RDU on columns 3..9
      u_risky[i] <- .u_rdu7(trials$o_risky[i, 3:9], trials$p_risky[i, 3:9],
                            gamma[i], beta[i])
      u_safe[i]  <- .u_rdu7(trials$o_safe[i, 3:9],  trials$p_safe[i, 3:9],
                            gamma[i], beta[i])
    } else {
      ## 2-outcome CPT on columns 1..2
      u_risky[i] <- .u_cpt2(trials$o_risky[i, 1], trials$o_risky[i, 2],
                            trials$p_risky[i, 1], gamma[i], beta[i])
      u_safe[i]  <- .u_cpt2(trials$o_safe[i, 1],  trials$o_safe[i, 2],
                            trials$p_safe[i, 1],  gamma[i], beta[i])
    }
  }
  drift <- theta * (u_risky^(1 / beta) - u_safe^(1 / beta))

  list(drift = as.numeric(drift),
       threshold = as.numeric(threshold),
       ndt = ndt,
       rel_sp = rep(0.5, n))
}

effective_cpt_cs_7o_sp_dr_skew <- function(raw, trials) {
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

  ## All trials have 7-outcome complex and 2-outcome simple.
  n <- nrow(trials$o_complex)
  u_complex <- numeric(n)
  for (i in seq_len(n)) {
    u_complex[i] <- .u_rdu7(trials$o_complex[i, ], trials$p_complex[i, ],
                            gamma_c, beta)
  }
  u_simple <- .u_cpt2(trials$o_simple[, 1], trials$o_simple[, 2],
                      trials$p_simple[, 1], gamma_s, beta)
  drift <- theta * (u_complex^(1 / beta) - u_simple^(1 / beta) + zeta)

  list(drift = as.numeric(drift),
       threshold = rep(as.numeric(threshold), n),
       ndt = rep(as.numeric(ndt), n),
       rel_sp = rep(as.numeric(rel_sp), n))
}

## Dispatch table: model_key or family → effective-params function.
## The effective functions are written for the FULL model (all deltas/extras).
## For simpler models, missing params are zero-padded via build_full_raw().
EFFECTIVE <- list(
  cpt_ccss_n_r_a_s     = effective_cpt_ccss_n_r_a_s,
  cpt_ccss_7o_n_r_a_s  = effective_cpt_ccss_7o_n_r_a_s,
  mv_ccss_n_r_a_s      = effective_mv_ccss_n_r_a_s,
  cpt_cs_sp_dr_skew    = effective_cpt_cs_sp_dr_skew,
  cpt_cs_7o_sp_dr_skew = effective_cpt_cs_7o_sp_dr_skew,
  mv_cs_sp_dr_skew     = effective_mv_cs_sp_dr_skew,
  cpt_ccss    = effective_cpt_ccss_n_r_a_s,
  cpt_ccss_7o = effective_cpt_ccss_7o_n_r_a_s,
  mv_ccss     = effective_mv_ccss_n_r_a_s,
  cpt_cs      = effective_cpt_cs_sp_dr_skew,
  cpt_cs_7o   = effective_cpt_cs_7o_sp_dr_skew,
  mv_cs       = effective_mv_cs_sp_dr_skew
)
