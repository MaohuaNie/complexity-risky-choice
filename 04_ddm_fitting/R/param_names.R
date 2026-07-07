## ============================================================================
## R/param_names.R — map participant_params indices to readable parameter names
##
## Aim:     Provide FULL_PARAM_NAMES (the full raw-parameter set per family),
##          EFFECTIVE_FN_NAME (which effective_* transform each family uses),
##          .active_full_names() (which named params a given model estimates),
##          and build_full_raw() (expand a model's K estimated params into a
##          full-length, zero-padded named raw vector). Downstream PPC code uses
##          these to feed the effective_* functions in transforms.R.
## Inputs:  A (family, model) pair and, for build_full_raw(), a raw-parameter
##          vector; no files read.
## Outputs: In-memory lookups/functions; nothing written to disk.
## Usage:   Not run directly; source()'d by the PPC-generation stage.
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM fitting stage).
## Pipeline order and dependencies are documented in ../README.md.
## ============================================================================
##
## Missing params (models with fewer than the full set, e.g. baseline = 5) are
## padded with 0: a delta of 0 means "no condition effect", which is the correct
## limiting case.

## Full param names for the MOST COMPLEX model in each family.
## Simpler models use a prefix of this list (indices 1..K).
FULL_PARAM_NAMES <- list(
  cpt_ccss = c("beta_raw", "theta_raw", "threshold_raw", "ndt_raw", "gamma_raw",
               "delta_beta", "delta_theta", "delta_threshold", "delta_gamma"),
  cpt_ccss_7o = c("beta_raw", "theta_raw", "threshold_raw", "ndt_raw", "gamma_raw",
                  "delta_beta", "delta_theta", "delta_threshold", "delta_gamma"),
  mv_ccss  = c("beta", "theta_raw", "threshold_raw", "ndt_raw", "eta",
               "delta_beta", "delta_theta", "delta_threshold", "delta_eta"),
  cpt_cs   = c("beta_raw", "theta_raw", "threshold_raw", "ndt_raw", "gamma_raw",
               "sp_raw", "zeta", "delta_gamma"),
  cpt_cs_7o = c("beta_raw", "theta_raw", "threshold_raw", "ndt_raw", "gamma_raw",
                "sp_raw", "zeta", "delta_gamma"),
  mv_cs    = c("beta", "theta_raw", "threshold_raw", "ndt_raw", "eta",
               "sp_raw", "zeta", "delta_eta")
)

## Which effective_* function to use for each family
EFFECTIVE_FN_NAME <- list(
  cpt_ccss    = "effective_cpt_ccss_n_r_a_s",
  cpt_ccss_7o = "effective_cpt_ccss_7o_n_r_a_s",
  mv_ccss     = "effective_mv_ccss_n_r_a_s",
  cpt_cs      = "effective_cpt_cs_sp_dr_skew",
  cpt_cs_7o   = "effective_cpt_cs_7o_sp_dr_skew",
  mv_cs       = "effective_mv_cs_sp_dr_skew"
)

## ---------------------------------------------------------------------------
## Active-parameter mapping: for each (family, model), which FULL_PARAM_NAMES
## entries correspond to the positions in the model's raw_k vector?
##
## Positions 1..5 are always the base params (first 5 of FULL_PARAM_NAMES).
## Positions 6+ are the "extras" (deltas for CCSS, mechanism params for CS),
## packed in canonical order — but only the ACTIVE ones are present.
##
## Canonical extra order matches how each Stan file assigns participant_params:
##   CCSS: [delta_beta, delta_theta, delta_threshold, delta_gamma]
##         (for mv_ccss: last is delta_eta)
##   CS:   [sp_raw,     zeta,        delta_gamma]
##         (for mv_cs:  last is delta_eta)
##
## Model-suffix letter/word → extra name:
##   CCSS: n=theta, r=beta, a=threshold, s=gamma-or-eta
##   CS:   sp=starting-point, dr=drift-adjust (zeta), skew=gamma-or-eta shift
## ---------------------------------------------------------------------------

.active_full_names <- function(family, model) {
  base <- FULL_PARAM_NAMES[[family]][1:5]
  if (is.null(model) || model == "baseline") return(base)

  is_mv   <- grepl("^mv_", family)
  is_ccss <- grepl("ccss", family)

  if (is_ccss) {
    letters_present <- strsplit(model, "_")[[1]]
    canonical <- c("r", "n", "a", "s")  # canonical order in Stan files
    active    <- canonical[canonical %in% letters_present]
    delta_map <- if (is_mv)
      c(r = "delta_beta", n = "delta_theta", a = "delta_threshold", s = "delta_eta")
    else
      c(r = "delta_beta", n = "delta_theta", a = "delta_threshold", s = "delta_gamma")
    extras <- unname(delta_map[active])
  } else {
    ## CS
    words_present <- strsplit(model, "_")[[1]]
    canonical <- c("sp", "dr", "skew")
    active    <- canonical[canonical %in% words_present]
    extra_map <- if (is_mv)
      c(sp = "sp_raw", dr = "zeta", skew = "delta_eta")
    else
      c(sp = "sp_raw", dr = "zeta", skew = "delta_gamma")
    extras <- unname(extra_map[active])
  }
  c(base, extras)
}

## Build a full-length named raw vector from a model's K estimated params.
## When `model` is supplied, extras are placed into their correct canonical
## slots (by name). When NULL, falls back to positional mapping — which is
## correct only for baseline or the FULL model of the family.
## Missing params are set to 0 (neutral — no shift effect).
build_full_raw <- function(family, raw_k, model = NULL) {
  full_names <- FULL_PARAM_NAMES[[family]]
  K_full  <- length(full_names)
  K_model <- length(raw_k)

  full_raw <- setNames(rep(0, K_full), full_names)

  if (!is.null(model)) {
    active_names <- .active_full_names(family, model)
    if (length(active_names) != K_model) {
      stop(sprintf(
        "build_full_raw: model '%s' in family '%s' expects %d raw params but got %d",
        model, family, length(active_names), K_model))
    }
    full_raw[active_names] <- raw_k
  } else {
    ## Legacy positional fallback
    full_raw[seq_len(K_model)] <- raw_k
  }
  full_raw
}
