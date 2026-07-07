## ============================================================================
## R/preprocess.R — turn a study data RDS into a Stan-ready data list
##
## Aim:     One preprocessor per family (dispatched via PREPROCESSOR). Each
##          reads the study RDS, drops catch trials, restricts to the relevant
##          test_part conditions, reorients outcomes/probabilities into the
##          direction Stan expects (risky-vs-safe for CCSS, complex-vs-simple
##          for CS), derives covariates (drift moments, per-option skewness),
##          runs sanity checks, and builds a contiguous participant index.
## Inputs:  A study data RDS (e.g. data/final_df_study2.rds) passed as data_rds.
## Outputs: A list per call: $stan_data (passed to model$sample()) and $id_map
##          (participant index <-> subject ID). Nothing written to disk here;
##          fit.R writes id_map.csv.
## Usage:   Not run directly; source()'d and called by R/fit.R.
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM fitting stage).
## Pipeline order and dependencies are documented in ../README.md.
## ============================================================================
##
## Each family entry point returns:
##   stan_data : list passed directly to cmdstan model$sample()
##   id_map    : data.frame(participant, subject) for re-joining results

suppressPackageStartupMessages({
  library(dplyr)
})

.reorder_ccss <- function(d) {
  ## Reorder outcomes within each option so the larger outcome comes first
  ## (matches Stan: `risky_o1 > risky_o2` etc.)
  d |>
    mutate(
      idx_A = if_else(O_A1 >= O_A2, 1L, -1L),
      idx_B = if_else(O_B1 >= O_B2, 1L, -1L),
      risky_o1 = if_else(idx_A == 1L, O_A1, O_A2),
      risky_o2 = if_else(idx_A == 1L, O_A2, O_A1),
      risky_p1 = if_else(idx_A == 1L, P_A1, P_A2),
      risky_p2 = if_else(idx_A == 1L, P_A2, P_A1),
      safe_o1  = if_else(idx_B == 1L, O_B1, O_B2),
      safe_o2  = if_else(idx_B == 1L, O_B2, O_B1),
      safe_p1  = if_else(idx_B == 1L, P_B1, P_B2),
      safe_p2  = if_else(idx_B == 1L, P_B2, P_B1),
      con      = if_else(test_part == "cc", 1L, -1L),
      cho      = if_else(true_response == "f", 1L, -1L)
    )
}

.reorder_cs <- function(d) {
  ## For CS: re-order within option, then assign complex/simple based on test_part
  ##   test_part == "cs" → A is complex, B is simple
  ##   test_part == "sc" → B is complex, A is simple
  d |>
    mutate(
      A_o_hi = pmax(O_A1, O_A2), A_o_lo = pmin(O_A1, O_A2),
      A_p_hi = if_else(O_A1 >= O_A2, P_A1, P_A2),
      A_p_lo = if_else(O_A1 >= O_A2, P_A2, P_A1),
      B_o_hi = pmax(O_B1, O_B2), B_o_lo = pmin(O_B1, O_B2),
      B_p_hi = if_else(O_B1 >= O_B2, P_B1, P_B2),
      B_p_lo = if_else(O_B1 >= O_B2, P_B2, P_B1),
      oa_complex = if_else(test_part == "cs", 1L, -1L),
      complex_o1 = if_else(oa_complex == 1L, A_o_hi, B_o_hi),
      complex_o2 = if_else(oa_complex == 1L, A_o_lo, B_o_lo),
      complex_p1 = if_else(oa_complex == 1L, A_p_hi, B_p_hi),
      complex_p2 = if_else(oa_complex == 1L, A_p_lo, B_p_lo),
      simple_o1  = if_else(oa_complex == 1L, B_o_hi, A_o_hi),
      simple_o2  = if_else(oa_complex == 1L, B_o_lo, A_o_lo),
      simple_p1  = if_else(oa_complex == 1L, B_p_hi, A_p_hi),
      simple_p2  = if_else(oa_complex == 1L, B_p_lo, A_p_lo),
      cho = if_else(true_response == "f", 1L, -1L),
      chose_complex = if_else(
        (oa_complex == 1L & cho ==  1L) | (oa_complex == -1L & cho == -1L),
        1L, -1L),
      accuracy_flipped = if_else(chose_complex == 1L, 0L, 1L)
    )
}

.lottery_skew <- function(p_hi, o_hi, o_lo) {
  p <- pmin(pmax(p_hi, 1e-6), 1 - 1e-6)
  s <- (1 - 2 * p) / sqrt(p * (1 - p))
  s * sign(o_hi - o_lo)
}

.common_checks <- function(d) {
  stopifnot(nrow(d) > 0)
  if (any(!is.finite(d$rt))) stop("Non-finite rt values.")
  if (any(d$rt <= 0)) stop("Non-positive rt values.")
}

.participant_index <- function(d) {
  subjects <- d |> distinct(subject) |> arrange(subject) |> pull(subject)
  id <- as.integer(match(d$subject, subjects))
  list(subjects = subjects, id = id, L = length(subjects))
}

## ---------------------------------------------------------------------------
## Two-outcome (study 2 / study 3) family entry points.
## ---------------------------------------------------------------------------

preprocess_cpt_ccss <- function(data_rds) {
  d <- readRDS(data_rds) |>
    filter(skew_level != "catch" & test_part %in% c("cc", "ss")) |>
    .reorder_ccss()
  .common_checks(d)
  idx <- .participant_index(d)

  stan_data <- list(
    N = nrow(d), L = idx$L,
    participant = idx$id,
    cho = as.integer(d$cho),
    rt  = as.numeric(d$rt),
    o_risky = cbind(d$risky_o1, d$risky_o2),
    o_safe  = cbind(d$safe_o1,  d$safe_o2),
    p_risky = cbind(d$risky_p1, d$risky_p2),
    p_safe  = cbind(d$safe_p1,  d$safe_p2),
    con = as.integer(d$con),
    starting_point = 0.5
  )
  list(stan_data = stan_data,
       id_map    = data.frame(participant = seq_len(idx$L),
                              subject     = idx$subjects))
}

preprocess_mv_ccss <- function(data_rds) {
  d <- readRDS(data_rds) |>
    filter(skew_level != "catch" & test_part %in% c("cc", "ss")) |>
    .reorder_ccss() |>
    mutate(evd = EV_diff, sdd = SD_diff, skew = Skew_diff)
  .common_checks(d)
  idx <- .participant_index(d)

  stan_data <- list(
    N = nrow(d), L = idx$L,
    participant = idx$id,
    cho = as.integer(d$cho),
    rt  = as.numeric(d$rt),
    evd = as.numeric(d$evd),
    sdd = as.numeric(d$sdd),
    skew = as.numeric(d$skew),
    con = as.integer(d$con),
    starting_point = 0.5
  )
  list(stan_data = stan_data,
       id_map    = data.frame(participant = seq_len(idx$L),
                              subject     = idx$subjects))
}

preprocess_cpt_cs <- function(data_rds) {
  d <- readRDS(data_rds) |>
    filter(skew_level != "catch" & test_part %in% c("cs", "sc")) |>
    .reorder_cs()
  .common_checks(d)
  idx <- .participant_index(d)

  stan_data <- list(
    N = nrow(d), L = idx$L,
    participant = idx$id,
    ## `cho` passed to Stan is chose_complex (±1 aligned with complex vs simple),
    ## NOT the raw keypress. This matches drift_t which is computed in the
    ## complex−simple direction via the reoriented complex_*/simple_* columns.
    ## Using raw cho here breaks sign alignment on sc trials (where oa_complex=-1)
    ## and makes `rel_sp_ll = accuracy_flipped + cho * rel_sp` step out of [0,1].
    cho = as.integer(d$chose_complex),
    accuracy_flipped = as.integer(d$accuracy_flipped),
    rt  = as.numeric(d$rt),
    o_complex = cbind(d$complex_o1, d$complex_o2),
    o_simple  = cbind(d$simple_o1,  d$simple_o2),
    p_complex = cbind(d$complex_p1, d$complex_p2),
    p_simple  = cbind(d$simple_p1,  d$simple_p2),
    starting_point = 0.5   # used by non-sp models; sp models use rel_sp instead
  )
  list(stan_data = stan_data,
       id_map    = data.frame(participant = seq_len(idx$L),
                              subject     = idx$subjects))
}

preprocess_mv_cs <- function(data_rds) {
  d <- readRDS(data_rds) |>
    filter(skew_level != "catch" & test_part %in% c("cs", "sc")) |>
    .reorder_cs()
  ## Per-option skewness, always computed from direction-mapped complex_*/simple_*
  ## columns built by .reorder_cs() — independent of any Skewness_c/s columns in
  ## the RDS (whose semantics we can't guarantee). This is both `skew_c` and
  ## `skew_s`, in complex-option / simple-option order.
  d <- d |> mutate(
    skew_c = .lottery_skew(complex_p1, complex_o1, complex_o2),
    skew_s = .lottery_skew(simple_p1,  simple_o1,  simple_o2)
  )
  ## Direction-adjust evd and sdd to complex - simple direction.
  ## The RDS stores EV_diff = EV_A - EV_B and SD_diff = SD_A - SD_B; for
  ## test_part == "sc" trials, A is the simple option, so we flip the sign.
  ## skew is already in complex - simple direction (from skew_c - skew_s below)
  ## so all three drift terms share a consistent direction.
  d <- d |> mutate(evd = oa_complex * EV_diff,
                   sdd = oa_complex * SD_diff)
  .common_checks(d)
  idx <- .participant_index(d)

  stan_data <- list(
    N = nrow(d), L = idx$L,
    participant = idx$id,
    ## `cho` = chose_complex (±1), not raw keypress — see preprocess_cpt_cs comment.
    cho = as.integer(d$chose_complex),
    accuracy_flipped = as.integer(d$accuracy_flipped),
    rt  = as.numeric(d$rt),
    evd = as.numeric(d$evd),
    sdd = as.numeric(d$sdd),
    ## Per-option skewness (used by models with delta_eta: skew, sp_skew, dr_skew, sp_dr_skew)
    skew_c = as.numeric(d$skew_c),
    skew_s = as.numeric(d$skew_s),
    ## Option-difference skewness (used by models without delta_eta: baseline, sp, dr, sp_dr).
    ## This definition makes the non-skew models the limiting case of the skew models when
    ## delta_eta = 0:  (eta + 0) * skew_c - (eta - 0) * skew_s  =  eta * (skew_c - skew_s).
    skew = as.numeric(d$skew_c - d$skew_s),
    starting_point = 0.5   # used by non-sp models; sp models use rel_sp instead
  )
  list(stan_data = stan_data,
       id_map    = data.frame(participant = seq_len(idx$L),
                              subject     = idx$subjects))
}

## ---------------------------------------------------------------------------
## Study-1 (number-of-outcomes) CPT variants: 7-outcome complex gambles.
##
## Stan array layout:
##   CCSS: o_risky/p_risky/o_safe/p_safe are N x 9
##         cols 1-2 = 2-outcome version (used on SS trials, con == -1)
##         cols 3-9 = 7-outcome version (used on CC trials, con == +1)
##   CS:   o_complex/p_complex are N x 7  (always 7-outcome)
##         o_simple/p_simple   are N x 2  (always 2-outcome)
##
## Data upstream in final_df_study1.rds already stores outcomes in descending
## order within each option (required for rank-dependent CPT weighting).
## ---------------------------------------------------------------------------

preprocess_cpt_ccss_7o <- function(data_rds) {
  d <- readRDS(data_rds) |>
    filter(skew_level != "catch" & test_part %in% c("cc", "ss")) |>
    mutate(
      cho = if_else(true_response == "f", 1L, -1L),
      con = if_else(test_part == "cc", 1L, -1L)
    )
  .common_checks(d)
  idx <- .participant_index(d)

  ## 9-column arrays: [O_A1, O_A2, complex_OA1..complex_OA7]
  o_risky <- cbind(d$O_A1, d$O_A2,
                   d$complex_OA1, d$complex_OA2, d$complex_OA3, d$complex_OA4,
                   d$complex_OA5, d$complex_OA6, d$complex_OA7)
  p_risky <- cbind(d$P_A1, d$P_A2,
                   d$complex_PA1, d$complex_PA2, d$complex_PA3, d$complex_PA4,
                   d$complex_PA5, d$complex_PA6, d$complex_PA7)
  o_safe  <- cbind(d$O_B1, d$O_B2,
                   d$complex_OB1, d$complex_OB2, d$complex_OB3, d$complex_OB4,
                   d$complex_OB5, d$complex_OB6, d$complex_OB7)
  p_safe  <- cbind(d$P_B1, d$P_B2,
                   d$complex_PB1, d$complex_PB2, d$complex_PB3, d$complex_PB4,
                   d$complex_PB5, d$complex_PB6, d$complex_PB7)

  stan_data <- list(
    N = nrow(d), L = idx$L,
    participant = idx$id,
    cho = as.integer(d$cho),
    rt  = as.numeric(d$rt),
    o_risky = o_risky, o_safe = o_safe,
    p_risky = p_risky, p_safe = p_safe,
    con = as.integer(d$con),
    starting_point = 0.5,
    ## Precomputed moments for PPC behavioural metrics (pass-through, not used by Stan)
    evd = if ("EV_diff"   %in% names(d)) as.numeric(d$EV_diff)   else NULL,
    sdd = if ("SD_diff"   %in% names(d)) as.numeric(d$SD_diff)   else NULL,
    skew= if ("Skew_diff" %in% names(d)) as.numeric(d$Skew_diff) else NULL
  )
  list(stan_data = stan_data,
       id_map    = data.frame(participant = seq_len(idx$L),
                              subject     = idx$subjects))
}

preprocess_cpt_cs_7o <- function(data_rds) {
  d <- readRDS(data_rds) |>
    filter(skew_level != "catch" & test_part %in% c("cs", "sc")) |>
    mutate(
      cho = if_else(true_response == "f", 1L, -1L),
      oa_complex = if_else(test_part == "cs", 1L, -1L),
      chose_complex = if_else(
        (oa_complex == 1L & cho ==  1L) | (oa_complex == -1L & cho == -1L),
        1L, -1L),
      accuracy_flipped = if_else(chose_complex == 1L, 0L, 1L)
    )
  .common_checks(d)
  idx <- .participant_index(d)

  ## Build complex (7-outcome) arrays dynamically based on which side is complex.
  ## For test_part == "cs": complex = A, simple = B.
  ## For test_part == "sc": complex = B, simple = A.
  pick <- function(cs, sc) if_else(d$oa_complex == 1L, cs, sc)

  o_complex <- cbind(
    pick(d$complex_OA1, d$complex_OB1),
    pick(d$complex_OA2, d$complex_OB2),
    pick(d$complex_OA3, d$complex_OB3),
    pick(d$complex_OA4, d$complex_OB4),
    pick(d$complex_OA5, d$complex_OB5),
    pick(d$complex_OA6, d$complex_OB6),
    pick(d$complex_OA7, d$complex_OB7)
  )
  p_complex <- cbind(
    pick(d$complex_PA1, d$complex_PB1),
    pick(d$complex_PA2, d$complex_PB2),
    pick(d$complex_PA3, d$complex_PB3),
    pick(d$complex_PA4, d$complex_PB4),
    pick(d$complex_PA5, d$complex_PB5),
    pick(d$complex_PA6, d$complex_PB6),
    pick(d$complex_PA7, d$complex_PB7)
  )
  o_simple <- cbind(pick(d$O_B1, d$O_A1), pick(d$O_B2, d$O_A2))
  p_simple <- cbind(pick(d$P_B1, d$P_A1), pick(d$P_B2, d$P_A2))

  stan_data <- list(
    N = nrow(d), L = idx$L,
    participant = idx$id,
    ## `cho` = chose_complex (±1), not raw keypress — see preprocess_cpt_cs comment.
    cho = as.integer(d$chose_complex),
    accuracy_flipped = as.integer(d$accuracy_flipped),
    rt  = as.numeric(d$rt),
    o_complex = o_complex, o_simple = o_simple,
    p_complex = p_complex, p_simple = p_simple,
    starting_point = 0.5
  )
  list(stan_data = stan_data,
       id_map    = data.frame(participant = seq_len(idx$L),
                              subject     = idx$subjects))
}

## Dispatch table: family → preprocess function
PREPROCESSOR <- list(
  cpt_ccss    = preprocess_cpt_ccss,
  cpt_ccss_7o = preprocess_cpt_ccss_7o,
  mv_ccss     = preprocess_mv_ccss,
  cpt_cs      = preprocess_cpt_cs,
  cpt_cs_7o   = preprocess_cpt_cs_7o,
  mv_cs       = preprocess_mv_cs
)
