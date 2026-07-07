## ============================================================================
## R/trials.R — build the trial template used to simulate every dataset
##
## Aim:     Load one preprocessed Study RDS, take a single subject's trials,
##          and recode them into the covariate layout each model family
##          expects (CCSS vs CS; CPT vs MV). The recovery reuses this real
##          stimulus structure so that cross-participant differences come only
##          from parameters.
## Inputs:  data_rds (a preprocessed Study RDS, e.g. data/final_df_study2.rds),
##          family, optional template_subject.
## Outputs: a data.frame of recoded trials (returned in memory; no files).
## Usage:   source("R/trials.R"); load_trials_template("data/final_df_study2.rds",
##          family = "cpt_ccss")  (called from run_recovery.R).
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM parameter
## recovery). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

## Load the study RDS and return a data.frame of trials for one "template"
## subject, filtered/recoded to the family's input layout.
##
## Returns columns (by family):
##   CCSS families:
##     o_risky1, o_risky2, p_risky1, p_risky2,
##     o_safe1,  o_safe2,  p_safe1,  p_safe2,        -- CPT only
##     evd, sdd, skew,                               -- MV only
##     con (+1 / -1)
##   CS families:
##     o_complex1, o_complex2, p_complex1, p_complex2,
##     o_simple1,  o_simple2,  p_simple1,  p_simple2,  -- CPT only
##     evd, sdd, skew_c, skew_s,                        -- MV only
##     accuracy_flipped
## ---------------------------------------------------------------------------

load_trials_template <- function(data_rds, family,
                                 template_subject = NULL) {
  if (!file.exists(data_rds)) stop("Data RDS not found: ", data_rds)
  d <- readRDS(data_rds)

  if (family %in% c("cpt_ccss", "mv_ccss")) {
    d <- d %>%
      filter(skew_level != "catch" & test_part %in% c("cc", "ss"))
  } else if (family %in% c("cpt_cs", "mv_cs")) {
    d <- d %>%
      filter(skew_level != "catch" & test_part %in% c("cs", "sc"))
  } else {
    stop("Unknown family: ", family)
  }

  ## pick a template subject (default: first by trial count)
  if (is.null(template_subject)) {
    counts <- d %>% count(subject, sort = TRUE)
    template_subject <- counts$subject[1]
  }
  t0 <- d %>% filter(subject == template_subject)
  if (nrow(t0) == 0) stop("Template subject has no trials: ", template_subject)

  ## CCSS: reorder within each option so first outcome is larger (matches Stan)
  if (family %in% c("cpt_ccss", "mv_ccss")) {
    t0 <- t0 %>%
      mutate(
        idx_A = if_else(O_A1 >= O_A2, 1L, -1L),
        idx_B = if_else(O_B1 >= O_B2, 1L, -1L),
        o_risky1 = if_else(idx_A == 1, O_A1, O_A2),
        o_risky2 = if_else(idx_A == 1, O_A2, O_A1),
        p_risky1 = if_else(idx_A == 1, P_A1, P_A2),
        p_risky2 = if_else(idx_A == 1, P_A2, P_A1),
        o_safe1  = if_else(idx_B == 1, O_B1, O_B2),
        o_safe2  = if_else(idx_B == 1, O_B2, O_B1),
        p_safe1  = if_else(idx_B == 1, P_B1, P_B2),
        p_safe2  = if_else(idx_B == 1, P_B2, P_B1),
        con      = if_else(test_part == "cc", 1L, -1L)
      )
    if (family == "mv_ccss") {
      ## MV uses precomputed differences (name in raw rds: EV_diff, SD_diff, Skew_diff)
      t0 <- t0 %>% mutate(evd = EV_diff, sdd = SD_diff, skew = Skew_diff)
    }
  }

  ## CS: build complex-vs-simple lotteries; match the Stan CS preprocessing
  if (family %in% c("cpt_cs", "mv_cs")) {
    t0 <- t0 %>%
      mutate(
        A_o_hi = pmax(O_A1, O_A2), A_o_lo = pmin(O_A1, O_A2),
        A_p_hi = if_else(O_A1 >= O_A2, P_A1, P_A2),
        A_p_lo = if_else(O_A1 >= O_A2, P_A2, P_A1),
        B_o_hi = pmax(O_B1, O_B2), B_o_lo = pmin(O_B1, O_B2),
        B_p_hi = if_else(O_B1 >= O_B2, P_B1, P_B2),
        B_p_lo = if_else(O_B1 >= O_B2, P_B2, P_B1),
        oa_complex = if_else(test_part == "cs", 1L, -1L),
        o_complex1 = if_else(oa_complex == 1L, A_o_hi, B_o_hi),
        o_complex2 = if_else(oa_complex == 1L, A_o_lo, B_o_lo),
        p_complex1 = if_else(oa_complex == 1L, A_p_hi, B_p_hi),
        p_complex2 = if_else(oa_complex == 1L, A_p_lo, B_p_lo),
        o_simple1  = if_else(oa_complex == 1L, B_o_hi, A_o_hi),
        o_simple2  = if_else(oa_complex == 1L, B_o_lo, A_o_lo),
        p_simple1  = if_else(oa_complex == 1L, B_p_hi, A_p_hi),
        p_simple2  = if_else(oa_complex == 1L, B_p_lo, A_p_lo)
      )
    if (family == "mv_cs") {
      ## MV CS needs per-option skewness (complex and simple). If the RDS
      ## already has Skewness_c / Skewness_s we use them; otherwise we
      ## compute the closed-form skewness for a 2-outcome lottery:
      ##   skew = (1 - 2 p_hi) / sqrt(p_hi * (1 - p_hi)) * sign(o_hi - o_lo)
      ## (standardised third central moment of a binary distribution, signed
      ## so that o_hi > o_lo is the reference orientation)
      lottery_skew <- function(p_hi, o_hi, o_lo) {
        p <- pmin(pmax(p_hi, 1e-6), 1 - 1e-6)
        s <- (1 - 2 * p) / sqrt(p * (1 - p))
        s * sign(o_hi - o_lo)
      }
      if (all(c("Skewness_c", "Skewness_s") %in% names(t0))) {
        t0 <- t0 %>% mutate(skew_c = Skewness_c, skew_s = Skewness_s)
      } else {
        t0 <- t0 %>% mutate(
          skew_c = lottery_skew(p_complex1, o_complex1, o_complex2),
          skew_s = lottery_skew(p_simple1,  o_simple1,  o_simple2)
        )
      }
      t0 <- t0 %>% mutate(evd = EV_diff, sdd = SD_diff)
    }
  }

  ## minimal shared columns
  t0$trial_idx <- seq_len(nrow(t0))
  t0
}
