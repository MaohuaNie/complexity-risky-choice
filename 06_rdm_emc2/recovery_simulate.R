## ============================================================================
## recovery_simulate.R — simulate one synthetic Study-2 dataset for RDM recovery
##
## Aim:     Generate a synthetic dataset from the hierarchical RDM with ALL
##          group-level parameters drawn from reasonable ranges (not anchored to
##          the real fit) over the real Study-2 stimuli. Step 1 of the parameter-
##          recovery pipeline that checks whether every RDM parameter is
##          recoverable under the study2_rdm.R design.
##            - 50 simulated participants (default)
##            - real Study-2 CS/SC stimuli (EVD/SDD/SkewD), one real subject's
##              trial sequence assigned to each simulated participant
##            - design: A ~ 1 (shared start-point range) -> 13 free parameters
## Inputs:  ../data/final_df_study2.rds (real stimuli; auto-discovered or --real_data)
## Outputs: <--out> .rds holding the simulated data plus known group- and
##          participant-level truths
## Usage:   Rscript recovery_simulate.R --seed 1 --out recovery/sim_seed01.rds
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (RDM / EMC2 robustness
## analysis). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

suppressPackageStartupMessages({
  library(optparse); library(EMC2); library(dplyr)
})

opt_list <- list(
  make_option("--seed",      type = "integer",   default = 1),
  make_option("--out",       type = "character", default = "sim.rds"),
  make_option("--n_subj",    type = "integer",   default = 50),
  make_option("--real_data", type = "character", default = NA_character_,
              help = "Path to real Study 2 data; auto-searched if NA."),
  make_option("--max_reject", type = "integer",  default = 30,
              help = "Max resample passes for pathological participants.")
)
opt <- parse_args(OptionParser(option_list = opt_list))
set.seed(opt$seed)

## ---------------------------------------------------------------
## 1. Generative ranges  (EDIT HERE to recalibrate)
##    Positive parameters are drawn on the NATURAL scale and then
##    log-transformed to EMC2's sampled scale. Drift slopes are drawn
##    directly on the additive log-drift scale.
## ---------------------------------------------------------------

draw_u <- function(lo, hi) runif(1, lo, hi)

# group-level means (natural scale unless noted)
g <- list(
  v_simple      = draw_u(1.0, 3.5),   # drift intercept, simple accumulator
  v_complex     = draw_u(1.0, 3.5),   # drift intercept, complex accumulator
  slope_EVD_c   = draw_u(0.0,  0.6),  # EVD slope, complex acc (positive: follow EV)
  slope_EVD_s   = draw_u(-0.6, 0.0),  # EVD slope, simple  acc (negative mirror)
  slope_SDD_c   = draw_u(-0.3, 0.3),  # risk slope, complex acc (free sign)
  slope_SDD_s   = draw_u(-0.3, 0.3),
  slope_SkewD_c = draw_u(-0.3, 0.3),  # skew slope, complex acc (free sign)
  slope_SkewD_s = draw_u(-0.3, 0.3),
  s_complex     = draw_u(0.5, 1.5),   # within-trial noise, complex (s_simple = 1)
  B_simple      = draw_u(1.0, 3.0),   # threshold, simple
  B_complex     = draw_u(1.0, 3.0),   # threshold, complex
  A             = draw_u(0.3, 1.5),   # shared start-point range
  t0            = draw_u(0.20, 0.60)  # non-decision time (s), within analysis-prior support
)

# group-level SDs, on the sampled (log/additive) scale
sd_intercept <- function() draw_u(0.10, 0.40)
sd_slope     <- function() draw_u(0.05, 0.25)
sd_t0        <- function() draw_u(0.05, 0.20)

## Map to EMC2 sampled-scale, named by the model's parameter labels.
true_mu <- c(
  "v_lRsimple"        = log(g$v_simple),
  "v_lRcomplex"       = log(g$v_complex),
  "v_lRsimple:EVD"    = g$slope_EVD_s,
  "v_lRcomplex:EVD"   = g$slope_EVD_c,
  "v_lRsimple:SDD"    = g$slope_SDD_s,
  "v_lRcomplex:SDD"   = g$slope_SDD_c,
  "v_lRsimple:SkewD"  = g$slope_SkewD_s,
  "v_lRcomplex:SkewD" = g$slope_SkewD_c,
  "s_lRcomplex"       = log(g$s_complex),
  "B_lRsimple"        = log(g$B_simple),
  "B_lRcomplex"       = log(g$B_complex),
  "A"                 = log(g$A),
  "t0"                = log(g$t0)
)
true_sd <- c(
  "v_lRsimple" = sd_intercept(), "v_lRcomplex" = sd_intercept(),
  "v_lRsimple:EVD" = sd_slope(), "v_lRcomplex:EVD" = sd_slope(),
  "v_lRsimple:SDD" = sd_slope(), "v_lRcomplex:SDD" = sd_slope(),
  "v_lRsimple:SkewD" = sd_slope(), "v_lRcomplex:SkewD" = sd_slope(),
  "s_lRcomplex" = sd_slope(),
  "B_lRsimple" = sd_intercept(), "B_lRcomplex" = sd_intercept(),
  "A" = sd_intercept(), "t0" = sd_t0()
)

cat(sprintf("seed %d — true group means (sampled scale):\n", opt$seed))
print(round(true_mu, 3))

## ---------------------------------------------------------------
## 2. Real Study-2 stimuli (EVD/SDD/SkewD), prepared as in study2_rdm.R
## ---------------------------------------------------------------

candidates <- c(opt$real_data, "../data/final_df_study2.rds",
                "final_df_study2.rds",
                "outcome-as-term/data/derived/final_df_study2.rds")
candidates <- candidates[!is.na(candidates)]
data_path  <- candidates[file.exists(candidates)][1]
if (is.na(data_path)) stop("Could not find real Study 2 data. Tried:\n  ",
                           paste(candidates, collapse = "\n  "))
cat("Using real data from:", data_path, "\n")

study2 <- readRDS(data_path)
dat_cs <- study2 %>%
  filter(test_part %in% c("cs", "sc")) %>%
  mutate(
    rt         = as.numeric(rt),
    oa_complex = ifelse(test_part == "cs", 1, -1),
    EVD   = EV_diff   * oa_complex,
    SDD   = SD_diff   * oa_complex,
    SkewD = Skew_diff * oa_complex,
    subject = as.character(subject)
  ) %>%
  filter(is.finite(rt), rt > 0) %>%
  mutate(across(c(EVD, SDD, SkewD), ~ as.numeric(scale(.))))   # standardise, as in the real fit

real_subjects <- unique(dat_cs$subject)

## Assign each simulated participant a randomly chosen real subject's trial
## sequence (preserves realistic per-participant trial counts + stimulus structure).
donor <- sample(real_subjects, opt$n_subj, replace = TRUE)
skeleton <- do.call(rbind, lapply(seq_len(opt$n_subj), function(i) {
  rows <- dat_cs[dat_cs$subject == donor[i], c("EVD", "SDD", "SkewD")]
  data.frame(
    subjects = factor(sprintf("sim%03d", i),
                      levels = sprintf("sim%03d", seq_len(opt$n_subj))),
    rt    = NA_real_,
    R     = factor("simple", levels = c("simple", "complex")),
    EVD   = rows$EVD, SDD = rows$SDD, SkewD = rows$SkewD
  )
}))
skeleton$S <- skeleton$R

## ---------------------------------------------------------------
## 3. Design — A ~ 1 (shared start-point range)
## ---------------------------------------------------------------

design_RDM_h <- design(
  data     = skeleton,
  model    = RDM,
  matchfun = function(d) d$S == d$lR,
  formula  = list(
    v  ~ 0 + lR + lR:(EVD + SDD + SkewD),
    s  ~ 0 + lR,
    B  ~ 0 + lR,
    A  ~ 1,
    t0 ~ 1
  ),
  constants = c(s_lRsimple = 0)
)

par_names <- names(sampled_pars(design_RDM_h))
stopifnot(setequal(par_names, names(true_mu)))
true_mu <- true_mu[par_names]; true_sd <- true_sd[par_names]

## ---------------------------------------------------------------
## 4. Per-participant draws + simulation, with pathology resampling
##    (reject participants with <2% or >98% complex choices)
## ---------------------------------------------------------------

draw_p_subj <- function(idx) {
  m <- matrix(NA, length(idx), length(true_mu),
              dimnames = list(levels(skeleton$subjects)[idx], names(true_mu)))
  for (k in seq_along(true_mu)) m[, k] <- rnorm(length(idx), true_mu[k], true_sd[k])
  m
}
subj_prop_complex <- function(dat) {
  as.data.frame(dat) %>%
    mutate(cc = as.integer(R == "complex")) %>%
    group_by(subjects) %>% summarise(p = mean(cc), n = n(), .groups = "drop")
}

p_subj <- draw_p_subj(seq_len(opt$n_subj))
sim_dat <- make_data(p_subj, design = design_RDM_h, n_trials = NULL,
                     data = skeleton[, c("subjects", "R", "S", "EVD", "SDD", "SkewD")])

for (pass in seq_len(opt$max_reject)) {
  pc <- subj_prop_complex(sim_dat)
  bad <- which(pc$p < 0.02 | pc$p > 0.98)
  if (length(bad) == 0) break
  cat(sprintf("  resample pass %d: %d pathological participant(s)\n", pass, length(bad)))
  p_subj[bad, ] <- draw_p_subj(bad)
  sim_dat <- make_data(p_subj, design = design_RDM_h, n_trials = NULL,
                       data = skeleton[, c("subjects", "R", "S", "EVD", "SDD", "SkewD")])
}
pc <- subj_prop_complex(sim_dat)
n_bad <- sum(pc$p < 0.02 | pc$p > 0.98)
if (n_bad > 0) cat(sprintf("  WARNING: %d participant(s) still pathological after %d passes\n",
                           n_bad, opt$max_reject))

cat(sprintf("\nSimulated complex-choice proportion: median %.2f, range [%.2f, %.2f]\n",
            median(pc$p), min(pc$p), max(pc$p)))
rt_all <- as.data.frame(sim_dat)$rt
cat(sprintf("Simulated RT (s): median %.2f, 5-95%% [%.2f, %.2f], max %.2f\n",
            median(rt_all), quantile(rt_all, .05), quantile(rt_all, .95), max(rt_all)))

## ---------------------------------------------------------------
## 5. Save
## ---------------------------------------------------------------

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(
  data        = sim_dat,
  true_mu     = true_mu,          # group-level truths (sampled scale)
  true_sd     = true_sd,
  per_subject = p_subj,           # per-participant truths (sampled scale)
  prop_complex = pc,              # per-participant complex-choice proportion
  seed        = opt$seed,
  n_subj      = opt$n_subj,
  par_names   = par_names
), file = opt$out)
cat("\nWrote:", opt$out, "\n")
