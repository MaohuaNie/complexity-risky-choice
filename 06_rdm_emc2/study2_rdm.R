## ============================================================================
## study2_rdm.R — fit the hierarchical Racing Diffusion Model to Study-2 CS/SC trials
##
## Aim:     Fit a hierarchical RDM (EMC2) to the complex-vs-simple gamble trials
##          (test parts "cs"/"sc") of Study 2, estimating per-accumulator drift,
##          within-trial noise, and threshold contrasts. This is the primary
##          race-model analysis whose robustness the recovery pipeline checks.
## Inputs:  ../data/final_df_study2.rds
## Outputs: emc_RDM_CS_study2.RData (fitted EMC2 object `emc`),
##          Table1_RDM_CS_study2.html (APA parameter table, natural scale)
## Usage:   Rscript study2_rdm.R
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (RDM / EMC2 robustness
## analysis). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

rm(list = ls())

library(EMC2)
library(dplyr)
library(gt)

set.seed(123456)

## ---------------------------------------------------------------
## 1. Data: prepare CS/SC subset
## ---------------------------------------------------------------
## Study 2 stores RT directly in `rt` (seconds) — no ms→s conversion.

data_path <- file.path("..", "data", "final_df_study2.rds")
stopifnot(file.exists(data_path))
study2 <- readRDS(data_path)

dat_cs <- study2 %>%
  filter(test_part %in% c("cs", "sc")) %>%
  mutate(
    rt  = as.numeric(rt),                        # already in seconds
    cho = ifelse(true_response == "f", 1, -1),   # key to ±1

    # complex option is on the A side ('f') in cs, on the B side ('j') in sc
    oa_complex = ifelse(test_part == "cs", 1, -1),

    choosing_complex = ifelse(
      (oa_complex ==  1 & cho ==  1) |
      (oa_complex == -1 & cho == -1),
      1, 0
    ),

    # complex − simple re-signing of predictors
    EVD   = EV_diff   * oa_complex,
    SDD   = SD_diff   * oa_complex,
    SkewD = Skew_diff * oa_complex,

    R = factor(
      ifelse(choosing_complex == 1, "complex", "simple"),
      levels = c("simple", "complex")
    ),

    subjects = factor(subject)
  ) %>%
  select(subjects, rt, R, EVD, SDD, SkewD) %>%
  filter(is.finite(rt), rt > 0) %>%
  droplevels() %>%
  mutate(across(c(EVD, SDD, SkewD), ~ as.numeric(scale(.))))

# S mirrors R for EMC2's race plumbing
dat_cs$S <- dat_cs$R


## ---------------------------------------------------------------
## 2. Hierarchical RDM design
## ---------------------------------------------------------------

design_RDM_h <- design(
  data     = dat_cs,
  model    = RDM,
  matchfun = function(d) d$S == d$lR,
  formula  = list(
    v  ~ 0 + lR + lR:(EVD + SDD + SkewD),
    s  ~ 0 + lR,
    B  ~ 0 + lR,   # per-accumulator threshold: response caution for simple vs complex
    A  ~ 1,        # shared start-point range (bias not separately estimated)
    t0 ~ 1
  ),
  constants = c(s_lRsimple = 0)  # anchor: s_simple = 1
)

mapped_pars(design_RDM_h)


## ---------------------------------------------------------------
## 3. Priors
## ---------------------------------------------------------------

pr_h <- prior(design_RDM_h, type = "standard")

# 3.1 Drift intercepts & slopes (predictors z-scored)
idx_v <- grep("^v_", names(pr_h$par$mean))
pr_h$par$mean[idx_v] <- 0
pr_h$par$sd  [idx_v] <- 0.5

# 3.2 Within-trial noise — neutral
pr_h$par$mean["s_lRcomplex"] <- 0
pr_h$par$sd  ["s_lRcomplex"] <- 0.35

# 3.3 Boundaries — neutral, no a priori bias
if (all(c("B_lRsimple", "B_lRcomplex") %in% names(pr_h$par$mean))) {
  pr_h$par$mean[c("B_lRsimple", "B_lRcomplex")] <- log(1.0)
  pr_h$par$sd  [c("B_lRsimple", "B_lRcomplex")] <- 0.4
}

# 3.4 Start-range & non-decision time
pr_h$par$mean["A"]  <- log(0.30); pr_h$par$sd["A"]  <- 0.4
pr_h$par$mean["t0"] <- log(0.30); pr_h$par$sd["t0"] <- 0.3


## ---------------------------------------------------------------
## 4. Build EMC object and fit
## ---------------------------------------------------------------
## The fit is expensive; run once, then reload from the .RData file.

fit_file <- "emc_RDM_CS_study2.RData"

if (!file.exists(fit_file)) {
  emc_h <- make_emc(dat_cs, design_RDM_h, prior = pr_h, type = "hierarchical")

  # On sciCORE, SLURM allocates SLURM_CPUS_PER_TASK cores. We keep 4 chains
  # and distribute the remaining cores across chains.
  n_ch  <- 4
  total <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "16"))
  cpc   <- max(1L, total %/% n_ch)

  cat(sprintf("Fitting with %d chains x %d cores/chain (total=%d)\n",
              n_ch, cpc, n_ch * cpc))

  emc <- fit(
    emc_h,
    burnin          = 1500,
    samples         = 3000,
    thin            = 1,
    n.chains        = n_ch,
    cores_per_chain = cpc,
    fileName        = fit_file
  )
}

load(fit_file)  # restores object `emc`


## ---------------------------------------------------------------
## 5. Summary table (natural scale)
## ---------------------------------------------------------------

summary_df <- summary(emc)

mu_tbl <- summary_df$mu %>%
  as.data.frame() %>%
  tibble::rownames_to_column("parameter")

log_params <- c("A",
                "B_lRsimple", "B_lRcomplex",
                "s_lRcomplex", "t0")

mu_tbl <- mu_tbl %>%
  mutate(
    across(
      c(`2.5%`, `50%`, `97.5%`),
      ~ ifelse(parameter %in% log_params, exp(.x), .x)
    )
  )

mu_tbl_clean <- mu_tbl %>%
  mutate(across(`2.5%`:`97.5%`, ~ round(.x, 3))) %>%
  rename(CI_low = `2.5%`, Median = `50%`, CI_high = `97.5%`)

tbl_apa <- mu_tbl_clean %>%
  select(parameter, Median, CI_low, CI_high) %>%
  rename(
    Parameter        = parameter,
    `95% CI (Lower)` = CI_low,
    `95% CI (Upper)` = CI_high
  ) %>%
  gt() %>%
  tab_header(
    title = md("**Table 1**<br>Hierarchical RDM parameter estimates — Study 2 (natural scale)")
  ) %>%
  fmt_number(
    columns  = c(Median, `95% CI (Lower)`, `95% CI (Upper)`),
    decimals = 3
  ) %>%
  tab_style(
    style     = list(cell_text(weight = "bold")),
    locations = cells_column_labels(everything())
  ) %>%
  tab_options(
    table.font.names                  = "Times New Roman",
    table.font.size                   = 12,
    table.width                       = pct(90),
    data_row.padding                  = px(4),
    table.border.top.style            = "solid",
    table.border.top.width            = px(1),
    table.border.bottom.style         = "solid",
    table.border.bottom.width         = px(1),
    column_labels.border.top.width    = px(1),
    column_labels.border.bottom.width = px(1),
    row.striping.include_table_body   = FALSE
  ) %>%
  tab_source_note(
    source_note =
      "Note. Parameters originally on the log scale (A, B_lRsimple, B_lRcomplex, s_lRcomplex, t0) were exponentiated. CI = credible interval."
  )

gtsave(tbl_apa, "Table1_RDM_CS_study2.html")
