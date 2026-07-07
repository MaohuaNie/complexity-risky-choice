## ============================================================================
## recovery_fit.R — refit one simulated dataset and extract true-vs-recovered
##
## Aim:     Fit one simulated dataset (from recovery_simulate.R) with the SAME
##          design (A ~ 1) and priors as study2_rdm.R, then extract recovered
##          estimates for ALL parameters at both the group and participant level
##          alongside the known truths. Step 2 of the RDM recovery pipeline.
## Inputs:  --sim <sim .rds from recovery_simulate.R>
## Outputs: --out <.RData fitted EMC2 object>,
##          --summary <.rds compact true-vs-recovered summary>
## Usage:   Rscript recovery_fit.R --sim recovery/sim_seed01.rds \
##              --out recovery/fit_seed01.RData --summary recovery/sum_seed01.rds
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (RDM / EMC2 robustness
## analysis). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

suppressPackageStartupMessages({
  library(optparse); library(EMC2); library(dplyr)
})

opt_list <- list(
  make_option("--sim",     type = "character", default = NULL),
  make_option("--out",     type = "character", default = "recovery_fit.RData"),
  make_option("--summary", type = "character", default = "recovery_summary.rds"),
  make_option("--burnin",  type = "integer",   default = 3000),
  make_option("--samples", type = "integer",   default = 6000),
  make_option("--chains",  type = "integer",   default = 4)
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$sim), file.exists(opt$sim))

sim <- readRDS(opt$sim)
dat <- sim$data
cat(sprintf("Loaded sim: %d subjects, %d trials\n",
            length(unique(dat$subjects)), nrow(dat)))

## ---------------------------------------------------------------
## Design (A ~ 1) + priors — identical to the real Study-2 analysis
## ---------------------------------------------------------------

design_RDM_h <- design(
  data     = dat,
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

pr_h <- prior(design_RDM_h, type = "standard")
idx_v <- grep("^v_", names(pr_h$par$mean))
pr_h$par$mean[idx_v] <- 0;               pr_h$par$sd[idx_v] <- 0.5
pr_h$par$mean["s_lRcomplex"] <- 0;       pr_h$par$sd["s_lRcomplex"] <- 0.35
pr_h$par$mean[c("B_lRsimple","B_lRcomplex")] <- log(1.0)
pr_h$par$sd  [c("B_lRsimple","B_lRcomplex")] <- 0.4
pr_h$par$mean["A"]  <- log(0.30);        pr_h$par$sd["A"]  <- 0.4
pr_h$par$mean["t0"] <- log(0.30);        pr_h$par$sd["t0"] <- 0.3

## ---------------------------------------------------------------
## Fit
## ---------------------------------------------------------------

emc_h <- make_emc(dat, design_RDM_h, prior = pr_h, type = "hierarchical")
n_ch  <- opt$chains
total <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "16"))
cpc   <- max(1L, total %/% n_ch)
cat(sprintf("Fitting %d chains x %d cores/chain; burnin %d, samples %d\n",
            n_ch, cpc, opt$burnin, opt$samples))

emc <- fit(emc_h, burnin = opt$burnin, samples = opt$samples, thin = 1,
           n.chains = n_ch, cores_per_chain = cpc, fileName = opt$out)
load(opt$out)

## ---------------------------------------------------------------
## Group-level recovered summary (median, 95% CrI, Rhat, ESS)
## ---------------------------------------------------------------

mu_summary <- summary(emc)$mu   # [param x {2.5%,50%,97.5%,Rhat,ESS}]
pars <- sim$par_names

group <- data.frame(
  parameter = pars,
  true      = as.numeric(sim$true_mu[pars]),
  rec_lo    = mu_summary[pars, "2.5%"],
  rec_md    = mu_summary[pars, "50%"],
  rec_hi    = mu_summary[pars, "97.5%"],
  rhat      = mu_summary[pars, "Rhat"],
  ess       = mu_summary[pars, "ESS"],
  row.names = NULL
)

## ---------------------------------------------------------------
## Participant-level recovered estimates (median + 95% CrI per subject)
## Robust extraction: try credint(selection="alpha"), else raw alpha draws.
## ---------------------------------------------------------------

subj_ids <- levels(dat$subjects)

extract_alpha <- function(emc) {
  ## summary(emc, selection = "alpha") returns a list keyed by SUBJECT
  ## (sim001, ...); each element is a matrix with rows = parameters and
  ## columns = {2.5%, 50%, 97.5%, Rhat, ESS}. (Confirmed on this EMC2 build.)
  ## capture.output() suppresses summary()'s verbose per-subject printing.
  sa <- NULL
  invisible(capture.output(sa <- summary(emc, selection = "alpha")))
  stopifnot(is.list(sa), length(sa) > 0)
  do.call(rbind, lapply(names(sa), function(sj) {
    m <- sa[[sj]]
    data.frame(subject   = sj,
               parameter = rownames(m),
               rec_lo    = m[, "2.5%"],
               rec_md    = m[, "50%"],
               rec_hi    = m[, "97.5%"],
               row.names = NULL)
  }))
}

subj_rec <- extract_alpha(emc)

## Attach the per-participant TRUTHS (sim$per_subject: subjects x params)
truth_long <- as.data.frame(as.table(sim$per_subject))
names(truth_long) <- c("subject", "parameter", "true")
truth_long$subject   <- as.character(truth_long$subject)
truth_long$parameter <- as.character(truth_long$parameter)
subj_rec$subject     <- as.character(subj_rec$subject)

subject <- merge(truth_long, subj_rec, by = c("subject", "parameter"))

cat(sprintf("\nParticipant-level rows recovered: %d (expect %d = %d subj x %d par)\n",
            nrow(subject), length(subj_ids) * length(pars),
            length(subj_ids), length(pars)))

## ---------------------------------------------------------------
## Save compact summary
## ---------------------------------------------------------------

recovery <- list(
  seed    = sim$seed,
  group   = group,       # group-level true vs recovered (+ Rhat/ESS)
  subject = subject,     # participant-level true vs recovered
  par_names = pars
)
saveRDS(recovery, opt$summary)
cat("\nWrote:", opt$out, "and", opt$summary, "\n")
cat("\nGroup-level recovery:\n"); print(group, digits = 3)
cat(sprintf("\nmax Rhat = %.3f, min ESS = %.0f\n", max(group$rhat), min(group$ess)))
