#!/usr/bin/env Rscript
## ============================================================================
## install_deps.R — one-shot installer for the recovery pipeline's R deps
##
## Aim:     Install every CRAN package the recovery pipeline needs, plus
##          cmdstanr and (if absent) the CmdStan backend. Run once before the
##          first `Rscript run_recovery.R`.
## Inputs:  none (queries installed.packages() and the CMDSTAN env var).
## Outputs: none on disk beyond installed packages / CmdStan toolchain.
## Usage:   Rscript install_deps.R
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM parameter
## recovery). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

cran_pkgs <- c("optparse", "dplyr", "tidyr", "tibble", "readr",
               "ggplot2", "posterior", "future", "future.apply", "rtdists")

to_install <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install)) {
  cat("Installing CRAN packages:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  cat("Installing cmdstanr from mc-stan.org...\n")
  install.packages("cmdstanr",
    repos = c("https://mc-stan.org/r-packages/", "https://cloud.r-project.org"))
}

## Install CmdStan backend if no CMDSTAN env var is set and no existing install.
if (Sys.getenv("CMDSTAN") == "" &&
    tryCatch(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)),
             error = function(e) TRUE)) {
  cat("Installing CmdStan via cmdstanr::install_cmdstan()...\n")
  cmdstanr::install_cmdstan(cores = max(1, parallel::detectCores() - 1))
} else {
  cat("CmdStan already available:",
      tryCatch(cmdstanr::cmdstan_path(), error = function(e) "(env)"), "\n")
}

cat("\nAll dependencies OK.\n")
