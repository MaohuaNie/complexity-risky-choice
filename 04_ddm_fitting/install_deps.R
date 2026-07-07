#!/usr/bin/env Rscript
## ============================================================================
## install_deps.R — one-shot dependency installer for the fitting pipeline
##
## Aim:     Install the CRAN packages the pipeline needs, install cmdstanr if
##          absent, and install the CmdStan backend if not already available.
##          Idempotent: already-installed components are skipped.
## Inputs:  Network access to CRAN / the Stan R-universe; optionally the CMDSTAN
##          environment variable if a CmdStan install already exists.
## Outputs: Installed R packages and (if needed) a CmdStan toolchain; nothing
##          written into the project tree.
## Usage:   Rscript install_deps.R
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM fitting stage).
## Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

cran_pkgs <- c("optparse", "dplyr", "tidyr", "tibble", "readr",
               "digest", "posterior")

to_install <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install)) {
  cat("Installing CRAN packages:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  cat("Installing cmdstanr from mc-stan.org ...\n")
  install.packages("cmdstanr",
    repos = c("https://stan-dev.r-universe.dev",
              "https://mc-stan.org/r-packages/",
              "https://cloud.r-project.org"))
}

if (Sys.getenv("CMDSTAN") == "" &&
    tryCatch(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)),
             error = function(e) TRUE)) {
  cat("Installing CmdStan backend ...\n")
  cmdstanr::install_cmdstan(cores = max(1, parallel::detectCores() - 1))
} else {
  cat("CmdStan already available:",
      tryCatch(cmdstanr::cmdstan_path(), error = function(e) "(env)"), "\n")
}

cat("\nAll dependencies OK.\n")
