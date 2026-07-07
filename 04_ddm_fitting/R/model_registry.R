## ============================================================================
## R/model_registry.R — catalogue of Stan models this pipeline fits
##
## Aim:     Build MODEL_REGISTRY, one row per (family, model) that can be fit,
##          and provide verify_registry() to check the referenced Stan files
##          exist. Six families are enumerated (cpt_ccss, cpt_ccss_7o, mv_ccss,
##          cpt_cs, cpt_cs_7o, mv_cs); CCSS families cover baseline + all 15
##          non-empty subsets of {n,r,a,s}, CS families cover baseline + subsets
##          of {sp,dr}. Each row records n_params (group-level parameter count),
##          which drives init generation.
## Inputs:  Nothing at load time beyond dplyr/tibble; verify_registry() reads
##          the stan/ directory listing.
## Outputs: MODEL_REGISTRY (in-memory tibble) and verify_registry(); no files
##          written.
## Usage:   Not run directly; source()'d by run_fit.R / run_all_fits.R.
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (DDM fitting stage).
## Pipeline order and dependencies are documented in ../README.md.
## ============================================================================
##
## Row fields:
##   family    : cpt_ccss | cpt_ccss_7o | mv_ccss | cpt_cs | cpt_cs_7o | mv_cs
##   model     : suffix identifying the variant (e.g. "baseline", "n_r_a_s")
##   stan_file : path under stan/ (conventionally {family}_{model}.stan)
##   n_params  : number of group-level parameters the Stan model declares
##               (= length of mu/sigma vectors = rows of `participant_params`)
##
## To add a model, drop its Stan file into stan/ and append a row here; to
## exclude one, remove/omit its row (run_all_fits.R will simply not fit it).

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

## --- CCSS families: baseline + all 15 non-empty subsets of {n,r,a,s} --------
ccss_delta_models <- local({
  letters <- c("n", "r", "a", "s")
  combos <- unlist(lapply(1:4, function(k)
    apply(combn(letters, k), 2, paste, collapse = "_")), use.names = FALSE)
  c("baseline", combos)
})

## --- CS families: baseline + 3 non-empty subsets of {sp, dr} ---------------
## `skew` is intentionally excluded from the fitting pipeline — the skew-variant
## Stan files (*_skew, *_sp_skew, *_dr_skew, *_sp_dr_skew) remain on disk and
## can be re-added by extending `letters` below if we later want to re-examine
## the delta_gamma / delta_eta mechanism.
cs_extra_models <- local({
  letters <- c("sp", "dr")
  combos <- unlist(lapply(seq_along(letters), function(k)
    apply(combn(letters, k), 2, paste, collapse = "_")), use.names = FALSE)
  c("baseline", combos)
})

## --- n_params per family/model (base + number of active extras) -------------
## CCSS (CPT, MV):  base = 5 (beta/theta/threshold/ndt, gamma or eta)
##                  + k deltas (one per letter in model name)
## CS   (CPT, MV):  base = 5 (same)
##                  + sp, dr each add 1 param
ccss_n_params <- function(model) {
  if (model == "baseline") 5L else 5L + length(strsplit(model, "_")[[1]])
}
cs_n_params <- function(model) {
  if (model == "baseline") 5L else 5L + length(strsplit(model, "_")[[1]])
}

MODEL_REGISTRY <- bind_rows(
  tibble(family = "cpt_ccss", model = ccss_delta_models,
         n_params = vapply(ccss_delta_models, ccss_n_params, integer(1))),
  tibble(family = "cpt_ccss_7o", model = ccss_delta_models,
         n_params = vapply(ccss_delta_models, ccss_n_params, integer(1))),
  tibble(family = "mv_ccss",  model = ccss_delta_models,
         n_params = vapply(ccss_delta_models, ccss_n_params, integer(1))),
  tibble(family = "cpt_cs",   model = cs_extra_models,
         n_params = vapply(cs_extra_models, cs_n_params, integer(1))),
  tibble(family = "cpt_cs_7o", model = cs_extra_models,
         n_params = vapply(cs_extra_models, cs_n_params, integer(1))),
  tibble(family = "mv_cs",    model = cs_extra_models,
         n_params = vapply(cs_extra_models, cs_n_params, integer(1)))
) |> mutate(
  stan_file = sprintf("stan/%s_%s.stan", family, model)
)

## ---------------------------------------------------------------------------
## Verify the Stan files listed above actually exist under stan/.
## Called once at startup; errors loudly on typos.
## ---------------------------------------------------------------------------
verify_registry <- function(stan_dir = "stan") {
  missing <- MODEL_REGISTRY$stan_file[!file.exists(MODEL_REGISTRY$stan_file)]
  if (length(missing) > 0) {
    stop("Stan files missing for registered models:\n  ",
         paste(missing, collapse = "\n  "))
  }
  invisible(TRUE)
}
