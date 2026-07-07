#!/usr/bin/env Rscript
## ============================================================================
## make_rtdist_one.R — empirical-vs-predicted RT-distribution PPC for ONE model
##
## Aim:     Render a single model's RT-distribution figure (observed histogram
##          vs predicted density, split by chosen option / complexity). Logic is
##          the same as plot_ppc.R::plot_rt_distribution but touches ONLY this
##          one plot, so it can be run for the best model without overwriting
##          the family-level figures.
## Inputs:  results/<study>/<family>/<model>/posterior_predictives.csv.
## Outputs: results/<study>/<family>/ppc_rt_distribution_<model>.pdf/.png
##          (or --out).
## Usage:   Rscript 04_ppc/make_rtdist_one.R --study study2 --family mv_ccss \
##            --model n_a_s --cond ccss
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (posterior predictive
## checks). Pipeline order and dependencies are documented in ../../README.md.
## ============================================================================
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(ggplot2); library(readr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--study",       type = "character"),
  make_option("--family",      type = "character"),
  make_option("--model",       type = "character"),
  make_option("--cond",        type = "character", help = "ccss | cs"),
  make_option("--results_dir", type = "character", default = "results"),
  make_option("--title",       type = "character", default = NULL),
  make_option("--out",         type = "character", default = NULL)
)))
stopifnot(opt$cond %in% c("ccss", "cs"))
is_ccss <- opt$cond == "ccss"

fam_dir  <- file.path(opt$results_dir, opt$study, opt$family)
ppc_path <- file.path(fam_dir, opt$model, "posterior_predictives.csv")
if (!file.exists(ppc_path)) stop("Not found: ", ppc_path,
                                 "\n  (run generate_ppc.R for this model first)")
cat("Reading", ppc_path, "...\n")
ppc <- read_csv(ppc_path, show_col_types = FALSE)

if (is_ccss) {
  obs  <- ppc %>% distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
    mutate(Condition = ifelse(observed_choice == 1, "Choosing A", "Choosing B"))
  pred <- ppc %>% mutate(Condition = ifelse(pred_choice == 1, "Choosing A", "Choosing B"))
} else {
  ## Recover which side was complex per trial, then map predicted choice to complex/simple
  ppc <- ppc %>% mutate(
    oa_complex         = as.integer(sign(chose_complex_obs * observed_choice)),
    chose_complex_pred = as.integer(oa_complex * pred_choice)
  )
  obs  <- ppc %>% distinct(participant, trial_in_participant, .keep_all = TRUE) %>%
    mutate(Condition = ifelse(chose_complex_obs == 1, "Choosing Complex", "Choosing Simple"))
  pred <- ppc %>% mutate(Condition = ifelse(chose_complex_pred == 1, "Choosing Complex", "Choosing Simple"))
}

theme_apa <- function(bs = 11) theme_classic(base_size = bs, base_family = "sans") +
  theme(text = element_text(colour = "black"), strip.background = element_blank(),
        panel.grid = element_blank())

ttl <- if (!is.null(opt$title)) opt$title else
  sprintf("Empirical vs. Predicted RT Distributions (%s / %s)", opt$family, opt$model)

p <- ggplot() +
  geom_histogram(data = obs, aes(x = observed_rt, y = after_stat(density), fill = Condition),
                 alpha = 0.4, bins = 50, position = "identity") +
  geom_density(data = pred, aes(x = pred_rt, colour = Condition), linewidth = 0.8) +
  scale_fill_manual(values = c("#FF9999", "#66CCCC")) +
  scale_colour_manual(values = c("#CC3333", "#009999")) +
  coord_cartesian(xlim = c(0, 20)) +
  labs(title = ttl, x = "RT", y = "Density") +
  theme_apa()

out_pdf <- if (!is.null(opt$out)) opt$out else
  file.path(fam_dir, paste0("ppc_rt_distribution_", opt$model, ".pdf"))
ggsave(out_pdf, p, width = 7, height = 5)
ggsave(sub("\\.pdf$", ".png", out_pdf), p, width = 7, height = 5, dpi = 300)
cat("Wrote", out_pdf, "\n")
