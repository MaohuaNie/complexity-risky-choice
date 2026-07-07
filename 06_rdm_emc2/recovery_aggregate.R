## ============================================================================
## recovery_aggregate.R — pool recovery replicates into a table and scatters
##
## Aim:     Read all per-replicate recovery summaries and produce, for every
##          parameter: individual-level Pearson r + 95% CrI coverage (pooled over
##          all simulated participants across datasets), group-level r + coverage
##          (across datasets), and a convergence summary (max Rhat, min ESS).
##          Final step of the RDM recovery pipeline.
## Inputs:  recovery/sum_seed*.rds (from recovery_fit.R); apa_theme.R if present
## Outputs: recovery/recovery_param_table.csv,
##          recovery/recovery_{individual,group}.{pdf,png}
## Usage:   Rscript recovery_aggregate.R
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (RDM / EMC2 robustness
## analysis). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(readr); library(tidyr)
})
if (file.exists("apa_theme.R")) source("apa_theme.R")

files <- list.files("recovery", "^sum_seed.*\\.rds$", full.names = TRUE)
stopifnot(length(files) > 0)
cat(sprintf("Aggregating %d recovery datasets.\n", length(files)))

sums <- lapply(files, readRDS)

## ---- collect group- and subject-level long tables ----
group_all <- bind_rows(lapply(sums, function(s)
  cbind(seed = s$seed, s$group)))
subject_all <- bind_rows(lapply(sums, function(s)
  cbind(seed = s$seed, s$subject)))

par_levels <- sums[[1]]$par_names

## ---- convergence across datasets ----
conv <- group_all %>% group_by(seed) %>%
  summarise(max_rhat = max(rhat), min_ess = min(ess), .groups = "drop")
cat("\n=== Convergence across datasets ===\n")
cat(sprintf("max Rhat: range [%.3f, %.3f]; datasets with maxRhat>1.05: %d/%d; >1.1: %d\n",
            min(conv$max_rhat), max(conv$max_rhat),
            sum(conv$max_rhat > 1.05), nrow(conv), sum(conv$max_rhat > 1.1)))
cat(sprintf("min ESS : range [%.0f, %.0f]; datasets with minESS<400: %d\n",
            min(conv$min_ess), max(conv$min_ess), sum(conv$min_ess < 400)))

## ---- per-parameter recovery metrics ----
metric <- function(df) {
  df %>% summarise(
    r        = cor(true, rec_md),
    coverage = mean(rec_lo <= true & true <= rec_hi),
    bias     = mean(rec_md - true),
    n        = n(),
    .groups  = "drop")
}

ind <- subject_all %>% group_by(parameter) %>% metric() %>%
  rename(r_ind = r, cov_ind = coverage, bias_ind = bias, n_ind = n)
grp <- group_all %>% group_by(parameter) %>% metric() %>%
  rename(r_grp = r, cov_grp = coverage, bias_grp = bias, n_grp = n)

tab <- left_join(ind, grp, by = "parameter") %>%
  mutate(parameter = factor(parameter, levels = par_levels)) %>%
  arrange(parameter) %>%
  transmute(parameter,
            r_individual    = round(r_ind, 2),
            cov_individual  = round(100 * cov_ind, 0),
            bias_individual = round(bias_ind, 3),
            r_group         = round(r_grp, 2),
            cov_group       = round(100 * cov_grp, 0))

write_csv(tab, "recovery/recovery_param_table.csv")
cat("\n=== Per-parameter recovery ===\n")
print(as.data.frame(tab), row.names = FALSE)

## ---- faceted scatter plots ----
scatter <- function(df, title) {
  df <- df %>% mutate(parameter = factor(parameter, levels = par_levels),
                      miss = !(rec_lo <= true & true <= rec_hi))
  ## per-parameter Pearson r, shown in each facet (APA: no leading zero)
  r_lab <- df %>% group_by(parameter) %>%
    summarise(r = cor(true, rec_md), .groups = "drop") %>%
    mutate(label = paste0("r = ", sub("^(-?)0\\.", "\\1.", sprintf("%.2f", r))))
  ggplot(df, aes(true, rec_md)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(ymin = rec_lo, ymax = rec_hi),
                   colour = "grey75", linewidth = 0.25, alpha = 0.5) +
    geom_point(aes(colour = miss), size = 0.7) +
    geom_text(data = r_lab, aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.12, vjust = 1.4,
              size = 3, fontface = "bold") +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "red"), guide = "none") +
    facet_wrap(~ parameter, scales = "free", ncol = 4) +
    labs(title = title, x = "True (sampled scale)", y = "Recovered (median + 95% CrI)") +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank())
}

p_ind <- scatter(subject_all, "Participant-level recovery (all datasets pooled)")
p_grp <- scatter(group_all,   "Group-level recovery (one point per dataset)")

ggsave("recovery/recovery_individual.pdf", p_ind, width = 10, height = 8)
ggsave("recovery/recovery_individual.png", p_ind, width = 10, height = 8, dpi = 200)
ggsave("recovery/recovery_group.pdf", p_grp, width = 10, height = 8)
ggsave("recovery/recovery_group.png", p_grp, width = 10, height = 8, dpi = 200)

cat("\nWrote: recovery/recovery_param_table.csv and recovery_{individual,group}.{pdf,png}\n")
