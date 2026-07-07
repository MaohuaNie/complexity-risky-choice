## ============================================================================
## apa_theme.R — shared APA-7-style ggplot2 theme + parameter labels
##
## Aim:     Provide theme_apa() and the PARAM_LABELS lookup so recovery plots
##          use a consistent APA look and descriptive parameter names. Sourced
##          by recovery_aggregate.R (helper, not run directly).
## Inputs:  none
## Outputs: none (defines theme_apa() and PARAM_LABELS in the calling session)
## Usage:   source("apa_theme.R")
## ----------------------------------------------------------------------------
## Part of the complexity-under-risk replication package (RDM / EMC2 robustness
## analysis). Pipeline order and dependencies are documented in ../README.md.
## ============================================================================

suppressMessages(library(ggplot2))

theme_apa <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title           = element_text(face = "plain",
                                          hjust = 0,
                                          size  = base_size * 1.05,
                                          margin = margin(b = 6)),
      plot.subtitle        = element_text(hjust = 0,
                                          size  = base_size * 0.95,
                                          colour = "grey25",
                                          margin = margin(b = 8)),
      axis.title.x         = element_text(margin = margin(t = 8)),
      axis.title.y         = element_text(margin = margin(r = 8)),
      axis.text            = element_text(colour = "black"),
      axis.ticks           = element_line(colour = "black", linewidth = 0.35),
      axis.line            = element_line(colour = "black", linewidth = 0.35),
      panel.grid           = element_blank(),
      strip.background     = element_blank(),
      strip.text           = element_text(face = "plain",
                                          size = base_size * 0.95,
                                          margin = margin(b = 4, t = 2)),
      legend.background    = element_blank(),
      legend.key           = element_blank(),
      legend.title         = element_text(size = base_size * 0.9),
      legend.text          = element_text(size = base_size * 0.85),
      plot.caption         = element_text(size  = base_size * 0.8,
                                          colour = "grey25",
                                          hjust = 0,
                                          margin = margin(t = 8))
    )
}

## Descriptive labels for the model's parameters — always use these
## instead of the raw names (s_lRcomplex, EVD, etc.) in plots/tables.
PARAM_LABELS <- c(
  noise_ratio        = "Within-trial noise ratio (complex \u00f7 simple)",
  log_noise_ratio    = "Log of within-trial noise ratio",
  threshold_complex  = "Decision threshold for complex option",
  threshold_simple   = "Decision threshold for simple option",
  drift_complex      = "Drift rate toward complex option",
  drift_simple       = "Drift rate toward simple option",
  EVD                = "Expected-value advantage of complex",
  SDD                = "Risk advantage of complex (SD difference)",
  SkewD              = "Skewness advantage of complex"
)
