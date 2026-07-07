# ============================================================================
# data_cleaning_study3.R — raw -> analysis-ready data for Study 3 (compound probability)
#
# Aim:     Build the analysis-ready dataset from the raw experiment export,
#          applying the documented participant/trial exclusions.
# Inputs:  probability_complexity/data/merged_data.rds (raw merged jsPsych export;
#            falls back to merged_data.csv if the RDS is absent)
#          probability_complexity/data/test_trial_study3_raw.csv (stimulus table)
# Outputs: probability_complexity/data/derived/result_ca_study3.rds
#          probability_complexity/data/derived/exclusion_report_study3.rds
#          probability_complexity/data/derived/final_df_study3.rds
# Exclusions (in order):
#   1. Participants with > 6 comprehension-question (MCQ) attempts.
#   2. Trial RT hard bounds: rt < 1 s or rt > 30 s.
#   3. Trial RT band: outside median +/- 3 SD per subject x category.
#   4. Participants who lost > 50% of their test trials to the above.
# Usage:   Run from the 02_data_preparation/ directory:
#            Rscript probability_complexity/scripts/data_cleaning_study3.R
#          (the script self-anchors with here::i_am()).
# ----------------------------------------------------------------------------
# Part of the complexity-under-risk replication package (data preparation).
# Pipeline order and dependencies are documented in ../../README.md.
# ============================================================================

set.seed(1234)

# ---- packages ----
need <- c(
  "dplyr","readr","tibble","janitor","ggplot2","here","stringr","tidyr",
  "fs","purrr"
)

missing <- setdiff(need, rownames(installed.packages()))
if (length(missing) > 0) {
  stop(
    "Missing required packages: ", paste(missing, collapse = ", "),
    "\nInstall them and rerun."
  )
}
invisible(lapply(need, library, character.only = TRUE))

# ---- project root ----
here::i_am("probability_complexity/scripts/data_cleaning_study3.R")

# ---- study root + paths (avoid repeating folder names) ----
STUDY_DIR <- here::here("probability_complexity")
DATA_DIR  <- file.path(STUDY_DIR, "data")

IN_RDS      <- file.path(DATA_DIR, "merged_data.rds")
IN_STIM_CSV <- file.path(DATA_DIR, "test_trial_study3_raw.csv")

OUT_DIR <- file.path(DATA_DIR, "derived")
fs::dir_create(OUT_DIR)



load_merged <- function(rds_path = IN_RDS, csv_path = IN_CSV) {
  if (file.exists(rds_path)) {
    readRDS(rds_path)
  } else if (file.exists(csv_path)) {
    readr::read_csv(csv_path, show_col_types = FALSE)
  } else stop("Neither merged RDS nor CSV found.")
}

raw <- load_merged()

# =========================================================
# 2) Minimal sanity checks
# =========================================================
n_total <- nrow(raw)
n_subj  <- dplyr::n_distinct(raw$subject)
cat("Rows:", n_total, "\nUnique subjects:", n_subj, "\n")

# In the preregistration, we planned to collect data from 150 participants.
# However, due to a technical issue between the hosting platform (JATOS) and
# Prolific, three participants completed the experiment but did not submit their
# completion codes. As a result, their data were not recorded on the server.
# The final dataset therefore contains 147 participants, as reflected in the
# current dataframe.

# =========================================================
# CA summary output
# =========================================================
result_ca <- raw %>%
  group_by(subject) %>%
  dplyr::summarize(
    accuracy_BNT = mean(accuracy_BNT, na.rm = TRUE),
    accuracy_HMT = mean(accuracy_HMT, na.rm = TRUE)
  ) %>%
  mutate(
    ca_average = rowMeans(across(c(accuracy_BNT, accuracy_HMT)), na.rm = TRUE)
  )


saveRDS(result_ca, file.path(OUT_DIR, "result_ca_study3.rds"))

# =========================================================
# Exclusion (PREREGISTERED)
# =========================================================

# --- 1) MCQ attempts: exclude > 6 --------------------------------------------
attempts <- raw %>%
  group_by(subject) %>%
  summarise(
    total_question_attempts = sum(trial_type_label == "question", na.rm = TRUE),
    .groups = "drop"
  )

subjects_mcq_keep <- attempts %>%
  filter(total_question_attempts <= 6) %>%
  pull(subject)

# --- 2) Base data: test, non-practice, MCQ-passing only ----------------------
dat0 <- raw %>%
  filter(trial_type_label == "test", test_part != "prac") %>%
  filter(subject %in% subjects_mcq_keep) %>%
  mutate(
    category = case_when(
      test_part == "ss" ~ "simple vs. simple",
      test_part == "cc" ~ "complex vs. complex",
      test_part %in% c("cs","sc") ~ "simple vs. complex",
      TRUE ~ as.character(test_part)
    ),
    rt = rt_trial / 1000
  )

# --- 3) Trial-level RT exclusions --------------------------------------------

# 3a) Hard bounds
dat1 <- dat0 %>%
  mutate(outlier_hard = rt < 1 | rt > 30)

# 3b) Median ± 3 SD per subject × category (handle sd==0)
stats_sc <- dat1 %>%
  filter(!outlier_hard) %>%
  group_by(subject, category) %>%
  summarise(
    med_rt = median(rt, na.rm = TRUE),
    sd_rt  = sd(rt, na.rm = TRUE),
    .groups = "drop"
  )

dat2 <- dat1 %>%
  left_join(stats_sc, by = c("subject","category")) %>%
  mutate(
    sd_rt = ifelse(is.na(sd_rt), 0, sd_rt),
    outlier_band = ifelse(
      sd_rt > 0,
      rt < (med_rt - 3 * sd_rt) | rt > (med_rt + 3 * sd_rt),
      FALSE
    ),
    to_drop = outlier_hard | outlier_band
  )

# --- 4) 50% removal rule per subject -----------------------------------------

# counts pre- and post-exclusion (only test/non-prac)
n_pre  <- dat0 %>% count(subject, name = "n_pre")
n_post <- dat2 %>% filter(!to_drop) %>% count(subject, name = "n_post")

drop50 <- n_pre %>%
  left_join(n_post, by = "subject") %>%
  mutate(
    n_post = coalesce(n_post, 0L),
    prop_removed = ifelse(n_pre > 0, (n_pre - n_post) / n_pre, 1),
    remove_subject = prop_removed > 0.5
  )

subjects_keep_final <- drop50 %>%
  filter(!remove_subject) %>%
  pull(subject)

data_exclusion <- dat2 %>%
  filter(!to_drop, subject %in% subjects_keep_final) %>%
  select(-outlier_hard, -outlier_band, -to_drop, -med_rt, -sd_rt)

# --- 5) concise exclusion report ---------------------------------------------
exclusion_report <- attempts %>%
  mutate(mcq_excluded = total_question_attempts > 6) %>%
  left_join(n_pre, by = "subject") %>%
  left_join(n_post, by = "subject") %>%
  left_join(drop50 %>% select(subject, prop_removed, remove_subject), by = "subject") %>%
  mutate(
    n_pre = coalesce(n_pre, 0L),
    n_post = coalesce(n_post, 0L),
    prop_removed = coalesce(prop_removed, 0),
    excluded_final = mcq_excluded | remove_subject
  ) %>%
  arrange(desc(excluded_final), subject)


saveRDS(exclusion_report, file.path(OUT_DIR, "exclusion_report_study3.rds"))

# =========================================================
# Stimuli table + trial index parsing
# =========================================================
data_exclusion <- data_exclusion %>%
  mutate(
    trial_index_base = suppressWarnings(as.integer(stringr::str_extract(optionA_Stimulus, "\\d+")))
  )

if (!file.exists(IN_STIM_CSV)) stop("Stimuli CSV not found: ", IN_STIM_CSV)

raw_stimuli_data <- readr::read_csv(IN_STIM_CSV, show_col_types = FALSE) %>%
  select(-1) %>%
  mutate(
    evd_bins = cut(
      EV_diff,
      breaks = c(-Inf, -15, -9, 1, 11, 21),
      labels = c("-21 to -19", "-11 to -9", "-1 to 1", "9 to 11", "19 to 21"),
      right = TRUE, include.lowest = TRUE
    )
  )

# =========================================================
# Merge stimuli info into cleaned behavioral data
# =========================================================
final_df <- data_exclusion %>%
  full_join(
    raw_stimuli_data,
    by = c("skew" = "skew", "trial_index_base" = "trial_index_raw")
  ) %>%
  filter(trial_type_label == "test" & test_part != "prac") %>%
  mutate(
    skew_level = case_when(
      skew == "lr"      ~ "left_vs_right",
      skew == "rl"      ~ "right_vs_left",
      skew == "ns"      ~ "no_skew_vs_no_skew",
      skew == "control" ~ "catch",
      TRUE              ~ NA_character_
    ),
    response = purrr::map_chr(response, ~ as.character(.x)[1]),
    true_response = case_when(
      risk_index == 1  ~ response,
      risk_index == -1 ~ if_else(response == "f", "j", "f"),
      TRUE             ~ response
    )
  ) %>%
  select(
    subject, rt, test_part, skew_level, true_response,
    P_A1, O_A1, P_A2, O_A2, P_B1, O_B1, P_B2, O_B2,
    EV_diff, SD_diff, Skew_diff, evd_bins,
    complex_PA1, complex_PA2, complex_PA3, complex_PA4, complex_PA5, complex_PA6,
    complex_PB1, complex_PB2, complex_PB3, complex_PB4, complex_PB5, complex_PB6
  )

# =========================================================
# Save outputs
# =========================================================

saveRDS(final_df, file.path(OUT_DIR, "final_df_study3.rds"))