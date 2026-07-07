# Stimulus tables (lottery definitions)

The exact lottery stimuli used in the three experiments — one file per study.
Each row is one test trial: a choice between **option A** and **option B**, each
defined by a **simple** (2-outcome) form and a **complex** form, where "complex"
implements that study's manipulation. The two forms of an option are matched on
expected value; the study manipulates *how the option is displayed/computed*, not
its economic value.

| File | Study | Complexity manipulation | Trials |
|---|---|---|---|
| `stimuli_study1_lottery_complexity.csv` | 1 | complex = **7 outcomes** vs. simple = 2 outcomes | 49 |
| `stimuli_study2_outcome_complexity.csv`   | 2 | complex = each outcome shown as a **5-term arithmetic expression** | 49 |
| `stimuli_study3_probability_complexity.csv`     | 3 | complex = **6-component compound probability** (multi-stage) | 212 |

## Columns common to all three

- `trial_index_raw` — trial number.
- `skew_level` / `skew` — the trial's skewness pairing (e.g. `left_vs_right`,
  `right_vs_left`, `no_skew_vs_no_skew`).
- **Simple (2-outcome) form** — for option A: `O_A1`, `P_A1`, `O_A2`, `P_A2`
  (outcomes and their probabilities); for option B: `O_B1`, `P_B1`, `O_B2`, `P_B2`.
- **Moments and moment differences** — `EV_A`/`EV_B` (expected value),
  `Var_A`/`Var_B`, `SD_A`/`SD_B`, `Skew_A`/`Skew_B`, and the A−B differences
  `EV_diff`, `SD_diff`, `Skew_diff`. These differences are the raw quantities that
  become the MVS drift inputs (EVD, SDD, SkewD) after the complex-vs-simple
  re-signing applied during preprocessing.

## Study-specific complex-form columns

- **Study 1** (`lottery_complexity`): `complex_OA1`–`complex_OA7` /
  `complex_PA1`–`complex_PA7` give the 7-outcome version of option A (and
  `complex_OB*`/`complex_PB*` for option B). The 7 outcomes preserve the
  distribution of the 2-outcome form.
- **Study 2** (`outcome_complexity`): `complex_outcome_A1_1`–`complex_outcome_A1_5`
  are the five arithmetic terms whose sum equals `O_A1`; `complex_outcome_A2_*`
  the five terms summing to `O_A2` (and the `B1_*`/`B2_*` sets for option B).
  The gamble is identical to the simple form; only the outcomes are displayed as
  sums the participant must compute.
- **Study 3** (`probability_complexity`): `complex_PA1`–`complex_PA6` /
  `complex_PB1`–`complex_PB6` are the compound-probability components that combine
  (multi-stage) to the simple-form outcome probabilities.

These tables are the stimulus set behind the online experiment
(`../../01_experiment/`); the per-participant trial data (with responses and RTs)
are in the derived `../final_df_study*.rds`.
