# Data preparation (raw → analysis-ready)

The full, transparent path from the **raw experiment export** to the
**analysis-ready** data used by every downstream stage, including all
participant and trial exclusions. One self-contained folder per study.

```
02_data_preparation/
├── lottery_complexity/   (Study 1)   ├── outcome_complexity/ (Study 2)   ├── probability_complexity/ (Study 3)
│   ├── data/
│   │   ├── merged_data.rds          ← RAW merged jsPsych export (all trials, all participants)
│   │   ├── test_trial_study*_raw.csv← the stimulus table (see ../../data/stimuli/)
│   │   └── derived/
│   │       ├── final_df_study*.rds        ← analysis-ready data (== ../../data/)
│   │       ├── result_ca_study*.rds       ← cognitive-ability scores
│   │       └── exclusion_report_study*.rds← per-participant exclusion record (who + why)
│   └── scripts/
│       └── data_cleaning_study*.R   ← the cleaning + exclusion pipeline
```

## Running it

Each `scripts/data_cleaning_study*.R` self-anchors with `here::i_am(...)`, so from
the `02_data_preparation/` directory:

```r
Rscript lottery_complexity/scripts/data_cleaning_study1.R
Rscript outcome_complexity/scripts/data_cleaning_study2.R
Rscript probability_complexity/scripts/data_cleaning_study3.R
```

Each reads its `data/merged_data.rds` + stimulus CSV and writes the three
`data/derived/` files. The `final_df_study*.rds` it produces are the same files
provided at the package root in `../data/` (which the analysis stages read).

Needs: `dplyr, tidyr, readr, tibble, janitor, stringr, purrr, fs, here, ggplot2`.

## Exclusions (documented in the scripts, recorded in `exclusion_report_*`)

Applied in each `data_cleaning_study*.R`:

1. **Comprehension (MCQ) failures** — exclude participants exceeding the
   allowed number of comprehension-question attempts (**> 8** in Studies 1–2,
   **> 6** in Study 3).
2. **RT hard bounds** — remove trials with RT **< 1 s or > 30 s**. (This is why
   the analysis-ready data have a minimum RT of exactly 1 s.)
3. **RT outliers** — remove trials outside the **per-subject × choice-category
   median ± 3 SD**.
4. **Excessive trial loss** — exclude any participant who lost **> 50%** of
   their trials to the RT screens (2–3).

Resulting sample sizes (these match the Ns reported in the paper):

| Study | Complexity type | Raw N | Comprehension | > 50% trials lost | **Final N** |
|---|---|---|---|---|---|
| 1 (`lottery_complexity`) | lottery complexity | 148 | 1 | 10 | **137** |
| 2 (`outcome_complexity`)   | outcome complexity | 147 | 1 | 8  | **138** |
| 3 (`probability_complexity`)     | probability complexity | 150 | 8 | 10 | **132** |

(Raw N is 148/147 rather than the recruited 150 for Studies 1–2 because 2 and 3
participants' data were lost in the Prolific→JATOS transfer — see the supplement.)

The exclusion criteria and counts are also reported in the main text / appendix;
this stage provides the raw data and code so the full pipeline is reproducible.

See ../README.md for the overall pipeline.
