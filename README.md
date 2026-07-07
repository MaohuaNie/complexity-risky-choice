# The Impact of Complexity on Risky Choice — Replication Package

Analysis code and data for the three-study project modelling how **complexity**
affects choice and response time in decisions under risk, using an
**evidence-accumulation (drift-diffusion) framework**. Starting from the raw
experiment export (with a fully documented cleaning/exclusion pipeline), the
package reproduces the model-free behavioural results, the hierarchical Bayesian
DDM fits (MVS main model; CPT robustness check), the posterior predictive checks,
the DDM parameter recovery, and the race-model (RDM) robustness analysis.

Three complexity manipulations, one per study (study labels follow the paper):
- **Study 1 — lottery complexity** (number-of-outcomes manipulation; folder `lottery_complexity`)
- **Study 2 — outcome complexity** (arithmetic-expression outcomes; folder `outcome_complexity`)
- **Study 3 — probability complexity** (compound probabilities; folder `probability_complexity`)

Each study has two conditions: **CS** (options differ in complexity) and **CCSS**
(options equally complex vs. equally simple).

## Repository layout

```
complexity_replication/
├── README.md              ← this file (start here)
├── ENVIRONMENT.md         ← R / CmdStan / EMC2 versions and how to install deps
├── data/                  ← derived data + raw stimulus tables (see "Data" below)
│   └── stimuli/           ← lottery definitions per study (CSV) + column dictionary
├── 01_experiment/         ← the online experiment (jsPsych/JATOS) that collected the data
├── 02_data_preparation/   ← RAW data → documented exclusions → analysis-ready (per study)
├── 03_behavioural/        ← model-free behavioural analysis + figures (an R notebook)
├── 04_ddm_fitting/        ← hierarchical Bayesian DDM fitting + posterior predictive checks
│   ├── R/                 ← pipeline helpers (design, priors, preprocess, transforms, fit)
│   ├── stan/              ← the 72 Stan models (+ a README explaining the naming)
│   ├── run_fit.R          ← fit ONE model on ONE study
│   ├── run_all_fits.R     ← driver for the full model set
│   ├── post/              ← posterior predictive checks (generate_ppc + plot_ppc*)
│   └── scripts/           ← SLURM submitters (fitting + PPC) — cluster
├── 05_ddm_recovery/       ← simulation-based DDM parameter recovery
├── 06_rdm_emc2/           ← Racing Diffusion Model robustness analysis (EMC2)
└── supplement/            ← LaTeX source for the online supplemental material
```

Fitting and its posterior predictive checks live in the **same** stage
(`04_ddm_fitting/`) because they share the `R/` helpers, the Stan models, and the
`results/` output tree. Run every command in that stage from the
`04_ddm_fitting/` directory (e.g. `Rscript post/generate_ppc.R …`).

Every script begins with a header block stating its **aim, inputs, outputs, and
how to run it**. Folder-level details are in each stage's own `README`.

## Data

`data/` contains the **derived, model-ready** data only:

| File | Study | Complexity type (manipulation) |
|---|---|---|
| `final_df_study1.rds` | Study 1 | lottery complexity (number of outcomes) |
| `final_df_study2.rds` | Study 2 | outcome complexity (arithmetic expressions) |
| `final_df_study3.rds` | Study 3 | probability complexity (compound probabilities) |
| `result_ca_study{1,2,3}.rds` | — | per-participant cognitive-ability scores |

These RDS files are already pre-processed (skew recoding, response flip, RT in
seconds, test-trial filtering). The **raw experiment export** and the exact
**cleaning/exclusion pipeline** that produces them are in `02_data_preparation/`
(one self-contained folder per study, with per-participant exclusion reports).

`data/stimuli/` holds the **raw stimulus tables** — the exact lottery definitions
used in each study (the 2-outcome simple form, the study-specific complex form,
and the computed moments), one CSV per study. See `data/stimuli/README.md` for
the column dictionary.

The fitting and recovery scripts default to a stage-local `data/` folder, so
`04_ddm_fitting/data` and `05_ddm_recovery/data` are **symlinks to the shared
`../data/`**. If your download tool does not preserve symlinks, recreate them
(`ln -s ../data 04_ddm_fitting/data`) or copy `data/` into those stages.

Note: the behavioural notebook (`03_behavioural/`) still uses the original
manipulation labels internally as its dataset keys and output-file names
(`number_of_outcome`, `outcome-as-term`, `compound_prob`), which map to
Studies 1, 2, and 3 respectively — i.e. to lottery, outcome, and probability
complexity (this mapping is also given in the notebook header).

## How to reproduce

First install dependencies — see **ENVIRONMENT.md**.

The stages are ordered; folders are numbered in run order.

0. **Online experiment** (`01_experiment/`) — the jsPsych/JATOS experiment used
   to collect the data. Documentation only (needs its JATOS assets to run); it
   is the origin of the raw data.

1. **Data preparation** (`02_data_preparation/`, laptop) — from each study
   folder, `scripts/data_cleaning_study*.R` reads the raw `merged_data.rds`,
   applies the documented exclusions (comprehension failures; RT outliers; >50%
   trial loss), and writes the analysis-ready `final_df_study*.rds`,
   `result_ca_study*.rds`, and an `exclusion_report_study*.rds`. The derived
   files are the ones provided at the package root in `data/`.

2. **Behavioural results** (`03_behavioural/`, laptop) — knit
   `behavioural_analysis.Rmd`. Point its `study_paths` at `../data/final_df_study*.rds`.
   Produces the combined choice + RT figures per study.

3. **DDM fitting** (`04_ddm_fitting/`, cluster) — the core, computationally heavy
   stage. `run_fit.R` fits one `(study, family, model)`; `run_all_fits.R` /
   `scripts/submit_fit_array.sh` run the full set as a SLURM array. The model set
   and per-study inclusion rules are defined in `R/model_registry.R`. Fitted
   output is written under `results/<study>/<family>/<model>/`. See
   `stan/README.md` for the model-naming key.

4. **Posterior predictive checks** (`04_ddm_fitting/post/`, cluster then laptop) —
   after fitting, `post/generate_ppc.R` simulates posterior-predictive
   choices/RTs (large `posterior_predictives.csv` per model); the
   `post/plot_ppc*.R` scripts and `post/make_rtdist_one.R` build the PPC and
   RT-distribution figures. Run these from the `04_ddm_fitting/` directory so
   they find `R/` and `results/`.

5. **DDM parameter recovery** (`05_ddm_recovery/`, cluster) — simulates synthetic
   datasets from known parameters over the real Study-2 trial structure, refits
   with the same Stan models, and compares recovered vs. true parameters.

6. **Race-model robustness** (`06_rdm_emc2/`, cluster) — `study2_rdm.R` fits the
   hierarchical RDM to Study-2 CS/SC trials with EMC2; the `recovery_*.R` scripts
   run its parameter recovery.

## Scope

This package covers the **runnable analysis pipeline** end to end — experiment →
raw data → documented exclusions → behavioural analysis → fitting → PPCs →
recovery → RDM. Downstream result-formatting utilities (LaTeX table builders,
model-ranking summaries) are intentionally not included; they only reformat
outputs of the stages above.
