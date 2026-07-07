# 05_ddm_recovery — simulation-based DDM parameter recovery

Self-contained stage that **validates the DDM fitting stage** (`../03_ddm_fitting/`).
For each hierarchical DDM model (MVS and CPT families) it:

1. draws true participant parameters from plausible ranges on the raw
   (Stan) scale (`R/priors.R`),
2. maps them to per-trial DDM parameters exactly as the Stan model does
   (`R/transforms.R`) over the **real Study-2 trial structure**
   (`R/trials.R`),
3. simulates choices + response times with `rtdists` and rejection-sampling
   (`R/simulate.R`),
4. refits each synthetic dataset with the same Stan model via `cmdstanr`
   (`R/fit.R`), and
5. compares recovered vs. true parameters, on the raw scale
   (`R/summarize.R`).

This stage runs independently of the empirical fits; it only needs a
preprocessed Study RDS as a trial template. See the top-level
[`../README.md`](../README.md) for the full pipeline order and how this stage
relates to fitting (03) and posterior predictive checks (04).

## Models

| model key           | Stan file (`stan/`)          | # params |
|---------------------|------------------------------|----------|
| `cpt_ccss_n_r_a_s`  | `cpt_ccss_n_r_a_s.stan`      | 9        |
| `mv_ccss_n_r_a_s`   | `mv_ccss_n_r_a_s.stan`       | 9        |
| `cpt_cs_sp_dr`      | `cpt_cs_sp_dr.stan`          | 7        |
| `mv_cs_sp_dr`       | `mv_cs_sp_dr.stan`           | 7        |
| `cpt_cs_sp_dr_skew` | `cpt_cs_sp_dr_skew.stan`     | 8        |
| `mv_cs_sp_dr_skew`  | `mv_cs_sp_dr_skew.stan`      | 8        |

The full model registry lives in `R/priors.R` (`PRIORS`); the Stan-file map is
in `R/fit.R` (`STAN_PATHS`).

## Layout

```
05_ddm_recovery/
  README.md
  install_deps.R            # one-shot dependency installer
  run_recovery.R            # main CLI: simulate -> fit -> compare, one model
  R/
    priors.R                # generative hierarchical priors + sampler
    transforms.R            # raw -> effective params (mirrors each Stan model)
    trials.R                # load + recode a trial template from a Study RDS
    simulate.R              # rtdists DDM simulation + Stan-data builder
    fit.R                   # cmdstanr compile + sample
    summarize.R             # merge truth + posterior, stats, recovery plot
  stan/                     # 6 Stan sources (copied from the fitting stage)
  scripts/
    submit_recovery.sh       # SLURM wrapper (single run)
    submit_recovery_long.sh  # SLURM wrapper with a larger time/QoS budget
    submit_recovery_array.sh # SLURM array over (study, model)
    aggregate_partial.R      # rebuild summaries from a partial/killed run
    check_one_dataset.R      # per-dataset recovery + diagnostics inspector
  data/                     # preprocessed Study RDS (trial template; add here)
  stan_cache/               # compiled Stan binaries (created on first run)
  logs/                     # SLURM stdout/stderr (created on first run)
  results/                  # per-run output
```

`data/` is not shipped in git; place a preprocessed Study RDS
(e.g. `final_df_study2.rds`) here before running.

## Install (once)

```bash
cd 05_ddm_recovery
Rscript install_deps.R
```

Installs `optparse dplyr tidyr tibble readr ggplot2 posterior future
future.apply rtdists`, plus `cmdstanr` and CmdStan if not already available.
On a cluster, prefer site R and CmdStan modules and edit the (commented)
`module load` lines in the `scripts/*.sh` submission scripts.

## Smoke test

```bash
Rscript run_recovery.R \
  --model cpt_ccss_n_r_a_s \
  --data  data/final_df_study2.rds \
  --n_datasets 2 --n_subjects 10 \
  --chains 4 --parallel_chains 4 --parallel_datasets 1 \
  --warmup 1000 --sampling 1000 \
  --out results/cpt_ccss_n_r_a_s_smoke
```

## Production run

```bash
sbatch scripts/submit_recovery.sh cpt_ccss_n_r_a_s 50 30 study2
# or all six models in one array submission:
sbatch --array=0-5 scripts/submit_recovery_array.sh study2 all 50 30
```

The last argument selects the trial template (`study1` / `study2` / `study3` →
`data/final_df_<study>.rds`). CS-family recoveries are slow; use
`submit_recovery_long.sh` (or the array script's larger budget) for them.

## Outputs per run

```
results/<model>_<study>_dsN_LS/
  sim_dataset_{d}.rds      # simulated trials + true raw params
  fit_dataset_{d}.rds      # cmdstanr draws + summary + diagnostics
  recovery_long.csv        # tidy true-vs-estimated, every dataset
  recovery_stats.csv       # r / mae / rmse per (dataset, level, param)
  recovery.pdf / .png      # participant-level scatter grid, pooled
  recovery_caption.txt     # APA-style figure caption
  run_meta.rds             # options, elapsed time, timestamp
```

If a run is killed before aggregation, rebuild the summaries from the finished
datasets with `Rscript scripts/aggregate_partial.R <results_dir>`, and inspect
any single dataset with `Rscript scripts/check_one_dataset.R <results_dir> <d> --plot`.

## Notes

- All comparisons are on the raw/unconstrained scale — the space Stan samples
  in — so true and estimated values are directly comparable.
- `mv_cs_*` models need per-option skewness. If `Skewness_c` / `Skewness_s`
  are absent from the RDS, `R/trials.R` derives them from the closed-form
  skewness of a two-outcome lottery.
- If `simulate_dataset()` fails to find acceptable params, the prior is too
  wide for that model — adjust the ranges in `R/priors.R` or relax the
  rejection checks in `R/simulate.R`. Raise `--adapt_delta` if you see
  divergences.
```
