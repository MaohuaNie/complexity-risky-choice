# RDM Parameter Recovery

Parameter-recovery pipeline for the hierarchical Racing Diffusion Model (RDM)
fit in `study2_rdm.R`. It checks whether *every* group- and participant-level
parameter of the model is recoverable under the exact design and priors used in
the real Study-2 analysis. See `../README.md` for how this folder fits into the
overall replication package.

## Why

The real fit (`emc_RDM_CS_study2.RData`) estimates per-accumulator contrasts on
the CS/SC trials — notably the within-trial noise ratio `s_complex / s_simple`
and the threshold contrast `B_complex / B_simple`, which can be poorly
identified when one response type is rare. Before interpreting the real
estimates we need evidence that these (and all other) parameters are recoverable
in this design. The pipeline simulates synthetic data with **known**
group-level parameters drawn from reasonable ranges (not anchored to the real
fit) over the real Study-2 stimuli, refits with the same EMC2 spec, and compares
recovered to true for every parameter.

## Files

| File | What it does |
|---|---|
| `recovery_simulate.R` | Simulates one dataset. Draws all 13 group-level parameters from reasonable ranges, assigns each of 50 simulated participants a real subject's trial sequence (real EVD/SDD/SkewD stimuli), draws per-participant parameters, and forward-simulates with `EMC2::make_data()`. Rejects/resamples pathological participants (<2% or >98% complex choices). Saves the data plus the known group- and participant-level truths. |
| `recovery_fit.R` | Refits one simulated dataset with the same EMC2 hierarchical RDM design (`A ~ 1`) and priors as `study2_rdm.R`; extracts recovered group- and participant-level estimates (median + 95% CrI, Rhat, ESS) alongside the truths. |
| `recovery_aggregate.R` | Reads all `recovery/sum_seed*.rds`; produces a per-parameter recovery table (Pearson r, 95% CrI coverage, bias, at both individual and group level) and faceted true-vs-recovered scatter plots, plus a convergence summary. |
| `apa_theme.R` | Shared APA-style ggplot2 theme + descriptive parameter labels, sourced by `recovery_aggregate.R`. |
| `submit_recovery.sh` | SLURM array (1–30); each task = one seed → simulate → fit → summary replicate. |
| `smoke_recovery.sh` | Tiny end-to-end pipeline check (n_subj = 10, 2 chains, 200 burnin/samples). **Not** a real recovery; numbers are nonsense. |
| `submit_study2_rdm.sh` | SLURM job for the real Study-2 RDM fit (`study2_rdm.R`), not part of recovery. |
| `RECOVERY_README.md` | This file. |

## Prerequisites

- `../data/final_df_study2.rds` — real Study 2 data (source of the stimuli). The
  simulator auto-discovers it (checks `../data/`, then a couple of legacy
  relative paths).
- EMC2 R package installed in the user's R library (see the `R_LIBS_USER` line
  in the SLURM scripts — cluster/user-specific).

## How to run

```bash
# Optional smoke check (confirms the pipeline runs end-to-end)
sbatch smoke_recovery.sh

# Full sweep — 30 replicates in parallel
sbatch submit_recovery.sh
squeue -u "$USER" -h -o '%.12i %.6T %.10M %j' | grep recov

# When all 30 are done
Rscript recovery_aggregate.R
```

Each fit uses 16 CPUs / 32G / ≤24h. Outputs land in `recovery/`.

## Results

`recovery_aggregate.R` writes:

- `recovery/recovery_param_table.csv` — per-parameter individual- and
  group-level Pearson r, 95% CrI coverage, and bias.
- `recovery/recovery_individual.{pdf,png}` — participant-level true-vs-recovered
  scatters (all datasets pooled), faceted by parameter.
- `recovery/recovery_group.{pdf,png}` — group-level scatters (one point per
  dataset), faceted by parameter.
- Console — convergence summary (max Rhat, min ESS across datasets).

Reading rubric: a parameter recovers well when its scatter falls on the diagonal
and the 95% CrIs cover the truth ~95% of the time. Flat or biased scatters (e.g.
for the `s` contrast) indicate that parameter is not reliably identified in this
design and should be interpreted with caution or dropped.
