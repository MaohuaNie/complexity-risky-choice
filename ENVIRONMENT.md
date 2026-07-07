# Software environment

All analyses were run in **R 4.4.2**. The DDM models were fit with **CmdStan
2.38.0** via the **cmdstanr** R package; the race model with the **EMC2** R
package. Heavy fitting was run on an HPC cluster (SLURM); the scripts under each
stage's `scripts/` folder are the job submitters and are cluster-specific
(edit the `module load` line, the account/QOS, and the project path for your
system). Everything except the full fits can be reproduced on a laptop.

## Install the dependencies

Each modelling stage ships a one-shot installer:

```bash
Rscript 04_ddm_fitting/install_deps.R      # DDM fitting (also installs CmdStan)
Rscript 05_ddm_recovery/install_deps.R     # DDM recovery (same package set)
```

The race-model stage needs **EMC2** (install once):
```r
install.packages("EMC2")   # or the version pinned in 06_rdm_emc2/RECOVERY_README.md
```

## R packages by stage

| Stage | Packages |
|---|---|
| **Data preparation** (`02_`) | dplyr, tidyr, readr, tibble, janitor, stringr, purrr, fs, here, ggplot2 |
| **Behavioural** (`03_`) | dplyr, tidyr, ggplot2, purrr, cowplot, lme4, performance |
| **DDM fitting** (`04_`) | cmdstanr, posterior, optparse, dplyr, tidyr, tibble, readr, digest |
| **PPC** (`04_`, in `post/`) | cmdstanr, posterior, loo, data.table, pbmcapply, parallel, ggplot2, patchwork, ragg, rtdists, dplyr, tidyr, readr, optparse, arm |
| **DDM recovery** (`05_`) | (same as fitting) + rtdists |
| **RDM / EMC2** (`06_`) | EMC2, dplyr, ggplot2, patchwork, readr |

## CmdStan / EMC2 notes

- **CmdStan** is installed by `install_deps.R` via `cmdstanr::install_cmdstan()`
  if the `CMDSTAN` environment variable is not already set. The Stan files use
  the array syntax required by CmdStan ≥ 2.33.
- **EMC2** on the cluster was loaded from a personal R library
  (`R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4`); locally, a normal
  `install.packages("EMC2")` is sufficient.

## Compute notes

- Full DDM fits (72 models × 3 studies) are days of cluster time; the fitting
  `scripts/submit_fit_array.sh` runs them as a SLURM array. Do not attempt the
  full set on a laptop.
- PPC generation writes very large `posterior_predictives.csv` files (~1 GB for
  CCSS models); keep them on scratch/cluster storage, not in the repo.
- Recovery and single-model PPC plots are laptop-feasible for one model at a time.
