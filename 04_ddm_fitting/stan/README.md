# Stan models

These 72 files are the drift-diffusion (DDM) model set fit in the study. They are
**systematic variants of a single hierarchical, non-centred Wiener model**, so
they share the same structure and differ only in which "shift" parameters are
active. Rather than reading all 72, read one representative model per family
(e.g. `mv_ccss_n_a_s.stan`) and use the naming key below.

## File-name convention: `<family>_<model>.stan`

**Family** = value representation × condition:

| Family | Value model | Condition | Notes |
|---|---|---|---|
| `mv_ccss`     | Mean–Variance–Skewness (MVS) | CCSS (options equally complex) | main-text model |
| `mv_cs`       | MVS | CS (options differ in complexity) | main-text model |
| `cpt_ccss`    | Cumulative Prospect Theory (CPT) | CCSS | robustness check (2-outcome studies) |
| `cpt_cs`      | CPT | CS | robustness check (2-outcome studies) |
| `cpt_ccss_7o` | CPT | CCSS | 7-outcome variant, **Study 1 only** |
| `cpt_cs_7o`   | CPT | CS | 7-outcome variant, **Study 1 only** |

**Model** = which shift parameters are active (the suffix lists the active ones):

- **CCSS models** — subsets of `{n, r, a, s}` (baseline + all 15 non-empty subsets = 16 files/family):
  - `n` = signal-to-noise shift `Δθ`
  - `r` = risk-preference shift `Δβ` (MVS) / utility-curvature shift `Δβ_CPT` (CPT)
  - `a` = decision-threshold shift `Δα`
  - `s` = skewness-preference shift `Δη` (MVS) / probability-weighting shift `Δγ` (CPT)
  - e.g. `mv_ccss_baseline` (no shifts), `mv_ccss_n_a_s` (Δθ+Δα+Δη), `mv_ccss_n_r_a_s` (all four)
- **CS models** — subsets of `{sp, dr, skew}` (8 files/family):
  - `sp` = starting-point bias `z`
  - `dr` = drift-rate adjustment `ζ`
  - `skew` = differential skewness/probability-weighting sensitivity — **present on disk but excluded from the fitted set** (see the manuscript's "Differential Skewness Sensitivity" section). Only `baseline`, `sp`, `dr`, `sp_dr` are fit.

The exact set of models fit per study is defined in `../R/model_registry.R` (do not
infer it from the file listing — the `*_skew` CS files and the `_7o` families are
conditionally included).

## Model structure (common to all files)

- **Data**: per-trial choice (`cho` ∈ {−1,+1}), response time (`rt`), the value
  predictors (`evd`, `sdd`, `skew`), the complexity dummy (`con` ∈ {−1,+1}) for
  CCSS, and the participant index. `starting_point` is fixed for CCSS models and
  estimated (`z`) in CS models.
- **Parameters**: a group-mean vector `mu`, a between-participant SD vector
  `sigma`, a Cholesky factor `L_corr` of the correlation matrix, and standard-normal
  deviations `z` (non-centred parameterisation). Positive parameters (θ, α, T₀)
  use a softplus `log(1+exp(·))` link; the starting point uses a probit link.
- **Likelihood**: the Wiener first-passage-time distribution over (choice, RT),
  with the drift rate built from the value representation and the active shifts.

All models use the array syntax required by CmdStan ≥ 2.33.
