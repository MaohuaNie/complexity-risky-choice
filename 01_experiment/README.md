# Online experiment (data-collection stage)

The browser-based experiments that collected the data, built with **jsPsych
7.3.3** and deployed on a **JATOS** server. This stage produces the raw trial
data that (after pre-processing) become the derived `../data/final_df_study*.rds`
used by every downstream analysis stage.

## Files

| File | Study | Manipulation |
|---|---|---|
| `study1_online_experiment_jatos.html` | Study 1 — lottery complexity | number of outcomes (7-outcome vs. 2-outcome lottery) |

Studies 2 and 3 use essentially the **same experiment code** — only the stimuli
(the lottery displays) differ, implementing the arithmetic-expression and
compound-probability manipulations. They are therefore not reproduced separately
here; this Study 1 file documents the construction for all three studies.

## How the experiment is built

The single HTML file defines the whole jsPsych timeline (see the header comment
at the top of the file for the ordered walkthrough): image preload → browser
check + fullscreen + consent → instructions and worked examples → comprehension
quizzes that loop until answered correctly → the main lottery-choice trials
(CS and CCSS conditions) with per-trial feedback and reaction-time recording →
questionnaires and cognitive-ability items → debrief. Trial data are submitted
to JATOS as CSV via `jatos.submitResultData()`.

## Running it

This file is not self-contained — it depends on assets that live alongside it in
the JATOS study-assets folder and are **not bundled** in this package:

- `jatos.js` (provided by the JATOS server),
- the local jsPsych plugins under `dist/`,
- the image stimuli under `ques_prac/` and the main-trial image folders.

To run or inspect it, import this HTML together with those assets as a **JATOS
study** (or serve the original `study_assets_root/` folder over http). It is
included here to **document how the experiment was constructed**, not as a
one-click reproducible build.

See ../README.md for the full pipeline.
