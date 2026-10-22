# Analysis code

The code that produced the Project Fulcrum outputs in the parent directory.

These are **verbatim exports of the cells that actually ran**, in execution
order, rather than a cleaned-up pipeline written afterwards. That is the point:
a tidied retelling is code that never produced the published numbers. The cells
originally ran in one persistent Python kernel, so each script is meant to be
run top-to-bottom in a single process.

| File | What it does |
|---|---|
| `00_fetch_inputs.py` | Downloads `antibody_agonism.csv` from this repository and checks it is the expected panel (500 rows, no missing values, `log10_KD == log10_koff - log10_kon`). Exits non-zero if not. |
| `01_analysis.py` | The analysis: held-out split, exploratory binning, model selection, the quasi-likelihood F-tests, mean-variance estimation, 5-fold CV, the final GLM, bootstrap CIs for the optima, and the four-panel figure. |
| `02_write_outputs.py` | Builds the report and trap inventory from the fitted model. Continues the kernel state left by `01_analysis.py`. |

## Running it

```bash
python 00_fetch_inputs.py
python 01_analysis.py
```

Requires `pandas`, `numpy`, `scipy`, `statsmodels`, `matplotlib`. `statsmodels`
is the only non-obvious dependency; **scikit-learn is not needed** - the
stratified 25% split is done with a seeded `numpy` generator, `rng(42)`.

## What was changed, and what was not

The analysis cells are unedited. Two additions were necessary to make them run
outside the platform they were written in, both confined to a clearly marked
shim at the top of `01_analysis.py`:

- `host.artifact_path()`, which resolved the input CSV on the platform, is
  replaced by a stub returning the local path to `antibody_agonism.csv`.
- `apply_figure_style()` and `META_GREY`, a house figure style, are replaced by
  a minimal matplotlib equivalent. **The regenerated figure differs
  cosmetically from the published PNG** - fonts, tick direction, panel
  proportions. Every number it reports is unaffected.

Cells that raised an exception are omitted, with one exception that is marked
inline: cell 10 failed partway through, on a `statsmodels` call that does not
exist for GLM results, but the lines before the failure had already defined
model objects that later cells use. Those lines are kept and the failure is
documented where it occurred. The quasi-likelihood F-test it was reaching for
is implemented by hand as `quasiF()` in the next cell.

Also omitted: cells that only listed platform artifacts or pushed files to
GitHub. They are plumbing, not analysis.

## Reproducibility note

Every random element is seeded, so the whole analysis reproduces exactly:

| Where | Seed | What it draws |
|---|---|---|
| cell 5 | `default_rng(42)` | the stratified 25% held-out split |
| cell 15 | `default_rng(7)` | 300 parametric-bootstrap replicates testing whether the site x KD interaction is an artifact of the 0 floor |
| cell 16 | `default_rng(seed)` | 5-fold CV fold assignment |
| cell 17 | `default_rng(3)` | 600 bootstrap refits for the CIs on the fitted optima |

The scripts were re-run from a clean directory against the public CSV before
being committed.
