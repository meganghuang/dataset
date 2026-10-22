# Project Fulcrum — the final model, written out

The recommended model from [`fulcrum_report.md`](fulcrum_report.md), stated explicitly: every term, every
fitted coefficient, and a worked arithmetic example. Fitted on the 375 training antibodies
(seed-42 per-site stratified holdout); the 125 held-out antibodies were scored once, at the end.

## Formula notation

```
agonism_pct ~ site + log10_KD + log10_KD^2 + log10_koff + log10_koff^2 + plate
```

quasi-Poisson family, **log link**, `site` and `plate` as factors. This notation names the terms only —
the coefficients are implicit, and the log lives in the family/link argument, not in the formula.

## The same model, explicit

Write `K` = `log10_KD` and `F` = `log10_koff`. For antibody *i* at site *s* on plate *p*:

```
agonism_pct_i ~ quasi-Poisson,  mean mu_i,  Var(y_i) = phi * mu_i,  phi = 11.568

log(mu_i) = eta_i

eta_i = -26.0345
        + 1.7962*[site=2] + 2.4726*[site=3] + 2.1378*[site=4] + 0.7822*[site=5]
        -4.6329*K -0.2319*K^2
        -3.0743*F -0.5190*F^2
        + 0.3935*[plate=2] + 0.4995*[plate=3] + 0.2776*[plate=4]
        + 0.4582*[plate=5] -0.0141*[plate=6]

mu_i = exp(eta_i)
```

`[site=3]` is an indicator: 1 if the antibody binds site 3, 0 otherwise. **Site 1 and plate 1 are the
reference levels** — they have no term because their effect is absorbed into the intercept, so
"site 3 = +2.4726" means *relative to site 1*. Fourteen parameters in total:
1 intercept, 4 site contrasts, 4 kinetic slopes, 5 plate contrasts.

Coefficients are on the log scale, so `exp(coef)` is a **multiplier** on activity.

## Coefficients

| Term | Coefficient | SE | 95% CI | exp(coef) | p |
|---|---|---|---|---|---|
| `Intercept (site 1, plate 1 reference)` | -26.0345 | 3.1399 | -32.1886, -19.8803 | 4.936e-12 | 1.1e-16 |
| `site 2 vs site 1` | +1.7962 | 0.2816 | +1.2441, +2.3482 | 6.026 | 1.8e-10 |
| `site 3 vs site 1` | +2.4726 | 0.2762 | +1.9313, +3.0139 | 11.85 | 3.4e-19 |
| `site 4 vs site 1` | +2.1378 | 0.2775 | +1.5939, +2.6817 | 8.481 | 1.3e-14 |
| `site 5 vs site 1` | +0.7822 | 0.3162 | +0.1624, +1.4019 | 2.186 | 0.013 |
| `log10_KD` | -4.6329 | 0.7420 | -6.0873, -3.1785 | 0.009727 | 4.3e-10 |
| `(log10_KD)^2` | -0.2319 | 0.0409 | -0.3122, -0.1517 | 0.793 | 1.5e-08 |
| `log10_koff` | -3.0743 | 0.4853 | -4.0256, -2.1231 | 0.04622 | 2.4e-10 |
| `(log10_koff)^2` | -0.5190 | 0.0733 | -0.6627, -0.3753 | 0.5951 | 1.4e-12 |
| `plate 2 vs plate 1` | +0.3935 | 0.1928 | +0.0157, +0.7714 | 1.482 | 0.041 |
| `plate 3 vs plate 1` | +0.4995 | 0.1974 | +0.1126, +0.8863 | 1.648 | 0.011 |
| `plate 4 vs plate 1` | +0.2776 | 0.1981 | -0.1107, +0.6659 | 1.32 | 0.16 |
| `plate 5 vs plate 1` | +0.4582 | 0.1926 | +0.0807, +0.8357 | 1.581 | 0.017 |
| `plate 6 vs plate 1` | -0.0141 | 0.2178 | -0.4410, +0.4128 | 0.986 | 0.95 |

Dispersion 11.568; pseudo-R2 0.887 on the training set.
`exp(intercept)` = 4.94e-12 is not interpretable on its own: K is around -10
and K^2 around 100, so the kinetic terms pull strongly the other way.
The same table with machine-readable precision is in
[`fulcrum_coefficients.csv`](fulcrum_coefficients.csv).

## Worked example — antibody S3_63

Site 3, plate 2, `log10_KD` = -9.127, `log10_koff` = -2.990:

| Term | Value |
|---|---|
| Intercept | -26.0345 |
| site 3 | +2.4726 |
| plate 2 | +0.3935 |
| -4.6329 x (-9.127) | +42.2847 |
| -0.2319 x (-9.127)^2 | -19.3190 |
| -3.0743 x (-2.990) | +9.1914 |
| -0.5190 x (-2.990)^2 | -4.6391 |
| **eta** | **4.3497** |

mu = exp(4.3497) = **77.5 %** of reference; observed 92.09 %.
At the fitted optimum (site 3, plate 1, K = -9.99, F = -2.96) the same equation gives
62.1 %.

## Why a log *link* and not a log-transformed response

The link applies to the fitted mean, not to the data:

- **log link (used here):** `log(E[y]) = eta`, i.e. `E[y] = exp(eta)`. `agonism_pct` stays on its raw
  percent scale and is never transformed. 113 of 500 antibodies read exactly 0; an observed zero is an
  ordinary draw from a distribution with a small positive mean, so every row is kept.
- **log transform + OLS:** `log(y) = eta + error`. The 113 zeros have to be dropped or fudged, which
  conditions on the outcome, and the fit estimates `E[log y]`, not `log E[y]`.

The variance model comes with the family: quasi-Poisson assumes Var(y) = phi * E[y], matching the measured
sd ~ mean^0.57 (95% CI 0.24-0.89) rather than the constant variance a log-transform fit would assume.

## Held-out performance

On the 125 untouched antibodies: R2 0.935, RMSE 4.39 % with the plate effect averaged out
(R2 0.930, RMSE 4.54 % predicting at each antibody's observed plate),
against R2 0.11 / RMSE 16.2 for a plain linear fit — the 0.94 and 4.1-4.4 quoted in the write-up.

## Reproducing the fit

```python
import pandas as pd, numpy as np, statsmodels.api as sm, statsmodels.formula.api as smf

df = pd.read_csv("antibody_agonism.csv")
rng = np.random.default_rng(42)
test_idx = []
for s, g in df.groupby("site"):
    test_idx += list(rng.choice(g.index.values, size=int(round(0.25*len(g))), replace=False))
tr = df.drop(index=test_idx).copy()
tr["siteC"]  = pd.Categorical(tr.site,  categories=[1, 2, 3, 4, 5])
tr["plateC"] = pd.Categorical(tr.plate, categories=[1, 2, 3, 4, 5, 6])

final = smf.glm(
    "agonism_pct ~ siteC + log10_KD + I(log10_KD**2) + log10_koff + I(log10_koff**2) + plateC",
    tr, family=sm.families.Poisson(sm.families.links.Log()),
).fit(scale="X2")      # scale="X2" is what makes it quasi-Poisson
```

Full pipeline in [`code/`](code/).

## What the formula does not say

Two terms were deliberately left out, each for a stated reason (see `fulcrum_report.md`): a site x KD
interaction, significant at p = 2e-5 but concentrated in the 44 %-censored site 1 and costing held-out
R2; and anything extrapolating past the panel's kon range, since the fitted kon optimum (10^7.03) sits
just beyond the largest measured value (10^6.997). `receptor_density` never reached significance
(p = 0.21-0.31) and is an assay-cell property, not a design lever.
