"""
Project Fulcrum - building the published report and trap inventory

VERBATIM EXPORT. These are the cells that actually produced the published
outputs, in execution order, not a tidied retelling. They originally ran in one
persistent Python kernel, so run this file top-to-bottom in a single process.

Cells that raised are omitted; one partial include is marked inline.
Environment: python 3 with pandas, numpy, scipy, statsmodels, matplotlib.
statsmodels is the only non-obvious requirement (scikit-learn is NOT needed -
the stratified split is done with numpy, seeded rng(42)).

Run 00_fetch_inputs.py first - it puts antibody_agonism.csv beside this file.
"""

# Continues the kernel state left by 01_analysis.py. To run standalone:
#     exec(open("01_analysis.py").read())
# first. The {{artifact:...}} markers these cells emit are platform-internal;
# the copies published in this repo have them rewritten to repo-relative paths.

# ---- cell 23 ----
fig_v="8219a1e8-1923-411e-9154-31c3a9778fb4"; csv_v="285a7a92-a56f-467c-aaf0-720c06dac2af"
ci=lambda n:(final.params[n]-1.96*final.bse[n], final.params[n]+1.96*final.bse[n])
rows_coef=[]
for n,lab in [('log10_KD','log10 KD'),('I(log10_KD ** 2)','(log10 KD)^2'),('log10_koff','log10 koff'),('I(log10_koff ** 2)','(log10 koff)^2')]:
    lo,hi=ci(n); rows_coef.append(f"| `{lab}` | {final.params[n]:+.3f} | {final.bse[n]:.3f} | {lo:+.3f}, {hi:+.3f} | {final.pvalues[n]:.1e} |")
site_rows=[]
for s in [2,3,4,5]:
    n=f'siteC[T.{s}]'; lo,hi=ci(n); c=final.params[n]
    site_rows.append(f"| site {s} vs site 1 | {c:+.3f} | x{np.exp(c):.2f} | x{np.exp(lo):.2f}-x{np.exp(hi):.2f} |")
plate_rows=[]
for p_ in range(2,7):
    n=f'plateC[T.{p_}]'; lo,hi=ci(n); c=final.params[n]
    plate_rows.append(f"| plate {p_} vs plate 1 | {c:+.3f} | x{np.exp(c):.2f} | x{np.exp(lo):.2f}-x{np.exp(hi):.2f} |")
nl=chr(10)
print("coef tables built:", len(rows_coef), len(site_rows), len(plate_rows))

# ---- cell 24 ----
report = f"""# Project Fulcrum - what switches the receptor on

**Analysis of 500 candidate agonist antibodies.** 375 used for exploration and fitting, **125 held out
and touched once**, at the end.

![Project Fulcrum results]({{{{!artifact:{fig_v}}}}})

*(a) Agonism by epitope site, individual antibodies in grey, site mean +/-95% CI in blue. (b) Agonism against
off-rate, all 500 antibodies; curve is the fitted model at site 3 with affinity held at its optimum.
(c) Fitted response surface over the kinetic space the team actually engineers, site 3, plate effect averaged
out; grey points are the panel's measured kinetics; star is the predicted optimum. (d) Held-out predictions
against observation for the 125 antibodies never used in fitting; dashed line is y = x.*

## Headline

Agonism is governed by **three things with interior optima, none of them linear**: epitope distance,
affinity, and off-rate. A plain linear fit misses all three and concludes that distance does not matter.

| Model | Held-out R2 | RMSE (% of reference) |
|---|---|---|
| Plain linear fit: `dist + kon + koff` | **0.11** | 16.2 |
| Recommended log-link model | **0.94** | 4.1-4.4 |

## Why the plain linear fit disappoints

Fitting `agonism ~ epitope_dist_nm + log10_kon + log10_koff` by OLS gives distance a coefficient of
**+0.39 %/nm (p = 0.069)** - effectively "distance has no effect". That is wrong, and it is wrong for a
specific reason: **activity peaks at site 3 and falls away on both sides** (panel a), so a straight line
through a symmetric hump has slope near zero. The same happens in off-rate. Three separate problems have
to be fixed before the data make sense:

1. **Non-monotone site effect.** Site means run 2.5 -> 14.3 -> 30.2 -> 20.5 -> 5.3 %. Site must enter as a
   factor, not as a linear distance. (`site` and `epitope_dist_nm` are redundant - same information.)
2. **Wrong scale.** Within each site the kinetic response has the same *shape* but a different *amplitude*:
   site 3 swings 1.6 -> 56 %, site 1 only 1.4 -> 3.7 %. Site acts as a **multiplicative gain**. On the raw
   scale this forces a thicket of site x kinetics interactions (R2 0.60 -> 0.79); a log link absorbs it
   almost entirely. Measured in censoring-free strata, the noise scales as sd ~ mean^0.57
   (95% CI 0.24-0.89) - variance proportional to the mean, Poisson-like, so a **quasi-Poisson GLM with a
   log link** fits both the multiplicative mean and the variance.
3. **113 of 500 values sit exactly on the 0 floor** (22.6%). A log *transform* would require dropping them,
   and dropping zeros means conditioning on the outcome - it flattens slopes hardest in exactly the
   low-activity sites. A log *link* models E[y] and keeps every row. This matters: the apparent
   site x kinetics interaction is strong on log-transformed positives (p ~ 9e-11) and much weaker once the
   zeros are retained.

## The collinearity question

`log10_KD` is exactly `log10_koff - log10_kon` (verified to 3.6e-15), so all three cannot enter together.
But the resolution is more interesting than "drop one". Added one at a time to the additive model, three
terms came up significant: `kon^2` (p = 6e-4), `koff^2`, and the `kon:koff` cross-term (p = 1.2e-7) - which
is exactly what a **quadratic in KD** expands to, since KD = koff - kon. Re-parameterising to
`KD + KD^2 + koff + koff^2` reproduces the full quadratic surface in (kon, koff) with one fewer parameter
and no loss of fit (F(1,365) = 0.021, p = 0.88). Decomposed, the surface is **-0.237*KD^2 - 0.536*koff^2**:
concave in both, i.e. an optimal **affinity** and, separately, an optimal **residence time**. No cubic
terms and no KD:koff term are needed (all p > 0.29).

## Recommended model

Quasi-Poisson GLM, log link, fitted on the 375 training antibodies:

```
agonism_pct ~ site + log10_KD + log10_KD^2 + log10_koff + log10_koff^2 + plate
```

Dispersion 11.6, pseudo-R2 0.887. Coefficients are on the log scale, so each is a **multiplier** on activity.

| Kinetic term | Coefficient | SE | 95% CI | p |
|---|---|---|---|---|
{nl.join(rows_coef)}

| Site (gain) | Coefficient | Multiplier | 95% CI |
|---|---|---|---|
{nl.join(site_rows)}

| Plate (nuisance) | Coefficient | Multiplier | 95% CI |
|---|---|---|---|
{nl.join(plate_rows)}

**Fitted optimum** (joint maximum of the kinetic surface; Hessian negative definite, so a true maximum;
95% CIs from 600 bootstrap refits):

| Quantity | Optimum | 95% CI | In panel range? |
|---|---|---|---|
| log10 KD | **-9.99** (KD ~ 100 pM) | -10.16 ... -9.84 | yes (-11.8 ... -5.1) |
| log10 koff | **-2.96** (~15 min residence) | -3.03 ... -2.89 | yes (-5.0 ... -1.0) |
| log10 kon | **7.03** | 6.84 ... 7.23 | **no - at the boundary** (max 7.00) |

## Answers to the six questions

1. **Across sites:** strongly **non-monotone**, peaking at site 3 (9 nm). Relative to site 1 the gain is
   x6.0 (site 2), **x11.9 (site 3)**, x8.5 (site 4), x2.2 (site 5). Site 1 (3 nm, membrane-proximal) is
   nearly dead - 44% of its antibodies read exactly zero.
2. **With kon and koff:** kon is monotone increasing and roughly linear in log space; **koff is not** - it
   has a clear interior optimum at 10^-2.96 s^-1. A straight line in either raw or log koff is not
   reasonable. The honest description is a quadratic surface over (log kon, log koff).
3. **Interactions:** mostly an artifact of mis-specification. site x koff looks significant (p = 1.6e-4)
   until the kon curvature is included, then drops to p = 0.099. A **site x KD** interaction does survive
   (p = 2e-5) and is not explained by the 0 floor - I tested this with 300 parametric-bootstrap replicates
   under a no-interaction model with the floor applied, which reproduced the zero fraction (0.221 vs 0.243
   observed) but produced spurious F of at most 3.55 against 4.63 observed. However, the fitted **optimal
   affinity is the same at every site where it is measurable** (vertex -9.98, -10.00, -9.71, -10.74 for
   sites 2-5); the interaction is in curvature *amplitude*, concentrated in site 1, whose shape is
   unidentifiable because it sits on the floor. **It does not improve held-out prediction and I left it
   out.**
4. **Plate matters, modestly.** Plates span a x1.67 range worst-to-best (plate 3 x1.65 vs plate 6 x0.99
   relative to plate 1; p = 0.0097). Because antibodies were spread evenly across plates - the site x plate
   design is balanced - **plate does not confound the site or kinetic conclusions**; it only inflates
   residual noise. Keeping it in sharpens the coefficients and improves ranking (Spearman rho 0.884 ->
   0.918) but barely changes RMSE.
5. **Best predictive model:** the one above. Held-out R2 = **0.94**, RMSE 4.1-4.4 % of reference on the 125
   untouched antibodies, matching the 5-fold CV estimate within the training set (R2 0.94) - no
   overfitting. Adding the interaction terms *lowers* held-out performance.
6. **Recommendations:** below.

## What I would tell the team

**Take forward:** the top of the ranked list is **site-3 antibodies with KD near 100 pM and koff near
10^-3 s^-1**. The top 15 by predicted activity are all site 3, with observed agonism of 64-92 % (median
74 %). Of the top 50, 32 are site 3, 16 site 4, 2 site 2 - and **none from sites 1 or 5**. Full ranking
attached.

**Design rules for the next round:**
- **Target the 9 nm epitope.** This is the largest single lever - x11.9 over the membrane-proximal site.
  Site 4 (12 nm) is a reasonable hedge at x8.5. Do not spend further effort on sites 1 or 5.
- **Tune to intermediate affinity, not maximum affinity.** The tightest binders in the panel are *worse*
  than 100 pM binders. Affinity maturation past ~100 pM is predicted to lose activity.
- **Tune off-rate to ~10^-3 s^-1** (~15 min residence). Both faster and slower are worse, independently of
  affinity - consistent with needing the receptor held together long enough to signal but not permanently.
- **Raise kon as far as you can.** Within this panel higher is always better, with no interior optimum
  reached.
- **Calibrate the plates.** A x1.67 uncontrolled swing is a large fraction of the differences being used to
  rank candidates. Put a reference agonist on every plate next round.

## Where this model cannot speak

- **The kon optimum is extrapolation.** The fitted optimum, log10 kon = 7.03, lies just past the panel's
  maximum of 7.00. "Push kon higher" is the model reading off its own boundary; whether kon has its own
  optimum beyond 10^7 M^-1 s^-1 **cannot be answered with these data**.
- **One dose, one replicate.** Every antibody was tested once at a single fixed dose, so nothing here
  separates potency from efficacy - a candidate that looks weak may simply have a right-shifted dose
  response. There is no within-antibody replication, so the dispersion (11.6) mixes assay noise and
  biological variation and cannot be decomposed.
- **`receptor_density` does nothing here** (p = 0.21-0.31 in every model). It is a property of the assay
  cells, not a design variable, and its range (0.49-1.97) is narrow - this is not evidence that receptor
  expression is irrelevant in vivo.
- **22.6% of the panel is at the 0 floor**, where the data carry no information about *how* inactive an
  antibody is. Site 1's shape in particular is unidentifiable. A censored-regression (Tobit-style)
  treatment would be the next refinement.
- **Correlation, not mechanism.** The optima are consistent with a receptor-crosslinking story, but this is
  an observational screen: the kinetics were measured, not assigned. The only way to test the design rules
  is to build antibodies at the predicted optimum and see whether they land near 80 %.
- **Site 3 is best within 3-15 nm.** The true geometric optimum may lie between the sampled distances;
  5 discrete sites cannot resolve it. A site at 10-11 nm is untested and might be better.

**Attached:** [fulcrum_ranked_candidates.csv]({{{{!artifact:{csv_v}}}}}) - all 500 antibodies ranked by
predicted agonism with the plate effect averaged out, flagged for held-out membership.
"""
open('fulcrum_report.md','w',encoding='utf-8').write(report)
print("report written:", len(report), "chars")

# ---- cell 25 ----
tricky = """================================================================================
PROJECT FULCRUM - THE TRICKY PARTS OF THIS DATASET
A catalogue of the traps in antibody_agonism.csv, what each one does to a naive
analysis, and how to detect it. Every number below was computed from the data.
================================================================================

SUMMARY: twelve traps. Three are structural (redundancy and exact collinearity
among the predictors), three are shape problems (non-monotone effects where a
linear term reads as null), three are error-structure problems (multiplicative
mean, variance proportional to the mean, 22.6% of the outcome on a floor), and
three are inference traps (an interaction manufactured by mis-specification, a
significant term that hurts prediction, an optimum that lies off the edge of the
measured range).

The combined cost of missing them: held-out R2 of 0.11 instead of 0.94.

--------------------------------------------------------------------------------
GROUP A - STRUCTURAL: predictors that are not independent
--------------------------------------------------------------------------------

TRAP 1. site and epitope_dist_nm are the same variable.
  Distance is fixed per site: site 1 = 3 nm, 2 = 6, 3 = 9, 4 = 12, 5 = 15, with
  exactly one distance value per site. Putting both in a model is a rank
  deficiency; the choice between them is a modelling decision, not a data
  question. Factor(site) spends 4 df and can represent any shape; linear
  epitope_dist_nm spends 1 df and can only represent a straight line.
  DETECT: df.groupby('site').epitope_dist_nm.nunique() -> all 1.

TRAP 2. log10_KD is exactly log10_koff minus log10_kon.
  Verified to floating-point: max |KD - (koff - kon)| = 3.6e-15. All three in one
  model is exact collinearity. The README warns about this, but the warning
  understates the interesting part - see Trap 4.
  DETECT: compute the difference and look at its maximum, not at a correlation
  matrix. The pairwise correlations are unremarkable (kon vs koff = -0.049,
  KD vs koff = 0.812, KD vs kon = -0.623) and would not by themselves reveal an
  exact identity.

TRAP 3. kon and koff are deliberately NOT tied to each other.
  r = -0.049. This is unlike real SPR panels, where affinity maturation usually
  moves both. It means the design CAN separate the on-rate effect from the
  off-rate effect - an opportunity, not a problem - but it also means KD carries
  genuinely different information from either rate, so "just use KD" discards
  real signal (pseudo-R2 0.733 for site+KD+KD^2 vs 0.865 once koff is added back).

--------------------------------------------------------------------------------
GROUP B - SHAPE: effects that are real but non-monotone, so a linear term reads
as null
--------------------------------------------------------------------------------

TRAP 4. The site effect is a hump, and a linear distance term misses it entirely.
  Mean agonism by site: 2.5, 14.3, 30.2, 20.5, 5.3 %. It rises to site 3 (9 nm)
  and falls away on both sides - roughly symmetric, so a straight line through it
  has almost no slope. Fitting agonism ~ epitope_dist_nm + kon + koff gives
  distance a coefficient of +0.39 %/nm with p = 0.069: the single largest driver
  in the dataset (x11.9 between best and worst site) is reported as not
  significant.
  THIS IS THE DATASET'S CENTRAL TRAP. A p-value near 0.07 on a null-looking
  coefficient invites exactly the wrong conclusion: "epitope position doesn't
  matter, focus on affinity".
  DETECT: plot the outcome against the ordinal predictor before fitting anything,
  or fit it as a factor and compare.

TRAP 5. koff has an interior optimum.
  Binned mean agonism against log10_koff: 8.1, 23.6, 26.8, 22.0, 6.5, 1.5 across
  sextiles from -4.65 to -1.30. Peak near 10^-3, collapsing at both ends. A
  linear koff term - in raw or in log scale - fits a monotone decline through a
  hump. The fitted optimum is log10_koff = -2.96 (95% CI -3.03 to -2.89), i.e.
  about 15 minutes of residence time.
  DETECT: quantile-binned means, or add a quadratic and test it
  (F(1,366) = 59.4, p = 1.2e-13 for koff^2 on top of KD + koff).

TRAP 6. Affinity has an interior optimum too, so "tighter is better" is false.
  The tightest binders in the panel are less active than 100 pM binders: binned
  means against log10_KD are 18.2, 28.4, 23.1, 12.2, 5.5, 1.2 - the peak is in
  the second sextile, not the first. Fitted optimum log10_KD = -9.99
  (95% CI -10.16 to -9.84). Any recommendation of the form "affinity-mature the
  lead further" is contradicted by the data.

--------------------------------------------------------------------------------
GROUP C - ERROR STRUCTURE: the outcome is not conditionally Gaussian with
constant variance
--------------------------------------------------------------------------------

TRAP 7. Site is a multiplicative gain, not an additive offset.
  Within each site the kinetic response has the same SHAPE but a different
  AMPLITUDE. Across koff quartiles, site 3 swings 1.6 -> 56.2 % while site 1
  swings only 1.4 -> 3.7 %. On the raw scale this must be absorbed by a thicket
  of site x kinetics interactions (R2 0.600 additive -> 0.785 with interactions,
  12 extra parameters). Under a log link the same structure is additive
  (pseudo-R2 0.815 -> 0.835 for the same 12 parameters).
  WHY IT MATTERS: on the wrong scale the interactions are real and necessary,
  and an analyst who adds them gets a defensible but uninterpretable model with
  19 parameters instead of 7. The interaction is an artifact of the scale, not
  biology.
  DETECT: compare within-group ratios rather than within-group differences;
  check whether the coefficient of variation is more stable than the SD.

TRAP 8. Variance grows with the mean.
  Estimated in the 17 strata (site x koff quartile x kon tercile) that contain NO
  censored values, so the estimate is not distorted by the floor:
  log(SD) ~ log(mean) has slope 0.567 (95% CI 0.241 to 0.892). That is
  SD proportional to sqrt(mean), i.e. VARIANCE proportional to the MEAN -
  Poisson-like, not constant (slope 0) and not constant-CV (slope 1). Residual
  dispersion in the final model is 11.6, so it is heavily over-dispersed relative
  to a true Poisson; quasi-Poisson (or Tweedie) is the right quasi-likelihood.
  OLS with constant variance will have badly heteroscedastic residuals and
  standard errors that are too small where activity is high.
  DETECT: estimate the mean-variance power, and do it only in strata with no
  censoring - censored strata bias the SD downward.

TRAP 9. 113 of 500 outcomes (22.6%) sit EXACTLY on the 0 floor.
  Not small values - exact zeros, from flooring. Distribution is wildly uneven:
  44% of site 1, 12-14% of sites 2-4, 29% of site 5; and 60% of the fastest-koff
  sextile. Two consequences:
    (a) A log TRANSFORM of the outcome is impossible without dropping or offsetting
        them. Dropping them is conditioning on the outcome - a selection effect -
        and it flattens relationships hardest in exactly the low-activity sites,
        which manufactures interactions. On log-transformed positives (284 of 375
        training rows) the site x kinetics interaction tests at p = 9.3e-11; with
        the zeros retained under a log LINK the same comparison is far weaker.
    (b) The floor destroys information about HOW inactive an antibody is. Site 1's
        dose-response shape is simply unidentifiable.
  THE FIX THAT WORKS: a log LINK models E[y] rather than log(y), so every row
  including the zeros is retained. This is the single most important modelling
  decision in the dataset.

--------------------------------------------------------------------------------
GROUP D - INFERENCE: ways the analysis misleads even after the model fits
--------------------------------------------------------------------------------

TRAP 10. A spurious interaction manufactured by omitted curvature.
  With the kinetics specified as kon + koff + koff^2, a site x koff interaction
  tests at F(8,359) = 3.98, p = 1.6e-4 - apparently strong evidence that the
  off-rate optimum differs by epitope. Once the kinetic surface is specified
  correctly (adding the kon curvature, i.e. KD + KD^2 + koff + koff^2), the same
  interaction falls to p = 0.069, and to p = 0.099 with plate in the model.
  The "interaction" was omitted main-effect curvature leaking into the
  interaction terms. Always get the functional form of the main effects right
  BEFORE testing interactions.

  Related, and the prettiest structure in the dataset: added one at a time, three
  terms all look significant - kon^2 (p = 6.2e-4), koff^2, and kon:koff
  (p = 1.2e-7). That triplet is exactly what a QUADRATIC IN KD expands to, since
  KD = koff - kon gives KD^2 = kon^2 - 2*kon*koff + koff^2. Re-parameterising to
  KD + KD^2 + koff + koff^2 reproduces the full quadratic surface with one fewer
  parameter and no loss of fit (F(1,365) = 0.021, p = 0.88). The fitted surface
  decomposes as -0.237*KD^2 - 0.536*koff^2: concave in both, i.e. an optimal
  affinity AND, separately, an optimal residence time.

TRAP 11. A statistically significant interaction that makes prediction WORSE.
  A site x KD interaction survives correct specification at p = 2.2e-5, and it is
  NOT an artifact of the floor - I tested that with 300 parametric-bootstrap
  replicates simulated under a no-interaction model with the floor applied, which
  reproduced the observed zero fraction (0.221 simulated vs 0.243 observed) but
  produced a spurious F of at most 3.55 (median 0.92) against 4.63 observed.
  So the term is real. It is still the wrong thing to put in the model:
    - the fitted optimal affinity is essentially IDENTICAL at every site where it
      is measurable (vertex -9.98, -10.00, -9.71, -10.74 for sites 2-5); only the
      curvature amplitude differs;
    - it is concentrated in site 1, which is 44% censored and whose fitted
      quadratic is convex (no interior optimum at all) - an unidentifiable shape;
    - excluding site 1 drops it to p = 0.032, and sites 2-4 only to p = 0.029;
    - and in 5-fold CV it REDUCES held-out R2, from 0.940 (additive) to 0.932.
  Significance and predictive value are different questions, and with n = 500 a
  real-but-tiny effect can be both significant and harmful to include.

TRAP 12. The fitted optimum in kon lies OFF THE EDGE of the measured range.
  Joint maximum of the kinetic surface: log10_kon = 7.027 (95% CI 6.84 to 7.23).
  The panel's maximum measured log10_kon is 6.997. The model is reading an
  optimum off its own boundary, so "push kon higher" is extrapolation, not a
  finding. Whether kon has a true interior optimum beyond 10^7 M^-1 s^-1 cannot
  be answered with these data. By contrast the koff and KD optima (-2.96, -9.99)
  are comfortably interior and trustworthy.
  DETECT: always compare a fitted vertex against the observed range of the
  variable before reporting it as a recommendation.

--------------------------------------------------------------------------------
TWO THINGS THAT LOOK LIKE TRAPS AND ARE NOT
--------------------------------------------------------------------------------

NOT A TRAP. receptor_density is a red herring.
  Marginally flat (binned means 11.8, 14.8, 14.6, 15.1, 15.5, 16.6) and never
  significant in any model: p = 0.236, 0.307, 0.206 depending on specification.
  It is a property of the assay cells rather than a design variable, and its range
  (0.49 to 1.97) is narrow. Including it costs a df and buys nothing. Note this is
  NOT evidence that receptor expression is irrelevant biologically - only that it
  does not vary enough here to matter.

NOT A CONFOUND. Plate is a real batch effect but a balanced one.
  Plate matters: a x1.67 spread worst-to-best (plate 3 at x1.65 vs plate 6 at
  x0.99 relative to plate 1), p = 0.0097. Plates were run on different days and
  were not calibrated. BUT the site x plate design is balanced (9 to 23 per cell,
  all 30 cells populated), and kinetics do not vary across plates (plate mean
  log10_kon 5.37 to 5.60, log10_koff -3.15 to -2.85). So plate does NOT bias the
  site or kinetic coefficients - it only inflates residual noise. This is the
  opposite of the usual batch-confounding situation and it is worth stating
  explicitly rather than assuming.
  THE REAL SUBTLETY: plate is unknowable for a future antibody on a future plate.
  For inference, keep it in. For prospective prediction, average the linear
  predictor over the six plate levels rather than silently predicting at the
  reference plate - predicting at plate 1 would systematically under-predict,
  since plate 1 is the second-worst of the six.

--------------------------------------------------------------------------------
WHAT A CLEAN ANALYSIS LOOKS LIKE
--------------------------------------------------------------------------------

  agonism_pct ~ site + log10_KD + log10_KD^2 + log10_koff + log10_koff^2 + plate

  quasi-Poisson family, log link, site and plate as factors.

  Held-out R2 = 0.94 on 125 antibodies never used in fitting (RMSE 4.1-4.4 % of
  reference), matching the within-training 5-fold CV estimate of 0.94.
  The naive linear fit gets held-out R2 = 0.11.

  The ordered decisions that get you there:
    1. hold out a test set before looking at anything;
    2. site as a factor, not linear distance            (Trap 4);
    3. log link, not log transform, because of the floor (Trap 9);
    4. quasi-Poisson because variance ~ mean             (Trap 8);
    5. quadratic in KD and in koff, not linear rates     (Traps 5, 6, 10);
    6. plate in for inference, averaged for prediction   (not a confound);
    7. drop receptor_density                             (red herring);
    8. leave the significant site x KD interaction OUT   (Trap 11);
    9. flag the kon optimum as extrapolation             (Trap 12).
================================================================================
"""
open('fulcrum_tricky_parts.txt','w',encoding='utf-8').write(tricky)
print("written:", len(tricky), "chars")
