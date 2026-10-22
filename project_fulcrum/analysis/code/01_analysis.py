"""
Project Fulcrum - model fitting, diagnostics, cross-validation, figure

VERBATIM EXPORT. These are the cells that actually produced the published
outputs, in execution order, not a tidied retelling. They originally ran in one
persistent Python kernel, so run this file top-to-bottom in a single process.

Cells that raised are omitted; one partial include is marked inline.
Environment: python 3 with pandas, numpy, scipy, statsmodels, matplotlib.
statsmodels is the only non-obvious requirement (scikit-learn is NOT needed -
the stratified split is done with numpy, seeded rng(42)).

Run 00_fetch_inputs.py first - it puts antibody_agonism.csv beside this file.
"""

# ============================================================================
# COMPATIBILITY SHIM - not part of the original analysis.
# The cells below ran inside an agent platform that supplied two things this
# standalone script does not have:
#   host.artifact_path(id) -> resolved path of the input CSV
#   apply_figure_style() / META_GREY / panel_letter() -> house figure style
# These stubs restore just enough for the cells to run unchanged. The figure
# will differ cosmetically from the published PNG (fonts, tick direction);
# every number it reports is unaffected.
# ============================================================================
import os as _os

class _Host:
    @staticmethod
    def artifact_path(_version_id):
        p = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "antibody_agonism.csv")
        if not _os.path.exists(p):
            raise FileNotFoundError("run 00_fetch_inputs.py first: " + p)
        return p
host = _Host()

META_GREY = "#8b949e"
def apply_figure_style():
    import matplotlib as _m
    _m.rcParams.update({"figure.dpi": 150, "savefig.dpi": 300,
                        "axes.spines.top": False, "axes.spines.right": False,
                        "xtick.direction": "out", "ytick.direction": "out",
                        "legend.frameon": False, "font.size": 9})

def panel_letter(ax, letter):
    ax.text(-0.17, 1.06, letter, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="bottom", ha="left")
# ======================= end shim ===========================================

# ---- cell 2 ----
import pandas as pd, numpy as np
df = pd.read_csv(host.artifact_path("a01d1688-9c09-4653-bdc4-49d9c6f339bf"))
print(df.shape, list(df.columns))
print(df.isna().sum().sum(), "missing")
print(df.describe().T.round(3).to_string())
print("\n--- site x plate counts ---")
print(pd.crosstab(df.site, df.plate).to_string())
print("\n--- epitope_dist per site ---")
print(df.groupby('site').epitope_dist_nm.agg(['min','max','nunique']).to_string())
print("\n--- zeros in outcome ---")
print((df.agonism_pct==0).sum(), "exact zeros;", (df.agonism_pct<=0.5).sum(), "<=0.5")

# ---- cell 3 ----
print("--- mean agonism by site ---")
print(df.groupby('site').agonism_pct.agg(['mean','median','std','count',lambda s:(s==0).mean()]).rename(columns={'<lambda_0>':'frac_zero'}).round(3).to_string())
print("\n--- mean agonism by plate ---")
print(df.groupby('plate').agonism_pct.agg(['mean','median','std','count',lambda s:(s==0).mean()]).rename(columns={'<lambda_0>':'frac_zero'}).round(3).to_string())
print("\n--- predictor correlations ---")
cols=['epitope_dist_nm','log10_kon','log10_koff','log10_KD','receptor_density']
print(df[cols].corr().round(3).to_string())
print("\n--- KD identity check ---")
print("max |KD-(koff-kon)| =", (df.log10_KD-(df.log10_koff-df.log10_kon)).abs().max())
print("\n--- receptor_density by plate ---")
print(df.groupby('plate').receptor_density.agg(['mean','std']).round(3).to_string())
print("\n--- kinetics by plate (confound check) ---")
print(df.groupby('plate')[['log10_kon','log10_koff']].mean().round(3).to_string())

# ---- cell 5 ----
import statsmodels.api as sm, statsmodels.formula.api as smf
print("statsmodels", sm.__version__)

rng = np.random.default_rng(42)
test_idx = []
for s, g in df.groupby('site'):
    k = int(round(0.25*len(g)))
    test_idx += list(rng.choice(g.index.values, size=k, replace=False))
te = df.loc[sorted(test_idx)].copy()
tr = df.drop(index=test_idx).copy()
print("train", tr.shape, "test", te.shape)
print(tr.site.value_counts().sort_index().to_dict(), te.site.value_counts().sort_index().to_dict())

# ---- cell 6 ----
def binned(d, col, q=6, by=None):
    d = d.copy()
    d['b'] = pd.qcut(d[col], q, duplicates='drop')
    g = ['b'] if by is None else [by,'b']
    return d.groupby(g, observed=True).agg(n=('agonism_pct','size'), mean=('agonism_pct','mean'),
                                           med=('agonism_pct','median'), fz=('agonism_pct', lambda s:(s==0).mean()),
                                           x=(col,'mean')).round(2)
for c in ['log10_kon','log10_koff','log10_KD','receptor_density']:
    print(f"=== agonism vs {c} (train) ===")
    print(binned(tr,c).to_string()); print()

# ---- cell 7 ----
print("=== koff bins within each site (train) ===")
print(binned(tr,'log10_koff',q=4,by='site').to_string())
print("\n=== kon bins within each site (train) ===")
print(binned(tr,'log10_kon',q=4,by='site').to_string())
print("\n=== mean/variance relation (train, by fitted-ish bins of site) ===")
print(tr.groupby('site').agonism_pct.agg(['mean','var','std']).assign(cv=lambda d:d['std']/d['mean']).round(3).to_string())

# ---- cell 8 ----
tr['siteC']=tr.site.astype('category'); tr['plateC']=tr.plate.astype('category')
te['siteC']=te.site.astype('category'); te['plateC']=te.plate.astype('category')
pos = tr[tr.agonism_pct>0].copy(); pos['logy']=np.log(pos.agonism_pct)
print("positives in train:", len(pos), "of", len(tr))

def rep(m,label):
    print(f"{label:52s} R2={m.rsquared:.3f} adjR2={m.rsquared_adj:.3f} AIC={m.aic:.1f} df={int(m.df_model)}")
    return m

print("\n--- RAW scale (train, all rows incl zeros) ---")
rep(smf.ols("agonism_pct ~ epitope_dist_nm + log10_kon + log10_koff", tr).fit(), "naive linear (dist,kon,koff)")
rep(smf.ols("agonism_pct ~ siteC + log10_kon + log10_koff", tr).fit(), "site factor + kon + koff")
rep(smf.ols("agonism_pct ~ siteC + log10_kon + log10_koff + I(log10_koff**2)", tr).fit(), "site + kon + koff + koff^2")
rep(smf.ols("agonism_pct ~ siteC*(log10_kon + log10_koff + I(log10_koff**2))", tr).fit(), "site x (kon + koff + koff^2)  [interactions]")

print("\n--- LOG scale (train, positives only) ---")
rep(smf.ols("logy ~ epitope_dist_nm + log10_kon + log10_koff", pos).fit(), "naive linear (dist,kon,koff)")
rep(smf.ols("logy ~ siteC + log10_kon + log10_koff", pos).fit(), "site factor + kon + koff")
m_log = rep(smf.ols("logy ~ siteC + log10_kon + log10_koff + I(log10_koff**2)", pos).fit(), "site + kon + koff + koff^2")
m_logx = rep(smf.ols("logy ~ siteC*(log10_kon + log10_koff + I(log10_koff**2))", pos).fit(), "site x (kon + koff + koff^2)  [interactions]")
print("\nF-test: do interactions add anything on the LOG scale?")
print(m_logx.compare_f_test(m_log))

# ---- cell 9 ----
t = tr.copy()
t['kb'] = pd.qcut(t.log10_koff, 4, labels=False)
t['nb'] = pd.qcut(t.log10_kon, 3, labels=False)
g = t.groupby(['site','kb','nb']).agonism_pct.agg(n='size', mean='mean', sd='std', fz=lambda s:(s==0).mean())
clean = g[(g.fz==0)&(g.n>=6)].copy()
print("strata with zero censoring and n>=6:", len(clean))
print(clean.round(2).to_string())
import scipy.stats as st
lm = st.linregress(np.log(clean['mean']), np.log(clean['sd']))
print(f"\nlog(sd) ~ log(mean):  slope p = {lm.slope:.3f}  (95% CI {lm.slope-1.96*lm.stderr:.3f}..{lm.slope+1.96*lm.stderr:.3f})  r2={lm.rvalue**2:.3f}")
print("  p=0 -> constant additive noise ; p=1 -> constant CV (multiplicative)")
print(f"  implied CV at p=1: {np.exp(lm.intercept):.3f}" if abs(lm.slope-1)<0.3 else "")

# ---- cell 10 ----
# NOTE: cell 10 raised AttributeError on its line 17 (GLMResults has no
# NOTE: .compare_f_test). Lines 1-16 had already run and defined qp/info/m0-m3,
# NOTE: which later cells use, so they are kept here. The quasi-likelihood
# NOTE: F-test it was reaching for is implemented as quasiF() in cell 11 below.
Log = sm.families.links.Log()
def qp(f, d=tr):
    m = smf.glm(f, d, family=sm.families.Poisson(Log)).fit(scale='X2')
    return m
def info(m, label):
    dev=m.deviance; df=m.df_resid
    print(f"{label:50s} scale={m.scale:7.3f} dev/df={dev/df:6.3f} pseudoR2={1-dev/m.null_deviance:.3f} df_m={int(m.df_model)}")
    return m

print("--- quasi-Poisson, log link (train, ALL 375 rows incl. zeros) ---")
m0=info(qp("agonism_pct ~ epitope_dist_nm + log10_kon + log10_koff"), "naive linear (dist,kon,koff)")
m1=info(qp("agonism_pct ~ siteC + log10_kon + log10_koff"), "site + kon + koff")
m2=info(qp("agonism_pct ~ siteC + log10_kon + log10_koff + I(log10_koff**2)"), "site + kon + koff + koff^2")
m3=info(qp("agonism_pct ~ siteC*(log10_kon + log10_koff + I(log10_koff**2))"), "site x (kon+koff+koff^2) [interactions]")

print("\nF-test (quasi-likelihood) interactions vs additive:")

# ---- cell 11 ----
import scipy.stats as st
def quasiF(big, small, label):
    ddf = small.df_resid - big.df_resid
    phi = big.pearson_chi2/big.df_resid
    F = ((small.deviance-big.deviance)/ddf)/phi
    p = st.f.sf(F, ddf, big.df_resid)
    print(f"{label:38s} F({ddf},{int(big.df_resid)})={F:7.3f}  p={p:.3g}")
    return p

quasiF(m3, m2, "site x (kon+koff+koff^2)")
m2a=qp("agonism_pct ~ siteC*log10_kon + log10_koff + I(log10_koff**2)")
m2b=qp("agonism_pct ~ siteC*(log10_koff + I(log10_koff**2)) + log10_kon")
quasiF(m2a, m2, "site x kon only")
quasiF(m2b, m2, "site x koff curve only")

print("\n--- other candidate terms on the additive log-link model ---")
m2c=qp("agonism_pct ~ siteC + log10_kon + log10_koff + I(log10_koff**2) + receptor_density")
quasiF(m2c, m2, "+ receptor_density")
m2d=qp("agonism_pct ~ siteC + log10_kon + log10_koff + I(log10_koff**2) + plateC")
quasiF(m2d, m2, "+ plate factor")
m2e=qp("agonism_pct ~ siteC + log10_kon + I(log10_kon**2) + log10_koff + I(log10_koff**2)")
quasiF(m2e, m2, "+ kon^2 (is kon curved too?)")
m2f=qp("agonism_pct ~ siteC + log10_kon + log10_koff + I(log10_koff**2) + log10_kon:log10_koff")
quasiF(m2f, m2, "+ kon:koff")

# ---- cell 12 ----
cands = {
 "A site+kon+koff+koff^2"          : "agonism_pct ~ siteC + log10_kon + log10_koff + I(log10_koff**2)",
 "B site+KD+KD^2"                  : "agonism_pct ~ siteC + log10_KD + I(log10_KD**2)",
 "C site+KD+KD^2+koff+koff^2"      : "agonism_pct ~ siteC + log10_KD + I(log10_KD**2) + log10_koff + I(log10_koff**2)",
 "D site+full quadratic(kon,koff)" : "agonism_pct ~ siteC + log10_kon + log10_koff + I(log10_kon**2) + I(log10_koff**2) + log10_kon:log10_koff",
 "E site+KD+KD^2+koff"             : "agonism_pct ~ siteC + log10_KD + I(log10_KD**2) + log10_koff",
}
fits={}
for k,f in cands.items():
    m=qp(f); fits[k]=m
    print(f"{k:34s} dev/df={m.deviance/m.df_resid:6.3f} pseudoR2={1-m.deviance/m.null_deviance:.4f} df_m={int(m.df_model)}")

print("\nNested quasi-F tests:")
quasiF(fits["D site+full quadratic(kon,koff)"], fits["B site+KD+KD^2"], "full surface vs KD-only")
quasiF(fits["C site+KD+KD^2+koff+koff^2"], fits["B site+KD+KD^2"], "+ koff curve on top of KD")
quasiF(fits["D site+full quadratic(kon,koff)"], fits["C site+KD+KD^2+koff+koff^2"], "surface vs KD+koff (C)")
quasiF(fits["C site+KD+KD^2+koff+koff^2"], fits["E site+KD+KD^2+koff"], "+ koff^2 on top of KD+koff")

print("\n--- coefficient pattern of the full surface D (does it match a KD quadratic?) ---")
d=fits["D site+full quadratic(kon,koff)"]
for n in ['I(log10_kon ** 2)','I(log10_koff ** 2)','log10_kon:log10_koff']:
    print(f"  {n:24s} {d.params[n]: .4f}")
print("  a KD^2 term predicts kon^2 = koff^2 = a and kon:koff = -2a")

# ---- cell 13 ----
C = "agonism_pct ~ siteC + log10_KD + I(log10_KD**2) + log10_koff + I(log10_koff**2)"
mC = qp(C)
tests = {
 "+ plate factor"              : C+" + plateC",
 "+ receptor_density"          : C+" + receptor_density",
 "site x KD curve"             : "agonism_pct ~ siteC*(log10_KD + I(log10_KD**2)) + log10_koff + I(log10_koff**2)",
 "site x koff curve"           : "agonism_pct ~ siteC*(log10_koff + I(log10_koff**2)) + log10_KD + I(log10_KD**2)",
 "site x all kinetics"         : "agonism_pct ~ siteC*(log10_KD + I(log10_KD**2) + log10_koff + I(log10_koff**2))",
 "+ cubic KD"                  : C+" + I(log10_KD**3)",
 "+ cubic koff"                : C+" + I(log10_koff**3)",
 "+ KD:koff"                   : C+" + log10_KD:log10_koff",
}
for k,f in tests.items():
    quasiF(qp(f), mC, k)

print("\n--- with plate in the model, retest the interactions ---")
Cp = C+" + plateC"
mCp = qp(Cp)
quasiF(qp("agonism_pct ~ siteC*(log10_koff + I(log10_koff**2)) + log10_KD + I(log10_KD**2) + plateC"), mCp, "site x koff | plate in model")
quasiF(qp("agonism_pct ~ siteC*(log10_KD + I(log10_KD**2)) + log10_koff + I(log10_koff**2) + plateC"), mCp, "site x KD   | plate in model")
quasiF(qp(Cp+" + receptor_density"), mCp, "receptor_density | plate in model")
print(f"\nmodel C+plate: dev/df={mCp.deviance/mCp.df_resid:.3f} pseudoR2={1-mCp.deviance/mCp.null_deviance:.4f}")

# ---- cell 14 ----
mKD = qp("agonism_pct ~ siteC*(log10_KD + I(log10_KD**2)) + log10_koff + I(log10_koff**2) + plateC")
p = mKD.params
print("per-site quadratic in log10_KD  ->  vertex (optimal log10 KD) and curvature")
for s in [1,2,3,4,5]:
    b1 = p['log10_KD'] + (p.get(f'siteC[T.{s}]:log10_KD',0) if s!=1 else 0)
    b2 = p['I(log10_KD ** 2)'] + (p.get(f'siteC[T.{s}]:I(log10_KD ** 2)',0) if s!=1 else 0)
    vtx = -b1/(2*b2) if b2!=0 else np.nan
    print(f"  site {s}: b1={b1: .3f} b2={b2: .4f}  vertex={vtx: .2f}  {'(concave=optimum)' if b2<0 else '(CONVEX=no interior optimum)'}")

print("\n--- empirical: mean agonism by KD quartile within site (train) ---")
print(binned(tr,'log10_KD',q=4,by='site').to_string())

print("\n--- is the interaction driven by the near-floor sites 1 and 5? ---")
sub = tr[~tr.site.isin([1])].copy(); sub['siteC']=sub.site.astype('category'); sub['plateC']=sub.plate.astype('category')
base = smf.glm(C+" + plateC", sub, family=sm.families.Poisson(Log)).fit(scale='X2')
full = smf.glm("agonism_pct ~ siteC*(log10_KD + I(log10_KD**2)) + log10_koff + I(log10_koff**2) + plateC", sub, family=sm.families.Poisson(Log)).fit(scale='X2')
quasiF(full, base, "site x KD, sites 2-5 only")
sub2 = tr[~tr.site.isin([1,5])].copy(); sub2['siteC']=sub2.site.astype('category'); sub2['plateC']=sub2.plate.astype('category')
base2 = smf.glm(C+" + plateC", sub2, family=sm.families.Poisson(Log)).fit(scale='X2')
full2 = smf.glm("agonism_pct ~ siteC*(log10_KD + I(log10_KD**2)) + log10_koff + I(log10_koff**2) + plateC", sub2, family=sm.families.Poisson(Log)).fit(scale='X2')
quasiF(full2, base2, "site x KD, sites 2-4 only")

# ---- cell 15 ----
mu_hat = mCp.fittedvalues.values
phi = mCp.pearson_chi2/mCp.df_resid
print(f"dispersion phi={phi:.3f}; observed train zero fraction={(tr.agonism_pct==0).mean():.3f}")
F_obs = ((mCp.deviance-mKD.deviance)/(mCp.df_resid-mKD.df_resid))/(mKD.pearson_chi2/mKD.df_resid)
print(f"observed site x KD quasi-F = {F_obs:.3f}")

sim = tr.copy()
rng = np.random.default_rng(7)
Fs=[]; zfs=[]
f_add = C+" + plateC"
f_int = "agonism_pct ~ siteC*(log10_KD + I(log10_KD**2)) + log10_koff + I(log10_koff**2) + plateC"
for i in range(300):
    y = np.maximum(0.0, mu_hat + rng.normal(0, np.sqrt(phi*mu_hat)))
    sim['agonism_pct'] = y
    zfs.append((y==0).mean())
    try:
        a = smf.glm(f_add, sim, family=sm.families.Poisson(Log)).fit(scale='X2')
        b = smf.glm(f_int, sim, family=sm.families.Poisson(Log)).fit(scale='X2')
        ddf = a.df_resid-b.df_resid
        Fs.append(((a.deviance-b.deviance)/ddf)/(b.pearson_chi2/b.df_resid))
    except Exception: pass
Fs=np.array(Fs)
print(f"\nsimulated under NO interaction + floor at 0  (n={len(Fs)} reps)")
print(f"  simulated zero fraction: mean={np.mean(zfs):.3f}  (observed {(tr.agonism_pct==0).mean():.3f})")
print(f"  spurious site x KD F: median={np.median(Fs):.2f}  90th pct={np.percentile(Fs,90):.2f}  95th={np.percentile(Fs,95):.2f}  max={Fs.max():.2f}")
print(f"  P(F_sim >= F_obs={F_obs:.2f}) = {(Fs>=F_obs).mean():.3f}")

# ---- cell 16 ----
from collections import defaultdict
def cv_eval(formula, kind, data=tr, K=5, seed=11, drop_plate_at_predict=False):
    rng=np.random.default_rng(seed); fold=np.empty(len(data),int)
    d=data.reset_index(drop=True)
    for s,g in d.groupby('site'):
        idx=g.index.values.copy(); rng.shuffle(idx)
        fold[idx]=np.arange(len(idx))%K
    yh=np.empty(len(d)); 
    for k in range(K):
        trn=d[fold!=k].copy(); val=d[fold==k].copy()
        for c in ('siteC','plateC'):
            trn[c]=pd.Categorical(trn[c.replace('C','')], categories=sorted(d[c.replace('C','')].unique()))
            val[c]=pd.Categorical(val[c.replace('C','')], categories=sorted(d[c.replace('C','')].unique()))
        if kind=='ols': m=smf.ols(formula,trn).fit()
        else: m=smf.glm(formula,trn,family=sm.families.Poisson(Log)).fit(scale='X2')
        v=val.copy()
        if drop_plate_at_predict: pass
        yh[fold==k]=np.clip(m.predict(v).values,0,None)
    y=d.agonism_pct.values
    return dict(RMSE=np.sqrt(np.mean((y-yh)**2)), MAE=np.mean(np.abs(y-yh)),
                R2=1-np.sum((y-yh)**2)/np.sum((y-y.mean())**2), rho=st.spearmanr(y,yh).statistic)

cands = [
 ("1 OLS naive linear (dist,kon,koff)", "agonism_pct ~ epitope_dist_nm + log10_kon + log10_koff", 'ols'),
 ("2 OLS site + kon + koff",            "agonism_pct ~ siteC + log10_kon + log10_koff", 'ols'),
 ("3 OLS site + quad kinetics + plate", "agonism_pct ~ siteC + log10_KD + I(log10_KD**2) + log10_koff + I(log10_koff**2) + plateC", 'ols'),
 ("4 OLS site x quad kinetics + plate", "agonism_pct ~ siteC*(log10_KD + I(log10_KD**2) + log10_koff + I(log10_koff**2)) + plateC", 'ols'),
 ("5 GLM log: additive C",              C, 'glm'),
 ("6 GLM log: additive C + plate",      C+" + plateC", 'glm'),
 ("7 GLM log: C + plate + site:KD",     "agonism_pct ~ siteC*(log10_KD + I(log10_KD**2)) + log10_koff + I(log10_koff**2) + plateC", 'glm'),
 ("8 GLM log: site x all kin + plate",  "agonism_pct ~ siteC*(log10_KD + I(log10_KD**2) + log10_koff + I(log10_koff**2)) + plateC", 'glm'),
]
rows=[]
for name,f,kind in cands:
    r=cv_eval(f,kind); r['model']=name; rows.append(r)
cvt=pd.DataFrame(rows).set_index('model')[['RMSE','MAE','R2','rho']].round(3)
print("5-fold CV within TRAIN (stratified by site), raw-scale metrics:")
print(cvt.to_string())

# ---- cell 17 ----
final = smf.glm(C+" + plateC", tr, family=sm.families.Poisson(Log)).fit(scale='X2')
print(final.summary().tables[1])
print(f"\ndispersion={final.pearson_chi2/final.df_resid:.2f}  pseudoR2={1-final.deviance/final.null_deviance:.3f}")

def optima(p):
    b1,b2 = p['log10_KD'], p['I(log10_KD ** 2)']
    c1,c2 = p['log10_koff'], p['I(log10_koff ** 2)']
    # log mu = .. + b1*(koff-kon) + b2*(koff-kon)^2 + c1*koff + c2*koff^2
    KDs = -b1/(2*b2)                      # d/dkon = 0
    koffs = -(c1 + b1 + 2*b2*KDs)/(2*c2)  # d/dkoff = 0 given KD=KDs
    return KDs, koffs, koffs-KDs          # KD*, koff*, kon*
KD_s,koff_s,kon_s = optima(final.params)
print(f"\njoint optimum of the kinetic surface:")
print(f"  log10_KD*   = {KD_s:7.3f}  (KD   = {10**KD_s:.3g} M)")
print(f"  log10_koff* = {koff_s:7.3f}  (koff = {10**koff_s:.3g} 1/s ; residence time {1/10**koff_s:.0f} s)")
print(f"  log10_kon*  = {kon_s:7.3f}  (kon  = {10**kon_s:.3g} 1/M/s)")

# bootstrap CIs for the optima
rng=np.random.default_rng(3); boot=[]
for i in range(600):
    bs = tr.sample(len(tr), replace=True, random_state=int(rng.integers(1e9))).copy()
    bs['siteC']=pd.Categorical(bs.site,categories=[1,2,3,4,5]); bs['plateC']=pd.Categorical(bs.plate,categories=[1,2,3,4,5,6])
    try:
        m=smf.glm(C+" + plateC", bs, family=sm.families.Poisson(Log)).fit(scale='X2')
        boot.append(optima(m.params))
    except Exception: pass
boot=np.array(boot)
for nm,i in [('log10_KD*',0),('log10_koff*',1),('log10_kon*',2)]:
    lo,hi=np.percentile(boot[:,i],[2.5,97.5])
    print(f"  {nm:12s} 95% CI  {lo:7.3f} .. {hi:7.3f}   (n={len(boot)} boots)")

# ---- cell 18 ----
plate_cols=[f'plateC[T.{i}]' for i in range(2,7)]
def pred_avg_plate(m, d):
    """linear predictor averaged over the 6 plates (prospective: unknown plate)"""
    lp=[]
    for pl in [1,2,3,4,5,6]:
        dd=d.copy(); dd['plateC']=pd.Categorical([pl]*len(dd),categories=[1,2,3,4,5,6])
        lp.append(np.log(np.clip(m.predict(dd).values,1e-12,None)))
    return np.exp(np.mean(lp,axis=0))

def metrics(y,yh):
    return dict(RMSE=np.sqrt(np.mean((y-yh)**2)), MAE=np.mean(np.abs(y-yh)),
                R2=1-np.sum((y-yh)**2)/np.sum((y-y.mean())**2), rho=st.spearmanr(y,yh).statistic)

yte=te.agonism_pct.values
naive=smf.ols("agonism_pct ~ epitope_dist_nm + log10_kon + log10_koff", tr).fit()
m_noplate=smf.glm(C, tr, family=sm.families.Poisson(Log)).fit(scale='X2')
rows={
 "naive OLS linear (dist,kon,koff)": metrics(yte, np.clip(naive.predict(te).values,0,None)),
 "GLM log-link additive, no plate" : metrics(yte, np.clip(m_noplate.predict(te).values,0,None)),
 "GLM log-link + plate (known)"    : metrics(yte, np.clip(final.predict(te).values,0,None)),
 "GLM log-link + plate (averaged)" : metrics(yte, pred_avg_plate(final,te)),
}
print("HELD-OUT performance on the 125 untouched antibodies:")
print(pd.DataFrame(rows).T.round(3).to_string())
print(f"\n(train CV R2 of chosen model was 0.938-0.940; test n=125, {(te.agonism_pct==0).sum()} zeros)")

# ---- cell 19 ----
import matplotlib as mpl, matplotlib.pyplot as plt
apply_figure_style()
full = df.copy(); full['siteC']=pd.Categorical(full.site,categories=[1,2,3,4,5]); full['plateC']=pd.Categorical(full.plate,categories=[1,2,3,4,5,6])
full['pred']=pred_avg_plate(final, full)

fig = plt.figure(figsize=(7.2,6.2))
gs = fig.add_gridspec(2,2, hspace=0.42, wspace=0.30)
ax1,ax2,ax3,ax4 = [fig.add_subplot(gs[i,j]) for i,j in [(0,0),(0,1),(1,0),(1,1)]]
focal='#1f6feb'; grey=META_GREY

# (a) site is a hump
sm_ = df.groupby('site').agonism_pct.agg(['mean','sem','count'])
ax1.errorbar(sm_.index, sm_['mean'], yerr=1.96*sm_['sem'], fmt='o-', color=focal, ms=5, lw=1.4, capsize=3)
for s in sm_.index: ax1.scatter(np.full((df.site==s).sum(),s)+np.random.default_rng(s).normal(0,.06,(df.site==s).sum()),
                                df.loc[df.site==s,'agonism_pct'], s=3, color=grey, alpha=.35, zorder=0)
ax1.set_xticks([1,2,3,4,5]); ax1.set_xlabel('Epitope site (3→15 nm from membrane)')
ax1.set_ylabel('Agonism (% of reference)')
ax1.set_title('Activity peaks at mid-distance, not at the extremes')
ax1.annotate('site 3\n(9 nm)', xy=(3,30.2), xytext=(3.45,52), fontsize=6,
             arrowprops=dict(arrowstyle='-', lw=.6, color='0.3'))

# (b) koff optimum, data + fitted curve at site 3
k=np.linspace(-5,-1,120)
base=pd.DataFrame({'log10_koff':k,'log10_KD':KD_s,'siteC':pd.Categorical([3]*120,categories=[1,2,3,4,5])})
base['plateC']=pd.Categorical([1]*120,categories=[1,2,3,4,5,6]); base['log10_KD']=KD_s
ax2.scatter(df.log10_koff, df.agonism_pct, s=4, color=grey, alpha=.4, zorder=0)
ax2.plot(k, pred_avg_plate(final, base), color=focal, lw=1.6)
ax2.axvline(koff_s, color='0.35', ls=':', lw=.9)
ax2.set_xlabel(r'log$_{10}$ $k_{off}$  (s$^{-1}$)'); ax2.set_ylabel('Agonism (% of reference)')
ax2.set_title('Off-rate has an interior optimum')
ax2.annotate(f'optimum\n{koff_s:.2f}', xy=(koff_s,70), xytext=(-4.85,72), fontsize=6)
ax2.text(.97,.60,'curve: site 3,\nKD at optimum', transform=ax2.transAxes, ha='right', fontsize=6, color=focal)

# (c) response surface in the engineered space
kon_g=np.linspace(4,7,70); koff_g=np.linspace(-5,-1,70); KX,KF=np.meshgrid(kon_g,koff_g)
gd=pd.DataFrame({'log10_kon':KX.ravel(),'log10_koff':KF.ravel()})
gd['log10_KD']=gd.log10_koff-gd.log10_kon
gd['siteC']=pd.Categorical([3]*len(gd),categories=[1,2,3,4,5]); gd['plateC']=pd.Categorical([1]*len(gd),categories=[1,2,3,4,5,6])
Z=pred_avg_plate(final,gd).reshape(KX.shape)
cf=ax3.contourf(KX,KF,Z,levels=12,cmap='viridis')
ax3.scatter(df.log10_kon, df.log10_koff, s=3, color='w', alpha=.35, lw=0)
ax3.plot(kon_s,koff_s,marker='*',ms=11,color='#e8590c',mec='k',mew=.4)
ax3.annotate('predicted\noptimum', xy=(kon_s,koff_s), xytext=(5.4,-1.6), fontsize=6, color='w',
             arrowprops=dict(arrowstyle='->', lw=.7, color='w'))
ax3.set_xlabel(r'log$_{10}$ $k_{on}$  (M$^{-1}$s$^{-1}$)'); ax3.set_ylabel(r'log$_{10}$ $k_{off}$  (s$^{-1}$)')
ax3.set_title('Optimum sits at the panel edge in $k_{on}$')
cb=fig.colorbar(cf,ax=ax3,pad=.02); cb.set_label('Predicted agonism (%)', fontsize=6); cb.ax.tick_params(labelsize=5)

# (d) held-out observed vs predicted
yh_te=pred_avg_plate(final,te)
ax4.scatter(yh_te, yte, s=11, color=focal, alpha=.75, lw=0)
lim=[-3,95]; ax4.plot(lim,lim,color='0.4',ls='--',lw=.8)
ax4.set_xlim(lim); ax4.set_ylim(lim)
ax4.set_xlabel('Predicted agonism (%)'); ax4.set_ylabel('Observed agonism (%)')
ax4.set_title('Held-out antibodies: $R^2$ = 0.94')
ax4.text(.05,.90,'125 antibodies\nnever used in fitting', transform=ax4.transAxes, fontsize=6, va='top')
for a,l in zip([ax1,ax2,ax3,ax4],'abcd'): panel_letter(a,l); a.margins(0.04)
ax1.margins(0.08)
fig.savefig('fulcrum_figure.png', dpi=300, bbox_inches='tight')
print("saved")

# ---- cell 20 ----
ax1.set_title('Activity peaks at mid-distance')
fig.savefig('fulcrum_figure.png', dpi=300, bbox_inches='tight')
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
ov=[(a.get_text(),b.get_text()) for i,(a,ba) in enumerate(texts) for b,bb in texts[i+1:] if ba.overlaps(bb)]
print("overlaps:", ov)

# ---- cell 21 ----
print("--- naive linear fit: what it wrongly concludes ---")
print(naive.summary().tables[1])

print("\n--- site gain multipliers (vs site 1), from final model ---")
for s in [2,3,4,5]:
    c=final.params[f'siteC[T.{s}]']; se=final.bse[f'siteC[T.{s}]']
    print(f"  site {s}: x{np.exp(c):5.2f}  (95% CI {np.exp(c-1.96*se):.2f}-{np.exp(c+1.96*se):.2f})")
print("\n--- plate multipliers (vs plate 1) ---")
for p_ in range(2,7):
    c=final.params[f'plateC[T.{p_}]']; se=final.bse[f'plateC[T.{p_}]']
    print(f"  plate {p_}: x{np.exp(c):5.2f}  (95% CI {np.exp(c-1.96*se):.2f}-{np.exp(c+1.96*se):.2f})")
print(f"  worst->best plate spread: x{np.exp(final.params[plate_cols].max()-final.params[plate_cols].min()):.2f}")

full['in_test']=full.ab_id.isin(te.ab_id)
full['resid_ratio']=full.agonism_pct/full.pred.clip(1e-9)
rank=full.sort_values('pred',ascending=False)[['ab_id','site','log10_kon','log10_koff','log10_KD','plate','agonism_pct','pred','in_test']].copy()
rank['rank']=np.arange(1,len(rank)+1)
print("\n=== TOP 15 by model-predicted agonism (plate effect averaged out) ===")
print(rank.head(15).round(3).to_string(index=False))
rank.round(4).to_csv('fulcrum_ranked_candidates.csv', index=False)
print(f"\nobserved agonism of top-15: median {rank.head(15).agonism_pct.median():.1f}%, range {rank.head(15).agonism_pct.min():.1f}-{rank.head(15).agonism_pct.max():.1f}%")
print("sites represented in top 50:", rank.head(50).site.value_counts().sort_index().to_dict())
