"""Task 4: redundancy clusters, VIFs, and two feature shortlists."""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor

from common import FEATURES, TARGET, TIMING, load, anomalies, ROOT

THRESH = 0.8
df = load()
A = anomalies(df, FEATURES + [TARGET])
A["year"] = df.year - 2013
FEATS = FEATURES + ["year"]
C = A[FEATS].corr()
r_t = A[FEATS].corrwith(A[TARGET])

# ---- clusters: connected components of the |r| > 0.8 graph ---------------
unseen, clusters = set(FEATS), []
while unseen:
    stack, comp = [unseen.pop()], set()
    while stack:
        f = stack.pop(); comp.add(f)
        nb = {g for g in FEATS if g != f and abs(C.loc[f, g]) > THRESH and g in unseen}
        unseen -= nb; stack += nb
    clusters.append(sorted(comp, key=lambda f: -abs(r_t[f])))
clusters.sort(key=lambda c: -len(c))
print(f"Redundancy clusters (|r| > {THRESH} between anomalies):")
rows = []
for i, cl in enumerate(clusters, 1):
    best = cl[0]
    best_pre = next((f for f in cl if TIMING[f] == "pre-peak"), None)
    if len(cl) > 1:
        print(f"  cluster {i}: " + ", ".join(f"{f} (r={r_t[f]:+.2f}, {TIMING[f]})" for f in cl))
        print(f"     keep (any): {best};  keep (forecast-safe): {best_pre}")
    for f in cl:
        rows.append({"cluster": i, "feature": f, "timing": TIMING[f], "r_target": r_t[f],
                     "cluster_size": len(cl), "keep_any": f == best, "keep_pre_peak": f == best_pre})
pd.DataFrame(rows).round(3).to_csv(ROOT / "feature_clusters.csv", index=False)
print("  all other features are singletons")
pairs = [(a, b, C.loc[a, b]) for i, a in enumerate(FEATS) for b in FEATS[i + 1:]
         if abs(C.loc[a, b]) > THRESH]
print("  pairs above threshold:", ", ".join(f"{a}~{b} ({r:+.2f})" for a, b, r in pairs))

# ---- shortlists ---------------------------------------------------------
SHORTLISTS = {
    "forecast_safe": {
        "tmax_Aug": "Strongest pre-peak signal (r=+0.39); hot Augusts -> later greendown. Stands in for vpd_AugSep (r=0.72 with it).",
        "mid_greenup_p50": "Late spring green-up -> late autumn (r=+0.32); nearly uncorrelated with the weather features, so adds unique information.",
        "precip_AugSep": "Wet late summer -> earlier (r=-0.31); keeps its effect jointly; chosen over pdsi_AugSep, which is weaker once tmax_Aug is in.",
        "cdd_Sep": "Only pre-peak cold-forcing term (r=-0.26); best forecast-safe member of the temperature cluster and the hook for a degree-day model.",
    },
    "explanatory": {
        "tmax_Aug": "Same as above.",
        "mid_greenup_p50": "Same as above.",
        "precip_AugSep": "Same as above.",
        "cdd_Sep1_Oct15": "Best member of the cold cluster (r=-0.38). The window runs past the mean peak (~Oct 5), so it is explanatory, not forecast-safe.",
    },
}


def vif_table(fs):
    X = sm.add_constant(A[fs])
    return pd.Series([variance_inflation_factor(X.values, i + 1) for i in range(len(fs))], index=fs)


Z = (A[FEATS + [TARGET]] - A[FEATS + [TARGET]].mean()) / A[FEATS + [TARGET]].std()
out = []
for name, picks in SHORTLISTS.items():
    fs = list(picks)
    vif = vif_table(fs)
    m = sm.OLS(Z[TARGET], sm.add_constant(Z[fs])).fit(cov_type="cluster", cov_kwds={"groups": df.year})
    print(f"\nShortlist: {name}  (in-sample R2 on anomalies = {m.rsquared:.3f}, adj {m.rsquared_adj:.3f}; "
          "screening only, not a validated model)")
    for f in fs:
        print(f"  {f:16s} r={r_t[f]:+.2f}  VIF={vif[f]:.2f}  std.beta={m.params[f]:+.2f} (p={m.pvalues[f]:.3f})  {picks[f]}")
        out.append({"shortlist": name, "feature": f, "timing": TIMING[f], "r_target": r_t[f],
                    "VIF": vif[f], "std_beta_joint": m.params[f], "p_joint_clustered": m.pvalues[f],
                    "joint_R2": m.rsquared, "reason": picks[f]})
pd.DataFrame(out).round(3).to_csv(ROOT / "feature_shortlists.csv", index=False)

# VIF on the full candidate set, for reference
print("\nVIF, all year-varying features together (shows the redundancy):")
print(vif_table(FEATS).round(1).sort_values(ascending=False).to_string())
print("saved analysis/feature_clusters.csv, analysis/feature_shortlists.csv")
