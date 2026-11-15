"""Task 5: are feature effects consistent across regions?"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from common import (REGIONS, REGION_COLORS, REGION_MARKERS, FEATURES, TARGET,
                    load, anomalies, save, plt, ROOT)

# top 3 by |r| + the other two shortlist members
FEATS = ["tmax_Aug", "cdd_Sep1_Oct15", "mid_greenup_p50", "precip_AugSep", "cdd_Sep"]
UNITS = {"tmax_Aug": "°C", "cdd_Sep1_Oct15": "°C·days", "mid_greenup_p50": "day",
         "precip_AugSep": "mm", "cdd_Sep": "°C·days"}
df = load()
A = anomalies(df, FEATURES + [TARGET])

rows, tests = [], []
for f in FEATS:
    for r in REGIONS:
        s = A[A.region == r]
        m = smf.ols(f"{TARGET} ~ {f}", s).fit()
        rows.append({"feature": f, "region": r, "slope": m.params[f], "se": m.bse[f],
                     "p": m.pvalues[f], "r": np.corrcoef(s[f], s[TARGET])[0, 1], "n": len(s)})
    pooled = smf.ols(f"{TARGET} ~ {f}", A).fit(cov_type="cluster", cov_kwds={"groups": A.year})
    rows.append({"feature": f, "region": "Pooled", "slope": pooled.params[f], "se": pooled.bse[f],
                 "p": pooled.pvalues[f], "r": np.corrcoef(A[f], A[TARGET])[0, 1], "n": len(A)})
    # heterogeneity: does letting the slope vary by region improve fit?
    m0 = smf.ols(f"{TARGET} ~ C(region) + {f}", A).fit()
    m1 = smf.ols(f"{TARGET} ~ C(region) * {f}", A).fit()
    F, p, _ = m1.compare_f_test(m0)
    tests.append({"feature": f, "F_slope_heterogeneity": F, "p": p})
S = pd.DataFrame(rows)
T = pd.DataFrame(tests)
S["slope_pm_se"] = S.apply(lambda x: f"{x.slope:+.3f} ± {x.se:.3f}", axis=1)
print(S.pivot(index="feature", columns="region", values="slope_pm_se")[REGIONS + ["Pooled"]].to_string())
print("\nTest of equal slopes across regions (F-test, interaction vs shared slope):")
print(T.round(3).to_string(index=False))
S.round(4).to_csv(ROOT / "feature_slopes_by_region.csv", index=False)
T.round(4).to_csv(ROOT / "feature_slope_heterogeneity.csv", index=False)

# ---- figure -------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(12, 7.2))
for ax, f in zip(axes.flat, FEATS):
    for r in REGIONS:
        s = A[A.region == r]
        c = REGION_COLORS[r]
        ax.scatter(s[f], s[TARGET], color=c, marker=REGION_MARKERS[r], s=26,
                   edgecolor="white", lw=0.5, label=r, zorder=3)
        b = S[(S.feature == f) & (S.region == r)].iloc[0]
        xs = np.array([s[f].min(), s[f].max()])
        ax.plot(xs, b.slope * xs + (s[TARGET].mean() - b.slope * s[f].mean()), color=c, lw=1.5)
    b = S[(S.feature == f) & (S.region == "Pooled")].iloc[0]
    xs = np.array([A[f].min(), A[f].max()])
    ax.plot(xs, b.slope * xs, color="#1a1a19", lw=2.2, ls="--", label="Pooled", zorder=4)
    h = A[A.helene_2024 == 1]
    ax.scatter(h[f], h[TARGET], s=90, facecolor="none", edgecolor="#c0262d", lw=1.5,
               zorder=5, label="Helene 2024")
    ax.axhline(0, color="#b5b5b0", lw=0.7); ax.axvline(0, color="#b5b5b0", lw=0.7)
    p = T.loc[T.feature == f, "p"].iloc[0]
    ax.set_title(f"{f}\npooled slope {b.slope:+.3f} days per {UNITS[f]}; equal-slope test p={p:.2f}", fontsize=9)
    ax.set_xlabel(f"{f} anomaly ({UNITS[f]})")
    ax.set_ylabel("Peak-date anomaly (days)")
leg_ax = axes.flat[-1]
leg_ax.axis("off")
hnd, lab = axes.flat[0].get_legend_handles_labels()
leg_ax.legend(hnd, lab, loc="center", fontsize=10, title="Region / line")
fig.suptitle("Within-region anomalies: per-region fits (solid) vs one shared slope (dashed)",
             fontsize=11)
fig.tight_layout()
save(fig, "05_feature_consistency.png")

# slope forest plot
fig, axes = plt.subplots(1, len(FEATS), figsize=(13, 2.8), sharey=True)
ylabels = REGIONS + ["Pooled"]
for ax, f in zip(axes, FEATS):
    s = S[S.feature == f].set_index("region").loc[ylabels]
    for i, (reg, row) in enumerate(s.iterrows()):
        col = REGION_COLORS.get(reg, "#1a1a19")
        ax.errorbar(row.slope, i, xerr=1.96 * row.se, fmt=REGION_MARKERS.get(reg, "o"),
                    color=col, ms=6, capsize=3, lw=1.5)
    ax.axvline(0, color="#8a8a85", lw=0.8)
    ax.set_title(f, fontsize=9)
    ax.set_xlabel(f"Slope (days per {UNITS[f]})")
    ax.grid(axis="y", visible=False)
axes[0].set_yticks(range(len(ylabels)), ylabels)
axes[0].invert_yaxis()
fig.suptitle("Per-region slopes with 95% CI (anomaly regressions)", fontsize=11, y=1.04)
save(fig, "05_slope_forest.png")
