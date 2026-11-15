"""Task 1: target exploration (peak_doy)."""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

from common import (REGIONS, REGION_COLORS, REGION_MARKERS, TARGET, load,
                    doy_to_label, date_axis, save, plt, ROOT)

df = load()
print(f"rows with target: {len(df)}  years {df.year.min()}-{df.year.max()}")

# ---- per-region summary --------------------------------------------------
g = df.groupby("region")[TARGET]
summ = pd.DataFrame({"n": g.size(), "mean": g.mean(), "sd": g.std(),
                     "min": g.min(), "max": g.max()}).loc[REGIONS]
for c in ["mean", "min", "max"]:
    summ[c + "_date"] = summ[c].map(doy_to_label)
summ["range_days"] = summ["max"] - summ["min"]
print("\nPer-region peak_doy summary\n", summ.round(1).to_string())

# ---- trends --------------------------------------------------------------
rows = []
for r in REGIONS:
    s = df[df.region == r]
    m = smf.ols(f"{TARGET} ~ year", s).fit()
    lo, hi = m.conf_int().loc["year"] * 10
    rows.append({"model": r, "slope_days_per_decade": m.params.year * 10,
                 "ci95_lo": lo, "ci95_hi": hi, "p": m.pvalues.year, "n": len(s)})
pooled = smf.ols(f"{TARGET} ~ C(region) + year", df)
for label, fit in [("pooled (OLS SE)", pooled.fit()),
                   ("pooled (year-clustered SE)",
                    pooled.fit(cov_type="cluster", cov_kwds={"groups": df.year}))]:
    lo, hi = fit.conf_int().loc["year"] * 10
    rows.append({"model": label, "slope_days_per_decade": fit.params.year * 10,
                 "ci95_lo": lo, "ci95_hi": hi, "p": fit.pvalues.year, "n": len(df)})
# same trend excluding Helene rows
fit_nh = smf.ols(f"{TARGET} ~ C(region) + year", df[df.helene_2024 == 0]).fit(
    cov_type="cluster", cov_kwds={"groups": df[df.helene_2024 == 0].year})
lo, hi = fit_nh.conf_int().loc["year"] * 10
rows.append({"model": "pooled, no Helene (clustered SE)",
             "slope_days_per_decade": fit_nh.params.year * 10, "ci95_lo": lo,
             "ci95_hi": hi, "p": fit_nh.pvalues.year, "n": int((df.helene_2024 == 0).sum())})
trend = pd.DataFrame(rows)
print("\nLinear trend of peak_doy on year\n", trend.round(3).to_string(index=False))

# ---- outliers ------------------------------------------------------------
z = (df[TARGET] - g.transform("mean")) / g.transform("std")
df["z"] = z
out = df.loc[z.abs() > 2.5, ["region", "year", TARGET, "z", "helene_2024"]]
print("\nOutliers |z| > 2.5:\n", out.round(2).to_string(index=False) if len(out) else "  none")
out2 = df.loc[z.abs() > 2.0, ["region", "year", TARGET, "z", "helene_2024"]]
print("\n(for context) |z| > 2.0:\n", out2.round(2).to_string(index=False))
print("\nHelene rows:\n", df.loc[df.helene_2024 == 1, ["region", "year", TARGET, "z"]].round(2).to_string(index=False))

# ---- MODIS vs VIIRS ------------------------------------------------------
v = df.dropna(subset=["peak_doy_viirs"])
mv = []
for r in REGIONS + ["All"]:
    s = v if r == "All" else v[v.region == r]
    d = s.peak_doy_viirs - s[TARGET]
    rr, pp = stats.pearsonr(s[TARGET], s.peak_doy_viirs)
    mv.append({"region": r, "n": len(s), "bias_viirs_minus_modis": d.mean(),
               "sd_diff": d.std(), "rmse": np.sqrt((d ** 2).mean()), "pearson_r": rr, "p": pp})
# anomaly agreement (within-region) is what matters for year-to-year signal
va = v.assign(m_a=v[TARGET] - v.groupby("region")[TARGET].transform("mean"),
              v_a=v.peak_doy_viirs - v.groupby("region").peak_doy_viirs.transform("mean"))
ra, pa = stats.pearsonr(va.m_a, va.v_a)
mv.append({"region": "All (within-region anomalies)", "n": len(va),
           "bias_viirs_minus_modis": np.nan, "sd_diff": (va.v_a - va.m_a).std(),
           "rmse": np.sqrt(((va.v_a - va.m_a) ** 2).mean()), "pearson_r": ra, "p": pa})
mv = pd.DataFrame(mv)
print("\nMODIS vs VIIRS (2013+)\n", mv.round(3).to_string(index=False))

# Do the two sensors trend differently over the overlap?
vv = v.assign(diff=v.peak_doy_viirs - v[TARGET])
drift = []
for col in [TARGET, "peak_doy_viirs", "diff"]:
    f = smf.ols(f"{col} ~ C(region) + year", vv).fit(cov_type="cluster", cov_kwds={"groups": vv.year})
    lo, hi = f.conf_int().loc["year"] * 10
    drift.append({"series": col, "slope_days_per_decade": f.params.year * 10,
                  "ci95_lo": lo, "ci95_hi": hi, "p": f.pvalues.year})
drift = pd.DataFrame(drift)
print("\nTrend over 2013-2025 overlap (region intercepts, year-clustered SE)\n",
      drift.round(3).to_string(index=False))
drift.round(4).to_csv(ROOT / "modis_viirs_drift.csv", index=False)

summ.round(2).to_csv(ROOT / "target_summary.csv")
trend.round(4).to_csv(ROOT / "target_trend.csv", index=False)
mv.round(4).to_csv(ROOT / "modis_viirs_agreement.csv", index=False)

# ---- figure: small multiples --------------------------------------------
fig, axes = plt.subplots(1, 4, figsize=(13, 3.6), sharey=True)
for ax, r in zip(axes, REGIONS):
    s = df[df.region == r]
    c = REGION_COLORS[r]
    ax.plot(s.year, s[TARGET], color=c, lw=2, marker=REGION_MARKERS[r], ms=4,
            label="MODIS (target)", zorder=3)
    sv = s.dropna(subset=["peak_doy_viirs"])
    ax.plot(sv.year, sv.peak_doy_viirs, color="#7a7a75", lw=1.2, ls="--",
            marker="o", ms=3, mfc="white", label="VIIRS (validation)", zorder=2)
    m = smf.ols(f"{TARGET} ~ year", s).fit()
    ax.plot(s.year, m.fittedvalues, color=c, lw=1, alpha=0.6, ls=":")
    ax.axhline(s[TARGET].mean(), color="#b5b5b0", lw=0.8)
    h = s[s.helene_2024 == 1]
    if len(h):
        ax.scatter(h.year, h[TARGET], s=110, facecolor="none", edgecolor="#c0262d",
                   lw=1.8, zorder=4, label="Helene 2024")
    ax.set_title(f"{r}  (trend {m.params.year*10:+.1f} d/decade)")
    ax.set_xlabel("Year")
axes[0].set_ylabel("Peak-color proxy date (mid-greendown)")
date_axis(axes[0])
handles = {}
for ax in axes:
    for hnd, lab in zip(*ax.get_legend_handles_labels()):
        handles.setdefault(lab, hnd)
fig.legend(handles.values(), handles.keys(), loc="upper center", ncol=3,
           bbox_to_anchor=(0.5, 1.06))
fig.suptitle("Peak-color date by year and region (dotted = linear trend, gray = region mean)",
             y=1.13, fontsize=11)
save(fig, "01_peak_by_year.png")

# ---- figure: MODIS vs VIIRS scatter --------------------------------------
fig, ax = plt.subplots(figsize=(4.6, 4.4))
for r in REGIONS:
    s = v[v.region == r]
    ax.scatter(s[TARGET], s.peak_doy_viirs, color=REGION_COLORS[r],
               marker=REGION_MARKERS[r], s=30, label=r, edgecolor="white", lw=0.5)
lims = [v[[TARGET, "peak_doy_viirs"]].min().min() - 2, v[[TARGET, "peak_doy_viirs"]].max().max() + 2]
ax.plot(lims, lims, color="#8a8a85", lw=0.8, ls="--", label="1:1")
ax.set_xlim(lims); ax.set_ylim(lims)
date_axis(ax, "x"); date_axis(ax, "y")
ax.set_xlabel("MODIS mid-greendown date")
ax.set_ylabel("VIIRS mid-greendown date")
ax.set_title("MODIS vs VIIRS, 2013–2025")
ax.legend(loc="upper left", fontsize=8)
save(fig, "01_modis_vs_viirs.png")
