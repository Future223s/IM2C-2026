"""Step 4: Does location change how a site responds to weather? (interactions for delta)

(a) spatial coherence of yearly anomalies vs distance (exponential decay fit, Mantel permutation test)
(b) per-location slopes of the peak anomaly on forecast-safe weather anomalies, equal-slope F-tests
(c) pooled interaction tests x_a * s_c (s = centered static moderator), year-clustered SE, Holm correction
(d) leave-one-year-out MAE: does any interaction improve a baseline delta model?
"""
import itertools

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.optimize import curve_fit
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import Ridge

from common import *  # noqa

panel, locs, train, ORDER = load()
L = locs.set_index("location").loc[ORDER]

# ---------------------------------------------------------------- (a) coherence vs distance
anom = train.assign(a=train.peak_doy - train.groupby("location").peak_doy.transform("mean"))
W = anom.pivot(index="year", columns="location", values="a")[ORDER]
R = W.corr()


def hav(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * 6371 * np.arcsin(np.sqrt(a))


D = pd.DataFrame([[hav(L.lat[a], L.lon[a], L.lat[b], L.lon[b]) for b in ORDER] for a in ORDER],
                 index=ORDER, columns=ORDER)
iu = np.triu_indices(len(ORDER), 1)
d_, r_ = D.values[iu], R.values[iu]
f_exp = lambda d, L0: np.exp(-d / L0)  # noqa: E731
(L0,), cov = curve_fit(f_exp, d_, r_, p0=[500])
L0_se = np.sqrt(cov[0, 0])
# Mantel test: permute location labels of the distance matrix
rng = np.random.default_rng(0)
obs = np.corrcoef(d_, r_)[0, 1]
perm = []
for _ in range(9999):
    p = rng.permutation(len(ORDER))
    perm.append(np.corrcoef(D.values[np.ix_(p, p)][iu], r_)[0, 1])
p_mantel = (np.sum(np.array(perm) <= obs) + 1) / (len(perm) + 1)
print(f"coherence: e-folding distance {L0:.0f} ± {L0_se:.0f} km; Mantel r = {obs:.2f}, p = {p_mantel:.4f}")

fig, ax = plt.subplots(figsize=(8, 4.4))
ax.scatter(d_, r_, s=34, color=BLUE, edgecolor="white", lw=1, zorder=3)
xs = np.linspace(0, d_.max() * 1.05, 100)
ax.plot(xs, f_exp(xs, L0), color=INK, lw=1.8)
for (i, j), dd, rr in zip(zip(*iu), d_, r_):
    pair = {ORDER[i], ORDER[j]}
    if pair in ({"Boone", "Asheville"}, {"Binghamton", "Scranton"}, {"Montpelier", "Asheville"}):
        ax.annotate(f"{ORDER[i]}–{ORDER[j]}", (dd, rr), xytext=(-4, -11), textcoords="offset points",
                    fontsize=7.5, color=INK2, ha="right" if dd > 1000 else "left")
ax.set_xlabel("Distance between locations (km)")
ax.set_ylabel("Correlation of yearly peak anomalies")
ax.set_ylim(0, 1)
ax.set_title("Nearby locations have early/late years together; the link fades over ~1,000 km")
ax.text(0.98, 0.95, f"fit r = exp(−d / {L0:.0f} km)\nMantel test: r = {obs:.2f}, p = {p_mantel:.4f}",
        transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=INK)
save(fig, "E13_anomaly_coherence_vs_distance",
     f"Correlation of yearly peak-date anomalies (2001–2025) for all 55 location pairs against distance. "
     f"Exponential decay fit gives an e-folding distance of {L0:.0f} km; Mantel permutation p = {p_mantel:.4f}.")

fig, ax = plt.subplots(figsize=(7.6, 6.6))
im = ax.imshow(R.values, cmap=DIV_CMAP, vmin=-1, vmax=1)
for i in range(len(ORDER)):
    for j in range(len(ORDER)):
        if i != j:
            ax.text(j, i, f"{R.values[i, j]:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if abs(R.values[i, j]) > 0.65 else INK)
ax.set_xticks(range(len(ORDER)), ORDER, rotation=45, ha="right")
ax.set_yticks(range(len(ORDER)), ORDER)
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.75, label="r of yearly anomalies")
ax.set_title("Year-to-year co-variation of peak date between locations")
save(fig, "E14_anomaly_correlation_matrix",
     "Correlation of yearly peak-date anomalies between every pair of locations (north → south). Blocks show "
     "a northern group, a central Appalachian group and the southern Blue Ridge.")

# ---------------------------------------------------------------- (b) per-location slopes
FEATS = {"tmax_pre60to31": "Late-summer heat, P−60…P−31 (°C)",
         "tmean_pre30": "Mean temp, 30 d before peak (°C)",
         "cdd_pre30": "Cold degree-days, 30 d before (°C·d)",
         "precip_pre60": "Precip, 60 d before peak (mm)",
         "gdd5_MarMay": "Spring GDD5, Mar–May (°C·d)",
         "mid_greenup_p50": "Spring green-up date (DOY)"}
A = train.copy()
for c in ["peak_doy"] + list(FEATS):
    A[c + "_a"] = A[c] - A.groupby("location")[c].transform("mean")
rows, eq = [], []
for f in FEATS:
    for loc, g in A.groupby("location"):
        m = smf.ols(f"peak_doy_a ~ {f}_a", data=g).fit()
        rows.append(dict(feature=f, location=loc, slope=m.params[f"{f}_a"], se=m.bse[f"{f}_a"]))
    full = smf.ols(f"peak_doy_a ~ {f}_a * C(location)", data=A).fit()
    red = smf.ols(f"peak_doy_a ~ {f}_a + C(location)", data=A).fit()
    Fv, pv, _ = full.compare_f_test(red)
    pooled = smf.ols(f"peak_doy_a ~ {f}_a", data=A).fit(cov_type="cluster", cov_kwds={"groups": A["year"]})
    eq.append(dict(feature=f, pooled=pooled.params[f"{f}_a"], pooled_se=pooled.bse[f"{f}_a"],
                   pooled_p=pooled.pvalues[f"{f}_a"], F_equal_slopes=Fv, p_equal_slopes=pv))
S = pd.DataFrame(rows)
EQ = pd.DataFrame(eq)
S.to_csv(TAB / "E14_per_location_slopes.csv", index=False)
EQ.to_csv(TAB / "E15_equal_slope_tests.csv", index=False)
print(EQ.round(4).to_string())

fig, axes = plt.subplots(2, 3, figsize=(13, 8), sharey=True)
yy = np.arange(len(ORDER))[::-1]
for ax, f in zip(axes.flat, FEATS):
    s = S[S.feature == f].set_index("location").loc[ORDER]
    e = EQ.set_index("feature").loc[f]
    ax.axvspan(e.pooled - 1.96 * e.pooled_se, e.pooled + 1.96 * e.pooled_se, color=ORANGE, alpha=0.12, lw=0)
    ax.axvline(e.pooled, color=ORANGE, lw=1.6)
    ax.errorbar(s.slope, yy, xerr=1.96 * s.se, fmt="o", color=BLUE, capsize=2.5, ms=6)
    ax.axvline(0, color=INK2, lw=0.9)
    ax.set_title(f"{FEATS[f]}\npooled {e.pooled:+.3g} (p={e.pooled_p:.3f});  equal slopes p={e.p_equal_slopes:.2f}",
                 fontsize=9)
    ax.grid(axis="y", visible=False)
for ax in axes[:, 0]:
    ax.set_yticks(yy, ORDER)
fig.suptitle("Days of peak shift per unit of each weather anomaly, by location (95% CI; orange = pooled slope)",
             x=0.01, ha="left", fontsize=12, fontweight="bold")
fig.text(0.01, -0.01, "Locations ordered north (top) → south. Pooled slope uses year-clustered SE. "
         "'Equal slopes' = F-test of location × feature interaction.", fontsize=8, color=INK2)
fig.tight_layout()
save(fig, "E15_per_location_slopes",
     "Sensitivity of each location's peak-date anomaly to six forecast-safe weather/phenology anomalies. Orange "
     "line and band: pooled slope and 95% CI. The F-test asks whether slopes differ between locations.")

# ---------------------------------------------------------------- (c) pooled interaction tests
MODS = {"lat": "Latitude", "elev_median_m": "Elevation", "norm_tmean_Sep": "Sep temp normal",
        "evt_northern_hardwood": "Northern hardwood share", "evt_oak_hickory": "Oak–hickory share"}
for s_ in MODS:
    A[s_ + "_c"] = A["location"].map(L[s_] - L[s_].mean()) / L[s_].std()     # standardized moderator
rows = []
for f in FEATS:
    for s_ in MODS:
        m = smf.ols(f"peak_doy_a ~ {f}_a + {f}_a:{s_}_c", data=A).fit(
            cov_type="cluster", cov_kwds={"groups": A["year"]})
        t = f"{f}_a:{s_}_c"
        rows.append(dict(feature=f, moderator=s_, main=m.params[f"{f}_a"], inter=m.params[t],
                         inter_se=m.bse[t], z=m.tvalues[t], p=m.pvalues[t]))
I = pd.DataFrame(rows)
I["p_holm"] = multipletests(I.p, method="holm")[1]
I["inter_rel"] = I.inter / I.main      # change in slope per 1 SD of moderator, relative to main slope
I.to_csv(TAB / "E16_interaction_tests.csv", index=False)
print(I.sort_values("p").head(10).round(4).to_string())

Z = I.pivot(index="feature", columns="moderator", values="z").loc[list(FEATS), list(MODS)]
Pm = I.pivot(index="feature", columns="moderator", values="p").loc[list(FEATS), list(MODS)]
Ph = I.pivot(index="feature", columns="moderator", values="p_holm").loc[list(FEATS), list(MODS)]
fig, ax = plt.subplots(figsize=(9, 5.4))
im = ax.imshow(Z.values, cmap=DIV_CMAP, vmin=-4, vmax=4)
for i in range(Z.shape[0]):
    for j in range(Z.shape[1]):
        star = "†" if Ph.values[i, j] < 0.05 else ("*" if Pm.values[i, j] < 0.05 else "")
        ax.text(j, i, f"{Z.values[i, j]:+.1f}{star}\np={Pm.values[i, j]:.3f}", ha="center", va="center",
                fontsize=7.5, color="white" if abs(Z.values[i, j]) > 2.6 else INK)
ax.set_xticks(range(len(MODS)), [MODS[m] for m in MODS], rotation=20, ha="right")
ax.set_yticks(range(len(FEATS)), [FEATS[f] for f in FEATS])
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.8, label="z of interaction (year-clustered)")
ax.set_title("Does a location's sensitivity to weather depend on its geography or forest type?")
ax.text(0, -0.24, "Positive z: the feature's effect on peak grows (more positive) where the moderator is higher. "
                  "* p < .05 unadjusted;  † survives Holm correction across all 30 tests.",
        transform=ax.transAxes, fontsize=8, color=INK2)
save(fig, "E16_interaction_heatmap",
     "Interaction tests: weather anomaly × standardized location moderator, pooled over 275 location-years with "
     "year-clustered SE. Cells show z and unadjusted p; † marks tests surviving Holm correction for 30 tests.")

# ---------------------------------------------------------------- (d) LOYO with interactions
BASE = ["tmax_pre60to31", "precip_pre60", "mid_greenup_p50"]


def loyo(features, inters=(), alpha=1.0):
    errs, errc = [], []
    for y in sorted(A.year.unique()):
        tr, te = A[A.year != y], A[A.year == y]
        mu = tr.groupby("location").peak_doy.mean()

        def design(part):
            xm = tr.groupby("location")[features].mean()
            X = (part[features].values - xm.loc[part.location].values)
            sd = (tr[features] - xm.loc[tr.location].values).std().values
            cols = [X / sd]
            for f, s_ in inters:
                k = features.index(f)
                cols.append((X[:, k] / sd[k] * part[s_ + "_c"].values)[:, None])
            return np.hstack(cols)
        ytr = tr.peak_doy.values - mu.loc[tr.location].values
        model = Ridge(alpha=alpha).fit(design(tr), ytr)
        pred = mu.loc[te.location].values + model.predict(design(te))
        errs += list(te.peak_doy.values - pred)
        errc += list(te.peak_doy.values - mu.loc[te.location].values)
    return np.mean(np.abs(errs)), np.mean(np.abs(errc))


base_mae, clim_mae = loyo(BASE)
rows = [dict(model="Climatology (μ only)", MAE=clim_mae), dict(model="Base δ: " + ", ".join(BASE), MAE=base_mae)]
for f in BASE:
    for s_ in MODS:
        mae, _ = loyo(BASE, [(f, s_)])
        rows.append(dict(model=f"+ {f} × {MODS[s_]}", MAE=mae))
LY = pd.DataFrame(rows)
LY["delta_vs_base"] = LY.MAE - base_mae
LY.to_csv(TAB / "E17_loyo_interactions.csv", index=False)
print(LY.round(3).to_string())

top = LY.iloc[2:].sort_values("delta_vs_base").head(8)
fig, ax = plt.subplots(figsize=(9, 4))
yy = np.arange(len(top))[::-1]
ax.barh(yy, top.delta_vs_base, color=[BLUE if v < 0 else "#d3d2cc" for v in top.delta_vs_base], height=0.6)
for yi, v in zip(yy, top.delta_vs_base):
    ax.text(v + (0.002 if v >= 0 else -0.002), yi, f"{v:+.3f} d", va="center",
            ha="left" if v >= 0 else "right", fontsize=8, color=INK)
ax.axvline(0, color=INK2, lw=1)
ax.set_yticks(yy, top.model)
ax.set_xlabel("Change in leave-one-year-out MAE vs. base δ model (days; negative = better)")
ax.grid(axis="y", visible=False)
ax.set_title(f"Adding one location interaction to the base δ model (base MAE {base_mae:.2f} d, "
             f"climatology {clim_mae:.2f} d)")
lim = max(0.05, top.delta_vs_base.abs().max() * 1.6)
ax.set_xlim(-lim, lim)
save(fig, "E17_loyo_interaction_gain",
     f"Best eight single interactions (weather feature × location moderator) added to a provisional base δ "
     f"model ({', '.join(BASE)}). Base LOYO MAE {base_mae:.2f} d vs climatology {clim_mae:.2f} d.")
pd.Series(dict(L0_km=L0, L0_se=L0_se, mantel_r=obs, mantel_p=p_mantel, base_mae=base_mae,
               clim_mae=clim_mae)).to_csv(TAB / "E18_step4_summary.csv")
