"""Step 2: Which location descriptors explain the location mean peak date (mu)? n = 11 locations.

For each candidate x: OLS mu ~ x (t-test, df = 9), permutation p (10,000 shuffles), location-bootstrap 95% CI of the
slope, and leave-one-location-out (LOLO) MAE vs the grand-mean benchmark.
Also builds 'thermal dates' from the daily weather climatology.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from common import *  # noqa

panel, locs, train, ORDER = load()
mu = locs.set_index("location")["peak_doy_mean"]

# ---------------------------------------------------------------- thermal dates from daily climatology
wx = pd.read_csv(PANEL / "weather_daily.csv", parse_dates=["date"])
wx = wx[wx.date.dt.year.between(2001, 2025)]
wx["doy"] = wx.date.dt.day_of_year
clim = wx.groupby(["location", "doy"])[["tmean", "tmin"]].mean()
TH = list(range(8, 19))
therm = {}
for loc in ORDER:
    c = clim.loc[loc]["tmean"].rolling(15, center=True, min_periods=8).mean()
    c = c[c.index >= 182]          # autumn: from Jul 1
    for T in TH:
        below = c[c < T]
        therm.setdefault(T, {})[loc] = below.index.min() if len(below) else np.nan
therm = pd.DataFrame(therm)          # index location, columns threshold
therm.to_csv(TAB / "E03_thermal_dates.csv")

# ---------------------------------------------------------------- candidate table
X = locs.set_index("location").copy()
for T in TH:
    X[f"thermal_date_{T}C"] = therm[T]
CANDIDATES = {
    "lat": "Latitude (°N)", "lon": "Longitude (°E)", "daylen_Oct1_h": "Day length Oct 1 (h)",
    "elev_median_m": "Forest elevation, median (m)", "elev_p90_m": "Forest elevation, 90th pct (m)",
    "slope_mean_deg": "Mean slope (°)", "northness_mean": "Mean northness",
    "region_area_km2": "Region area (km²)",
    "norm_tmean_Sep": "Sep mean temp, normal (°C)", "norm_tmin_Oct": "Oct min temp, normal (°C)",
    "norm_precip_AugSep": "Aug–Sep precip, normal (mm)", "norm_gdd5_MarMay": "Spring GDD5, normal (°C·d)",
    "frac_hardwood": "Hardwood share of area", "frac_evergreen": "Evergreen share of area",
    "frac_developed": "Developed share of area",
    "evt_northern_hardwood": "Northern hardwood share (EVT)", "evt_oak_hickory": "Oak–hickory share (EVT)",
    "evt_cove_mesic": "Cove/mesic share (EVT)",
    "thermal_date_12C": "Date climatological Tmean < 12 °C",
    "thermal_date_15C": "Date climatological Tmean < 15 °C",
}
X[list(CANDIDATES)].assign(peak_doy_mean=mu).to_csv(TAB / "E04_location_predictors.csv")

rng = np.random.default_rng(42)
y = mu.loc[ORDER].values
n = len(y)
lolo_grand = np.mean([abs(y[i] - np.delete(y, i).mean()) for i in range(n)])


def lolo_mae(x, y):
    err = []
    for i in range(len(y)):
        m = np.ones(len(y), bool); m[i] = False
        b = np.polyfit(x[m], y[m], 1)
        err.append(y[i] - np.polyval(b, x[i]))
    return np.mean(np.abs(err)), np.array(err)


def univariate(x, y, nperm=10000, nboot=5000):
    fit = sm.OLS(y, sm.add_constant(x)).fit()
    r = np.corrcoef(x, y)[0, 1]
    perm = np.array([np.corrcoef(x, rng.permutation(y))[0, 1] for _ in range(nperm)])
    p_perm = (np.sum(np.abs(perm) >= abs(r)) + 1) / (nperm + 1)
    boots = []
    for _ in range(nboot):
        idx = rng.integers(0, len(y), len(y))
        if np.ptp(x[idx]) > 0:
            boots.append(np.polyfit(x[idx], y[idx], 1)[0])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    mae, _ = lolo_mae(x, y)
    return dict(slope=fit.params[1], se=fit.bse[1], p_t=fit.pvalues[1], p_perm=p_perm, r=r, r2=r ** 2,
                boot_lo=lo, boot_hi=hi, lolo_mae=mae)


rows = []
for c, label in CANDIDATES.items():
    x = X.loc[ORDER, c].values.astype(float)
    rows.append(dict(predictor=c, label=label, **univariate(x, y)))
U = pd.DataFrame(rows)
U["lolo_gain_d"] = lolo_grand - U["lolo_mae"]
U = U.sort_values("lolo_mae").reset_index(drop=True)
U.to_csv(TAB / "E05_univariate_mu_models.csv", index=False)
print(f"grand-mean LOLO MAE = {lolo_grand:.2f} d")
print(U[["predictor", "slope", "se", "r", "r2", "p_t", "p_perm", "lolo_mae"]].round(3).to_string())

# ---------------------------------------------------------------- Fig E05: univariate ranking
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 6.6), sharey=True, gridspec_kw=dict(width_ratios=[1, 1]))
yy = np.arange(len(U))[::-1]
a1.barh(yy, U["r"], color=[RED if v > 0 else BLUE for v in U["r"]], height=0.62)
for yi, (r, p) in zip(yy, zip(U["r"], U["p_perm"])):
    star = "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""
    a1.text(r + (0.03 if r >= 0 else -0.03), yi, f"{r:+.2f}{star}", va="center",
            ha="left" if r >= 0 else "right", fontsize=8, color=INK)
a1.axvline(0, color=INK2, lw=1)
a1.set_xlim(-1.15, 1.15)
a1.set_yticks(yy, U["label"])
a1.set_xlabel("Correlation with mean peak date (n = 11)")
a1.set_title("Correlation (red = later peak where higher)")
a1.grid(axis="y", visible=False)
a2.barh(yy, U["lolo_mae"], color=[BLUE if v < lolo_grand else "#d3d2cc" for v in U["lolo_mae"]], height=0.62)
a2.axvline(lolo_grand, color=ORANGE, lw=2, ls="--")
a2.text(lolo_grand + 0.05, len(U) - 0.4, f"No predictor\n(grand mean) {lolo_grand:.2f} d", color=INK2,
        fontsize=8, va="top")
for yi, v in zip(yy, U["lolo_mae"]):
    a2.text(v + 0.04, yi, f"{v:.2f}", va="center", fontsize=8, color=INK)
a2.set_xlabel("Leave-one-location-out MAE of mean peak date (days)")
a2.set_title("Out-of-sample error (lower = better)")
a2.grid(axis="y", visible=False)
fig.text(0.01, -0.01, "Stars: permutation p < .05/.01/.001 (10,000 shuffles). Blue bars beat the grand mean. "
         "Sorted by out-of-sample error.", fontsize=8, color=INK2)
fig.suptitle("Which single location descriptor best explains a site's typical peak date?", x=0.01, ha="left",
             fontsize=12, fontweight="bold")
fig.tight_layout()
save(fig, "E05_univariate_predictors",
     "Every candidate location descriptor, ranked by leave-one-location-out error when used alone to predict the "
     "site's mean peak date. Left: correlation with permutation significance. Right: out-of-sample MAE vs. the "
     "no-predictor benchmark (orange).")

# ---------------------------------------------------------------- Fig E06: scatter small multiples
show = ["lat", "elev_median_m", "lon", "norm_tmean_Sep", "evt_northern_hardwood", "thermal_date_12C"]
fig, axes = plt.subplots(2, 3, figsize=(12, 7.4), sharey=True)
yt = [268, 272, 276, 280, 284]
for ax, c in zip(axes.flat, show):
    x = X.loc[ORDER, c].values.astype(float)
    b = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    # 95% confidence band of the fit
    fit = sm.OLS(y, sm.add_constant(x)).fit()
    pr = fit.get_prediction(sm.add_constant(xs)).summary_frame(alpha=0.05)
    ax.fill_between(xs, pr["mean_ci_lower"], pr["mean_ci_upper"], color=BLUE, alpha=0.12, lw=0)
    ax.plot(xs, pr["mean"], color=BLUE, lw=2)
    ax.scatter(x, y, s=60, color=INK, edgecolor="white", lw=1.5, zorder=3)
    offs = {"Shenandoah": (6, -3, "left"), "Montpelier": (5, 3, "left"), "Asheville": (5, 3, "left"),
            "Valley Head": (-5, -11, "right"), "Frostburg": (5, -11, "left")}
    for xi, yi, nme in zip(x, y, ORDER):
        if nme in offs:
            dx_, dy_, ha_ = offs[nme]
            ax.annotate(nme, (xi, yi), xytext=(dx_, dy_), textcoords="offset points", fontsize=7.5,
                        color=INK2, ha=ha_)
    ax.set_ylim(265, 291)
    row = U.set_index("predictor").loc[c]
    ax.set_title(CANDIDATES[c], fontsize=10)
    ax.text(0.03, 0.97, f"r = {row.r:+.2f}  (perm p = {row.p_perm:.3f})\n"
                        f"slope {row.slope:+.3g} d per unit,  LOLO MAE {row.lolo_mae:.2f} d",
            transform=ax.transAxes, va="top", fontsize=8, color=INK)
    if c == "thermal_date_12C":
        xt = [int(v) for v in np.linspace(x.min(), x.max(), 4)]
        ax.set_xticks(xt, [doy_label(v) for v in xt])
for ax in axes[:, 0]:
    ax.set_yticks(yt, [doy_label(v) for v in yt])
    ax.set_ylabel("Mean peak date")
fig.suptitle("Location mean peak date vs. six location descriptors (line = OLS fit, band = 95% CI)",
             x=0.01, ha="left", fontsize=12, fontweight="bold")
fig.tight_layout()
save(fig, "E06_mu_vs_descriptors",
     "Mean peak date of each location against latitude, forest elevation, longitude, September normal temperature, "
     "northern-hardwood share and the climatological date mean temperature falls below 12 °C. n = 11; Shenandoah "
     "sits above every fit.")

# ---------------------------------------------------------------- Fig E07: thermal threshold scan
rows = []
for T in TH:
    x = therm.loc[ORDER, T].values.astype(float)
    r = np.corrcoef(x, y)[0, 1]
    rows.append(dict(T=T, r=r, lolo=lolo_mae(x, y)[0], slope=np.polyfit(x, y, 1)[0]))
TS = pd.DataFrame(rows)
TS.to_csv(TAB / "E06_thermal_threshold_scan.csv", index=False)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.8))
a1.plot(TS["T"], TS["r"], marker="o", color=BLUE)
a1.set_xlabel("Threshold T (°C)")
a1.set_ylabel("r with mean peak date")
a1.set_ylim(0, 1)
a1.set_title("Correlation vs. threshold (flat = not cherry-picked)")
a2.plot(TS["T"], TS["lolo"], marker="o", color=BLUE)
a2.axhline(lolo_grand, color=ORANGE, ls="--", lw=2)
a2.text(TH[0], lolo_grand + 0.05, "grand mean", color=INK2, fontsize=8)
a2.set_xlabel("Threshold T (°C)")
a2.set_ylabel("LOLO MAE (days)")
a2.set_title("Out-of-sample error vs. threshold")
fig.suptitle("Does peak follow the date a site's average temperature drops below T °C?",
             x=0.01, ha="left", fontsize=12, fontweight="bold")
fig.tight_layout()
save(fig, "E07_thermal_threshold_scan",
     "For thresholds 8–18 °C, the date each location's 2001–2025 climatological (15-day smoothed) mean "
     "temperature first falls below T, related to the location's mean peak date. Left: correlation. Right: "
     "leave-one-location-out error.")

# ---------------------------------------------------------------- Fig E08: collinearity among predictors
core = ["lat", "lon", "elev_median_m", "slope_mean_deg", "norm_tmean_Sep", "norm_tmin_Oct",
        "norm_precip_AugSep", "norm_gdd5_MarMay", "frac_hardwood", "frac_developed",
        "evt_northern_hardwood", "evt_oak_hickory", "evt_cove_mesic", "thermal_date_12C"]
C = X.loc[ORDER, core].assign(peak_doy_mean=y).corr()
lab = [CANDIDATES.get(c, "Mean peak date") for c in C.columns]
fig, ax = plt.subplots(figsize=(9.5, 8.2))
im = ax.imshow(C.values, cmap=DIV_CMAP, vmin=-1, vmax=1)
for i in range(len(C)):
    for j in range(len(C)):
        v = C.values[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6.8,
                color="white" if abs(v) > 0.65 else INK)
ax.set_xticks(range(len(C)), lab, rotation=55, ha="right", fontsize=8)
ax.set_yticks(range(len(C)), lab, fontsize=8)
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.7, label="Pearson r across the 11 locations")
ax.set_title("Location descriptors are strongly inter-correlated (n = 11)")
save(fig, "E08_predictor_correlations",
     "Correlation between location descriptors across the 11 sites. Latitude, temperature normals, species mix and "
     "thermal date form one block, so with n = 11 they cannot be separated statistically.")

pd.Series(dict(lolo_grand=lolo_grand)).to_csv(TAB / "E07_benchmarks.csv")
