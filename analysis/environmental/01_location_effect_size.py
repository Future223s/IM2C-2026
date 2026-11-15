"""Step 1: How big is the location effect? Map, distributions, variance decomposition, paired tests."""
import itertools

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

from common import *  # noqa

panel, locs, train, ORDER = load()
d = train[["location", "year", "peak_doy"]].copy()

# ---------------------------------------------------------------- Fig 1: map
fig, ax = plt.subplots(figsize=(7.2, 7.6))
sc = ax.scatter(locs.lon, locs.lat, c=locs.peak_doy_mean, cmap=SEQ_CMAP, s=170,
                edgecolor="white", linewidth=2, zorder=3, vmin=270, vmax=284)
# 40 km radius circles (approximate, unclipped) for scale
for _, r in locs.iterrows():
    k = 40 / 111.0
    t = np.linspace(0, 2 * np.pi, 100)
    ax.plot(r.lon + k / np.cos(np.radians(r.lat)) * np.cos(t), r.lat + k * np.sin(t),
            color=MUTED, lw=0.8, ls="--", zorder=1)
ov = {"Shenandoah": (10, -4, "left"), "Valley Head": (-10, -4, "right"), "Boone": (10, -4, "left"),
      "Asheville": (-10, -4, "right"), "Frostburg": (-10, 0, "right"), "Roanoke": (10, -4, "left"),
      "State College": (-10, 0, "right"), "Old Forge": (-10, 0, "right"),
      "Binghamton": (-10, 0, "right"), "Scranton": (10, -4, "left"), "Montpelier": (10, -4, "left")}
for _, r in locs.iterrows():
    dx, dy, ha = ov[r.location]
    ax.annotate(f"{r.location}\n{doy_label(r.peak_doy_mean)}", (r.lon, r.lat), xytext=(dx, dy),
                textcoords="offset points", ha=ha, va="center", fontsize=8.5, color=INK)
cb = fig.colorbar(sc, ax=ax, shrink=0.55, pad=0.02)
cb.set_label("Mean peak (MODIS mid-greendown), 2001–2025")
cb.set_ticks([271, 274, 277, 280, 283])
cb.set_ticklabels([doy_label(v) for v in [271, 274, 277, 280, 283]])
ax.set_xlabel("Longitude (°)")
ax.set_ylabel("Latitude (°)")
ax.set_aspect(1 / np.cos(np.radians(40)))
ax.set_title("Typical peak date at the 11 study locations")
ax.text(0, -0.09, "Dashed circles: 40 km study radius (before clipping to the ecoregion).",
        transform=ax.transAxes, fontsize=8, color=INK2)
save(fig, "E01_map_typical_peak",
     "Map of the 11 locations, colored by mean MODIS mid-greendown date 2001–2025. Northern sites peak ~Sep 28; "
     "Shenandoah is latest (~Oct 10).")

# ---------------------------------------------------------------- Fig 2: distributions
fig, ax = plt.subplots(figsize=(10, 4.6))
rng = np.random.default_rng(1)
for i, loc in enumerate(ORDER):
    v = d[d.location == loc].peak_doy
    ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=14, color=BLUE, alpha=0.35, lw=0)
    m, se = v.mean(), v.std(ddof=1) / np.sqrt(len(v))
    ax.errorbar(i, m, yerr=1.96 * se, fmt="o", color=INK, ms=7, capsize=4, lw=1.5, zorder=4)
ax.set_xticks(range(len(ORDER)), ORDER, rotation=30, ha="right")
yt = [262, 268, 274, 280, 286, 292]
ax.set_yticks(yt, [doy_label(v) for v in yt])
ax.set_ylabel("Peak date (MODIS mid-greendown)")
ax.set_title("Each location's 25 peak dates (dots) and mean ± 95% CI, north → south")
ax.grid(axis="x", visible=False)
save(fig, "E02_peak_by_location",
     "Peak date in every year 2001–2025 (dots) with the location mean and 95% CI (black). Locations are ordered "
     "north to south. Year-to-year spread (~±4 d) is similar to the spread between locations.")

# ---------------------------------------------------------------- Variance decomposition
m2 = smf.ols("peak_doy ~ C(location) + C(year)", data=d).fit()
aov = sm.stats.anova_lm(m2, typ=2)
ss = aov["sum_sq"]
share = (ss / ss.sum()).rename({"C(location)": "Location", "C(year)": "Year",
                                "Residual": "Residual"})
aov.to_csv(TAB / "E01_anova_two_way.csv")
print(aov)

fig, ax = plt.subplots(figsize=(8.5, 1.9))
left = 0
cols = [BLUE, ORANGE, "#d3d2cc"]
for (name, s), c in zip(share.items(), cols):
    ax.barh(0, s, left=left, color=c, height=0.5, edgecolor=SURFACE, linewidth=2)
    ax.text(left + s / 2, 0, f"{name}\n{s:.0%}", ha="center", va="center", fontsize=9,
            color="white" if c != "#d3d2cc" else INK)
    left += s
ax.set_xlim(0, 1)
ax.set_yticks([])
ax.set_xticks([0, .25, .5, .75, 1], ["0%", "25%", "50%", "75%", "100%"])
ax.grid(False)
F_loc, p_loc = aov.loc["C(location)", "F"], aov.loc["C(location)", "PR(>F)"]
F_yr, p_yr = aov.loc["C(year)", "F"], aov.loc["C(year)", "PR(>F)"]
ax.set_title("Where does the variation in peak date come from? (two-way ANOVA, 275 location-years)")
ax.text(0, -0.62, "Location = fixed differences between sites.  Year = weather shared by all sites that year.  "
        "Residual = local weather + noise.\n"
        f"Location: F = {F_loc:.1f}, p = {p_loc:.1e}.   Year: F = {F_yr:.1f}, p = {p_yr:.1e}.",
        transform=ax.transAxes, fontsize=8.5, color=INK2)
save(fig, "E03_variance_decomposition",
     f"Share of the total sum of squares in peak date explained by location (fixed site differences), year "
     f"(weather shared across the Appalachians) and the remainder. Location F = {F_loc:.1f}, year F = {F_yr:.1f}.")

# ---------------------------------------------------------------- Pairwise paired-by-year tests
w = d.pivot(index="year", columns="location", values="peak_doy")[ORDER]
rows = []
for a, b in itertools.combinations(ORDER, 2):
    diff = w[a] - w[b]
    t, p = stats.ttest_1samp(diff, 0)
    rows.append(dict(a=a, b=b, mean_diff=diff.mean(), sd_diff=diff.std(ddof=1), t=t, p=p))
pw = pd.DataFrame(rows)
pw["p_holm"] = multipletests(pw.p, method="holm")[1]
pw.to_csv(TAB / "E02_pairwise_paired_tests.csv", index=False)
n_sig = (pw.p_holm < 0.05).sum()
print(f"pairwise: {n_sig}/{len(pw)} significant after Holm")

M = pd.DataFrame(np.nan, index=ORDER, columns=ORDER)
S = pd.DataFrame("", index=ORDER, columns=ORDER)
for _, r in pw.iterrows():
    M.loc[r.a, r.b], M.loc[r.b, r.a] = r.mean_diff, -r.mean_diff
    star = "***" if r.p_holm < 0.001 else "**" if r.p_holm < 0.01 else "*" if r.p_holm < 0.05 else ""
    S.loc[r.a, r.b] = S.loc[r.b, r.a] = star
fig, ax = plt.subplots(figsize=(8.4, 7.2))
lim = np.nanmax(np.abs(M.values))
im = ax.imshow(M.values, cmap=DIV_CMAP, norm=TwoSlopeNorm(0, -lim, lim))
for i in range(len(ORDER)):
    for j in range(len(ORDER)):
        if i != j:
            v = M.values[i, j]
            ax.text(j, i, f"{v:+.1f}\n{S.values[i, j]}", ha="center", va="center", fontsize=7,
                    color="white" if abs(v) > 0.6 * lim else INK)
ax.set_xticks(range(len(ORDER)), ORDER, rotation=45, ha="right")
ax.set_yticks(range(len(ORDER)), ORDER)
ax.grid(False)
cb = fig.colorbar(im, ax=ax, shrink=0.7)
cb.set_label("Row minus column, mean days (red = row later)")
ax.set_title("How many days later is each row location than each column location?")
ax.text(0, -0.25, f"Paired by year (cancels shared weather); stars = Holm-adjusted p < .05/.01/.001. "
                  f"{n_sig} of {len(pw)} pairs differ significantly.", transform=ax.transAxes, fontsize=8,
        color=INK2)
save(fig, "E04_pairwise_differences",
     "Mean difference in peak date between every pair of locations, paired by year. Pairing removes the weather "
     "all sites share, so differences of ~1.5 d are detectable. Stars: Holm-adjusted significance.")

summary = dict(F_location=F_loc, p_location=p_loc, F_year=F_yr, p_year=p_yr,
               share_location=share.iloc[0], share_year=share.iloc[1], share_resid=share.iloc[2],
               pairs_sig=int(n_sig), pairs_total=len(pw))
pd.Series(summary).to_csv(TAB / "E00_step1_summary.csv")
print(pd.Series(summary).round(4))
