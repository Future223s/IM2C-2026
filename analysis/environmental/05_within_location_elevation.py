"""Step 5: Elevation (and aspect) WITHIN each location, from pixels.csv (4,000 forest pixels/location, 2001-2025).

Pixel anomaly = pixel peak_doy - median of its location-year (removes that year's weather and the location mean).
  (a) binned anomaly vs elevation per location
  (b) per-location linear slope, 95% CI from a block bootstrap over years (1,000 reps)
  (c) two-piece (hinge) breakpoint search per region group, bootstrap over years for the break
  (d) multiple regression: elevation + northness + forest fraction (location-year fixed effects)
Bootstraps use per-year sufficient statistics, so every replicate is exact OLS on the resampled years.
"""
import zipfile

import numpy as np
import pandas as pd

from common import *  # noqa

panel, locs, train, ORDER = load()
px_path = PANEL / "pixels.csv"
if not px_path.exists():                       # repo copy ships it zipped
    with zipfile.ZipFile(PANEL.parent / "panel_v1_pixels.zip") as z:
        z.extract("panel_v1/pixels.csv", PANEL.parent)
px = pd.read_csv(px_path)
px["anom"] = px.peak_doy - px.groupby(["location", "year"]).peak_doy.transform("median")
px = px[px.anom.abs() < 60]          # drop a handful of gross outliers (failed retrievals)
GROUP = {**{l: "North (VT, NY, PA)" for l in ORDER[:5]},
         **{l: "Central (MD, WV, VA)" for l in ORDER[5:9]},
         **{l: "South (NC Blue Ridge)" for l in ORDER[9:]}}
px["group"] = px.location.map(GROUP)
print(f"{len(px):,} pixel-years, {px.pixel_id.nunique():,} pixels")
rng = np.random.default_rng(7)
L_lat = locs.set_index("location").lat
NB = 1000

summ = px.drop_duplicates("pixel_id").groupby("location").elev_m.describe(percentiles=[.05, .5, .95]).loc[ORDER]
summ.to_csv(TAB / "E19_pixel_elevation_summary.csv")


def demean(df, cols, by=("location", "year")):
    return df[cols] - df.groupby(list(by))[cols].transform("mean")


# ---------------------------------------------------------------- (a) binned medians per location
fig, axes = plt.subplots(3, 4, figsize=(13, 8.6), sharey=True)
for ax, loc in zip(axes.flat, ORDER):
    d = px[px.location == loc]
    bins = np.arange(0, 2100, 100)
    d = d.assign(b=pd.cut(d.elev_m, bins, labels=bins[:-1] + 50))
    g = d.groupby("b", observed=True).anom.agg(["median", "count", lambda s: s.quantile(.25),
                                                 lambda s: s.quantile(.75)])
    g.columns = ["med", "n", "q1", "q3"]
    g = g[g.n >= 50]
    x = g.index.astype(float)
    ax.fill_between(x, g.q1, g.q3, color=BLUE, alpha=0.15, lw=0)
    ax.plot(x, g.med, color=BLUE, marker="o", ms=4)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.axvline(1350, color=ORANGE, lw=1, ls="--")
    ax.set_title(loc, fontsize=10)
    ax.set_xlim(0, 1900)
    ax.set_ylim(-14, 10)
for ax in axes.flat[len(ORDER):]:
    ax.axis("off")
axes.flat[len(ORDER)].text(0.05, 0.6, "Blue: median pixel anomaly per 100 m\nelevation bin (band = IQR),\nbins with "
                           "≥ 50 pixel-years.\nOrange dashed: 1,350 m break from\nthe earlier 4-region analysis.",
                           fontsize=8.5, color=INK2, transform=axes.flat[len(ORDER)].transAxes, va="top")
for ax in axes[:, 0]:
    ax.set_ylabel("Peak anomaly vs.\nlocation-year median (d)")
for ax in axes[-1, :3]:
    ax.set_xlabel("Pixel elevation (m)")
fig.suptitle("Within each location, how does peak date change with elevation?", x=0.01, ha="left",
             fontsize=12, fontweight="bold")
fig.tight_layout()
save(fig, "E18_within_location_binned",
     "Median peak-date anomaly of forest pixels by 100 m elevation bin within each location (anomaly relative to "
     "the location-year median, so weather and the location mean are removed).")


# ---------------------------------------------------------------- (b) per-location slopes, year block bootstrap
def year_stats(d, xcols):
    """Per-year sufficient statistics after location-year demeaning: X'X, X'y for each year."""
    X = demean(d, xcols).values
    yv = demean(d, ["anom"]).values[:, 0]
    out = {}
    for yr, idx in d.groupby("year").indices.items():
        Xi, yi = X[idx], yv[idx]
        out[yr] = (Xi.T @ Xi, Xi.T @ yi)
    return out


def boot_coef(stats, nb=NB):
    yrs = np.array(sorted(stats))
    XtX = np.array([stats[y][0] for y in yrs]); Xty = np.array([stats[y][1] for y in yrs])
    est = np.linalg.solve(XtX.sum(0), Xty.sum(0))
    reps = []
    for _ in range(nb):
        k = rng.integers(0, len(yrs), len(yrs))
        reps.append(np.linalg.solve(XtX[k].sum(0), Xty[k].sum(0)))
    reps = np.array(reps)
    return est, np.percentile(reps, 2.5, axis=0), np.percentile(reps, 97.5, axis=0)


rows = []
for loc in ORDER:
    d = px[px.location == loc].copy()
    d["e100"] = d.elev_m / 100
    est, lo, hi = boot_coef(year_stats(d, ["e100"]))
    est3, lo3, hi3 = boot_coef(year_stats(d, ["e100", "northness", "forest_frac"]))
    rows.append(dict(location=loc, slope_d_per_100m=est[0], lo=lo[0], hi=hi[0],
                     slope_ctrl=est3[0], lo_ctrl=lo3[0], hi_ctrl=hi3[0],
                     northness=est3[1], northness_lo=lo3[1], northness_hi=hi3[1],
                     forest_frac=est3[2], elev_p5=summ.loc[loc, "5%"], elev_p95=summ.loc[loc, "95%"],
                     elev_median=summ.loc[loc, "50%"]))
SL = pd.DataFrame(rows)
SL.to_csv(TAB / "E20_within_location_slopes.csv", index=False)
print(SL.round(3).to_string())

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw=dict(width_ratios=[1.25, 1]))
yy = np.arange(len(ORDER))[::-1]
a1.errorbar(SL.slope_d_per_100m, yy + 0.15, xerr=[SL.slope_d_per_100m - SL.lo, SL.hi - SL.slope_d_per_100m],
            fmt="o", color=BLUE, capsize=3, ms=6, label="Elevation only")
a1.errorbar(SL.slope_ctrl, yy - 0.15, xerr=[SL.slope_ctrl - SL.lo_ctrl, SL.hi_ctrl - SL.slope_ctrl],
            fmt="s", color=ORANGE, capsize=3, ms=5, label="+ northness + forest fraction")
a1.axvline(0, color=INK2, lw=1)
a1.set_yticks(yy, [f"{l}  ({int(p5)}–{int(p95)} m)" for l, p5, p95 in zip(SL.location, SL.elev_p5, SL.elev_p95)])
a1.set_xlabel("Days of peak shift per +100 m within the location")
a1.legend(loc="lower left", fontsize=8)
a1.set_title("Linear within-location elevation slope (95% CI, year bootstrap)", fontsize=10)
a1.grid(axis="y", visible=False)
a2.errorbar(SL.northness, yy, xerr=[SL.northness - SL.northness_lo, SL.northness_hi - SL.northness], fmt="o",
            color=BLUE, capsize=3, ms=6)
a2.axvline(0, color=INK2, lw=1)
a2.set_yticks(yy, SL.location)
a2.set_xlabel("Days per unit northness (+1 = north-facing)")
a2.set_title("Aspect effect (same model)", fontsize=10)
a2.grid(axis="y", visible=False)
fig.suptitle("Elevation and aspect effects inside each location (5th–95th pct elevation in brackets)", x=0.01,
             ha="left", fontsize=12, fontweight="bold")
fig.tight_layout()
save(fig, "E19_within_location_slopes",
     "Per-location linear slope of pixel peak anomaly on elevation (d per 100 m), alone and controlling for "
     "northness and forest fraction, with 95% CIs from a block bootstrap over years. Right: the aspect effect.")


# ---------------------------------------------------------------- (c) two-piece breakpoint per group
HGRID = np.arange(400, 1651, 50)


def hinge_stats(d, h):
    d = d.assign(e100=d.elev_m / 100, hinge=np.maximum(0, d.elev_m - h) / 100)
    return year_stats(d, ["e100", "hinge"])


def sse_from_stats(stats, yty, yrs):
    XtX = sum(stats[y][0] for y in yrs); Xty = sum(stats[y][1] for y in yrs)
    b = np.linalg.solve(XtX, Xty)
    return yty - b @ Xty, b


rows, prof = [], []
for grp in ["North (VT, NY, PA)", "Central (MD, WV, VA)", "South (NC Blue Ridge)"]:
    d = px[px.group == grp].copy()
    yv = demean(d, ["anom"]).values[:, 0]
    yty_year = pd.Series(yv ** 2).groupby(d.year.values).sum()
    yrs = sorted(d.year.unique())
    hmax = np.percentile(d.elev_m, 99)
    grid = [h for h in HGRID if (d.elev_m > h).sum() >= 1000]   # >= 1,000 pixel-years above the break
    S = {h: hinge_stats(d, h) for h in grid}
    lin = year_stats(d.assign(e100=d.elev_m / 100), ["e100"])
    yty = yty_year.sum()
    sse_lin = yty - np.linalg.solve(sum(lin[y][0] for y in yrs), sum(lin[y][1] for y in yrs)) @ \
        sum(lin[y][1] for y in yrs)
    sse = {h: sse_from_stats(S[h], yty, yrs)[0] for h in grid}
    for h in grid:
        prof.append(dict(group=grp, h=h, sse_rel=sse[h] / sse_lin))
    hbest = min(sse, key=sse.get)
    _, b = sse_from_stats(S[hbest], yty, yrs)
    # bootstrap years: re-select break and refit
    bh, below, above = [], [], []
    yarr = np.array(yrs)
    for _ in range(500):
        k = yarr[rng.integers(0, len(yarr), len(yarr))]
        yty_b = yty_year.loc[k].sum()
        s_b = {h: sse_from_stats(S[h], yty_b, k) for h in grid}
        hb = min(s_b, key=lambda h: s_b[h][0])
        bh.append(hb); below.append(s_b[hb][1][0]); above.append(s_b[hb][1][0] + s_b[hb][1][1])
    n_above = int((d.elev_m > hbest).sum())
    rows.append(dict(group=grp, break_m=hbest, break_lo=np.percentile(bh, 2.5), break_hi=np.percentile(bh, 97.5),
                     slope_below=b[0], below_lo=np.percentile(below, 2.5), below_hi=np.percentile(below, 97.5),
                     slope_above=b[0] + b[1], above_lo=np.percentile(above, 2.5), above_hi=np.percentile(above, 97.5),
                     pixel_years_above=n_above, sse_reduction_pct=100 * (1 - sse[hbest] / sse_lin),
                     elev_p99=hmax))
BR = pd.DataFrame(rows)
PR = pd.DataFrame(prof)
BR.to_csv(TAB / "E21_breakpoints_by_group.csv", index=False)
print(BR.round(2).to_string())

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4))
cols = {"North (VT, NY, PA)": BLUE, "Central (MD, WV, VA)": ORANGE, "South (NC Blue Ridge)": AQUA}
for grp, g in PR.groupby("group"):
    a1.plot(g.h, (1 - g.sse_rel) * 100, marker="o", ms=4, color=cols[grp], label=grp)
a1.axvline(1350, color=MUTED, ls="--", lw=1)
a1.text(1360, a1.get_ylim()[1] * 0.95 if a1.get_ylim()[1] > 0 else 0, "1,350 m (old)", fontsize=8,
        color=INK2, va="top")
a1.set_xlabel("Candidate break elevation (m)")
a1.set_ylabel("% reduction in residual SS vs. straight line")
a1.set_title("Where is the break? (higher = better fit)", fontsize=10)
a1.legend(fontsize=8, loc="upper left")
yy = np.arange(len(BR))[::-1]
for k, (yi, r) in enumerate(zip(yy, BR.itertuples())):
    a2.errorbar(r.slope_below, yi + 0.12, xerr=[[r.slope_below - r.below_lo], [r.below_hi - r.slope_below]],
                fmt="o", color=BLUE, capsize=3, label="below break" if k == 0 else None)
    a2.errorbar(r.slope_above, yi - 0.12, xerr=[[r.slope_above - r.above_lo], [r.above_hi - r.slope_above]],
                fmt="s", color=ORANGE, capsize=3, label="above break" if k == 0 else None)
    a2.annotate(f"break {r.break_m:.0f} m ({r.break_lo:.0f}–{r.break_hi:.0f})", (-3.4, yi), xytext=(0, 13),
                textcoords="offset points", fontsize=7.8, color=INK2, ha="left")
a2.axvline(0, color=INK2, lw=1)
a2.set_yticks(yy, BR.group)
a2.set_ylim(-0.6, len(BR) - 0.3)
a2.set_xlabel("Days per +100 m")
a2.set_title("Two-piece slopes with 95% CI (year bootstrap, break re-chosen)", fontsize=10)
a2.legend(fontsize=8, loc="lower left")
a2.grid(axis="y", visible=False)
fig.suptitle("Is there an elevation where peak timing starts dropping fast? Region groups compared", x=0.01,
             ha="left", fontsize=12, fontweight="bold")
fig.tight_layout()
save(fig, "E20_breakpoint_by_region",
     "Left: improvement of a two-piece (hinge) elevation model over a straight line for each candidate break, by "
     "region group. Right: slope below and above the best break with bootstrap CIs (break re-selected in every "
     "replicate).")

# ---------------------------------------------------------------- (c2) per-location breaks
rows = []
for loc in ORDER:
    d = px[px.location == loc].copy()
    grid = [h for h in HGRID if (d.elev_m > h).sum() >= 1000 and (d.elev_m < h).sum() >= 1000]
    if len(grid) < 3:
        continue
    yv = demean(d, ["anom"]).values[:, 0]
    yty_year = pd.Series(yv ** 2).groupby(d.year.values).sum()
    yrs = np.array(sorted(d.year.unique()))
    S = {h: hinge_stats(d, h) for h in grid}
    lin = year_stats(d.assign(e100=d.elev_m / 100), ["e100"])
    XtX = sum(lin[y][0] for y in yrs); Xty = sum(lin[y][1] for y in yrs)
    sse_lin = yty_year.sum() - np.linalg.solve(XtX, Xty) @ Xty
    fits = {h: sse_from_stats(S[h], yty_year.sum(), yrs) for h in grid}
    hb = min(fits, key=lambda h: fits[h][0])
    b = fits[hb][1]
    bh, ab = [], []
    for _ in range(300):
        k = yrs[rng.integers(0, len(yrs), len(yrs))]
        fb = {h: sse_from_stats(S[h], yty_year.loc[k].sum(), k) for h in grid}
        h2 = min(fb, key=lambda h: fb[h][0]); bh.append(h2); ab.append(fb[h2][1].sum())
    rows.append(dict(location=loc, lat=L_lat[loc], break_m=hb, break_lo=np.percentile(bh, 2.5),
                     break_hi=np.percentile(bh, 97.5), slope_below=b[0], slope_above=b.sum(),
                     above_lo=np.percentile(ab, 2.5), above_hi=np.percentile(ab, 97.5),
                     sse_gain_pct=100 * (1 - fits[hb][0] / sse_lin), px_years_above=int((d.elev_m > hb).sum()),
                     norm_tmean_Sep=locs.set_index("location").loc[loc, "norm_tmean_Sep"]))
LB = pd.DataFrame(rows)
LB.to_csv(TAB / "E23_breakpoints_by_location.csv", index=False)
print(LB.round(2).to_string())
strong = LB[(LB.sse_gain_pct >= 1.0) & (LB.above_hi < 0)]
r_lat = np.corrcoef(strong.lat, strong.break_m)[0, 1]
fig, ax = plt.subplots(figsize=(10, 5))
LBo = LB.set_index("location").reindex([l for l in ORDER if l in set(LB.location)])
yy = np.arange(len(LBo))[::-1]
for yi, (loc, r) in zip(yy, LBo.iterrows()):
    ok = loc in set(strong.location)
    edge = r.break_m in (HGRID.min(), )
    ax.errorbar(r.break_m, yi, xerr=[[r.break_m - r.break_lo], [r.break_hi - r.break_m]], fmt="o",
                color=BLUE if ok else MUTED, mfc=BLUE if ok else SURFACE, capsize=3, ms=8, mew=1.5)
    note = (f"{r.slope_above:+.1f} d/100 m above  (fit gain {r.sse_gain_pct:.1f}%)"
            + ("  — break at search edge" if edge else ""))
    ax.text(1720, yi, note, va="center", fontsize=8, color=INK if ok else INK2)
ax.set_yticks(yy, LBo.index)
ax.set_xlim(300, 2700)
ax.set_xticks(range(400, 1700, 200))
ax.set_xlabel("Best break elevation (m), 95% year-bootstrap range")
ax.grid(axis="y", visible=False)
ax.set_title("Per-location breakpoint of the two-piece elevation model (north → south)")
ax.text(0, -0.17, f"Filled = clear break (fit gain ≥ 1% and slope above significantly < 0). Clear breaks vs latitude: "
                  f"r = {r_lat:+.2f}, n = {len(strong)} — no reliable trend.", transform=ax.transAxes, fontsize=8,
        color=INK2)
save(fig, "E22_breakpoint_by_location",
     "Best break elevation per location with bootstrap range and the slope above it. Clear steep upper segments "
     "exist at Montpelier (>600 m) and Asheville (>1,200 m); elsewhere breaks are weak or sit at the search edge.")

# ---------------------------------------------------------------- (e) scales compared
geo = pd.read_csv(TAB / "E13_geography_model_summary.csv", index_col=0).iloc[:, 0]
cmp_ = pd.DataFrame([
    dict(scale="Between locations (μ ~ lat + elev, n = 11)", slope=geo["b_elev_per100"]),
    dict(scale="Within locations, pooled linear (all pixels)",
         slope=np.average(SL.slope_d_per_100m, weights=1 / ((SL.hi - SL.lo) / 3.92) ** 2)),
    *[dict(scale=f"Within, {r.group}: below {r.break_m:.0f} m", slope=r.slope_below) for r in BR.itertuples()],
    *[dict(scale=f"Within, {r.group}: above {r.break_m:.0f} m", slope=r.slope_above) for r in BR.itertuples()],
    dict(scale="Free-air lapse rate (6.5 °C/km) × Sep temp sensitivity", slope=None),
])
lapse = -0.65 * pd.read_csv(TAB / "E09_mu_model_coefficients.csv").query("term == 'norm_tmean_Sep'").coef.iloc[0]
cmp_.loc[cmp_.slope.isna(), "slope"] = lapse
cmp_.to_csv(TAB / "E22_elevation_scales.csv", index=False)
fig, ax = plt.subplots(figsize=(9.5, 4.2))
yy = np.arange(len(cmp_))[::-1]
ax.barh(yy, cmp_.slope, color=[ORANGE if "Between" in s else MUTED if "lapse" in s else BLUE for s in cmp_.scale],
        height=0.6)
for yi, v in zip(yy, cmp_.slope):
    ax.text(v + (0.08 if v >= 0 else -0.08), yi, f"{v:+.2f}", va="center", ha="left" if v >= 0 else "right",
            fontsize=8, color=INK)
ax.axvline(0, color=INK2, lw=1)
ax.set_yticks(yy, cmp_.scale, fontsize=8.5)
ax.set_xlabel("Days of peak shift per +100 m")
ax.set_title("Elevation effect at different scales")
ax.grid(axis="y", visible=False)
xl = ax.get_xlim(); ax.set_xlim(xl[0] - 0.6, xl[1] + 0.4)
save(fig, "E21_elevation_scales",
     "Elevation effect estimated between locations (orange), within locations (blue) and from the September "
     "temperature sensitivity of μ times the standard lapse rate (gray).")
