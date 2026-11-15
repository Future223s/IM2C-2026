"""Figures E1-E3 for the elevation side analysis. Reads pixels.csv + results.json."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import REGION_COLORS, REGION_MARKERS, date_axis, plt  # noqa: E402
from elevation_adjust import _two_piece, median_elev  # noqa: E402

FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
MOUNTAIN = ["Asheville", "Boone", "Shenandoah"]
INK, MUTED = "#2b2b29", "#5f5f5a"

res = json.load(open(HERE / "results.json"))
tp = res["two_piece_best"]
H, B, A = tp["hinge_m"], tp["below100"]["est"], tp["above100"]["est"]
LIN = res["headline_mountain"]["elev100"]

df = pd.read_csv(HERE / "pixels.csv")
df["anom"] = df.peak_doy - df.groupby(["region", "year"]).peak_doy.transform("median")
df["bin"] = (df.elev_m // 100 * 100).astype(int)
mtn = df[df.region.isin(MOUNTAIN)].copy()


def save(fig, name):
    fig.savefig(FIG / name, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("saved", (FIG / name).relative_to(HERE.parent.parent))


def binned(d, min_n=50):
    g = d.groupby("bin").anom
    b = pd.DataFrame({"n": g.size(), "med": g.median(), "q25": g.quantile(.25),
                      "q75": g.quantile(.75)}).reset_index()
    b = b[b.n >= min_n]
    b["x"] = b.bin + 50
    return b


# ---------- Fig E1 ----------
fig, ax = plt.subplots(figsize=(7.2, 4.6))
hb = ax.hexbin(mtn.elev_m, mtn.anom, gridsize=(60, 40), extent=(150, 1850, -30, 25),
               cmap="Greys", bins="log", mincnt=1, linewidths=0, alpha=0.45)
for i, r in enumerate(MOUNTAIN):
    b = binned(mtn[mtn.region == r])
    ax.plot(b.x + (i - 1) * 12, b.med, color=REGION_COLORS[r], marker=REGION_MARKERS[r],
            ms=4, lw=1.2, alpha=0.9, label=f"{r} (bin median)")
pb = binned(mtn)
ax.errorbar(pb.x, pb.med, yerr=[pb.med - pb.q25, pb.q75 - pb.med], fmt="o", color=INK,
            ms=5, lw=1, capsize=0, label="All mountain pixels: median and IQR", zorder=5)
# fitted lines; intercept = pixel-weighted mean of region intercepts → line through pooled mean
xs = np.linspace(mtn.elev_m.min(), mtn.elev_m.max(), 200)
f2 = np.array([_two_piece(x, H, B, A) for x in xs])
c2 = mtn.anom.mean() - np.mean([_two_piece(x, H, B, A) for x in mtn.elev_m.sample(20000, random_state=0)])
ax.plot(xs, f2 + c2, color=INK, lw=2, label=f"Two-piece fit (break {H:.0f} m)", zorder=6)
cl = mtn.anom.mean() - LIN["est"] * mtn.elev_m.mean() / 100
ax.plot(xs, cl + LIN["est"] * xs / 100, color=MUTED, lw=1.4, ls="--", label="Single-slope fit", zorder=6)
ax.axvline(H, color=MUTED, lw=0.8, ls=":")
ax.axhline(0, color="#8a8a85", lw=0.8)
txt = (f"Below {H:.0f} m: {B:+.1f} d/100 m (95% CI {tp['below100']['lo']:+.1f} to {tp['below100']['hi']:+.1f})\n"
       f"Above {H:.0f} m: {-A:.1f} days earlier per 100 m (95% CI {-tp['above100']['hi']:.1f}–{-tp['above100']['lo']:.1f})\n"
       f"Single slope: {-LIN['est']:.2f} days earlier per 100 m (95% CI {-LIN['hi']:.2f}–{-LIN['lo']:.2f})")
ax.text(0.02, 0.03, txt, transform=ax.transAxes, fontsize=7.5, color=INK, va="bottom",
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#c8c8c3", lw=0.6))
ax.set_ylim(-30, 25); ax.set_xlim(150, 1850)
ax.set_xlabel("Pixel mean elevation (m)")
ax.set_ylabel("Peak-color anomaly (days vs region-year median)")
ax.set_title("Peak color is flat with elevation up to ~1,350 m, then arrives much earlier",
             loc="left")
ax.legend(loc="upper right", fontsize=7)
fig.text(0.0, -0.02, "MODIS MCD12Q2 MidGreendown, 2001–2025, deciduous/mixed forest pixels (≥50%), "
         "Asheville, Boone and Shenandoah. Grey hexbins: all pixel-years (log density).",
         fontsize=6.5, color=MUTED)
save(fig, "E1_anomaly_vs_elevation.png")

# ---------- Fig E2 ----------
reg = pd.read_csv(HERE / "slopes_by_region.csv").set_index("region")
yr = pd.read_csv(HERE / "slopes_by_year.csv")
yr2 = pd.read_csv(HERE / "slopes_by_year_two_piece.csv")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.9), gridspec_kw={"width_ratios": [1, 1.15]})

rows = [("Mountain pooled, single slope", LIN, INK, "o"),
        (f"Mountain pooled, below {H:.0f} m", tp["below100"], INK, "o"),
        (f"Mountain pooled, above {H:.0f} m", tp["above100"], INK, "o")]
for r in ["Asheville", "Boone", "Shenandoah", "Raleigh"]:
    rows.append((f"{r} (single slope)", {"est": reg.loc[r, "slope"], "lo": reg.loc[r, "lo"],
                 "hi": reg.loc[r, "hi"]}, REGION_COLORS[r], REGION_MARKERS[r]))
for r in ["Asheville", "Boone"]:
    t = res[f"two_piece_{r}"]["above"]
    rows.append((f"{r}, above {H:.0f} m", t, REGION_COLORS[r], REGION_MARKERS[r]))
XL = (-6, 2.2)
for i, (lab, e, c, m) in enumerate(rows):
    y = -i
    est = e["est"]
    if est > XL[1]:   # Raleigh: off-scale; annotate instead
        a1.annotate(f"{est:+.1f} ({e['lo']:+.1f} to {e['hi']:+.1f}) →\nrange only ~90 m",
                    (XL[1], y), xytext=(XL[1] - 0.1, y), ha="right", va="center", fontsize=6.5, color=MUTED)
        continue
    a1.plot([e["lo"], e["hi"]], [y, y], color=c, lw=2, solid_capstyle="round")
    a1.plot(est, y, marker=m, color=c, ms=7, mec="white", mew=1)
a1.set_yticks([-i for i in range(len(rows))], [r[0] for r in rows], fontsize=7.5)
a1.axvline(0, color="#8a8a85", lw=0.8)
a1.set_xlim(*XL)
a1.set_xlabel("Days per 100 m (negative = earlier uphill)")
a1.set_title("Slope estimates with 95% year-bootstrap CI", loc="left")
a1.grid(axis="y", visible=False)

cats = [("Single slope", yr[yr.region == "Mountain pooled"].slope, INK),
        (f"Below {H:.0f} m", yr2.below, INK), (f"Above {H:.0f} m", yr2.above, INK)] + \
       [(r, yr[yr.region == r].slope, REGION_COLORS[r]) for r in ["Asheville", "Boone", "Shenandoah"]]
rng = np.random.default_rng(0)
for i, (lab, s, c) in enumerate(cats):
    a2.boxplot(s, positions=[i], widths=0.5, showfliers=False, patch_artist=True,
               boxprops=dict(fc="none", ec="#8a8a85", lw=0.8),
               medianprops=dict(color=INK, lw=1.4), whiskerprops=dict(color="#8a8a85", lw=0.8),
               capprops=dict(color="#8a8a85", lw=0.8))
    a2.scatter(i + rng.uniform(-0.15, 0.15, len(s)), s, s=14, color=c, alpha=0.75, lw=0, zorder=3)
a2.axhline(0, color="#8a8a85", lw=0.8)
a2.axvline(2.5, color="#c8c8c3", lw=0.8)
a2.set_xticks(range(len(cats)), [c[0] for c in cats], fontsize=7.5, rotation=20, ha="right")
a2.text(1, 1.02, "Mountain pooled", transform=a2.get_xaxis_transform(), ha="center", fontsize=7.5, color=MUTED)
a2.text(4, 1.02, "Single slope by region", transform=a2.get_xaxis_transform(), ha="center", fontsize=7.5, color=MUTED)
a2.set_ylabel("Slope fitted to one year (days per 100 m)")
a2.set_title("Year-by-year slopes, 2001–2025 (one dot per year)", loc="left", pad=16)
a2.grid(axis="x", visible=False)
fig.tight_layout()
save(fig, "E2_slopes_by_region_and_year.png")

# ---------- Fig E3 ----------
med = median_elev()
fig, ax = plt.subplots(figsize=(6.4, 4.2))
for r in ["Asheville", "Boone"]:
    d = df[df.region == r]
    base = d.groupby("year").peak_doy.median().mean()     # typical (2001–2025) region peak DOY
    b = binned(d)
    c = REGION_COLORS[r]
    ax.fill_between(b.x, base + b.q25, base + b.q75, color=c, alpha=0.12, lw=0)
    ax.plot(b.x, base + b.med, color=c, marker=REGION_MARKERS[r], ms=4.5, lw=1.4, label=f"{r}: bin median (IQR band)")
    xs = np.linspace(400, 1820, 200)
    fit = base + np.array([_two_piece(x, H, B, A) - _two_piece(med[r], H, B, A) for x in xs])
    ax.plot(xs, fit, color=c, lw=1, ls="--", alpha=0.9)
    ax.plot(med[r], base, marker="o", ms=9, mfc="white", mec=c, mew=1.6, zorder=5)
    ax.annotate(f"{r} median forest\nelevation {med[r]:.0f} m", (med[r], base),
                xytext=(med[r] - 330, base + 8) if r == "Asheville" else (med[r] + 80, base - 7),
                fontsize=7, color=INK, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
ax.axvline(H, color=MUTED, lw=0.8, ls=":")
ax.set_ylim(258, 296)
date_axis(ax, "y", step=7)
ax.set_xlabel("Pixel mean elevation (m)")
ax.set_ylabel("Typical peak-color date (MODIS MidGreendown)")
ax.set_title("Typical peak date vs elevation, Asheville and Boone, 2001–2025", loc="left")
ax.plot([], [], color=MUTED, ls="--", lw=1, label="Pooled two-piece rule")
ax.legend(loc="lower left", fontsize=7)
fig.text(0.0, -0.02, "Typical date = mean over years of the region-year median + binned median anomaly. "
         "MidGreendown is the mid-point of EVI decline, a proxy for peak color.", fontsize=6.5, color=MUTED)
save(fig, "E3_peak_date_vs_elevation.png")
