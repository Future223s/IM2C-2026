"""Shared helpers for the pre-modeling analysis scripts."""
from pathlib import Path
import datetime as dt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT.parent / "data" / "data"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(exist_ok=True)

REGIONS = ["Raleigh", "Asheville", "Boone", "Shenandoah"]
# Fixed categorical order (blue, orange, aqua, yellow) + marker shape as secondary encoding
REGION_COLORS = {"Raleigh": "#2a78d6", "Asheville": "#eb6834",
                 "Boone": "#1baf7a", "Shenandoah": "#eda100"}
REGION_MARKERS = {"Raleigh": "o", "Asheville": "s", "Boone": "^", "Shenandoah": "D"}

TARGET = "peak_doy"
PRE_PEAK = ["tmax_Aug", "tmin_Sep", "cdd_Sep", "precip_AugSep", "vpd_AugSep",
            "pdsi_AugSep", "srad_Sep", "mid_greenup_p50"]
PARTLY = ["cdd_Sep1_Oct15", "cold_nights_Sep1_Oct15"]
POST = ["tmin_SepOct", "first_frost_doy", "precip_Oct", "srad_SepOct"]
FEATURES = PRE_PEAK + PARTLY + POST  # year-varying candidates (year handled separately)
TIMING = {**{f: "pre-peak" for f in PRE_PEAK}, **{f: "partly" for f in PARTLY},
          **{f: "post" for f in POST}, "year": "pre-peak"}

plt.rcParams.update({
    "figure.dpi": 100, "savefig.dpi": 200, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e6e6e3", "grid.linewidth": 0.6,
    "axes.edgecolor": "#8a8a85", "axes.labelcolor": "#2b2b29",
    "xtick.color": "#5f5f5a", "ytick.color": "#5f5f5a", "font.size": 9,
    "axes.titlesize": 10, "legend.frameon": False,
})


def load(include_forecast=False):
    """model_df with the region column renamed to `region`, Charleston dropped."""
    df = pd.read_csv(DATA_DIR / "model_df.csv").rename(columns={"name": "region"})
    df = df[df.region != "Charleston"].copy()
    if not include_forecast:
        df = df[df[TARGET].notna()].copy()
    return df.reset_index(drop=True)


def anomalies(df, cols):
    """Subtract each region's mean (over the rows in df) from each column."""
    out = df[["region", "year", "helene_2024"]].copy()
    for c in cols:
        out[c] = df[c] - df.groupby("region")[c].transform("mean")
    return out


def doy_to_label(doy, year=2001):
    d = dt.date(year, 1, 1) + dt.timedelta(days=float(doy) - 1)
    return d.strftime("%b %-d")


def date_axis(ax, axis="y", step=5):
    """Relabel a DOY axis with calendar dates (non-leap year)."""
    lo, hi = (ax.get_ylim() if axis == "y" else ax.get_xlim())
    ticks = np.arange(np.ceil(lo / step) * step, hi + 0.01, step)
    labels = [doy_to_label(t) for t in ticks]
    if axis == "y":
        ax.set_yticks(ticks, labels)
    else:
        ax.set_xticks(ticks, labels)


def save(fig, name):
    path = FIG_DIR / name
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    print(f"saved {path.relative_to(ROOT.parent)}")
