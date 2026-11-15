"""Shared loading + figure style for the environmental (location-effect) analysis.

Data: data/panel_v1 (11 locations x 2001-2026). Run scripts from this folder with the fall-foliage env.
"""
import datetime as dt
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
PANEL = ROOT.parents[1] / "data" / "panel_v1"
FIG = ROOT / "figures"
TAB = ROOT / "tables"
FIG.mkdir(exist_ok=True)
TAB.mkdir(exist_ok=True)

# --- palette (dataviz reference palette, light mode) ---
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
RED = "#e34948"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SEQ_CMAP = LinearSegmentedColormap.from_list("seq_blue", SEQ)
DIV_CMAP = LinearSegmentedColormap.from_list(
    "div_br", ["#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f3a7a6", "#e66767", "#b8302f"])

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight", "savefig.facecolor": SURFACE,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "font.size": 9.5, "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.edgecolor": GRID, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
    "legend.frameon": False, "lines.linewidth": 2, "lines.markersize": 7,
})


def load():
    panel = pd.read_csv(PANEL / "panel.csv")
    locs = pd.read_csv(PANEL / "locations.csv")
    order = list(locs["location"])
    train = panel[panel["is_forecast_row"] == 0].copy()
    return panel, locs, train, order


def doy_label(doy, year=2001):
    return (dt.date(year, 1, 1) + dt.timedelta(days=float(doy) - 1)).strftime("%b %-d")


def save(fig, name, caption):
    """Save a figure and register its caption (figures/CAPTIONS.md is rebuilt from these)."""
    fig.savefig(FIG / f"{name}.png")
    plt.close(fig)
    cap = FIG / "captions.tsv"
    rows = {}
    if cap.exists():
        for line in cap.read_text().splitlines():
            k, v = line.split("\t", 1)
            rows[k] = v
    rows[name] = caption.replace("\n", " ")
    cap.write_text("\n".join(f"{k}\t{v}" for k, v in sorted(rows.items())) + "\n")
    print(f"  saved figures/{name}.png")


def year_cluster_ols(formula, data):
    import statsmodels.formula.api as smf
    return smf.ols(formula, data=data).fit(cov_type="cluster", cov_kwds={"groups": data["year"]})


def label_points(ax, x, y, names, dx=4, dy=3, overrides=None, size=8.5):
    """Direct labels next to scatter points; overrides: name -> (dx, dy, ha)."""
    overrides = overrides or {}
    for xi, yi, n in zip(x, y, names):
        ddx, ddy, ha = overrides.get(n, (dx, dy, "left"))
        ax.annotate(n, (xi, yi), xytext=(ddx, ddy), textcoords="offset points",
                    fontsize=size, color=INK2, ha=ha, va="bottom")
