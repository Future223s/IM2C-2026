"""Elevation side analysis: days of peak-color shift per 100 m. Reads only pixels.csv.

Method: within each region x year, subtract the region-year median peak DOY from every pixel
(removes that year's weather). Regress the anomaly on elevation (per 100 m) with region
intercepts, so the slope is a purely within-region relationship. CIs come from a block
bootstrap over years (1,000 reps); pixel-clustered SEs are reported as a cross-check.

Outputs: results.json, slopes_by_region.csv, slopes_by_year.csv, binned_medians.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

HERE = Path(__file__).resolve().parent
MOUNTAIN = ["Asheville", "Boone", "Shenandoah"]
REGIONS = ["Raleigh", "Asheville", "Boone", "Shenandoah"]
N_BOOT = 1000
HINGE_M = 1000
HINGE_GRID = list(range(900, 1601, 50))
rng = np.random.default_rng(42)
YEARS = np.arange(2001, 2026)
W = np.stack([np.bincount(rng.integers(0, len(YEARS), len(YEARS)), minlength=len(YEARS))
              for _ in range(N_BOOT)])   # N_BOOT x n_years resample counts (same for all models)


def load():
    df = pd.read_csv(HERE / "pixels.csv")
    df["anom"] = df.peak_doy - df.groupby(["region", "year"]).peak_doy.transform("median")
    df["elev100"] = df.elev_m / 100
    df["hinge100"] = np.clip(df.elev_m - HINGE_M, 0, None) / 100
    return df


def design(d, extra):
    """Region dummies (no global intercept) + extra columns."""
    regs = sorted(d.region.unique())
    X = np.column_stack([(d.region == r).to_numpy(float) for r in regs]
                        + [d[c].to_numpy(float) for c in extra])
    return X, regs


def year_stats(d, extra):
    """Per-year sufficient statistics X'X, X'y, y'y (so the year bootstrap is cheap)."""
    X, regs = design(d, extra)
    y = d.anom.to_numpy(float)
    yi = np.searchsorted(YEARS, d.year.to_numpy())
    k = X.shape[1]
    XtX, Xty, yty = np.zeros((len(YEARS), k, k)), np.zeros((len(YEARS), k)), np.zeros(len(YEARS))
    for j in range(len(YEARS)):
        m = yi == j
        XtX[j], Xty[j], yty[j] = X[m].T @ X[m], X[m].T @ y[m], y[m] @ y[m]
    return XtX, Xty, yty, regs


def solve_w(stats, w):
    XtX, Xty, yty, _ = stats
    A, b = np.tensordot(w, XtX, 1), w @ Xty
    beta = np.linalg.solve(A, b)
    return beta, w @ yty - beta @ b          # coefficients, SSE


def fit_boot(d, extra=("elev100",)):
    """OLS with region intercepts; 95% CI from a block bootstrap over years."""
    d = d.dropna(subset=list(extra))
    stats = year_stats(d, extra)
    regs = stats[3]
    k = stats[0].shape[1]
    beta, _ = solve_w(stats, np.ones(len(YEARS)))
    boots = np.array([solve_w(stats, w)[0] for w in W])
    names = [f"region[{r}]" for r in regs] + list(extra)
    out = {}
    for i, nm in enumerate(names[len(regs):], start=len(regs)):
        lo, hi = np.percentile(boots[:, i], [2.5, 97.5])
        out[nm] = {"est": beta[i], "lo": lo, "hi": hi, "boot_sd": boots[:, i].std()}
    out["_n"] = int(len(d))
    out["_n_pixels"] = int(d.pixel_id.nunique())
    out["_n_years"] = int(d.year.nunique())
    return out


def with_hinge(d, h):
    return d.assign(below100=np.minimum(d.elev_m, h) / 100,
                    above100=np.clip(d.elev_m - h, 0, None) / 100)


def profile_hinge(d, controls=()):
    """Pick the breakpoint by least squares over HINGE_GRID; the bootstrap re-picks it in
    every replicate, so the CIs include breakpoint uncertainty."""
    extra = ("below100", "above100") + tuple(controls)
    d = d.dropna(subset=list(controls))
    stats = {h: year_stats(with_hinge(d, h), extra) for h in HINGE_GRID}
    nreg = len(stats[HINGE_GRID[0]][3])

    def best(w):
        fits = {h: solve_w(st, w) for h, st in stats.items()}
        h = min(fits, key=lambda h: fits[h][1])
        return h, fits[h][0], {hh: f[1] for hh, f in fits.items()}

    h0, b0, sse = best(np.ones(len(YEARS)))
    bh, bb = [], []
    for w in W:
        h, b, _ = best(w)
        bh.append(h); bb.append(b)
    bb = np.array(bb)
    out = {"hinge_m": h0, "hinge_lo": float(np.percentile(bh, 2.5)),
           "hinge_hi": float(np.percentile(bh, 97.5)),
           "sse_profile": {str(h): v / len(d) for h, v in sse.items()}}
    for i, nm in enumerate(extra):
        lo, hi = np.percentile(bb[:, nreg + i], [2.5, 97.5])
        out[nm] = {"est": b0[nreg + i], "lo": lo, "hi": hi}
    return out


def cluster_se(d, formula="anom ~ C(region) + elev100"):
    m = smf.ols(formula, d).fit(cov_type="cluster",
                                cov_kwds={"groups": pd.factorize(d.pixel_id)[0]})
    ci = m.conf_int().loc["elev100"]
    return {"est": m.params["elev100"], "lo": ci[0], "hi": ci[1], "se": m.bse["elev100"]}


def simple_slope(d):
    x = d.elev100 - d.groupby("region").elev100.transform("mean")
    return float((x * d.anom).sum() / (x ** 2).sum())


def fmt(r, k="elev100"):
    e = r[k] if k in r else r
    unit = "d/100 m" if k.endswith("100") else "d per unit"
    return f"{e['est']:+.2f} {unit} (95% CI {e['lo']:+.2f} to {e['hi']:+.2f})"


def main():
    df = load()
    mtn = df[df.region.isin(MOUNTAIN)]
    res = {}

    # 1. headline: mountain regions pooled
    res["headline_mountain"] = fit_boot(mtn)
    res["headline_mountain_clusterSE"] = cluster_se(mtn)
    res["all_regions"] = fit_boot(df)
    print("Headline (mountain, pooled):", fmt(res["headline_mountain"]))
    print("  pixel-clustered CI:        ", fmt(res["headline_mountain_clusterSE"]))
    print("All 4 regions:               ", fmt(res["all_regions"]))

    # 2. per region
    rows = []
    for r in REGIONS:
        d = df[df.region == r]
        b = fit_boot(d)["elev100"]
        px = d.drop_duplicates("pixel_id").elev_m
        rows.append({"region": r, "slope": b["est"], "lo": b["lo"], "hi": b["hi"],
                     "elev_p5": px.quantile(.05), "elev_p50": px.median(),
                     "elev_p95": px.quantile(.95), "n_pixels": d.pixel_id.nunique(),
                     "n_pixel_years": len(d)})
    by_region = pd.DataFrame(rows)
    by_region.to_csv(HERE / "slopes_by_region.csv", index=False)
    print("\nPer region:\n", by_region.round(2).to_string(index=False))

    # 3. per year
    rows = []
    for (y, r), d in df.groupby(["year", "region"]):
        rows.append({"year": y, "region": r, "slope": simple_slope(d)})
    for y, d in mtn.groupby("year"):
        rows.append({"year": y, "region": "Mountain pooled", "slope": simple_slope(d)})
    by_year = pd.DataFrame(rows)
    by_year.to_csv(HERE / "slopes_by_year.csv", index=False)
    s = by_year[by_year.region == "Mountain pooled"].slope
    res["yearly_mountain"] = {"median": s.median(), "min": s.min(), "max": s.max(),
                              "iqr_lo": s.quantile(.25), "iqr_hi": s.quantile(.75),
                              "sd": s.std(), "n_negative": int((s < 0).sum()), "n": int(len(s)),
                              "min_year": int(by_year.loc[s.idxmin(), "year"]),
                              "max_year": int(by_year.loc[s.idxmax(), "year"])}
    print("\nYearly mountain slopes:", {k: round(float(v), 2) for k, v in res["yearly_mountain"].items()})
    print(by_year.pivot(index="year", columns="region", values="slope").round(2).to_string())

    # 4. linearity: binned medians + hinge model
    df["bin"] = (df.elev_m // 100 * 100).astype(int)
    g = df.groupby(["region", "bin"]).anom
    binned = pd.DataFrame({"n": g.size(), "median": g.median(),
                           "q25": g.quantile(.25), "q75": g.quantile(.75),
                           "mean": g.mean()}).reset_index()
    binned = binned[binned.n >= 50]
    binned.to_csv(HERE / "binned_medians.csv", index=False)
    print("\nBinned medians (>=50 pixel-years):")
    print(binned.pivot(index="bin", columns="region", values="median").round(1).to_string())

    # (a) the handoff's fixed 1,000 m hinge
    hm = fit_boot(with_hinge(mtn, HINGE_M), extra=("below100", "above100"))
    hc = fit_boot(mtn.assign(hinge100=np.clip(mtn.elev_m - HINGE_M, 0, None) / 100),
                  extra=("elev100", "hinge100"))
    res["two_piece_1000"] = {"below": hm["below100"], "above": hm["above100"],
                             "change": hc["hinge100"]}
    print(f"\nTwo-piece, fixed hinge {HINGE_M} m: below {fmt(hm, 'below100')}, "
          f"above {fmt(hm, 'above100')}; change {fmt(hc, 'hinge100')}")
    # (b) breakpoint chosen by least squares (re-chosen in each bootstrap replicate)
    for key, d, ctrl in [("two_piece_best", mtn, ()),
                         ("two_piece_best_controls", mtn, ("northness", "forest_frac")),
                         ("two_piece_best_AshBoone", mtn[mtn.region != "Shenandoah"], ())]:
        pr = profile_hinge(d, ctrl)
        res[key] = pr
        print(f"{key}: hinge {pr['hinge_m']} m (95% {pr['hinge_lo']:.0f}-{pr['hinge_hi']:.0f}); "
              f"below {fmt(pr, 'below100')}; above {fmt(pr, 'above100')}")
    hb = res["two_piece_best"]["hinge_m"]
    rows = []
    for y, d in with_hinge(mtn, hb).groupby("year"):
        X, _ = design(d, ("below100", "above100"))
        b = np.linalg.lstsq(X, d.anom.to_numpy(float), rcond=None)[0]
        rows.append({"year": y, "below": b[-2], "above": b[-1]})
    yr2 = pd.DataFrame(rows)
    yr2.to_csv(HERE / "slopes_by_year_two_piece.csv", index=False)
    res["yearly_two_piece"] = {"above_median": yr2.above.median(), "above_min": yr2.above.min(),
                               "above_max": yr2.above.max(), "above_n_negative": int((yr2.above < 0).sum()),
                               "below_median": yr2.below.median(), "below_min": yr2.below.min(),
                               "below_max": yr2.below.max()}
    print("  yearly two-piece:", {k: round(float(v), 2) for k, v in res["yearly_two_piece"].items()})
    px_m = mtn.drop_duplicates("pixel_id")
    res["n_pixels_above_hinge"] = px_m[px_m.elev_m > hb].groupby("region").size().to_dict()
    print("  pixels above hinge:", res["n_pixels_above_hinge"])
    for r in ["Asheville", "Boone"]:
        t = fit_boot(with_hinge(df[df.region == r], hb), extra=("below100", "above100"))
        res[f"two_piece_{r}"] = {"below": t["below100"], "above": t["above100"],
                                 "n_pixels_above": int((df[df.region == r]
                                                        .drop_duplicates("pixel_id").elev_m > hb).sum())}
        print(f"  {r} at {hb} m: below {fmt(t, 'below100')}, above {fmt(t, 'above100')}")

    # 5. aspect
    res["with_northness"] = fit_boot(mtn, extra=("elev100", "northness"))
    print("\nWith northness: elev", fmt(res["with_northness"]),
          "| northness", fmt(res["with_northness"], "northness"))
    res["with_controls"] = fit_boot(mtn, extra=("elev100", "northness", "forest_frac"))
    print("With northness + forest_frac: elev", fmt(res["with_controls"]),
          "| forest_frac", fmt(res["with_controls"], "forest_frac"))
    px = mtn.drop_duplicates("pixel_id")
    res["corr_elev_forest_frac"] = float(px.elev_m.corr(px.forest_frac))
    res["corr_elev_northness"] = float(px.elev_m.corr(px.northness))

    # 6. Helene sensitivity
    no_hel = mtn[~((mtn.year == 2024) & mtn.region.isin(["Asheville", "Boone"]))]
    res["no_helene"] = fit_boot(no_hel)
    print("Excluding 2024 Asheville/Boone:", fmt(res["no_helene"]))
    t = fit_boot(with_hinge(no_hel, hb), extra=("below100", "above100"))
    full = fit_boot(with_hinge(mtn, hb), extra=("below100", "above100"))
    res["two_piece_fixed_best"] = {"below": full["below100"], "above": full["above100"]}
    res["no_helene_two_piece"] = {"below": t["below100"], "above": t["above100"]}
    print(f"  two-piece at {hb} m, full: below {fmt(full, 'below100')}, above {fmt(full, 'above100')}")
    print(f"  two-piece at {hb} m, no Helene: below {fmt(t, 'below100')}, above {fmt(t, 'above100')}")

    json.dump(res, open(HERE / "results.json", "w"), indent=2, default=float)
    print("\nwrote results.json, slopes_by_region.csv, slopes_by_year.csv, binned_medians.csv")


if __name__ == "__main__":
    main()
