"""Stage 7 validation for data/panel_v1 -> panel_v1/BUILD_LOG.md (BUILD_DATASET.md §7).

Each check is PASS, FLAG (report to Arin, not blocking) or FAIL (blocking: stop and investigate).
"""
import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd

from locations import ORDER, LOCATIONS, CARRIED_OVER

REG_PANEL = ["peak_doy", "peak_doy_viirs", "tmax_Aug", "precip_AugSep", "cdd_Sep", "pdsi_AugSep",
             "first_frost_doy"]
REG_STATIC = ["elev_median_m", "frac_hardwood"]
TOL = 1e-6
# Columns that are percentiles. These were switched from ee.Reducer.percentile (approximate, depends on the
# other regions in the request) to exact pixel-center histograms, so they cannot match model_df.csv to 1e-6.
# Agreed with Arin 2026-11-15: exact match for everything else, these tolerances for the percentiles.
PCT_TOL = {"peak_doy": 1.5, "peak_doy_viirs": 1.5, "elev_median_m": 5.0}


def md_table(df, floatfmt="{:.2f}"):
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, (float, np.floating)):
                cells.append("" if pd.isna(v) else floatfmt.format(v))
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def validate(out_dir, raw_dir, notes_dir, old_model_df):
    out_dir, raw_dir, notes_dir = Path(out_dir), Path(raw_dir), Path(notes_dir)
    panel = pd.read_csv(out_dir / "panel.csv")
    locs = pd.read_csv(out_dir / "locations.csv")
    pdsi = pd.read_csv(out_dir / "pdsi.csv")
    old = pd.read_csv(old_model_df).rename(columns={"name": "location"})
    status = {}
    L = []   # log lines

    tr = panel[panel["year"].between(2001, 2025)]

    # ---------------------------------------------------------------- 1. regression
    L.append("## 1. Regression test (Boone, Asheville, Shenandoah vs `data/data/model_df.csv`)\n")
    new = tr[tr["location"].isin(CARRIED_OVER)].merge(
        locs[["location"] + REG_STATIC], on="location")
    o = old[old["location"].isin(CARRIED_OVER) & old["year"].between(2001, 2025)]
    m = new.merge(o[["location", "year"] + REG_PANEL + REG_STATIC], on=["location", "year"],
                  suffixes=("_new", "_old"), how="outer", indicator=True)
    rows, reg_ok = [], True
    if (m["_merge"] != "both").any():
        reg_ok = False
        L.append(f"Row mismatch: {m.loc[m['_merge'] != 'both', ['location', 'year', '_merge']].values.tolist()}\n")
    for c in REG_PANEL + REG_STATIC:
        a, b = m[f"{c}_new"], m[f"{c}_old"]
        nan_mismatch = int((a.isna() != b.isna()).sum())
        diff = (a - b).abs()
        mx = float(diff.max()) if diff.notna().any() else 0.0
        worst = ""
        if mx >= TOL:
            i = diff.idxmax()
            worst = f"{m.at[i, 'location']} {m.at[i, 'year']}: {b[i]:.4f} → {a[i]:.4f}"
        tol = PCT_TOL.get(c, TOL)
        n_bad = int((diff >= TOL).sum())
        ok = mx < tol and nan_mismatch == 0
        reg_ok &= ok
        rows.append(dict(column=c, n_compared=int(diff.notna().sum()), max_abs_diff=mx,
                         mean_abs_diff=float(diff.mean()) if diff.notna().any() else 0.0,
                         tolerance=tol, n_diff_ge_1e6=n_bad, nan_mismatch=nan_mismatch, worst=worst,
                         result="PASS" if ok else "FAIL"))
    reg = pd.DataFrame(rows)
    L.append(md_table(reg, "{:.3g}"))
    status["1 regression"] = "PASS" if reg_ok else "FAIL"
    reg.to_csv(raw_dir / "regression_detail.csv", index=False)
    if not reg_ok:
        # per-row detail of every differing value for investigation
        det = []
        for c in REG_PANEL + REG_STATIC:
            d = m[(m[f"{c}_new"] - m[f"{c}_old"]).abs() >= PCT_TOL.get(c, TOL)]
            for _, r in d.iterrows():
                det.append(dict(column=c, location=r["location"], year=r["year"],
                                old=r[f"{c}_old"], new=r[f"{c}_new"],
                                diff=r[f"{c}_new"] - r[f"{c}_old"]))
        pd.DataFrame(det).to_csv(raw_dir / "regression_diffs.csv", index=False)
        L.append("\nEvery differing value: `raw/regression_diffs.csv`.")

    # ---------------------------------------------------------------- 2. shape
    L.append("\n## 2. Shape\n")
    exp_rows = len(ORDER) * 26
    dups = int(panel.duplicated(["location", "year"]).sum())
    gaps = tr[tr["peak_doy"].isna()][["location", "year"]].values.tolist()
    f26 = panel[panel["year"] == 2026]["peak_doy"].notna().sum()
    L.append(f"- rows: {len(panel)} (expected {exp_rows}); duplicate (location, year): {dups}")
    L.append(f"- 2001–2025 rows with null `peak_doy`: {len(gaps)} {gaps if gaps else ''}")
    L.append(f"- 2026 rows with non-null `peak_doy`: {f26}")
    shape_ok = len(panel) == exp_rows and dups == 0 and f26 == 0
    status["2 shape"] = ("PASS" if shape_ok and not gaps else "FLAG" if shape_ok else "FAIL")

    # ---------------------------------------------------------------- 3. ranges
    L.append("\n## 3. Ranges\n")
    frac_cols = [c for c in locs.columns if c.startswith(("frac_", "evt_"))] + ["hardwood_share_of_forest"]
    checks = [
        ("panel.peak_doy", tr["peak_doy"], 230, 350),
        ("panel.tmax_Aug", panel["tmax_Aug"], 15, 35),
        ("panel.precip_AugSep", panel["precip_AugSep"], 0, np.inf),
        ("panel.pdsi_AugSep", panel["pdsi_AugSep"], -8, 8),
        ("pdsi.pdsi (daily file)", pdsi["pdsi"], -8, 8),
        ("panel.qa_retention", panel["qa_retention"], 0, 1),
        ("locations.elev_median_m", locs["elev_median_m"], 50, 1500),
    ] + [(f"locations.{c}", locs[c], 0, 1) for c in frac_cols]
    rr, range_ok = [], True
    for name, s, lo, hi in checks:
        s = s.dropna()
        bad = int(((s < lo) | (s > hi)).sum())
        range_ok &= bad == 0
        rr.append(dict(check=name, n=len(s), min=s.min() if len(s) else np.nan,
                       max=s.max() if len(s) else np.nan, allowed=f"[{lo}, {hi}]",
                       out_of_range=bad))
    L.append(md_table(pd.DataFrame(rr), "{:.3f}"))
    status["3 ranges"] = "PASS" if range_ok else "FAIL"

    # ---------------------------------------------------------------- 4. geography
    L.append("\n## 4. Geography sense-check (north → south)\n")
    g = locs[["location", "lat", "elev_median_m", "peak_doy_mean", "peak_doy_sd"]].copy()
    g["peak_date_typical"] = g["peak_doy_mean"].map(
        lambda d: (dt.date(2001, 1, 1) + dt.timedelta(days=float(d) - 1)).strftime("%b %-d"))
    L.append(md_table(g, "{:.1f}"))
    r_lat = np.corrcoef(g["lat"], g["peak_doy_mean"])[0, 1]
    north = g[g["location"].isin(["Montpelier", "Old Forge"])]["peak_doy_mean"].mean()
    south = g[~g["location"].isin(["Montpelier", "Old Forge"])]["peak_doy_mean"].mean()
    L.append(f"\n- r(lat, peak_doy_mean) = {r_lat:+.2f}. Montpelier + Old Forge mean {north:.1f} vs others {south:.1f}.")
    inverted = north > south
    status["4 geography"] = "FLAG" if inverted else "PASS"
    if inverted:
        L.append("- **FLAG: northern locations are not earlier on average. Not fixed (see spec).**")

    # ---------------------------------------------------------------- 5. missingness
    L.append("\n## 5. Missingness (non-null share)\n")
    target_like = ["peak_doy", "peak_date", "peak_anom", "peak_doy_p25", "peak_doy_p75", "peak_spread_days",
                   "mid_greendown_count", "forest_px_count", "qa_retention", "mid_greenup_p50",
                   "senescence_p50", "dormancy_p50"]
    rows, flagged = [], []
    for c in panel.columns:
        if c == "peak_doy_viirs":
            base, scope = panel[panel["year"].between(2013, 2025)], "2013–2025"
        elif c in target_like:
            base, scope = tr, "2001–2025"
        else:
            base, scope = panel, "all rows"
        share = base[c].notna().mean()
        rows.append(dict(file="panel", column=c, scope=scope, non_null=int(base[c].notna().sum()),
                         n=len(base), share=share))
        if share < 0.95:
            flagged.append(f"panel.{c} ({share:.0%} of {scope})")
    for c in locs.columns:
        share = locs[c].notna().mean()
        rows.append(dict(file="locations", column=c, scope="all", non_null=int(locs[c].notna().sum()),
                         n=len(locs), share=share))
        if share < 0.95:
            flagged.append(f"locations.{c} ({share:.0%})")
    miss = pd.DataFrame(rows)
    miss.to_csv(raw_dir / "missingness.csv", index=False)
    low = miss[miss["share"] < 1]
    L.append("Columns below 100% (full table: `raw/missingness.csv`):\n")
    L.append(md_table(low, "{:.3f}") if len(low) else "- none: every column is 100% populated within its scope.")
    L.append(f"\n- Columns under 95%: {flagged if flagged else 'none'}")
    status["5 missingness"] = "FLAG" if flagged else "PASS"

    # ---------------------------------------------------------------- 6. MODIS vs VIIRS
    L.append("\n## 6. MODIS vs VIIRS (2013–2025)\n")
    mv = panel[panel["year"].between(2013, 2025)].dropna(subset=["peak_doy", "peak_doy_viirs"])
    rows = []
    for loc in ORDER:
        d = mv[mv["location"] == loc]
        if len(d) < 3:
            continue
        diff = d["peak_doy_viirs"] - d["peak_doy"]
        rows.append(dict(location=loc, n=len(d), bias_d=diff.mean(), rmse_d=np.sqrt((diff ** 2).mean()),
                         r=np.corrcoef(d["peak_doy"], d["peak_doy_viirs"])[0, 1]))
    diff = mv["peak_doy_viirs"] - mv["peak_doy"]
    an = mv.assign(a=mv["peak_doy"] - mv.groupby("location")["peak_doy"].transform("mean"),
                   b=mv["peak_doy_viirs"] - mv.groupby("location")["peak_doy_viirs"].transform("mean"))
    rows.append(dict(location="All (pooled)", n=len(mv), bias_d=diff.mean(),
                     rmse_d=np.sqrt((diff ** 2).mean()), r=np.corrcoef(an["a"], an["b"])[0, 1]))
    L.append(md_table(pd.DataFrame(rows)))
    L.append("\n(Pooled r is between within-location anomalies.) FINDINGS §1 (4 regions): bias ≈ +8 d, "
             "anomaly r = 0.50.")
    gap = mv.assign(gap=diff).groupby("year")["gap"].mean()
    slope = np.polyfit(gap.index, gap.values, 1)[0] * 10
    L.append(f"- VIIRS − MODIS gap trend, mean over locations: {slope:+.1f} d/decade (FINDINGS: +7.7).")
    status["6 MODIS vs VIIRS"] = "PASS"

    # ---------------------------------------------------------------- 7. pixel counts
    L.append("\n## 7. MODIS forest pixel counts\n")
    px = locs[["location", "forest_px_count", "region_area_km2", "frac_hardwood"]].copy()
    px["flag"] = np.where(px["forest_px_count"] < 2000, "**< 2,000**", "")
    L.append(md_table(px, "{:,.3f}").replace(".000 |", " |"))
    weak = px[px["forest_px_count"] < 2000]["location"].tolist()
    status["7 pixel counts"] = "FLAG" if weak else "PASS"
    L.append(f"\n- Locations under 2,000 forest pixels: {weak if weak else 'none'} (Raleigh had 2,905).")

    # ---------------------------------------------------------------- 8. climatology baseline
    L.append("\n## 8. Climatology baseline (leave-one-year-out)\n")
    d = tr.dropna(subset=["peak_doy"])[["location", "year", "peak_doy"]].copy()
    s = d.groupby("location")["peak_doy"].transform("sum")
    n = d.groupby("location")["peak_doy"].transform("count")
    d["pred"] = (s - d["peak_doy"]) / (n - 1)
    d["err"] = d["peak_doy"] - d["pred"]
    rows = []
    for loc in ORDER:
        e = d[d["location"] == loc]["err"]
        rows.append(dict(location=loc, n=len(e), MAE=e.abs().mean(), RMSE=np.sqrt((e ** 2).mean())))
    rows.append(dict(location="**Overall**", n=len(d), MAE=d["err"].abs().mean(),
                     RMSE=np.sqrt((d["err"] ** 2).mean())))
    base = pd.DataFrame(rows)
    base.to_csv(raw_dir / "climatology_baseline.csv", index=False)
    L.append(md_table(base))
    L.append("\n4-region reference (FINDINGS §2): MAE 3.22 / RMSE 3.95.")
    status["8 baseline"] = "PASS"

    # ---------------------------------------------------------------- write log
    head = [f"# BUILD_LOG — panel_v1\n",
            f"Built {dt.datetime.now():%Y-%m-%d %H:%M} by `data/build_panel.py`. "
            f"Spec: `docs/BUILD_DATASET.md`.\n",
            "## Check summary\n", "| check | result |", "|---|---|"]
    head += [f"| {k} | **{v}** |" for k, v in status.items()]
    build = ["\n# Build notes\n"]
    for stage in ["regions", "pheno", "weather", "static", "features", "assemble", "deviations"]:
        f = notes_dir / f"{stage}.md"
        if f.exists():
            build.append(f"## {stage}\n")
            build.append(f.read_text())
    text = "\n".join(head) + "\n" + "\n".join(build) + "\n# Validation (Stage 7)\n\n" + "\n".join(L) + "\n"
    (out_dir / "BUILD_LOG.md").write_text(text)
    print("\n".join(f"{k}: {v}" for k, v in status.items()))
    return all(v != "FAIL" for v in status.values())
