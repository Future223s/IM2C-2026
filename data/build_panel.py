"""Build the 11-location Appalachian peak-foliage panel (v1) -> data/panel_v1/.

Spec: docs/BUILD_DATASET.md (build steps) and docs/DATASET_GUIDE.md (output columns).
Extends the Earth Engine pipeline in data/notebooks/DataPreprocessing.ipynb; helpers live in
ee_helpers.py. Never writes to data/data/ or analysis/.

Usage (fall-foliage env):
    python build_panel.py                      # all stages, using cached raw pulls
    python build_panel.py --stage regions      # one stage: regions|pheno|weather|static|features|assemble|bundle|validate|pixels
    python build_panel.py --stage weather --refresh   # re-pull from Earth Engine even if cached
"""
import argparse
import datetime as dt
import json
import re
import shutil
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from locations import (LOCATIONS, ORDER, BUFFER_M, FOREST_THRESHOLD, PHENO_YEARS, WX_YEARS,  # noqa: E402
                       WX_START, WX_END, HELENE, PROJECT, CARRIED_OVER)

DATA = Path(__file__).resolve().parent
OUT = DATA / "panel_v1"
RAW = OUT / "raw"
CHUNKS = RAW / "chunks"
NOTES = RAW / "notes"
OLD_MODEL_DF = DATA / "data" / "model_df.csv"
GUIDE_SRC = DATA.parent / "docs" / "DATASET_GUIDE.md"
for d in (OUT, RAW, CHUNKS, NOTES):
    d.mkdir(parents=True, exist_ok=True)

REFRESH = False
_EE_READY = False


def ee_init():
    global _EE_READY, ee, H
    if not _EE_READY:
        import ee  # noqa: F401
        import ee_helpers as H  # noqa: F401
        ee.Initialize(project=PROJECT)
        globals()["ee"], globals()["H"] = ee, H
        _EE_READY = True


def regions_fc():
    ee_init()
    return H.build_study_regions(LOCATIONS, BUFFER_M)


def cached(path, fn, parse_dates=None):
    """Read path if it exists (and --refresh not given), else run fn() and save its DataFrame."""
    path = Path(path)
    if path.exists() and not REFRESH:
        return pd.read_csv(path, parse_dates=parse_dates)
    df = fn()
    df.to_csv(path, index=False)
    return pd.read_csv(path, parse_dates=parse_dates)   # round-trip so cached / fresh are identical


def note(stage, text):
    """Append a line to this stage's build notes (collected into BUILD_LOG.md by `validate`)."""
    print(text)
    with open(NOTES / f"{stage}.md", "a") as f:
        f.write(text.rstrip() + "\n")


def reset_notes(stage):
    (NOTES / f"{stage}.md").unlink(missing_ok=True)


def to_location(df):
    return df.rename(columns={"name": "location"})


def sort_loc(df, *more):
    df = df.copy()
    df["_o"] = df["location"].map({n: i for i, n in enumerate(ORDER)})
    return df.sort_values(["_o", *more]).drop(columns="_o").reset_index(drop=True)


# =============================================================== Stage 1: regions
def stage_regions():
    reset_notes("regions")
    ee_init()
    regions = regions_fc()

    def pull():
        fc = regions.map(lambda f: f.set("area_km2", f.geometry().area(100).divide(1e6)))
        df = H.fc_to_df(fc)[["name", "ecoregion", "area_km2"]]
        return df

    reg = cached(RAW / "regions.csv", pull)
    note("regions", "### Regions (40 km buffer clipped to EPA L3 ecoregion)\n")
    note("regions", "| location | ecoregion | region_area_km2 | flag |\n|---|---|---|---|")
    for _, r in sort_loc(to_location(reg)).iterrows():
        flag = "**< 1,500 km²**" if r.area_km2 < 1500 else ""
        note("regions", f"| {r.location} | {r.ecoregion} | {r.area_km2:,.0f} | {flag} |")

    def pull_overlaps():
        rows = []
        feats = {n: regions.filter(ee.Filter.eq("name", n)).first().geometry() for n in ORDER}
        pts = {n: np.array(LOCATIONS[n][:2]) for n in ORDER}
        for i, a in enumerate(ORDER):
            for b in ORDER[i + 1:]:
                # only pairs whose centers are < 2 buffers apart can overlap (cheap pre-filter)
                dlat = (pts[a][1] - pts[b][1]) * 111.0
                dlon = (pts[a][0] - pts[b][0]) * 111.0 * np.cos(np.radians((pts[a][1] + pts[b][1]) / 2))
                dist = float(np.hypot(dlat, dlon))
                if dist >= 2 * BUFFER_M / 1000 + 2:
                    rows.append(dict(a=a, b=b, center_dist_km=dist, overlap_km2=0.0))
                    continue
                ov = feats[a].intersection(feats[b], ee.ErrorMargin(100)).area(100).divide(1e6).getInfo()
                rows.append(dict(a=a, b=b, center_dist_km=dist, overlap_km2=ov))
        return pd.DataFrame(rows)

    ov = cached(RAW / "region_overlaps.csv", pull_overlaps)
    close = ov.sort_values("center_dist_km").head(4)
    note("regions", "\nClosest pairs (center distance, overlap of clipped regions):\n")
    for _, r in close.iterrows():
        note("regions", f"- {r.a} – {r.b}: {r.center_dist_km:.0f} km apart, overlap {r.overlap_km2:,.1f} km²")
    n_ov = (ov.overlap_km2 > 0.01).sum()
    note("regions", f"\nPairs with any overlap: {n_ov}")
    return reg


# =============================================================== Stage 3a: phenology
def phen_from_hist(hist):
    """Exact p25/p50/p75 + count per (region, year, product, band) from the pixel histograms."""
    from ee_helpers import hist_percentile, COMMON
    rows = {}
    for (name, year, product, band), d in hist.groupby(["name", "year", "product", "band"]):
        r = rows.setdefault((name, year, product), dict(name=name, year=year, product=product))
        if band == "forest_px":
            r["forest_px_count"] = d["n"].sum()
            continue
        for q in (25, 50, 75):
            r[f"{band}_p{q}"] = hist_percentile(d["value"].values, d["n"].values, q)
        r[f"{band}_count"] = d["n"].sum()
    phen = pd.DataFrame(list(rows.values()))
    reg = pd.read_csv(RAW / "regions.csv")[["name", "ecoregion"]]
    phen = phen.merge(reg, on="name", how="left")
    keep = (["name", "ecoregion", "year", "product", "forest_px_count"]
            + [f"{b}_{s}" for b in COMMON for s in ("p25", "p50", "p75", "count")])
    return phen.reindex(columns=keep)


def stage_pheno():
    reset_notes("pheno")
    path = RAW / "phenology.csv"
    if path.exists() and not REFRESH:
        phen = pd.read_csv(path)
    else:
        ee_init()
        regions = regions_fc()
        hists = []
        for spec, product in ((H.MODIS, "MODIS"), (H.VIIRS, "VIIRS")):
            years, enc = H.phenology_setup(spec, product, PHENO_YEARS, regions)
            note("pheno", f"- {product}: years {years[0]}–{years[-1]} ({len(years)}), encoding `{enc}`")
            for year in years:
                f = CHUNKS / f"phenohist_{product}_{year}.csv"
                hists.append(cached(f, lambda: H.pull_phenology_hist_year(
                    spec, product, year, enc, regions, FOREST_THRESHOLD)))
                print(f"  {product} {year} done", end="\r")
            print()
        hist = pd.concat(hists, ignore_index=True)
        hist.to_csv(RAW / "phenology_histograms.csv", index=False)
        phen = phen_from_hist(hist)
        phen["qa_retention"] = phen["mid_greendown_count"] / phen["forest_px_count"]
        phen.to_csv(path, index=False)
        phen = pd.read_csv(path)

    note("pheno", "- Percentiles are EXACT (unweighted per-region pixel histograms of integer-day dates; "
                  "numpy 'linear' percentile definition), not `ee.Reducer.percentile` (see deviations).")
    mod = phen[phen["product"] == "MODIS"]
    last = mod["year"].max()
    note("pheno", f"- MODIS MCD12Q2 most recent year: **{last}**")
    have = mod[(mod["year"] == 2025) & mod["mid_greendown_p50"].notna()]["name"]
    miss = sorted(set(ORDER) - set(have))
    note("pheno", f"- MODIS 2025 target missing for: {miss if miss else 'none'}")
    odd = phen[(phen["mid_greendown_p50"] < 230) | (phen["mid_greendown_p50"] > 350)]
    note("pheno", f"- Rows with mid_greendown_p50 outside DOY 230–350: {len(odd)}")
    for _, r in odd.iterrows():
        note("pheno", f"  - {r['name']} {r.year} {r['product']}: {r.mid_greendown_p50:.1f}")
    return phen


# =============================================================== Stage 3b: weather
def stage_weather():
    reset_notes("weather")
    wpath, ppath = RAW / "weather_daily.csv", RAW / "pdsi.csv"
    need_ee = REFRESH or not (wpath.exists() and ppath.exists())
    if need_ee:
        ee_init()
        regions = regions_fc()

    if wpath.exists() and not REFRESH:
        wx = pd.read_csv(wpath, parse_dates=["date"])
    else:
        frames = []
        for year in WX_YEARS:
            f = CHUNKS / f"gridmet_{year}.csv"
            frames.append(cached(f, lambda: H.pull_daily_year(
                "IDAHO_EPSCOR/GRIDMET", H.GRIDMET_BANDS, year, WX_START, WX_END, 4638.3, regions)))
            print(f"  GRIDMET {year} done", end="\r")
        print()
        wx = pd.concat(frames, ignore_index=True)
        wx["date"] = pd.to_datetime(wx["date"])
        wx = wx.rename(columns={"tmmn": "tmin", "tmmx": "tmax"})
        wx["tmin"] -= 273.15                   # Kelvin -> °C
        wx["tmax"] -= 273.15
        wx["tmean"] = (wx["tmin"] + wx["tmax"]) / 2
        wx.to_csv(wpath, index=False)
        wx = pd.read_csv(wpath, parse_dates=["date"])

    if ppath.exists() and not REFRESH:
        pdsi = pd.read_csv(ppath, parse_dates=["date"])
    else:
        frames = []
        for year in WX_YEARS:
            f = CHUNKS / f"pdsi_{year}.csv"
            frames.append(cached(f, lambda: H.pull_daily_year(
                "GRIDMET/DROUGHT", ["pdsi"], year, WX_START, WX_END, 4638.3, regions)))
            print(f"  PDSI {year} done", end="\r")
        print()
        pdsi = pd.concat(frames, ignore_index=True)
        pdsi["date"] = pd.to_datetime(pdsi["date"])
        pdsi.to_csv(ppath, index=False)
        pdsi = pd.read_csv(ppath, parse_dates=["date"])

    note("weather", f"- gridMET daily: {len(wx):,} rows, {wx.date.min().date()} → **{wx.date.max().date()}** "
                    f"(last available date)")
    note("weather", f"- GRIDMET/DROUGHT pdsi: {len(pdsi):,} rows, {pdsi.date.min().date()} → {pdsi.date.max().date()}")
    days = wx.assign(year=wx.date.dt.year).groupby(["name", "year"]).date.nunique()
    short = days[(days < 275) & (days.index.get_level_values("year") < 2026)]
    note("weather", f"- Location-years (2001–2025) with < 275 days: {len(short)}")
    na = wx[["tmin", "tmax", "pr", "srad", "vpd"]].isna().sum()
    note("weather", f"- Missing daily values: {dict(na[na > 0]) or 'none'}")
    return wx, pdsi


# =============================================================== Stage 3c: static
EVT_GROUPS = ["evt_northern_hardwood", "evt_oak_hickory", "evt_cove_mesic",
              "evt_spruce_fir", "evt_pine", "evt_other"]


# Life form is in the EVT name: tree classes end in Forest / Woodland / Plantation / Treed /
# Flatwoods, plus treed pine barrens. Shrubland / Herbaceous / Barrens / Savanna classes are not forest.
_FOREST_WORDS = ("forest", "woodland", "plantation", "treed", "flatwoods", "pine barrens",
                 "pine-oak barrens")
# Explicit northern-hardwood systems (checked first, BUILD_DATASET.md 3c.4 keywords)
_NH_PHRASES = ("northern hardwood", "beech-maple", "maple-beech-birch", "maple-basswood",
               "laurentian-acadian hardwood", "acadian low-elevation hardwood")
# Remaining groups: the keyword that appears EARLIEST in the name wins (NatureServe names list the
# dominant type first, e.g. "Dry Oak-Pine" -> oak, "Northern Pine-Oak" -> pine).
_KEYWORDS = {
    "evt_spruce_fir":  [r"\bspruce", r"\bfir\b", r"\bhemlock"],
    "evt_pine":        [r"\bpine"],
    "evt_oak_hickory": [r"\boak\b", r"\bhickory"],
    "evt_cove_mesic":  [r"\bcove\b", r"\bmesic\b", r"\bmesophytic", r"\brich\b"],
}


def evt_is_forest(name):
    n = name.lower()
    return any(w in n for w in _FOREST_WORDS)


def evt_match(name):
    """(group, matched keyword) for a forest EVT class; (None, '') for non-forest."""
    if not isinstance(name, str) or not evt_is_forest(name):
        return None, ""
    n = name.lower()
    for ph in _NH_PHRASES:
        if ph in n:
            return "evt_northern_hardwood", ph
    n = n.replace("dry-mesic", "dry_mesic")      # "dry-mesic oak" is an oak type, not cove/mesic
    best = None
    for grp, pats in _KEYWORDS.items():
        for p in pats:
            m = re.search(p, n)
            if m and (best is None or m.start() < best[0]):
                best = (m.start(), grp, m.group(0))
    if best:
        return best[1], best[2]
    return "evt_other", ""


def evt_group(name):
    return evt_match(name)[0]


def stage_static():
    reset_notes("static")
    ee_init()
    regions = regions_fc()
    nlcd = H.nlcd_landcover()
    forest_mask = H.forest_mask_30m()
    dem = ee.Image("USGS/SRTMGL1_003").select("elevation")

    def per_region(fn, tag):
        frames = []
        for name in ORDER:
            f = CHUNKS / f"static_{tag}_{name.replace(' ', '_')}.csv"
            one = regions.filter(ee.Filter.eq("name", name))
            frames.append(cached(f, lambda: fn(one)))
            print(f"  {tag}: {name} done      ", end="\r")
        print()
        return pd.concat(frames, ignore_index=True)

    # --- 1. elevation of forest pixels (cell 17, unchanged) ---
    elev_reducer = (ee.Reducer.percentile([10, 50, 90], ["p10", "p50", "p90"])
                    .combine(ee.Reducer.mean(), sharedInputs=True))

    def elev_pull(coll):
        elev_df = H.fc_to_df(dem.updateMask(forest_mask).reduceRegions(
            collection=coll, reducer=elev_reducer, scale=30, tileScale=8))
        rename = {}
        for k, new in [("p10", "elev_p10_m"), ("p50", "elev_median_m"), ("p90", "elev_p90_m"),
                       ("mean", "elev_mean_m")]:
            for cand in (k, f"elevation_{k}"):
                if cand in elev_df.columns:
                    rename[cand] = new
        return elev_df.rename(columns=rename)[["name"] + list(rename.values())]

    elev_df = per_region(elev_pull, "elev")

    # exact percentiles: SRTM is integer metres, so the unweighted forest-pixel histogram is complete
    def elev_hist_pull(coll):
        fc = dem.updateMask(forest_mask).reduceRegions(
            collection=coll, reducer=ee.Reducer.frequencyHistogram().unweighted(), scale=30, tileScale=8)
        rows = []
        for f in fc.getInfo()["features"]:
            p = f["properties"]
            for v, n in p["histogram"].items():
                if v != "null":
                    rows.append(dict(name=p["name"], value=float(v), n=n))
        return pd.DataFrame(rows)

    eh = per_region(elev_hist_pull, "elevhist")
    eh.to_csv(RAW / "elevation_histograms.csv", index=False)
    exact = []
    for name, d in eh.groupby("name"):
        exact.append({"name": name, **{col: H.hist_percentile(d["value"].values, d["n"].values, q)
                                       for col, q in (("elev_p10_m", 10), ("elev_median_m", 50),
                                                      ("elev_p90_m", 90))}})
    # elev_mean_m keeps the notebook's (exact) mean reducer; the three percentiles are replaced
    elev_df = elev_df[["name", "elev_mean_m"]].merge(pd.DataFrame(exact), on="name")
    note("static", "- Elevation p10/p50/p90 are EXACT (unweighted 30 m forest-pixel histogram of integer SRTM "
                   "metres); `elev_mean_m` uses the notebook's mean reducer.")

    # --- 2. terrain + frac above 1,350 m (forest pixels) ---
    northness = ee.Terrain.aspect(dem).multiply(np.pi / 180).cos().rename("northness_mean")
    terrain = (dem.gt(1350).rename("frac_forest_above_1350m")
               .addBands(ee.Terrain.slope(dem).rename("slope_mean_deg"))
               .addBands(northness)
               .updateMask(forest_mask))

    def terrain_pull(coll):
        df = H.fc_to_df(terrain.reduceRegions(collection=coll, reducer=ee.Reducer.mean(),
                                              scale=30, tileScale=8))
        return df[["name", "frac_forest_above_1350m", "slope_mean_deg", "northness_mean"]]

    terr_df = per_region(terrain_pull, "terrain")

    # --- 3. NLCD cover fractions (cell 17, unchanged) ---
    cover = (nlcd.eq(41).rename("frac_deciduous")
             .addBands(nlcd.eq(42).rename("frac_evergreen"))
             .addBands(nlcd.eq(43).rename("frac_mixed"))
             .addBands(nlcd.eq(90).rename("frac_woody_wetland"))
             .addBands(nlcd.gte(21).And(nlcd.lte(24)).rename("frac_developed")))

    def cover_pull(coll):
        df = H.fc_to_df(cover.reduceRegions(collection=coll, reducer=ee.Reducer.mean(),
                                            scale=30, tileScale=8))
        return df[["name", "frac_deciduous", "frac_evergreen", "frac_mixed",
                   "frac_woody_wetland", "frac_developed"]]

    cover_df = per_region(cover_pull, "cover")
    cover_df["frac_hardwood"] = cover_df["frac_deciduous"] + cover_df["frac_mixed"]
    cover_df["hardwood_share_of_forest"] = cover_df["frac_hardwood"] / (
        cover_df["frac_hardwood"] + cover_df["frac_evergreen"])

    # --- 4. LANDFIRE EVT: per-region class histogram (30 m), grouped in pandas ---
    EVT_ID = "LANDFIRE/Vegetation/EVT/v1_4_0"
    evt_img = (ee.ImageCollection(EVT_ID).filter(ee.Filter.eq("system:index", "CONUS")).first())

    def evt_pull(coll):
        fc = evt_img.select("EVT").reduceRegions(collection=coll,
                                                 reducer=ee.Reducer.frequencyHistogram(),
                                                 scale=30, tileScale=8)
        rows = []
        for f in fc.getInfo()["features"]:
            p = f["properties"]
            for v, cnt in p["histogram"].items():
                rows.append(dict(name=p["name"], evt_value=int(float(v)), px=cnt))
        return pd.DataFrame(rows)

    hist = per_region(evt_pull, "evt")
    classes_path = RAW / "evt_classes_all.csv"
    if not classes_path.exists() or REFRESH:
        pd.DataFrame({"evt_value": evt_img.get("EVT_class_values").getInfo(),
                      "evt_name": evt_img.get("EVT_class_names").getInfo()}).to_csv(classes_path, index=False)
    classes = pd.read_csv(classes_path)
    hist = hist.merge(classes, on="evt_value", how="left")
    hist["is_forest"] = hist["evt_name"].fillna("").map(evt_is_forest)
    hist["group"] = hist["evt_name"].fillna("").map(evt_group)
    hist.to_csv(RAW / "evt_histogram.csv", index=False)

    cmap = (hist.groupby(["evt_value", "evt_name"], dropna=False)
                .agg(px_total=("px", "sum"), is_forest=("is_forest", "first"), group=("group", "first"))
                .reset_index().sort_values("px_total", ascending=False))
    cmap["matched_keyword"] = cmap["evt_name"].map(lambda s: evt_match(s)[1])
    cmap.to_csv(RAW / "evt_class_map.csv", index=False)

    fh = hist[hist["is_forest"]]
    evt = (fh.pivot_table(index="name", columns="group", values="px", aggfunc="sum", fill_value=0))
    evt = evt.reindex(columns=EVT_GROUPS, fill_value=0)
    evt = evt.div(evt.sum(axis=1), axis=0).reset_index()
    note("static", f"- LANDFIRE EVT asset: `{EVT_ID}` (image `CONUS`, LF 2016 Remap; the only EVT version "
                   f"in the Earth Engine catalog on {dt.date.today()}). Share of EVT forest pixels.")
    note("static", f"- EVT classes present in the 11 regions: {len(cmap)}, of which forest: "
                   f"{int(cmap.is_forest.sum())}. Full mapping in `raw/evt_class_map.csv`.")

    # --- 5. area ---
    reg = pd.read_csv(RAW / "regions.csv").rename(columns={"area_km2": "region_area_km2"})

    static = (reg[["name", "ecoregion", "region_area_km2"]]
              .merge(elev_df, on="name").merge(terr_df, on="name")
              .merge(cover_df, on="name").merge(evt, on="name", how="left"))
    static.to_csv(RAW / "static.csv", index=False)
    note("static", f"- static.csv: {len(static)} locations × {static.shape[1]} columns")
    return static


# =============================================================== Stage 4: features
def win(df, start, end):
    """Rows between two MM-DD dates (inclusive)."""
    md = df["date"].dt.strftime("%m-%d")
    return df[(md >= start) & (md <= end)]


cdd = lambda d: (20 - d["tmean"]).clip(lower=0)      # noqa: E731  cold degree-days below 20 °C
gdd5 = lambda d: (d["tmean"] - 5).clip(lower=0)      # noqa: E731  growing degree-days above 5 °C

# fixed windows (MM-DD, inclusive) -> used for both features and feature_timing.csv
FIXED = {
    "tmean_MarApr":           ("03-01", "04-30"),
    "gdd5_MarMay":            ("03-01", "05-31"),
    "tmax_JunJul":            ("06-01", "07-31"),
    "precip_JunJul":          ("06-01", "07-31"),
    "tmax_Aug":               ("08-01", "08-31"),
    "tmin_Sep":               ("09-01", "09-30"),
    "tmin_SepOct":            ("09-01", "10-31"),
    "cdd_Sep":                ("09-01", "09-30"),
    "cdd_Sep1_Oct15":         ("09-01", "10-15"),
    "cold_nights_Sep1_Oct15": ("09-01", "10-15"),
    "precip_AugSep":          ("08-01", "09-30"),
    "precip_Oct":             ("10-01", "10-31"),
    "vpd_AugSep":             ("08-01", "09-30"),
    "srad_Sep":               ("09-01", "09-30"),
    "srad_SepOct":            ("09-01", "10-31"),
    "pdsi_AugSep":            ("08-01", "09-30"),
    "first_frost_doy":        ("08-01", "11-30"),
}
# relative windows: (start offset, end offset) relative to P, inclusive
RELATIVE = {
    "tmean_pre30":       (-30, -1),
    "tmin_pre30":        (-30, -1),
    "cdd_pre30":         (-30, -1),
    "cold_nights_pre30": (-30, -1),
    "tmax_pre60to31":    (-60, -31),
    "precip_pre60":      (-60, -1),
    "vpd_pre60":         (-60, -1),
    "srad_pre30":        (-30, -1),
}


def load_raw():
    wx = pd.read_csv(RAW / "weather_daily.csv", parse_dates=["date"])
    wx["year"] = wx["date"].dt.year
    wx["doy"] = wx["date"].dt.day_of_year
    wx["md"] = wx["date"].dt.strftime("%m-%d")
    pdsi = pd.read_csv(RAW / "pdsi.csv", parse_dates=["date"])
    pdsi["year"] = pdsi["date"].dt.year
    phen = pd.read_csv(RAW / "phenology.csv")
    static = pd.read_csv(RAW / "static.csv")
    return wx, pdsi, phen, static


def peak_summary(phen):
    m = phen[(phen["product"] == "MODIS") & phen["year"].between(2001, 2025)].dropna(subset=["mid_greendown_p50"])
    rows = []
    for name, d in m.groupby("name"):
        slope = np.polyfit(d["year"], d["mid_greendown_p50"], 1)[0] if len(d) >= 3 else np.nan
        rows.append(dict(name=name, peak_doy_mean=d["mid_greendown_p50"].mean(),
                         peak_doy_sd=d["mid_greendown_p50"].std(ddof=1),
                         peak_doy_trend_d_per_decade=slope * 10, n_target_years=len(d)))
    return pd.DataFrame(rows)


def stage_features():
    reset_notes("features")
    wx, pdsi, phen, static = load_raw()
    g = ["name", "year"]

    def agg(start, end, col, how):
        return win(wx, start, end).groupby(g)[col].agg(how)

    # --- 1. fixed windows (cell 18, unchanged) + 2. new spring/summer ---
    feat = pd.DataFrame({
        # spring / early summer (new)
        "tmean_MarApr":        agg("03-01", "04-30", "tmean", "mean"),
        "gdd5_MarMay":         win(wx, "03-01", "05-31").assign(c=gdd5).groupby(g)["c"].sum(),
        "tmax_JunJul":         agg("06-01", "07-31", "tmax", "mean"),
        "precip_JunJul":       agg("06-01", "07-31", "pr", "sum"),
        # temperature
        "tmax_Aug":            agg("08-01", "08-31", "tmax", "mean"),
        "tmin_Sep":            agg("09-01", "09-30", "tmin", "mean"),
        "tmin_SepOct":         agg("09-01", "10-31", "tmin", "mean"),
        "cdd_Sep":             win(wx, "09-01", "09-30").assign(c=cdd).groupby(g)["c"].sum(),
        "cdd_Sep1_Oct15":      win(wx, "09-01", "10-15").assign(c=cdd).groupby(g)["c"].sum(),
        "cold_nights_Sep1_Oct15": win(wx, "09-01", "10-15").assign(c=lambda d: d["tmin"] < 5)
                                    .groupby(g)["c"].sum(),
        # water
        "precip_AugSep":       agg("08-01", "09-30", "pr", "sum"),
        "precip_Oct":          agg("10-01", "10-31", "pr", "sum"),
        "vpd_AugSep":          agg("08-01", "09-30", "vpd", "mean"),
        # sunlight
        "srad_Sep":            agg("09-01", "09-30", "srad", "mean"),
        "srad_SepOct":         agg("09-01", "10-31", "srad", "mean"),
    })
    feat["pdsi_AugSep"] = win(pdsi, "08-01", "09-30").groupby(g)["pdsi"].mean()

    # first fall frost: first day on/after Aug 1 with tmin <= 0 °C; none by Nov 30 -> 335, censored
    frost = (wx[(wx["md"] >= "08-01") & (wx["tmin"] <= 0)].groupby(g)["doy"].min())
    feat["first_frost_doy"] = frost
    feat = feat.reset_index()
    feat["frost_censored"] = feat["first_frost_doy"].isna().astype(int)
    feat["first_frost_doy"] = feat["first_frost_doy"].fillna(335)

    # --- 3. relative windows, anchored at P = round(peak_doy_mean), same P every year ---
    ps = peak_summary(phen)
    P = dict(zip(ps["name"], ps["peak_doy_mean"].round().astype(int)))
    wxp = wx.assign(P=wx["name"].map(P))
    wxp["rel"] = wxp["doy"] - wxp["P"]

    def rel(lo, hi):
        return wxp[(wxp["rel"] >= lo) & (wxp["rel"] <= hi)]

    relf = pd.DataFrame({
        "tmean_pre30":       rel(-30, -1).groupby(g)["tmean"].mean(),
        "tmin_pre30":        rel(-30, -1).groupby(g)["tmin"].mean(),
        "cdd_pre30":         rel(-30, -1).assign(c=cdd).groupby(g)["c"].sum(),
        "cold_nights_pre30": rel(-30, -1).assign(c=lambda d: d["tmin"] < 5).groupby(g)["c"].sum(),
        "tmax_pre60to31":    rel(-60, -31).groupby(g)["tmax"].mean(),
        "precip_pre60":      rel(-60, -1).groupby(g)["pr"].sum(),
        "vpd_pre60":         rel(-60, -1).groupby(g)["vpd"].mean(),
        "srad_pre30":        rel(-30, -1).groupby(g)["srad"].mean(),
    }).reset_index()
    ndays = rel(-60, -1).groupby(g).size()
    incomplete = ndays[ndays < 60]
    note("features", f"- Relative windows anchored at P = round(peak_doy_mean): "
                     + ", ".join(f"{n} {P[n]}" for n in ORDER))
    note("features", f"- Location-years with an incomplete 60-day relative window: {len(incomplete)}"
                     + (f" {list(incomplete.index)}" if len(incomplete) else ""))
    feat = feat.merge(relf, on=g, how="left")

    # --- 5. normals (2001-2025) and 4. day length -> location-level ---
    tr = feat[feat["year"].between(2001, 2025)]
    wtr = wx[wx["year"].between(2001, 2025)]
    norms = pd.DataFrame({
        "norm_tmean_Sep": win(wtr, "09-01", "09-30").groupby(g)["tmean"].mean().groupby("name").mean(),
        "norm_tmin_Oct": win(wtr, "10-01", "10-31").groupby(g)["tmin"].mean().groupby("name").mean(),
        "norm_precip_AugSep": tr.groupby("name")["precip_AugSep"].mean(),
        "norm_gdd5_MarMay": tr.groupby("name")["gdd5_MarMay"].mean(),
    }).reset_index()

    lat = {n: LOCATIONS[n][1] for n in ORDER}
    doys = np.arange(172, 366)
    dl = pd.DataFrame({"name": ORDER,
                       "daylen_Oct1_h": [H_day_length(lat[n], 274) for n in ORDER],
                       "doy_daylen_below_11h": [int(doys[np.argmax(H_day_length(lat[n], doys) < 11)])
                                                for n in ORDER]})

    feat.to_csv(RAW / "features_yearly.csv", index=False)
    loc_feats = ps.merge(norms, on="name").merge(dl, on="name")
    loc_feats.to_csv(RAW / "features_location.csv", index=False)
    note("features", f"- yearly features: {feat.shape[0]} rows × {feat.shape[1]} cols; "
                     f"location features: {loc_feats.shape}")
    return feat, loc_feats


def H_day_length(lat, doy):
    from ee_helpers import day_length   # pure numpy; no Earth Engine needed
    return day_length(lat, doy)


# =============================================================== Stage 5: assemble
PANEL_COLS = (
    ["location", "year", "is_forecast_row", "helene_2024", "ecoregion",
     # target (§3)
     "peak_doy", "peak_date", "peak_doy_mean", "peak_anom", "peak_doy_p25", "peak_doy_p75",
     "peak_spread_days", "peak_doy_viirs", "mid_greendown_count", "forest_px_count", "qa_retention"]
    # weather (§4b)
    + ["tmean_MarApr", "gdd5_MarMay", "tmax_JunJul", "precip_JunJul", "tmax_Aug", "tmin_Sep",
       "tmin_SepOct", "cdd_Sep", "cdd_Sep1_Oct15", "cold_nights_Sep1_Oct15", "precip_AugSep",
       "precip_Oct", "vpd_AugSep", "srad_Sep", "srad_SepOct", "pdsi_AugSep", "first_frost_doy",
       "frost_censored"]
    + list(RELATIVE)
    # phenology covariates (§4c)
    + ["mid_greenup_p50", "senescence_p50", "dormancy_p50"]
)

LOC_COLS = (["location", "lon", "lat", "state", "ecoregion", "region_area_km2", "forest_px_count",
             "elev_p10_m", "elev_median_m", "elev_p90_m", "elev_mean_m", "frac_forest_above_1350m",
             "slope_mean_deg", "northness_mean",
             "frac_deciduous", "frac_evergreen", "frac_mixed", "frac_woody_wetland", "frac_developed",
             "frac_hardwood", "hardwood_share_of_forest"]
            + EVT_GROUPS
            + ["daylen_Oct1_h", "doy_daylen_below_11h",
               "norm_tmean_Sep", "norm_tmin_Oct", "norm_precip_AugSep", "norm_gdd5_MarMay",
               "peak_doy_mean", "peak_doy_sd", "peak_doy_trend_d_per_decade", "n_target_years"])


def stage_assemble():
    reset_notes("assemble")
    wx, pdsi, phen, static = load_raw()
    feat = pd.read_csv(RAW / "features_yearly.csv")
    locf = pd.read_csv(RAW / "features_location.csv")

    modis = phen[phen["product"] == "MODIS"][[
        "name", "year", "mid_greendown_p50", "mid_greendown_p25", "mid_greendown_p75",
        "mid_greendown_count", "forest_px_count", "qa_retention",
        "mid_greenup_p50", "senescence_p50", "dormancy_p50"]].rename(columns={
            "mid_greendown_p50": "peak_doy", "mid_greendown_p25": "peak_doy_p25",
            "mid_greendown_p75": "peak_doy_p75"})
    modis["peak_spread_days"] = modis["peak_doy_p75"] - modis["peak_doy_p25"]
    viirs = (phen[phen["product"] == "VIIRS"][["name", "year", "mid_greendown_p50"]]
             .rename(columns={"mid_greendown_p50": "peak_doy_viirs"}))

    df = (feat.merge(modis, on=["name", "year"], how="left")
              .merge(viirs, on=["name", "year"], how="left")
              .merge(static[["name", "ecoregion"]], on="name", how="left")
              .merge(locf[["name", "peak_doy_mean"]], on="name", how="left"))
    df = to_location(df)
    df["helene_2024"] = [int((l, y) in HELENE) for l, y in zip(df["location"], df["year"])]
    df["is_forecast_row"] = (df["year"] == 2026).astype(int)
    df["peak_date"] = (pd.to_datetime(df["year"].astype(str) + "-01-01")
                       + pd.to_timedelta(df["peak_doy"] - 1, unit="D")).dt.strftime("%Y-%m-%d")
    df["peak_anom"] = df["peak_doy"] - df["peak_doy_mean"]
    missing = [c for c in PANEL_COLS if c not in df.columns]
    assert not missing, missing
    panel = sort_loc(df[PANEL_COLS], "year")
    panel.to_csv(OUT / "panel.csv", index=False)
    note("assemble", f"- panel.csv: {panel.shape[0]} rows × {panel.shape[1]} cols")

    # locations.csv
    fpx = (modis.groupby("name")["forest_px_count"].median().rename("forest_px_count").reset_index())
    meta = pd.DataFrame([{"name": n, "lon": LOCATIONS[n][0], "lat": LOCATIONS[n][1],
                          "state": LOCATIONS[n][2]} for n in ORDER])
    loc = (meta.merge(static, on="name").merge(fpx, on="name").merge(locf, on="name"))
    loc = sort_loc(to_location(loc))[LOC_COLS]
    loc.to_csv(OUT / "locations.csv", index=False)
    note("assemble", f"- locations.csv: {loc.shape[0]} rows × {loc.shape[1]} cols")

    # weather_daily.csv / pdsi.csv
    w = to_location(wx[["name", "date", "tmin", "tmax", "tmean", "pr", "srad", "vpd"]])
    w["date"] = w["date"].dt.strftime("%Y-%m-%d")
    sort_loc(w, "date").to_csv(OUT / "weather_daily.csv", index=False)
    p = to_location(pdsi[["name", "date", "pdsi"]])
    p["date"] = p["date"].dt.strftime("%Y-%m-%d")
    sort_loc(p, "date").to_csv(OUT / "pdsi.csv", index=False)
    note("assemble", f"- weather_daily.csv: {len(w):,} rows; pdsi.csv: {len(p):,} rows")

    # feature_timing.csv
    ref = dt.date(2001, 1, 1)   # non-leap reference year for MM-DD -> DOY
    md2doy = lambda md: (dt.date(2001, int(md[:2]), int(md[3:])) - ref).days + 1  # noqa: E731
    pmean = dict(zip(locf["name"], locf["peak_doy_mean"]))
    rows = []
    for n in ORDER:
        P = int(round(pmean[n]))
        specs = [(f, md2doy(s), md2doy(e), "fixed") for f, (s, e) in FIXED.items()]
        specs += [(f, P + lo, P + hi, "relative") for f, (lo, hi) in RELATIVE.items()]
        for f, s, e, kind in specs:
            dbp = pmean[n] - e
            if e < P:                 # window ends before the typical peak day P
                timing = "pre-peak"
            elif s < pmean[n]:
                timing = "partly"
            else:
                timing = "post"
            rows.append(dict(feature=f, location=n, window_type=kind, window_start_doy=s,
                             window_end_doy=e, peak_doy_mean=pmean[n], days_before_peak=dbp,
                             timing=timing))
    ft = pd.DataFrame(rows)
    ft.to_csv(OUT / "feature_timing.csv", index=False)
    note("assemble", f"- feature_timing.csv: {len(ft)} rows")

    shutil.copy(RAW / "evt_class_map.csv", OUT / "evt_class_map.csv")   # shared with the species teammate
    if not (OUT / "DATASET_GUIDE.md").exists():
        shutil.copy(GUIDE_SRC, OUT / "DATASET_GUIDE.md")
        note("assemble", "- DATASET_GUIDE.md copied from docs/")
    else:
        note("assemble", "- panel_v1/DATASET_GUIDE.md kept (it holds the build's deviations from docs/ version)")
    return panel, loc


# =============================================================== Stage 6: bundle
def stage_bundle():
    zpath = DATA / "panel_v1.zip"
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(OUT.iterdir()):
            if f.is_file() and f.name != "pixels.csv":
                z.write(f, f"panel_v1/{f.name}")
    print(f"wrote {zpath} ({zpath.stat().st_size / 1e6:.1f} MB)")
    px = OUT / "pixels.csv"
    if px.exists():
        with zipfile.ZipFile(DATA / "panel_v1_pixels.zip", "w", zipfile.ZIP_DEFLATED) as z:
            z.write(px, "panel_v1/pixels.csv")
        print("wrote panel_v1_pixels.zip")


# =============================================================== Stage 7: validate
def stage_validate():
    from validate_panel import validate
    ok = validate(OUT, RAW, NOTES, OLD_MODEL_DF)
    return ok


# =============================================================== Stage 3d: pixels (optional)
def stage_pixels():
    from pull_pixels import pull_pixels
    ee_init()
    pull_pixels(OUT / "pixels.csv", CHUNKS)


STAGES = {"regions": stage_regions, "pheno": stage_pheno, "weather": stage_weather,
          "static": stage_static, "features": stage_features, "assemble": stage_assemble,
          "validate": stage_validate, "bundle": stage_bundle, "pixels": stage_pixels}
DEFAULT_ORDER = ["regions", "pheno", "weather", "static", "features", "assemble", "validate", "bundle"]


def main():
    global REFRESH
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=list(STAGES), action="append",
                    help="run only this stage (repeatable); default: all except pixels")
    ap.add_argument("--refresh", action="store_true", help="re-pull from Earth Engine, ignoring caches")
    args = ap.parse_args()
    REFRESH = args.refresh
    for s in (args.stage or DEFAULT_ORDER):
        print(f"\n===== {s} =====")
        result = STAGES[s]()
        if s == "validate" and result is False:
            print("Validation FAILED - stopping before bundle. See panel_v1/BUILD_LOG.md")
            sys.exit(1)


if __name__ == "__main__":
    main()
