"""Earth Engine helpers ported from data/notebooks/DataPreprocessing.ipynb (cells 0-4, 18).

Ported verbatim. The only changes:
  * the study-region FeatureCollection, years and window dates are function arguments
    instead of notebook globals;
  * `pull_daily` can split a year into chunks (fallback for the 5,000-element getInfo limit);
  * locations are (lon, lat, state) tuples, so only coords[:2] are used for the point.
Raw outputs keep the notebook's `name` column; build_panel.py renames it to `location`.
"""
import datetime as dt

import ee
import numpy as np
import pandas as pd

ECOREGIONS_ID = "EPA/Ecoregions/2013/L3"
NLCD_ID = "USGS/NLCD_RELEASES/2021_REL/NLCD/2021"


# ---------------------------------------------------------------- cell 1: regions + forest mask
def build_study_regions(locations, buffer_m):
    """40 km buffer around each point, clipped to the EPA L3 ecoregion containing the point."""
    ecoregions = ee.FeatureCollection(ECOREGIONS_ID)
    feats = []
    for name, coords in locations.items():
        center = ee.Geometry.Point(list(coords[:2]))
        home_eco = ee.Feature(ecoregions.filterBounds(center).first())
        clipped = center.buffer(buffer_m).intersection(home_eco.geometry(), ee.ErrorMargin(100))
        feats.append(ee.Feature(clipped, {"name": name, "ecoregion": home_eco.get("us_l3name")}))
    return ee.FeatureCollection(feats)


def nlcd_landcover():
    return ee.Image(NLCD_ID).select("landcover")


def forest_mask_30m():
    nlcd = nlcd_landcover()
    return nlcd.eq(41).Or(nlcd.eq(43))       # deciduous + mixed forest (30 m)


def forest_on_grid(proj, threshold):
    # Fraction of each coarse pixel that is forest, thresholded -> 1/0 on the product's own grid
    return (forest_mask_30m().unmask(0)
            .reduceResolution(reducer=ee.Reducer.mean(), maxPixels=1024)
            .reproject(crs=proj)
            .gte(threshold))


def fc_to_df(fc):
    # FeatureCollection -> DataFrame without geometries (avoids the geemap.ee_to_df pandas bug)
    fc = fc.map(lambda f: ee.Feature(None, f.toDictionary()))
    return pd.DataFrame([f["properties"] for f in fc.getInfo()["features"]])


# ---------------------------------------------------------------- cells 2-3: phenology
MODIS = dict(
    id="MODIS/061/MCD12Q2",
    bands={"MidGreenup_1": "mid_greenup", "Senescence_1": "senescence",
           "MidGreendown_1": "mid_greendown", "Dormancy_1": "dormancy"},
    qa=lambda img: img.select("QA_Overall_1").lte(1),
)
VIIRS = dict(
    id="NASA/VIIRS/002/VNP22Q2",
    bands={"Date_Mid_Greenup_Phase_1": "mid_greenup",
           "Onset_Greenness_Decrease_1": "senescence",
           "Date_Mid_Senescence_Phase_1": "mid_greendown",
           "Onset_Greenness_Minimum_1": "dormancy"},
    qa=lambda img: (img.select("GLSP_QC_1").bitwiseAnd(7).lte(1)                     # bits 0-2: good / other quality
                    .And(img.select("GLSP_QC_1").rightShift(5).bitwiseAnd(7).eq(1))),  # bits 5-7: land
)
COMMON = ["mid_greenup", "senescence", "mid_greendown", "dormancy"]


def year_image(col, year):
    return ee.Image(col.filter(ee.Filter.calendarRange(year, year, "year")).first())


def detect_encoding(img, band, year, regions):
    # Decide how dates are stored by looking at the median raw value over the study area
    v = img.select(band).reduceRegion(ee.Reducer.median(), regions.geometry().bounds(),
                                      scale=5000, bestEffort=True).get(band).getInfo()
    if v is None:
        raise ValueError(f"No valid {band} pixels in {year} - try another year for detection")
    if v <= 800:
        enc = "day_of_year"
    elif v > 10_000:
        enc = "days_since_1970"
    else:
        enc = "days_since_2000"
    print(f"  raw median {band} ({year}) = {v:.0f}  ->  encoding: {enc}")
    return enc


def doy_offset(enc, year):
    # value to subtract from the raw date so it becomes day-of-year (Jan 1 = 1)
    if enc == "day_of_year":
        return 0
    epoch = dt.date(1970, 1, 1) if enc == "days_since_1970" else dt.date(2000, 1, 1)
    return (dt.date(year, 1, 1) - epoch).days - 1


def available_years(spec, years):
    col = ee.ImageCollection(spec["id"])
    available = set(col.aggregate_array("system:time_start")
                       .map(lambda t: ee.Date(t).get("year")).getInfo())
    return [y for y in years if y in available]


def pull_phenology_year(spec, product, year, enc, regions, threshold):
    """One year of the notebook's pull_phenology loop body (split out so each year can be cached)."""
    col = ee.ImageCollection(spec["id"])
    src = list(spec["bands"])
    reducer = (ee.Reducer.percentile([25, 50, 75])
               .combine(ee.Reducer.count(), sharedInputs=True))
    img = year_image(col, year)
    proj = img.select(src[0]).projection()
    forest = forest_on_grid(proj, threshold)
    doy = (img.select(src, COMMON)
           .subtract(doy_offset(enc, year))
           .updateMask(forest)
           .updateMask(spec["qa"](img)))
    stack = doy.addBands(ee.Image.constant(1).rename("forest_px")
                         .reproject(proj).updateMask(forest))
    fc = stack.reduceRegions(collection=regions, reducer=reducer,
                             crs=proj, scale=proj.nominalScale().getInfo(), tileScale=8)
    df = fc_to_df(fc)
    df["year"], df["product"] = year, product
    keep = (["name", "ecoregion", "year", "product", "forest_px_count"]
            + [f"{b}_{s}" for b in COMMON for s in ("p25", "p50", "p75", "count")])
    return df[[c for c in keep if c in df.columns]]


def pull_phenology_hist_year(spec, product, year, enc, regions, threshold):
    """Same masked DOY stack as pull_phenology_year, but returns the EXACT per-region pixel histogram
    (unweighted: one count per pixel whose center is in the region) instead of approximate percentiles.

    ee.Reducer.percentile is a bucketed approximation whose answer for one region depends on the other
    regions in the request (verified 2026-11-15, see BUILD_LOG). MODIS/VIIRS dates are integer days, so
    the histogram holds every pixel value and percentiles can be computed exactly in pandas.
    Output: long table name, year, product, band, value, n  (band 'forest_px' has value 1).
    """
    col = ee.ImageCollection(spec["id"])
    src = list(spec["bands"])
    img = year_image(col, year)
    proj = img.select(src[0]).projection()
    forest = forest_on_grid(proj, threshold)
    doy = (img.select(src, COMMON)
           .subtract(doy_offset(enc, year))
           .updateMask(forest)
           .updateMask(spec["qa"](img)))
    stack = doy.addBands(ee.Image.constant(1).rename("forest_px")
                         .reproject(proj).updateMask(forest))
    fc = stack.reduceRegions(collection=regions, reducer=ee.Reducer.frequencyHistogram().unweighted(),
                             crs=proj, scale=proj.nominalScale().getInfo(), tileScale=8)
    fc = fc.map(lambda f: ee.Feature(None, f.toDictionary()))
    rows = []
    for f in fc.getInfo()["features"]:
        p = f["properties"]
        for band in COMMON + ["forest_px"]:
            for v, n in (p.get(band) or {}).items():
                if v == "null":          # masked pixels are reported under a 'null' key; skip
                    continue
                rows.append(dict(name=p["name"], year=year, product=product, band=band,
                                 value=float(v), n=n))
    return pd.DataFrame(rows)


def hist_percentile(values, counts, q):
    """Exact percentile of the pixel values described by a histogram; same definition as
    numpy.percentile(method='linear') applied to the expanded pixel list."""
    order = np.argsort(values)
    v, c = np.asarray(values, float)[order], np.asarray(counts, float)[order]
    N = c.sum()
    if N == 0:
        return np.nan
    h = (N - 1) * q / 100.0
    cum = np.cumsum(c)                       # cum[i] = number of pixels with value <= v[i]
    lo_i = np.searchsorted(cum, np.floor(h) + 1)   # index of the floor(h)-th pixel (0-based)
    hi_i = np.searchsorted(cum, np.ceil(h) + 1)
    return v[lo_i] + (h - np.floor(h)) * (v[hi_i] - v[lo_i])


def phenology_setup(spec, product, years, regions):
    """The checks at the top of the notebook's pull_phenology: available years, bands, encoding."""
    col = ee.ImageCollection(spec["id"])
    years = available_years(spec, years)
    print(f"{product}: years available {years[0]}-{years[-1]} ({len(years)} years)")
    names = year_image(col, years[0]).bandNames().getInfo()
    missing = [b for b in spec["bands"] if b not in names]
    if missing:
        raise KeyError(f"{product} is missing bands {missing}. Available: {names}")
    src = list(spec["bands"])
    enc = detect_encoding(year_image(col, years[-1]), src[2], years[-1], regions)
    return years, enc


# ---------------------------------------------------------------- cell 4: daily weather
GRIDMET_BANDS = ["tmmn", "tmmx", "pr", "srad", "vpd"]


def pull_daily_range(col_id, bands, start, end, scale, regions):
    """Regional means for every image in [start, end). Body of the notebook's pull_daily loop."""
    col = ee.ImageCollection(col_id).select(bands).filterDate(start, end)

    def per_image(img):
        stats = img.reduceRegions(collection=regions, reducer=ee.Reducer.mean(), scale=scale)
        return stats.map(lambda f: ee.Feature(None, f.toDictionary())
                         .set("date", img.date().format("YYYY-MM-dd")))

    fc = col.map(per_image).flatten()
    df = pd.DataFrame([f["properties"] for f in fc.getInfo()["features"]])
    if df.empty:
        return df
    if len(bands) == 1 and "mean" in df.columns:
        df = df.rename(columns={"mean": bands[0]})
    df = df[["name", "date"] + [b for b in bands if b in df.columns]]
    return df


def pull_daily_year(col_id, bands, year, wx_start, wx_end, scale, regions):
    """One year, Mar 1 -> Dec 1 (exclusive). Falls back to two halves if getInfo is too large."""
    try:
        return pull_daily_range(col_id, bands, f"{year}-{wx_start}", f"{year}-{wx_end}",
                                scale, regions)
    except ee.EEException as e:
        print(f"    {year}: single request failed ({str(e)[:80]}...) -> splitting into halves")
        mid = f"{year}-07-16"
        a = pull_daily_range(col_id, bands, f"{year}-{wx_start}", mid, scale, regions)
        b = pull_daily_range(col_id, bands, mid, f"{year}-{wx_end}", scale, regions)
        return pd.concat([a, b], ignore_index=True)


# ---------------------------------------------------------------- cell 18: day length
def day_length(lat, doy):
    """Day length (hours) from latitude + day of year - Forsythe et al. (1995) CBM model."""
    theta = 0.2163108 + 2 * np.arctan(0.9671396 * np.tan(0.00860 * (doy - 186)))
    phi = np.arcsin(0.39795 * np.cos(theta))
    p = 0.8333
    x = (np.sin(np.radians(p)) + np.sin(np.radians(lat)) * np.sin(phi)) / (
        np.cos(np.radians(lat)) * np.cos(phi))
    return 24 - (24 / np.pi) * np.arccos(np.clip(x, -1, 1))
