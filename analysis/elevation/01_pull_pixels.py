"""Pull pixel-level MODIS peak-color dates + elevation for the 4 study regions.

Run once (slow, hits Earth Engine). Writes pixels.csv; all analysis reads only that file.
"""
import datetime as dt
from pathlib import Path

import ee
import pandas as pd

ee.Initialize(project="fall-foliage-project")

OUT = Path(__file__).resolve().parent / "pixels.csv"
N_PIXELS = 4000          # max forest pixels kept per region (getInfo limit is 5,000)
N_DRAW = 25000           # raw draw; most points land on non-forest and get filtered out
SEED = 42

locations = {   # [lon, lat]
    "Raleigh":    [-78.69649427302112, 35.7664191909458],
    "Boone":      [-81.68516978460896, 36.21834246106182],
    "Asheville":  [-82.55831963941482, 35.633546529345516],
    "Shenandoah": [-78.40685152551191, 38.499092552797826],
}

mcd = ee.ImageCollection("MODIS/061/MCD12Q2")
modis_proj = mcd.first().select("MidGreendown_1").projection()
modis_scale = modis_proj.nominalScale().getInfo()
years = [int(s[:4]) for s in mcd.aggregate_array("system:index").getInfo()]
print(f"MODIS years {years[0]}-{years[-1]}, scale {modis_scale:.1f} m")


def to_modis(img):
    """Area-average a fine-resolution image onto the MODIS grid."""
    return (img.reduceResolution(ee.Reducer.mean(), maxPixels=1024)
               .reproject(crs=modis_proj))


# --- static layers on the MODIS grid ---
nlcd = ee.Image("USGS/NLCD_RELEASES/2021_REL/NLCD/2021").select("landcover")
forest = nlcd.eq(41).Or(nlcd.eq(43))
forest_frac = to_modis(forest.unmask(0).toFloat()).rename("forest_frac")

dem = ee.Image("USGS/SRTMGL1_003").select("elevation")
elev = to_modis(dem.toFloat()).rename("elev")
aspect = ee.Terrain.aspect(dem)
northness = to_modis(aspect.multiply(3.141592653589793 / 180).cos()).rename("northness")

# --- peak-color (MidGreendown) DOY per year ---
gd_bands = []
for y in years:
    img = mcd.filter(ee.Filter.calendarRange(y, y, "year")).first()
    offset = (dt.date(y, 1, 1) - dt.date(1970, 1, 1)).days
    gd = (img.select("MidGreendown_1").subtract(offset).add(1)
             .updateMask(img.select("QA_Overall_1").lte(1))
             .rename(f"gd_{y}"))
    gd_bands.append(gd)

stack = (ee.Image.cat([elev, northness, forest_frac, ee.Image.pixelLonLat()] + gd_bands)
         .updateMask(forest_frac.gte(0.5)))

ecoregions = ee.FeatureCollection("EPA/Ecoregions/2013/L3")

frames = []
for name, (lon, lat) in locations.items():
    pt = ee.Geometry.Point([lon, lat])
    eco = ecoregions.filterBounds(pt).first()
    eco_name = eco.get("us_l3name").getInfo()
    geom = pt.buffer(40_000).intersection(eco.geometry(), ee.ErrorMargin(100))

    fc = stack.sample(region=geom, projection=modis_proj, scale=modis_scale,
                      numPixels=N_DRAW, seed=SEED, geometries=False,
                      dropNulls=False, tileScale=4)
    # dropNulls=False also returns fully masked (non-forest) pixels; keep forest pixels only
    # then keep a random subset of at most N_PIXELS so getInfo stays under 5,000
    fc = (fc.filter(ee.Filter.notNull(["elev"]))
            .randomColumn("rnd", SEED).sort("rnd").limit(N_PIXELS))
    feats = fc.getInfo()["features"]
    df = pd.DataFrame([f["properties"] for f in feats]).drop(columns="rnd")
    df.insert(0, "region", name)
    df.insert(1, "pixel_id", [f"{name[:3].upper()}_{i:05d}" for i in range(len(df))])
    print(f"{name} ({eco_name}): {len(df)} forest pixels sampled")
    frames.append(df)

wide = pd.concat(frames, ignore_index=True)
gd_cols = [f"gd_{y}" for y in years]
for c in gd_cols:
    if c not in wide:
        wide[c] = pd.NA
long = wide.melt(id_vars=["region", "pixel_id", "longitude", "latitude", "elev",
                          "northness", "forest_frac"],
                 value_vars=gd_cols, var_name="year", value_name="peak_doy")
long["year"] = long.year.str[3:].astype(int)
long = (long.rename(columns={"longitude": "lon", "latitude": "lat", "elev": "elev_m"})
            .dropna(subset=["peak_doy"])
            [["region", "pixel_id", "lon", "lat", "elev_m", "northness", "forest_frac",
              "year", "peak_doy"]]
            .sort_values(["region", "pixel_id", "year"]))
long.to_csv(OUT, index=False)

print(f"\nwrote {OUT.name}: {len(long):,} pixel-year rows")
px = long.drop_duplicates("pixel_id")
summ = px.groupby("region").elev_m.describe(percentiles=[.05, .5, .95])[
    ["count", "min", "5%", "50%", "95%", "max"]].round(0)
summ["pixel_years"] = long.groupby("region").size()
print(summ)
q = long.peak_doy.quantile([.01, .05, .5, .95, .99])
print("\npeak_doy quantiles:\n", q.round(1))
share = long.peak_doy.between(260, 330).mean()
print(f"share of peak_doy in 260-330: {share:.1%}")
if share < 0.9:
    raise SystemExit("DOY conversion looks wrong: <90% of dates in 260-330")
