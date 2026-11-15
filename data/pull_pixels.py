"""Stage 3d (optional): pixel-level MODIS peak dates + elevation for all LOCATIONS -> panel_v1/pixels.csv.

Generalized from analysis/elevation/01_pull_pixels.py (same sampling: up to 4,000 random forest
pixels per location on the MODIS grid, SRTM mean elevation, northness, forest fraction, seed 42).
Only change: reads locations.LOCATIONS, caches each location, and uses `location` instead of `region`.
"""
import datetime as dt
from pathlib import Path

import ee
import pandas as pd

from locations import LOCATIONS, ORDER, BUFFER_M, FOREST_THRESHOLD

N_PIXELS = 4000          # max forest pixels kept per location (getInfo limit is 5,000)
N_DRAW = 25000           # raw draw; most points land on non-forest and get filtered out
SEED = 42


def pull_pixels(out_path, chunk_dir):
    mcd = ee.ImageCollection("MODIS/061/MCD12Q2")
    modis_proj = mcd.first().select("MidGreendown_1").projection()
    modis_scale = modis_proj.nominalScale().getInfo()
    years = [int(s[:4]) for s in mcd.aggregate_array("system:index").getInfo()]
    print(f"MODIS years {years[0]}-{years[-1]}, scale {modis_scale:.1f} m")

    def to_modis(img):
        return (img.reduceResolution(ee.Reducer.mean(), maxPixels=1024)
                   .reproject(crs=modis_proj))

    nlcd = ee.Image("USGS/NLCD_RELEASES/2021_REL/NLCD/2021").select("landcover")
    forest = nlcd.eq(41).Or(nlcd.eq(43))
    forest_frac = to_modis(forest.unmask(0).toFloat()).rename("forest_frac")
    dem = ee.Image("USGS/SRTMGL1_003").select("elevation")
    elev = to_modis(dem.toFloat()).rename("elev")
    northness = to_modis(ee.Terrain.aspect(dem).multiply(3.141592653589793 / 180).cos()).rename("northness")

    gd_bands = []
    for y in years:
        img = mcd.filter(ee.Filter.calendarRange(y, y, "year")).first()
        offset = (dt.date(y, 1, 1) - dt.date(1970, 1, 1)).days
        gd_bands.append(img.select("MidGreendown_1").subtract(offset).add(1)
                           .updateMask(img.select("QA_Overall_1").lte(1)).rename(f"gd_{y}"))

    stack = (ee.Image.cat([elev, northness, forest_frac, ee.Image.pixelLonLat()] + gd_bands)
             .updateMask(forest_frac.gte(FOREST_THRESHOLD)))
    ecoregions = ee.FeatureCollection("EPA/Ecoregions/2013/L3")

    frames = []
    for name in ORDER:
        f = Path(chunk_dir) / f"pixels_{name.replace(' ', '_')}.csv"
        if f.exists():
            frames.append(pd.read_csv(f))
            continue
        lon, lat = LOCATIONS[name][:2]
        pt = ee.Geometry.Point([lon, lat])
        eco = ecoregions.filterBounds(pt).first()
        geom = pt.buffer(BUFFER_M).intersection(eco.geometry(), ee.ErrorMargin(100))
        fc = stack.sample(region=geom, projection=modis_proj, scale=modis_scale,
                          numPixels=N_DRAW, seed=SEED, geometries=False,
                          dropNulls=False, tileScale=4)
        fc = (fc.filter(ee.Filter.notNull(["elev"]))
                .randomColumn("rnd", SEED).sort("rnd").limit(N_PIXELS))
        df = pd.DataFrame([x["properties"] for x in fc.getInfo()["features"]]).drop(columns="rnd")
        df.insert(0, "location", name)
        df.insert(1, "pixel_id", [f"{name.replace(' ', '')[:3].upper()}_{i:05d}" for i in range(len(df))])
        df.to_csv(f, index=False)
        print(f"{name}: {len(df)} forest pixels sampled")
        frames.append(df)

    wide = pd.concat(frames, ignore_index=True)
    gd_cols = [f"gd_{y}" for y in years]
    for c in gd_cols:
        if c not in wide:
            wide[c] = pd.NA
    long = wide.melt(id_vars=["location", "pixel_id", "longitude", "latitude", "elev",
                              "northness", "forest_frac"],
                     value_vars=gd_cols, var_name="year", value_name="peak_doy")
    long["year"] = long.year.str[3:].astype(int)
    long = (long.rename(columns={"longitude": "lon", "latitude": "lat", "elev": "elev_m"})
                .dropna(subset=["peak_doy"])
                [["location", "pixel_id", "lon", "lat", "elev_m", "northness", "forest_frac",
                  "year", "peak_doy"]])
    long["_o"] = long["location"].map({n: i for i, n in enumerate(ORDER)})
    long = long.sort_values(["_o", "pixel_id", "year"]).drop(columns="_o")
    long.to_csv(out_path, index=False)
    print(f"wrote {out_path}: {len(long):,} pixel-year rows")
    share = long.peak_doy.between(260, 330).mean()
    print(f"share of peak_doy in 260-330: {share:.1%}")
    return long
