# Build instructions: Appalachian peak-foliage panel (v1)

**Audience:** Arin's Claude (coding agent) working in the `IMMC/` repo on Arin's Mac.
**Goal:** produce the standardized dataset described in `DATASET_GUIDE.md` (read it first; it is the spec for every output file and column) for **11 locations**, by **extending the existing Earth Engine pipeline**, not rewriting it.
**Ground rules:**
- **Check what exists before writing anything** (Stage 0). Reuse the functions in `data/notebooks/DataPreprocessing.ipynb` wherever possible. That notebook is the source of truth for how each Earth Engine dataset was processed yesterday.
- **Do not modify or overwrite** anything in `data/data/` or `analysis/`. The 4-region analysis scripts read those files. All new output goes to `data/panel_v1/`.
- Stop and report to Arin if a check in Stage 7 fails. Don't silently patch the data.

---

## Stage 0 — Reconnaissance (do this first, report a short summary)

1. List the repo. Expected layout (as of 2026-11-15):
   ```
   data/regions.py                    # stale: old 5 points, geemap only
   data/cleaning.py                   # stub
   data/notebooks/DataPreprocessing.ipynb   # THE pipeline (cells 0–19)
   data/notebooks/maps.ipynb
   data/data/{model_df,phenology_regions,weather_daily,weather_features,pdsi,static_features}.csv
   analysis/common.py, 01–05_*.py, FINDINGS.md
   analysis/elevation/{01_pull_pixels.py, 02_elevation_analysis.py, elevation_adjust.py, ELEVATION_FINDINGS.md, ...}
   ```
   If anything has changed (e.g. a `panel_v1/` already exists, or new scripts are present), read it and adapt. Don't duplicate work.
2. Read the notebook cell by cell. Map each piece to the stage below where it will be reused:
   | Notebook cell | Contents | Reused in |
   |---|---|---|
   | 0 | config: `PROJECT="fall-foliage-project"`, `UNIV_BUFFER=40_000`, `FOREST_THRESHOLD=0.5`, year ranges | Stage 1 |
   | 1 | `study_regions` (buffer ∩ EPA L3 ecoregion), NLCD forest mask, `forest_on_grid`, `fc_to_df` | Stage 2 |
   | 2–3 | `MODIS`/`VIIRS` specs, `detect_encoding`, `doy_offset`, `pull_phenology` | Stage 3a |
   | 4 | `pull_daily` for gridMET and GRIDMET/DROUGHT | Stage 3b |
   | 14, 18 | weather window features, frost, day length (`day_length`) | Stage 4 |
   | 17 | static: forest-pixel elevation percentiles (SRTM), NLCD cover fractions | Stage 3c |
   | 19 | assembly into `model_df` (MODIS target, VIIRS overlay, Helene flag) | Stage 5 |
3. Read `analysis/FINDINGS.md`, `analysis/elevation/ELEVATION_FINDINGS.md` and `analysis/elevation/01_pull_pixels.py`.
4. Confirm the Python env (`fall-foliage` conda env per FINDINGS) has `ee`, `geemap`, `pandas`, `numpy`. Confirm `ee.Initialize(project="fall-foliage-project")` works. If authentication fails, ask Arin to run `earthengine authenticate`. Don't try to work around it.
5. Record the existing values for **Boone, Asheville, Shenandoah** from `data/data/model_df.csv`. These three locations carry over unchanged, so they are the regression test in Stage 7.

## Stage 1 — Location config

Create `data/locations.py` (leave `regions.py` alone) with one dict that every later step imports:

```python
LOCATIONS = {  # name: (lon, lat, state)
    "Montpelier":    (-72.576, 44.260, "VT"),
    "Old Forge":     (-74.972, 43.710, "NY"),
    "Binghamton":    (-75.918, 42.099, "NY"),
    "Scranton":      (-75.662, 41.409, "PA"),
    "State College": (-77.860, 40.793, "PA"),
    "Frostburg":     (-78.928, 39.658, "MD"),
    "Valley Head":   (-80.037, 38.546, "WV"),
    "Shenandoah":    (-78.40685152551191, 38.499092552797826, "VA"),  # keep exact old point
    "Roanoke":       (-79.941, 37.271, "VA"),
    "Boone":         (-81.68516978460896, 36.21834246106182, "NC"),   # keep exact old point
    "Asheville":     (-82.55831963941482, 35.633546529345516, "NC"),  # keep exact old point
}
ORDER = list(LOCATIONS)     # north → south; use this order in every table and figure
BUFFER_M = 40_000
FOREST_THRESHOLD = 0.5
PHENO_YEARS = range(2001, 2026)
WX_YEARS = range(2001, 2027)
WX_START, WX_END = "03-01", "12-01"   # was 07-01; extended so spring/summer features exist
HELENE = {("Asheville", 2024), ("Boone", 2024)}
```
- The 8 new coordinates are town centers taken from memory. **Verify each one** (e.g. against a geocoder or Wikipedia) and correct it if it is off by more than ~0.02°. Do **not** change the three carried-over points.
- Print the EPA L3 ecoregion each point falls in, and each clipped region's area in km². If a clipped region is under ~1,500 km² (the ecoregion boundary cuts the circle a lot), flag it to Arin.
- Check pairwise overlap of the clipped regions. Binghamton and Scranton are about 80 km apart, so their circles nearly touch. Report any overlapping area.

## Stage 2 — Build script skeleton

Create `data/build_panel.py` (a script, not a notebook, so it can be rerun). Structure:
- Import config from `locations.py`. Port the notebook helpers (`fc_to_df`, `forest_on_grid`, `year_image`, `detect_encoding`, `doy_offset`, `pull_phenology`, `pull_daily`, `day_length`) **verbatim** into `data/ee_helpers.py`, changing only what is needed to take the location dict as an argument.
- **Cache every raw Earth Engine pull** to `data/panel_v1/raw/*.csv`, and skip the pull if the cache file exists (add a `--refresh` flag). Earth Engine pulls are slow, and you will rerun the feature code many times.
- Stages 3–6 are separate functions. `main()` runs them in order, with CLI flags `--stage` and `--refresh`.

## Stage 3 — Earth Engine pulls

### 3a. Phenology → `raw/phenology.csv`
Reuse `pull_phenology(MODIS, ...)` and `pull_phenology(VIIRS, ...)` unchanged, including the QA masks, the encoding detection and the forest mask. Output one row per (location, year, product) with p25/p50/p75/count for the 4 bands, plus `forest_px_count`.
- Confirm the most recent available MODIS MCD12Q2 year (it was 2025 in yesterday's pull). If 2025 is missing for any location, report it.
- Sanity check from notebook cell 3: flag any `mid_greendown_p50` outside DOY 230–350.

### 3b. Daily weather → `raw/weather_daily.csv`, `raw/pdsi.csv`
Reuse `pull_daily("IDAHO_EPSCOR/GRIDMET", ["tmmn","tmmx","pr","srad","vpd"], ...)` with the **new window Mar 1 → Dec 1**. Convert K → °C and add `tmean` exactly as in cell 4. Then do the same for `GRIDMET/DROUGHT` `pdsi`.
- Size check: 11 locations × 275 days ≈ 3,000 features per year, which is under the 5,000-element `getInfo` limit. If it errors anyway, split each year into two halves.
- Weather is averaged over the clipped region. The notebook does **not** apply the forest mask to weather; keep it that way for consistency with the old data.
- 2026 will stop about 2 days before today. Record the last available date in the build log.

### 3c. Static → `raw/static.csv`
Run everything over the clipped regions:
1. **Elevation (forest pixels):** cell 17 code, unchanged (SRTM `USGS/SRTMGL1_003`, p10/p50/p90/mean, scale 30, mask = deciduous + mixed). Add `frac_forest_above_1350m` = mean of `(elev > 1350)` over forest pixels.
2. **Terrain (forest pixels):** `ee.Terrain.slope(dem)` → `slope_mean_deg`. `ee.Terrain.aspect(dem)` → `northness = cos(radians(aspect))` → `northness_mean`.
3. **NLCD cover fractions:** cell 17 code, unchanged, including `frac_hardwood` and `hardwood_share_of_forest`.
4. **LANDFIRE EVT (species / forest type):** search the Earth Engine catalog for the newest CONUS LANDFIRE Existing Vegetation Type asset (expected under `LANDFIRE/Vegetation/EVT/...`; confirm the exact ID and version). Get its class table (class value → EVT name; it is usually in the image properties, or else the LANDFIRE CSV). Group classes by name keywords into:
   `evt_northern_hardwood` (northern hardwood, maple–beech–birch, Laurentian-Acadian hardwood), `evt_oak_hickory` (oak, oak–hickory, dry-mesic oak), `evt_cove_mesic` (cove, mesic, rich), `evt_spruce_fir` (spruce, fir, hemlock), `evt_pine` (pine), `evt_other` (all remaining forest classes).
   Compute each group's share of **forest** area, at 30 m. Save the full keyword → class mapping to `raw/evt_class_map.csv` so the species teammate can audit and regroup it. If no LANDFIRE EVT asset exists in the catalog, leave the `evt_*` columns empty and say so in the build log. Don't substitute another dataset without asking Arin.
5. **Area:** `region_area_km2` from the clipped geometry.

### 3d. (Optional, run last) Pixel table → `pixels.csv`
Generalize `analysis/elevation/01_pull_pixels.py` to read `LOCATIONS` (same sampling: up to 4,000 random forest pixels per location, MODIS grid, SRTM mean elevation, northness, forest fraction). This runs long. Do it only after Stages 4–7 pass, and ask Arin first if the estimated run time is over ~30 min.

## Stage 4 — Feature engineering (pandas, from the raw caches)

1. **Fixed-window features:** reproduce cell 18 exactly (`tmax_Aug`, `tmin_Sep`, `tmin_SepOct`, `cdd_Sep`, `cdd_Sep1_Oct15`, `cold_nights_Sep1_Oct15`, `precip_AugSep`, `precip_Oct`, `vpd_AugSep`, `srad_Sep`, `srad_SepOct`, `pdsi_AugSep`, `first_frost_doy`, `frost_censored`). Keep the same window functions and the same frost rule (first tmin ≤ 0 °C on or after Aug 1; if none, 335 with censored = 1).
2. **New spring/summer features:** `tmean_MarApr`, `gdd5_MarMay` (Σ max(0, tmean − 5)), `tmax_JunJul`, `precip_JunJul`.
3. **Relative windows:** first compute `peak_doy_mean` per location (MODIS `mid_greendown_p50`, 2001–2025). Let P = round(peak_doy_mean). Using DOY (so leap years shift by at most one day, which is acceptable), compute `tmean_pre30`, `tmin_pre30`, `cdd_pre30`, `cold_nights_pre30`, `tmax_pre60to31`, `precip_pre60`, `vpd_pre60`, `srad_pre30` with the windows defined in DATASET_GUIDE §4b. Use the **same P for every year** of a location, including 2026.
4. **Day length:** `daylen_Oct1_h` and `doy_daylen_below_11h` via `day_length()` from cell 18.
5. **Normals** (2001–2025 means per location): `norm_tmean_Sep`, `norm_tmin_Oct`, `norm_precip_AugSep`, `norm_gdd5_MarMay`.
6. **Peak summary** per location: `peak_doy_mean`, `peak_doy_sd`, `peak_doy_trend_d_per_decade` (OLS slope × 10), `n_target_years`.

## Stage 5 — Assemble outputs into `data/panel_v1/`

- `panel.csv`: start from the weather feature table (so that 2026 rows exist, as in cell 19). Left-join the MODIS target, renaming `mid_greendown_p50 → peak_doy`, `_p25/_p75 → peak_doy_p25/_p75`. Add `peak_spread_days`, `mid_greendown_count`, `forest_px_count`, `qa_retention`, `mid_greenup_p50`, `senescence_p50` and `dormancy_p50`; then VIIRS `peak_doy_viirs`, `ecoregion`, `helene_2024`, `is_forecast_row`, `peak_date`, `peak_doy_mean` and `peak_anom`. Column order should follow DATASET_GUIDE §4. Sort rows by `ORDER`, then by year.
- `locations.csv`: one row per location with every column in DATASET_GUIDE §4d.
- `weather_daily.csv`, `pdsi.csv`: cleaned copies of the raw files (column `location` instead of `name`, ISO dates).
- `feature_timing.csv`: for each pair of (fixed-window or relative feature, location), record `window_end_doy` and `peak_doy_mean`, and compute `days_before_peak = peak_doy_mean − window_end_doy`. Set `timing` to `pre-peak` if days_before_peak ≥ 1, to `partly` if the window starts before peak but ends on or after it, and to `post` otherwise.
- Copy `DATASET_GUIDE.md` into `panel_v1/`. If any column ended up different from the guide (renamed, dropped, or empty like `evt_*`), **update the guide** so it matches the data exactly.
- Use the column name `location` everywhere. The old files use `name`; don't carry that over.

## Stage 6 — Shareable bundle

Zip `panel_v1/` without `raw/` and without `pixels.csv` into `data/panel_v1.zip`. That zip is what Arin sends to teammates. Make `pixels.csv` a separate zip if it exists.

## Stage 7 — Validation checks (all must pass; write results to `panel_v1/BUILD_LOG.md`)

1. **Regression test:** for Boone, Asheville and Shenandoah, 2001–2025, compare the new `peak_doy`, `peak_doy_viirs`, `tmax_Aug`, `precip_AugSep`, `cdd_Sep`, `pdsi_AugSep`, `first_frost_doy`, `elev_median_m` and `frac_hardwood` against `data/data/model_df.csv`. Expect identical values (max abs diff < 1e-6). The code path and regions are the same, so any difference means something changed. **Investigate and report it. Don't proceed.**
2. **Shape:** `panel.csv` has 11 × 26 = 286 rows and no duplicate (location, year). `peak_doy` is non-null for every 2001–2025 row (report any gaps) and null for every 2026 row.
3. **Ranges:** `peak_doy` is between 230 and 350. `tmax_Aug` is between 15 and 35 °C. `precip_AugSep` ≥ 0. `pdsi` is between −8 and +8. Fractions are between 0 and 1. `elev_median_m` is between 50 and 1,500.
4. **Geography sense-check:** report `peak_doy_mean` by location, north → south. Expect northern and higher locations (Montpelier, Old Forge) to be earlier on average. If the ordering is clearly inverted, flag it (it could be a forest-mask or encoding problem), but don't "fix" it.
5. **Missingness table:** non-null count per column. Flag any column under 95% (except `peak_doy_viirs` before 2013, and `evt_*` if no asset was found).
6. **MODIS vs VIIRS:** per location bias, RMSE and r for 2013–2025. Compare with FINDINGS §1 (about +8 d bias).
7. **Pixel counts:** `forest_px_count` per location. Flag any location under ~2,000 MODIS forest pixels. Raleigh had 2,905 and was the weakest-signal region.
8. **Quick baseline:** leave-one-year-out climatology MAE/RMSE per location and overall. This is the number the delta model has to beat. Put it in the log.

## Stage 8 — Report back to Arin

A short message covering:
- where the files are;
- the regression-test result;
- the final coordinates (and any you corrected);
- any locations with weak data;
- the climatology baseline table;
- any deviations from DATASET_GUIDE.

Don't start modeling. That follows `plan.md`.
