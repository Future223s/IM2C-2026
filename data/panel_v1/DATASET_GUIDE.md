# Appalachian Peak Foliage Dataset — User Guide (v1)

**Who this is for:** each teammate (and their AI assistant). Load this file as context before working with the data.
**Status:** built 2026-11-15 by `data/build_panel.py`; all Stage 7 checks pass (see `BUILD_LOG.md` in the build folder). This copy has been updated to match the data exactly; differences from the planning version are marked *(build note)*.
**Project:** IMMC "Peak Peeping": predict when peak fall foliage happens this year (2026) at 11 Appalachian locations, and say how accurate the model would be 10 years out.
**Model form the team has agreed on:** `predicted peak = location mean + delta`. The *mean* is a location's typical peak date. The *delta* is how far this year lands from it. Your job is to work out which factors explain the delta, or the differences between locations, for your piece of the problem.

| Teammate | Focus | Main files and columns |
|---|---|---|
| Environmental (Arin) | Why locations differ: elevation, terrain, latitude/day length, climate normals | `locations.csv` (all of it), `panel.csv` → `peak_doy_mean`, `peak_anom` |
| Temporal | Year-to-year weather and phenology that push a given year early or late | `panel.csv` weather and phenology columns, `weather_daily.csv` |
| Species | How tree species composition shifts timing, and how it changes sensitivity to weather | `locations.csv` → `frac_*`, `evt_*` columns; `panel.csv` → `peak_anom` |

---

## 1. Files

All files are in `IMMC/data/panel_v1/`. Every table uses the same location names and spellings.

| File | Grain | Rows (approx.) | What it is |
|---|---|---|---|
| `panel.csv` | location × year, 2001–2026 | 11 × 26 = 286 | **The main table.** Target, weather features and phenology for each location-year. |
| `locations.csv` | location | 11 | Static features: coordinates, elevation, terrain, land cover, species mix, day length, climate normals, typical peak. |
| `weather_daily.csv` | location × day, Mar 1–Nov 30 | 78,441 | Daily gridMET means over each location's **whole clipped region** (weather is *not* forest-masked, same as the earlier analysis). Columns: `location, date, tmin, tmax, tmean, pr, srad, vpd`. Use it to build your own windows. |
| `pdsi.csv` | location × 5-day step | 15,686 | Palmer Drought Severity Index (gridMET drought, 5-day cadence), Mar–Nov. Columns: `location, date, pdsi`. |
| `feature_timing.csv` | feature × location | 25 × 11 = 275 | Whether each weather feature's window ends before that location's typical peak (see §5). Phenology covariates are not listed; their timing is given in §4c. |
| `pixels.csv` *(separate zip, large)* | forest pixel × year | ~1.1 M | Pixel-level peak date for up to 4,000 random forest pixels per location. Columns: `location, pixel_id, lon, lat, elev_m, northness, forest_frac, year, peak_doy` (MODIS grid, SRTM mean elevation per pixel). For within-location elevation work. |
| `evt_class_map.csv` | LANDFIRE EVT class | 113 | Every vegetation class present in the 11 regions, its pixel count, and the `evt_*` group it was assigned to. |
| `HOW_TO_USE_THIS_DATASET.md` | — | — | Instructions for each teammate's analysis (read with this file). |
| `BUILD_LOG.md` | — | — | Build notes and validation results. |
| `DATASET_GUIDE.md` | — | — | This file. |

**Keys:** `location` (string) and `year` (int). Join `locations.csv` onto `panel.csv` with `location`.

## 2. Locations

| location | state | lon | lat | Notes |
|---|---|---|---|---|
| Montpelier | VT | −72.576 | 44.260 | Northernmost. Northern hardwoods (maple, beech, birch). |
| Old Forge | NY | −74.972 | 43.710 | Adirondacks. Lots of lake and wetland area, which the forest mask removes. |
| Binghamton | NY | −75.918 | 42.099 | Appalachian Plateau. **Buffer nearly touches Scranton's** (~80 km apart). |
| Scranton | PA | −75.662 | 41.409 | *(build note)* Ridge and Valley ecoregion cuts the circle hard: only **1,064 km²** (3,215 forest pixels), the smallest region. |
| State College | PA | −77.860 | 40.793 | Ridge and Valley. |
| Frostburg | MD | −78.928 | 39.658 | Allegheny Plateau. |
| Valley Head | WV | −80.037 | 38.546 | High Alleghenies. *(build note)* The point sits in the **Ridge and Valley** ecoregion (Tygart Valley), so the region is 1,816 km² of valley-and-ridge, not the Allegheny plateau. |
| Shenandoah | VA | −78.407 | 38.499 | Shenandoah NP point carried over from the earlier analysis. Blue Ridge region is only 1,215 km² (narrow ridge). |
| Roanoke | VA | −79.941 | 37.271 | Mix of city and Blue Ridge. |
| Boone | NC | −81.685 | 36.218 | Carried over. |
| Asheville | NC | −82.558 | 35.634 | Carried over. Hit by Hurricane Helene in 2024. |

**Region definition:** each location is a **40 km radius circle clipped to the EPA Level III ecoregion that contains the center point**. All satellite and weather values are summarized over the **deciduous + mixed forest** pixels inside that region (NLCD 2021 classes 41 and 43; a coarse pixel counts as forest if at least 50% of it is forest). So a "location" means the forested landscape around the town, not the town itself.

Raleigh and Charleston from the earlier 4-region analysis are **not** included.

The 8 new coordinates were checked against OpenStreetMap (all within 0.01°). No two regions overlap (Binghamton–Scranton are 79 km apart; their clipped regions do not touch).

## 3. The target: `peak_doy`

- **What it is:** the median, across forest pixels, of **MODIS MCD12Q2 v061 `MidGreendown_1`**. *(build note)* All percentiles are **exact** pixel percentiles (numpy 'linear' definition) over the forest pixels whose **center** lies in the region; pixel dates are whole days, so every median is a whole or half day. This is the day of year when the vegetation index (EVI) is halfway through its autumn decline. It is our **proxy for peak color**.
- **Units:** day of year (Jan 1 = 1). DOY 274 = Oct 1 and DOY 305 = Nov 1 in non-leap years; subtract 1 in leap years. `peak_date` gives the same value as a calendar date.
- **Years with a target:** 2001–2025. **2026 rows have no target** (`is_forecast_row = 1`), because the 2026 MODIS product has not been released. 2026 is the year we are predicting.
- **Important:** mid-greendown is a satellite signal, not "peak color" as a person on the ground would report it. In the earlier analysis it came about **1–3 weeks before** news-reported peak color. The plan calibrates this offset separately (see `plan.md`). Model `peak_doy` as it is.
- **Related target columns:**
  - `peak_doy_p25`, `peak_doy_p75`: 25th and 75th percentile across pixels. `peak_spread_days` = p75 − p25, a measure of how uneven the location's timing is.
  - `peak_doy_viirs`: the same quantity from **VIIRS VNP22Q2**, 2013+ only. VIIRS runs about **8 days later** than MODIS, and **the gap grew about 7.7 days per decade from 2013 to 2025**. Use it only as a robustness check, never as a substitute target.
  - `mid_greendown_count`, `forest_px_count`, `qa_retention`: pixel counts. `qa_retention` is ≥ 0.998 everywhere, so it is useless as a weight. Use `peak_spread_days` or `forest_px_count` instead. `forest_px_count` is the same every year for a location (fixed NLCD 2021 mask).
- **Derived columns:**
  - `peak_doy_mean` (in `locations.csv`, and repeated on every row of `panel.csv`): the location's mean `peak_doy` over 2001–2025.
  - `peak_anom` = `peak_doy − peak_doy_mean`. **This is the "delta" for exploratory analysis.** It uses all years, so it leaks a little; for honest model scoring, compute the mean leaving out the test year (see `plan.md`).

**What we already know (from the 4-region NC/VA analysis, 2001–2025):** year-to-year SD is about **4 days** at every location. There is **no robust trend**.

**The number to beat (this dataset, 11 locations):** leave-one-year-out climatology ("predict the location's mean of the other years") scores **MAE 3.34 d / RMSE 4.13 d** overall. Per location MAE ranges from 2.51 (Old Forge) to 4.13 (Frostburg); see `BUILD_LOG.md` §8. (The 4-region data gave 3.22 / 3.95.)

## 4. Column dictionary

### 4a. `panel.csv` — identifiers and flags
| Column | Meaning |
|---|---|
| `location`, `year` | Keys |
| `is_forecast_row` | 1 for 2026 (no target), else 0 |
| `helene_2024` | 1 for Asheville and Boone in 2024 (Hurricane Helene defoliation). In the earlier analysis it was not an outlier, but keep it as a sensitivity check. |
| `ecoregion` | EPA L3 name |

### 4b. `panel.csv` — weather features (gridMET, mean over the clipped region; not forest-masked)
Naming rule: `{variable}_{window}`. Fixed windows use calendar months. Units: temperature °C, precipitation mm (total over the window), solar radiation W/m² (mean), VPD kPa (mean), degree-days °C·day.

| Column | Window | Definition |
|---|---|---|
| `tmean_MarApr` | Mar 1–Apr 30 | Mean daily temp (spring warmth, drives green-up) |
| `gdd5_MarMay` | Mar 1–May 31 | Growing degree-days above 5 °C |
| `tmax_JunJul`, `precip_JunJul` | Jun 1–Jul 31 | Early-summer heat and moisture |
| `tmax_Aug` | Aug | Mean daily max. **Strongest single signal in the old data: about +1.1 days later per +1 °C.** |
| `tmin_Sep` | Sep | Mean daily min |
| `tmin_SepOct` | Sep 1–Oct 31 | Mean daily min (post-peak for most locations) |
| `cdd_Sep` | Sep | Cold degree-days = Σ max(0, 20 − tmean) |
| `cdd_Sep1_Oct15` | Sep 1–Oct 15 | Same, longer window. Overlaps the peak. |
| `cold_nights_Sep1_Oct15` | Sep 1–Oct 15 | Number of days with tmin < 5 °C |
| `precip_AugSep` | Aug 1–Sep 30 | Total precipitation. Wetter → earlier in the old data. |
| `precip_Oct` | Oct | Total precipitation |
| `vpd_AugSep` | Aug 1–Sep 30 | Mean vapor pressure deficit (dryness). r = 0.72 with `tmax_Aug`. |
| `srad_Sep`, `srad_SepOct` | Sep / Sep–Oct | Mean shortwave radiation |
| `pdsi_AugSep` | Aug 1–Sep 30 | Mean Palmer Drought Severity Index (negative = drought) |
| `first_frost_doy` | Aug 1 onward | First day with tmin ≤ 0 °C. If there was no frost by Nov 30, set to 335 with `frost_censored = 1`. |
| `frost_censored` | — | 1 if no frost observed (see above) |
| **Relative windows** | | **Windows anchored to each location's own `peak_doy_mean`** (P = round(peak_doy_mean), the same P for every year including 2026, counted in DOY), so they mean the same thing in Vermont and North Carolina. P: Montpelier 272, Old Forge 271, Binghamton 274, Scranton 278, State College 277, Frostburg 275, Valley Head 275, Shenandoah 283, Roanoke 279, Boone 278, Asheville 281. |
| `tmean_pre30` | P−30 … P−1 | Mean temp in the 30 days before typical peak |
| `tmin_pre30` | P−30 … P−1 | Mean min temp |
| `cdd_pre30` | P−30 … P−1 | Cold degree-days (base 20 °C) |
| `cold_nights_pre30` | P−30 … P−1 | Days with tmin < 5 °C |
| `tmax_pre60to31` | P−60 … P−31 | Late-summer heat (the analogue of `tmax_Aug`) |
| `precip_pre60` | P−60 … P−1 | Total precipitation |
| `vpd_pre60` | P−60 … P−1 | Mean VPD |
| `srad_pre30` | P−30 … P−1 | Mean solar radiation |

### 4c. `panel.csv` — satellite phenology covariates (MODIS, forest-pixel medians, DOY)
| Column | Meaning | Usable for forecasting? |
|---|---|---|
| `mid_greenup_p50` | Spring mid-green-up. A late spring went with a late fall in the old data (r = +0.32). | Yes, known by June |
| `senescence_p50` | Start of autumn decline | Borderline (late Aug–Sep) |
| `dormancy_p50` | End of autumn decline | **No.** This is after peak, so using it would leak the answer. |

### 4d. `locations.csv` — static features
| Group | Columns | Notes |
|---|---|---|
| Geography | `lon`, `lat`, `state`, `ecoregion`, `region_area_km2`, `forest_px_count` | |
| Elevation (SRTM 30 m, forest pixels only; percentiles exact) | `elev_p10_m`, `elev_median_m`, `elev_p90_m`, `elev_mean_m`, `frac_forest_above_1350m` | 1,350 m is where peak timing started **dropping about 2.9 d per 100 m** in the earlier elevation analysis. It was flat below that. |
| Terrain | `slope_mean_deg`, `northness_mean` | Northness = cos(aspect). North-facing pixels peaked about 4 d earlier per unit of northness. |
| Land cover (NLCD 2021, share of region area) | `frac_deciduous`, `frac_evergreen`, `frac_mixed`, `frac_woody_wetland`, `frac_developed`, `frac_hardwood`, `hardwood_share_of_forest` | `frac_hardwood` = deciduous + mixed |
| Species / forest type (LANDFIRE EVT, share of forest area) | `evt_northern_hardwood`, `evt_oak_hickory`, `evt_cove_mesic`, `evt_spruce_fir`, `evt_pine`, `evt_other` | *(build note)* Asset `LANDFIRE/Vegetation/EVT/v1_4_0` (LF 2016 Remap, 30 m; the only EVT version in Earth Engine). Denominator = **all EVT tree classes** in the region (conifer included), *not* the NLCD deciduous+mixed mask. Grouping by class-name keywords: explicit northern-hardwood names first, otherwise the earliest-named genus wins ("Dry Oak-Pine" → oak, "Pine-Oak" → pine; "dry-mesic oak" → oak). Ruderal, floodplain, swamp, plantation and urban forest classes → `evt_other`. Six groups sum to 1. The class → group table is `evt_class_map.csv` (columns `evt_value, evt_name, px_total, is_forest, group, matched_keyword`) for auditing/regrouping. Note State College's `evt_other` = 0.37, mostly "Ruderal Forest – Northern and Central Hardwood and Conifer" (young regrowth forest). The species teammate may replace or extend these with FIA or USFS forest-type data. |
| Day length | `daylen_Oct1_h`, `doy_daylen_below_11h` | Calculated from latitude (Forsythe et al. 1995) |
| Climate normals (gridMET means 2001–2025) | `norm_tmean_Sep`, `norm_tmin_Oct`, `norm_precip_AugSep`, `norm_gdd5_MarMay` | For explaining **between-location** differences |
| Typical peak | `peak_doy_mean`, `peak_doy_sd`, `peak_doy_trend_d_per_decade`, `n_target_years` | Computed from 2001–2025 |

## 5. Forecast-safe vs. explanatory features

A feature is **forecast-safe** for a location only if its window ends **before that location's typical peak**. Because peak is earlier in Vermont than in North Carolina, a feature like `cdd_Sep1_Oct15` can be safe in one place and not another. `feature_timing.csv` lists, for each (feature, location):
`window_type` (fixed/relative), `window_start_doy`, `window_end_doy` (non-leap DOY), `peak_doy_mean`, `days_before_peak` (= peak_doy_mean − window_end_doy), `timing` ∈ {`pre-peak`, `partly`, `post`}.
*(build note)* Rule: `pre-peak` if the window ends before P = round(peak_doy_mean); `partly` if it starts before the mean peak but ends on/after P; otherwise `post`. (Using P instead of the raw "≥ 1 day" rule keeps the `*_pre30` windows pre-peak by construction.)

Result in this build: all spring/summer features, `tmax_Aug` and every relative window are pre-peak everywhere. **Every September feature (`tmin_Sep`, `cdd_Sep`, `srad_Sep`, `precip_AugSep`, `vpd_AugSep`, `pdsi_AugSep`) is `partly` for Montpelier and Old Forge** (their typical peak is Sep 28) and pre-peak elsewhere. `*_Sep1_Oct15`, `*_SepOct`, `first_frost_doy` are `partly` everywhere; `precip_Oct` is `post` in the two northern sites.

- **Forecasting this year's peak:** use only `pre-peak` features. This is why the relative windows (`*_pre30`, `*_pre60`) are there: they are pre-peak everywhere by construction.
- **Explaining why a year was early or late:** any feature is fine, but say which kind of claim you are making.

## 6. How to analyze it (shared rules so our results combine)

1. **Work in anomalies.** For year-to-year questions, subtract each location's mean from both the target and the feature (`x − mean_by_location(x)`), then relate the two. This is what the delta model will use.
2. **Effective sample size is about 26 years, not 286 rows.** Nearby locations share weather in the same year, so their rows are not independent. Report significance with **year-clustered standard errors** or on year-averaged data. With 11 locations spread from VT to NC, years are less tightly linked than they were in the 4-region data, but they are still correlated.
3. **Static features cannot explain a delta on their own.** For a fixed location, elevation or species mix is the same every year, so it is absorbed into the mean. Static factors matter in two ways only:
   - (a) they explain **why location means differ** (cross-sectional, n = 11 locations), or
   - (b) they **change how strongly a location responds to weather** (interaction: e.g. does the `tmax_Aug` effect grow with `evt_northern_hardwood`?).
   Frame environmental and species results as (a) or (b).
4. **Validate by leaving years out:** hold out one year at a time (all locations together), refit, predict, and report MAE/RMSE against the climatology baseline.
5. **Don't use** `dormancy_p50`, `peak_doy_p25/p75`, `peak_doy_viirs` or `mid_greendown_count` as predictors. They are measured at or after the target.
6. **Units:** days are DOY, temperatures °C, precipitation mm. Positive `peak_anom` = later than usual.

## 7. What to hand back (the "contract")

When you finish, send a short Markdown summary with this table so the results can go straight into the combined model:

| feature | role (delta main effect / interaction with ___ / explains mean) | effect size (days per unit, ± SE) | forecast-safe? | evidence (r, p with year-clustered SE, held-out-year MAE change) | notes |
|---|---|---|---|---|---|

Also list: features you tried that **did not** help, any data problems you found, and any new columns you created (with the exact formula so they can be rebuilt).

## 8. Quick-start loading (Python)

```python
import pandas as pd
P = "IMMC/data/panel_v1/"
panel = pd.read_csv(P + "panel.csv")
locs  = pd.read_csv(P + "locations.csv")
train = panel[panel.is_forecast_row == 0]

# within-location anomaly of any column
anom = lambda df, c: df[c] - df.groupby("location")[c].transform("mean")
train = train.assign(tmax_Aug_anom=anom(train, "tmax_Aug"))
```

## 9. Known caveats

- **MODIS vs. VIIRS trend disagreement** (see §3). Don't claim "peak is getting later/earlier" from one sensor.
- **2026 weather** runs to **2026-11-11** (gridMET) and 2026-11-06 (PDSI). `first_frost_doy = 335` with `frost_censored = 1` in 2026 (Shenandoah, Roanoke) means "no frost by Nov 11".
- **Pixel dates are integers**, so every `peak_doy` is a whole or half day. This is harmless.
- **Values differ slightly from the old 4-region `model_df.csv`** for Boone, Asheville and Shenandoah: satellite percentiles by up to ~1 d (mean 0.3 d), elevation percentiles by ~3 m, and pixel counts are 2–5% lower. The old values used an approximate Earth Engine percentile; these are exact. Weather columns are identical. Details in `BUILD_LOG.md`.
- **Typical peak does not follow latitude alone:** Montpelier/Old Forge ≈ Sep 28, but Frostburg and Valley Head (Oct 2) are earlier than Scranton/State College (Oct 4), and Shenandoah (Oct 10) is the latest of all, as in the old data.
- **Region overlap:** the Binghamton and Scranton regions sit about 80 km apart, so their 40 km buffers almost touch, and their weather is highly correlated.
- **The proxy offset:** mid-greendown ≠ visual peak. The comparison to news reports happens at the end of the plan, not in your analysis.

**Sources (Earth Engine IDs):** `MODIS/061/MCD12Q2`, `NASA/VIIRS/002/VNP22Q2`, `IDAHO_EPSCOR/GRIDMET`, `GRIDMET/DROUGHT`, `USGS/NLCD_RELEASES/2021_REL/NLCD/2021`, `USGS/SRTMGL1_003`, `EPA/Ecoregions/2013/L3`, LANDFIRE EVT (see the build doc for the exact asset).
