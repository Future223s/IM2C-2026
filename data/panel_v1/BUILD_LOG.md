# BUILD_LOG — panel_v1

Built 2026-11-15 14:51 by `data/build_panel.py`. Spec: `docs/BUILD_DATASET.md`.

## Check summary

| check | result |
|---|---|
| 1 regression | **PASS** |
| 2 shape | **PASS** |
| 3 ranges | **PASS** |
| 4 geography | **PASS** |
| 5 missingness | **PASS** |
| 6 MODIS vs VIIRS | **PASS** |
| 7 pixel counts | **PASS** |
| 8 baseline | **PASS** |

# Build notes

## regions

### Regions (40 km buffer clipped to EPA L3 ecoregion)
| location | ecoregion | region_area_km2 | flag |
|---|---|---|---|
| Montpelier | Northeastern Highlands | 4,941 |  |
| Old Forge | Northeastern Highlands | 4,594 |  |
| Binghamton | Northern Allegheny Plateau | 4,966 |  |
| Scranton | Ridge and Valley | 1,064 | **< 1,500 km²** |
| State College | Ridge and Valley | 3,816 |  |
| Frostburg | Central Appalachians | 2,085 |  |
| Valley Head | Ridge and Valley | 1,816 |  |
| Shenandoah | Blue Ridge | 1,215 | **< 1,500 km²** |
| Roanoke | Ridge and Valley | 2,308 |  |
| Boone | Blue Ridge | 4,591 |  |
| Asheville | Blue Ridge | 4,963 |  |

Closest pairs (center distance, overlap of clipped regions):
- Binghamton – Scranton: 79 km apart, overlap 0.0 km²
- Boone – Asheville: 102 km apart, overlap 0.0 km²
- Frostburg – Shenandoah: 136 km apart, overlap 0.0 km²
- Valley Head – Shenandoah: 142 km apart, overlap 0.0 km²

Pairs with any overlap: 0

## pheno

- MODIS: years 2001–2025 (25), encoding `days_since_1970`
- VIIRS: years 2013–2025 (13), encoding `days_since_2000`
- Percentiles are EXACT (unweighted per-region pixel histograms of integer-day dates; numpy 'linear' percentile definition), not `ee.Reducer.percentile` (see deviations).
- MODIS MCD12Q2 most recent year: **2025**
- MODIS 2025 target missing for: none
- Rows with mid_greendown_p50 outside DOY 230–350: 0

## weather

- gridMET daily: 78,441 rows, 2001-03-01 → **2026-11-11** (last available date)
- GRIDMET/DROUGHT pdsi: 15,686 rows, 2001-03-01 → 2026-11-06
- Location-years (2001–2025) with < 275 days: 0
- Missing daily values: none

## static

- Elevation p10/p50/p90 are EXACT (unweighted 30 m forest-pixel histogram of integer SRTM metres); `elev_mean_m` uses the notebook's mean reducer.
- LANDFIRE EVT asset: `LANDFIRE/Vegetation/EVT/v1_4_0` (image `CONUS`, LF 2016 Remap; the only EVT version in the Earth Engine catalog on 2026-11-15). Share of EVT forest pixels.
- EVT classes present in the 11 regions: 113, of which forest: 60. Full mapping in `raw/evt_class_map.csv`.
- static.csv: 11 locations × 23 columns

## features

- Relative windows anchored at P = round(peak_doy_mean): Montpelier 272, Old Forge 271, Binghamton 274, Scranton 278, State College 277, Frostburg 275, Valley Head 275, Shenandoah 283, Roanoke 279, Boone 278, Asheville 281
- Location-years with an incomplete 60-day relative window: 0
- yearly features: 286 rows × 28 cols; location features: (11, 11)

## assemble

- panel.csv: 286 rows × 45 cols
- locations.csv: 11 rows × 37 cols
- weather_daily.csv: 78,441 rows; pdsi.csv: 15,686 rows
- feature_timing.csv: 275 rows
- panel_v1/DATASET_GUIDE.md kept (it holds the build's deviations from docs/ version)
- pixels.csv (Stage 3d): 1,080,330 pixel-year rows, 4,000 forest pixels per location, 98.9% of peak_doy in DOY 260–330. Shipped separately as panel_v1_pixels.zip.

## deviations

### Percentiles: exact instead of approximate (decided by Arin, 2026-11-15)

- First build: the regression test failed on `peak_doy` / `peak_doy_viirs` (up to 1.2 d) and Boone `elev_median_m` (0.001 m);
  weather, frost and NLCD columns matched to ~1e-12.
- Cause (verified): `ee.Reducer.percentile` is a bucketed approximation whose result for one region depends on the *other*
  regions in the same `reduceRegions` request. The ported code with the notebook's original 5 regions reproduced the old
  2001 medians exactly (Boone 275.354009, Asheville 278.252745, Shenandoah 283.411840); with the 11 regions they became
  275.312 / 277.553 / 284.000, and Boone alone 275.450. The old `model_df.csv` values therefore carry ~±1 d of this noise.
- Fix: every phenology and elevation percentile is now computed **exactly** from an unweighted per-region pixel histogram
  (dates and SRTM are integers), using numpy's 'linear' percentile definition. A location's value no longer depends on the
  other locations. Raw histograms: `raw/phenology_histograms.csv`, `raw/elevation_histograms.csv`.
- Side effect: "unweighted" counts a pixel only if its **center** is inside the clipped region. The old weighted reducers also
  counted edge pixels that partly overlap the region. `forest_px_count` and `*_count` are therefore 2–5% lower than in
  `model_df.csv` (Boone −359, Asheville −299, Shenandoah −298; the same every year).
- Regression test now: exact (< 1e-6) for weather, frost and `frac_hardwood`; tolerance 1.5 d for `peak_doy` /
  `peak_doy_viirs` and 5 m for `elev_median_m`.
- `raw/phenology_approx_percentile.csv` and `raw/static_approx_percentile.csv` keep the first (approximate) build for reference.

# Validation (Stage 7)

## 1. Regression test (Boone, Asheville, Shenandoah vs `data/data/model_df.csv`)

| column | n_compared | max_abs_diff | mean_abs_diff | tolerance | n_diff_ge_1e6 | nan_mismatch | worst | result |
|---|---|---|---|---|---|---|---|---|
| peak_doy | 75 | 0.963 | 0.294 | 1.5 | 57 | 0 | Asheville 2024: 276.0375 → 277.0000 | PASS |
| peak_doy_viirs | 39 | 0.777 | 0.209 | 1.5 | 21 | 0 | Boone 2015: 282.2227 → 283.0000 | PASS |
| tmax_Aug | 75 | 1.6e-13 | 4.2e-14 | 1e-06 | 0 | 0 |  | PASS |
| precip_AugSep | 75 | 2.84e-13 | 3.83e-14 | 1e-06 | 0 | 0 |  | PASS |
| cdd_Sep | 75 | 3.87e-12 | 9.57e-13 | 1e-06 | 0 | 0 |  | PASS |
| pdsi_AugSep | 75 | 3.55e-15 | 3.71e-16 | 1e-06 | 0 | 0 |  | PASS |
| first_frost_doy | 75 | 0 | 0 | 1e-06 | 0 | 0 |  | PASS |
| elev_median_m | 75 | 3.5 | 2.83 | 5 | 75 | 0 | Shenandoah 2001: 563.5039 → 560.0000 | PASS |
| frac_hardwood | 75 | 2.22e-16 | 1.48e-16 | 1e-06 | 0 | 0 |  | PASS |

## 2. Shape

- rows: 286 (expected 286); duplicate (location, year): 0
- 2001–2025 rows with null `peak_doy`: 0 
- 2026 rows with non-null `peak_doy`: 0

## 3. Ranges

| check | n | min | max | allowed | out_of_range |
|---|---|---|---|---|---|
| panel.peak_doy | 275 | 262.000 | 291.000 | [230, 350] | 0 |
| panel.tmax_Aug | 286 | 21.818 | 32.638 | [15, 35] | 0 |
| panel.precip_AugSep | 286 | 80.077 | 549.844 | [0, inf] | 0 |
| panel.pdsi_AugSep | 286 | -4.860 | 5.874 | [-8, 8] | 0 |
| pdsi.pdsi (daily file) | 15686 | -6.113 | 7.618 | [-8, 8] | 0 |
| panel.qa_retention | 275 | 0.998 | 1.000 | [0, 1] | 0 |
| locations.elev_median_m | 11 | 392.000 | 955.000 | [50, 1500] | 0 |
| locations.frac_forest_above_1350m | 11 | 0.000 | 0.051 | [0, 1] | 0 |
| locations.frac_deciduous | 11 | 0.351 | 0.791 | [0, 1] | 0 |
| locations.frac_evergreen | 11 | 0.004 | 0.169 | [0, 1] | 0 |
| locations.frac_mixed | 11 | 0.043 | 0.269 | [0, 1] | 0 |
| locations.frac_woody_wetland | 11 | 0.000 | 0.113 | [0, 1] | 0 |
| locations.frac_developed | 11 | 0.014 | 0.278 | [0, 1] | 0 |
| locations.frac_hardwood | 11 | 0.594 | 0.919 | [0, 1] | 0 |
| locations.evt_northern_hardwood | 11 | 0.036 | 0.740 | [0, 1] | 0 |
| locations.evt_oak_hickory | 11 | 0.001 | 0.587 | [0, 1] | 0 |
| locations.evt_cove_mesic | 11 | 0.000 | 0.498 | [0, 1] | 0 |
| locations.evt_spruce_fir | 11 | 0.002 | 0.145 | [0, 1] | 0 |
| locations.evt_pine | 11 | 0.003 | 0.110 | [0, 1] | 0 |
| locations.evt_other | 11 | 0.040 | 0.365 | [0, 1] | 0 |
| locations.hardwood_share_of_forest | 11 | 0.790 | 0.993 | [0, 1] | 0 |

## 4. Geography sense-check (north → south)

| location | lat | elev_median_m | peak_doy_mean | peak_doy_sd | peak_date_typical |
|---|---|---|---|---|---|
| Montpelier | 44.3 | 445.0 | 271.6 | 3.5 | Sep 28 |
| Old Forge | 43.7 | 575.0 | 271.2 | 3.3 | Sep 28 |
| Binghamton | 42.1 | 435.0 | 274.1 | 4.1 | Oct 1 |
| Scranton | 41.4 | 465.0 | 277.8 | 4.1 | Oct 4 |
| State College | 40.8 | 392.0 | 277.5 | 4.0 | Oct 4 |
| Frostburg | 39.7 | 716.0 | 275.1 | 4.7 | Oct 2 |
| Valley Head | 38.5 | 955.0 | 275.0 | 4.2 | Oct 2 |
| Shenandoah | 38.5 | 560.0 | 283.4 | 4.0 | Oct 10 |
| Roanoke | 37.3 | 532.0 | 278.8 | 4.6 | Oct 5 |
| Boone | 36.2 | 929.0 | 278.5 | 3.8 | Oct 5 |
| Asheville | 35.6 | 881.0 | 281.0 | 4.0 | Oct 7 |

- r(lat, peak_doy_mean) = -0.77. Montpelier + Old Forge mean 271.4 vs others 277.9.

## 5. Missingness (non-null share)

Columns below 100% (full table: `raw/missingness.csv`):

- none: every column is 100% populated within its scope.

- Columns under 95%: none

## 6. MODIS vs VIIRS (2013–2025)

| location | n | bias_d | rmse_d | r |
|---|---|---|---|---|
| Montpelier | 13 | 9.31 | 9.95 | 0.58 |
| Old Forge | 13 | 9.46 | 10.01 | 0.57 |
| Binghamton | 13 | 9.23 | 9.84 | 0.66 |
| Scranton | 13 | 10.15 | 10.56 | 0.76 |
| State College | 13 | 10.23 | 10.72 | 0.67 |
| Frostburg | 13 | 10.38 | 10.75 | 0.81 |
| Valley Head | 13 | 11.23 | 11.59 | 0.77 |
| Shenandoah | 13 | 9.77 | 10.19 | 0.75 |
| Roanoke | 13 | 9.08 | 9.73 | 0.62 |
| Boone | 13 | 9.92 | 10.42 | 0.59 |
| Asheville | 13 | 8.54 | 9.09 | 0.52 |
| All (pooled) | 143 | 9.76 | 10.28 | 0.67 |

(Pooled r is between within-location anomalies.) FINDINGS §1 (4 regions): bias ≈ +8 d, anomaly r = 0.50.
- VIIRS − MODIS gap trend, mean over locations: +7.3 d/decade (FINDINGS: +7.7).

## 7. MODIS forest pixel counts

| location | forest_px_count | region_area_km2 | frac_hardwood | flag |
|---|---|---|---|---|
| Montpelier | 16,052 | 4,940.843 | 0.664 |  |
| Old Forge | 14,576 | 4,594.392 | 0.635 |  |
| Binghamton | 14,947 | 4,966.483 | 0.594 |  |
| Scranton | 3,215 | 1,063.936 | 0.619 |  |
| State College | 11,604 | 3,816.120 | 0.622 |  |
| Frostburg | 6,893 | 2,085.228 | 0.664 |  |
| Valley Head | 7,469 | 1,816.471 | 0.830 |  |
| Shenandoah | 5,476 | 1,215.426 | 0.919 |  |
| Roanoke | 7,055 | 2,307.619 | 0.629 |  |
| Boone | 18,346 | 4,591.085 | 0.775 |  |
| Asheville | 16,973 | 4,963.086 | 0.700 |  |

- Locations under 2,000 forest pixels: none (Raleigh had 2,905).

## 8. Climatology baseline (leave-one-year-out)

| location | n | MAE | RMSE |
|---|---|---|---|
| Montpelier | 25 | 2.78 | 3.60 |
| Old Forge | 25 | 2.51 | 3.37 |
| Binghamton | 25 | 3.27 | 4.16 |
| Scranton | 25 | 3.18 | 4.23 |
| State College | 25 | 3.39 | 4.10 |
| Frostburg | 25 | 4.13 | 4.77 |
| Valley Head | 25 | 3.55 | 4.28 |
| Shenandoah | 25 | 3.45 | 4.04 |
| Roanoke | 25 | 3.90 | 4.71 |
| Boone | 25 | 3.27 | 3.85 |
| Asheville | 25 | 3.28 | 4.11 |
| **Overall** | 275 | 3.34 | 4.13 |

4-region reference (FINDINGS §2): MAE 3.22 / RMSE 3.95.
