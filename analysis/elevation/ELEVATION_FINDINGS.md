# Elevation side analysis: peak fall color vs elevation

**Data.** MODIS MCD12Q2 v061 `MidGreendown_1` (the same proxy for peak color as the main model), 2001–2025, on the
MODIS 463 m grid. Pixels are deciduous or mixed forest (NLCD 2021 classes 41 and 43, at least 50% of the pixel). Regions are a
40 km buffer intersected with the EPA L3 ecoregion. Each region has up to 4,000 random forest pixels (Raleigh has 2,766).
Elevation is the SRTM 30 m mean over each pixel. In total there are 369,138 pixel-years (`pixels.csv`).

**Method.** Within each region × year, every pixel's peak DOY is compared with that region-year's median. This removes the year's
weather. The **anomaly** is then regressed on elevation with region intercepts. 95% CIs come from a block
bootstrap over years (1,000 reps). Mountain regions = Asheville, Boone and Shenandoah.
Scripts: `01_pull_pixels.py` (Earth Engine, run once) → `02_elevation_analysis.py` → `03_figures.py`;
`elevation_adjust.py` contains the adjustment rule.

## Headline

**The effect is not linear. Peak color is flat with elevation up to about 1,350 m. Above that it arrives
about 2.9 days earlier per 100 m (95% CI 2.6–3.3).**

| Model (mountain regions pooled) | Slope, days per 100 m (95% CI) |
|---|---|
| Single linear slope (the requested headline) | **−0.18** (−0.29 to −0.07) |
| Two-piece, below 1,350 m | **−0.07** (−0.19 to +0.05) |
| Two-piece, above 1,350 m | **−2.93** (−3.32 to −2.57) |

- The single slope is significant but small. It averages a flat stretch and a steep one, so **the single slope
  should not be used as the adjustment rule.** Pixel-clustered SEs give a narrower CI (−0.22 to −0.14).
  The year bootstrap is the more conservative of the two.
- The breakpoint was chosen by least squares on a 50 m grid (900–1,600 m). **All 1,000 bootstrap
  replicates chose 1,350 m**, and the CIs above include that re-selection. Using Asheville and Boone only, the break is
  1,300 m (95% 1,250–1,300 m) and the upper slope is the same (−2.93).
- The handoff's fixed 1,000 m hinge fits poorly: below −0.15, above −0.25, and the change in slope
  is not significant (−0.10, CI −0.25 to +0.06). The steepening is real, but it starts higher than 1,000 m.
- The steep segment rests on **325 sampled pixels above 1,350 m** (206 Asheville, 119 Boone). Shenandoah
  tops out at 1,174 m, so it only contributes to the flat segment.

## Per region (single slope, days per 100 m)

| Region | Slope (95% CI) | Elevation p5–p95 | Above 1,350 m (two-piece) |
|---|---|---|---|
| Asheville | +0.47 (+0.31 to +0.63) | 584–1,358 m | **−5.00** (−5.55 to −4.53) |
| Boone | −0.30 (−0.43 to −0.15) | 430–1,283 m | **−2.51** (−2.92 to −2.04) |
| Shenandoah | −0.88 (−1.03 to −0.74) | 307–952 m | n/a (max 1,174 m) |
| Raleigh | +3.67 (+2.45 to +4.66) | 67–155 m | n/a |

- Below the break, the regions **disagree in sign**. The low-elevation bumps differ by region: Asheville's 500–700 m pixels (the
  escarpment edge) are about 5 days early, while Boone's 500–700 m pixels are about 4 days late. These look like local landscape and forest-type effects,
  not an elevation gradient. Only Shenandoah has a consistent modest decline below 1,350 m.
- **Raleigh should not be interpreted.** Its elevation range is only about 90 m, and the slope more than doubles when lon/lat
  are added (+8). The slope is picking up river-bottom, lake-edge and urban-edge differences, not elevation.

## Stability across years (`slopes_by_year*.csv`, Fig E2 right)

- Single slope, mountain pooled: median −0.25, range −0.73 (2011) to +0.31 (2019). It is negative in 17 of 25 years.
  The sign is **not** stable from year to year.
- Two-piece **above 1,350 m: negative in all 25 years**, median −2.8, range −5.6 to −0.8. The steep upper segment is
  reliable every year, which is what makes the adjustment rule trustworthy.
- Two-piece below 1,350 m: median −0.13, range −0.66 to +0.45, centered near zero.

## Linearity (`binned_medians.csv`, Fig E1)

Binned medians (100 m bins, at least 50 pixel-years) stay within ±2 days of zero from 200 to 1,350 m in the pooled data. They then fall steadily:
−3 d (1,400s), −6 (1,500s), −9 (1,600s), −13 (1,700s) and −15 (1,800s). This is consistent with the shift from oak and cove
hardwoods to northern hardwoods (beech, birch, buckeye) near that elevation, which color earlier. That interpretation is plausible but was not tested here.

## Confounders and sensitivity

- **Aspect:** northness has a strong effect of its own. North-facing pixels peak **4.3 days earlier per unit of
  northness** (CI 4.0–4.7). Because northness is averaged over a 463 m pixel, most values fall roughly between −0.5 and +0.5.
  It is nearly uncorrelated with elevation (r = 0.01), so **the elevation slope barely moves** (−0.18 → −0.15).
- **Forest fraction** (not in the handoff, added as a check): pixels with less forest peak earlier (+13 d per unit, so a 50%-forest
  pixel peaks about 6.6 days earlier than a 100% one). With northness and forest fraction both added, the single slope is −0.31 and the
  two-piece model is −0.22 below and **−2.61 above** 1,350 m (CI −2.95 to −2.26). The break is unchanged. The conclusion holds.
- **Helene (dropping 2024 for Asheville and Boone):** the single slope stays at −0.18 (−0.30 to −0.07). The two-piece model goes from −0.07 / −2.93
  to −0.07 / −2.92. **Helene has no measurable effect.**

## Adjustment rule (`elevation_adjust.py`, `offset_table.csv`)

Use `adjust_for_elevation_two_piece(forecast_doy, target_elev_m, region_median_elev_m, hinge_m=1350,
slope_below=-0.07, slope_above=-2.93)`. The reference is the region's `elev_median_m` from `static_features.csv`
(Asheville 883 m, Boone 931 m). The linear `adjust_for_elevation` from the handoff is also included, for comparison only.

Offset in days relative to the regional forecast (negative = earlier):

| Target elevation | Asheville (ref 883 m) | Boone (ref 931 m) |
|---|---|---|
| 500 m | +0.3 | +0.3 |
| 650 m (Asheville city) | +0.2 | — |
| 800 m | +0.1 | +0.1 |
| 1,000 m | −0.1 | 0.0 |
| 1,200 m | −0.2 | −0.2 |
| 1,500 m | **−4.7** | **−4.7** |
| 1,800 m (Grandfather Mountain, Boone) | −13.5 | **−13.5** |
| 2,037 m (Mount Mitchell, Asheville) | **−20.4** ⚠ extrapolated | — |

Example: "If Asheville's forecast is Oct 20 for its median forest elevation of 883 m, a 1,500 m ridge peaks about 5 days earlier (Oct 15),
and Mount Mitchell peaks about 3 weeks earlier (around Sep 30)." The region-specific fits are steeper for Asheville (−18 d at 1,800 m) and similar for Boone
(−12 d). See the `offset_region_specific_d` column. The pooled rule is the more conservative choice.

## Figures (`figures/`)

- `E1_anomaly_vs_elevation.png`: Peak-color anomaly vs pixel elevation for the three mountain regions. Pooled bin medians with IQR, per-region medians, and two-piece and single-slope fits. Peak color is flat up to about 1,350 m, then arrives about 2.9 days earlier per 100 m.
- `E2_slopes_by_region_and_year.png`: Slope estimates with 95% CIs by model and region (left), and slopes fitted separately for each year (right). The upper-elevation slope is negative every year.
- `E3_peak_date_vs_elevation.png`: Typical peak calendar date vs elevation for Asheville and Boone, with the pooled two-piece rule. Suited to the editor letter.

## Caveats

- Each pixel is about 463 m across and averages over terrain. Elevation is the **pixel mean**, so ridgetops and summits are smoothed toward
  their slopes. Sampled pixel means reach only 1,819 m, so **anything above about 1,800 m (Mount Mitchell) is extrapolated**. Summit spruce-fir is also
  excluded by the deciduous/mixed forest mask.
- The slope is an **empirical within-region relationship, not a causal lapse rate**. Forest composition changes with
  elevation, and that change is part of the effect.
- Above 1,350 m the data come only from Asheville and Boone (325 pixels). Treat the steep segment as specific to the southern Blue Ridge.
- MidGreendown is the midpoint of EVI decline, which is a proxy for peak color, not a direct observation of it. The offsets
  carry over to whatever the main model forecasts only if that forecast uses the same MODIS target.
- QA filtering (`QA_Overall_1 ≤ 1`) removed almost nothing. Nearly all forest pixels have QA 0.
