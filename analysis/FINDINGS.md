# Pre-modeling findings: peak fall foliage

Data: `data/data/model_df.csv`. 4 regions × 2001–2025 = **100 rows with a target**. The handoff
expected about 96. 2025 does have a target, so only the 4 rows for 2026 are forecast rows.
All numbers below come from the scripts in `analysis/`. Run them in order (`python 01_target_eda.py`, …)
with the `fall-foliage` env. Each script also saves a CSV of its table to `analysis/`.

## 1. Target summary (`01_target_eda.py`)

| Region | Mean | SD (d) | Earliest | Latest | Range (d) |
|---|---|---|---|---|---|
| Raleigh | Oct 5 (278.5) | 3.9 | Sep 26 (2025) | Oct 14 | 18 |
| Asheville | Oct 7 (280.9) | 4.0 | Oct 2 | Oct 17 | 15 |
| Boone | Oct 5 (278.4) | 3.6 | Oct 1 | Oct 14 | 13 |
| Shenandoah | Oct 10 (283.4) | 3.9 | Oct 5 | Oct 18 | 13 |

- **Year-to-year SD is about 4 days in every region.** Regions differ by up to 5 days on average.
- **Trend:** the pooled slope with region intercepts is **−1.3 d/decade**. The naive OLS p is 0.017, but
  with year-clustered SEs (which allow for the shared weather) it is **p = 0.20** (95% CI −3.2 to +0.7). Only Raleigh has a
  clear trend on its own: **−2.9 d/decade, p = 0.006**. Asheville, Boone and Shenandoah are between −0.4 and −1.1, all with p > 0.3.
  Dropping the Helene rows changes the pooled slope to −1.1 (p = 0.27). Conclusion: **there is no robust trend in the MODIS target.**
- **Outliers:** none have |z| > 2.5. The most extreme points (|z| > 2) are Raleigh 2025 (Sep 26, z = −2.4),
  Boone 2007, Asheville 2007, Raleigh 2002, Asheville 2005 and Shenandoah 2021.
  **Helene 2024 is not an outlier.** Asheville (z = −1.2) and Boone (z = −1.1) are early, but inside the normal range.
- **MODIS–VIIRS agreement (2013–2025, 13 years per region):**

  | Region | Bias (VIIRS − MODIS) | RMSE | r |
  |---|---|---|---|
  | Raleigh | +4.6 d | 6.1 | 0.21 |
  | Asheville | +8.7 d | 9.2 | 0.53 |
  | Boone | +9.9 d | 10.5 | 0.51 |
  | Shenandoah | +9.7 d | 10.2 | 0.69 |
  | Within-region anomalies, pooled | — | 3.5 | 0.50 |

  VIIRS runs about 8 days later than MODIS on average. Their year-to-year anomalies agree only moderately (r = 0.50).
  **⚠ The two sensors also disagree about the trend.** Over 2013–2025 MODIS has a trend of −0.9 d/decade (p = 0.68).
  VIIRS has **+6.8 d/decade (p = 0.001)**, so the gap between them grows by **+7.7 d/decade
  (p < 0.001)**. This is visible in `01_peak_by_year.png`. The "is peak getting later?" answer
  therefore depends on which sensor is used. Mention this in the report's uncertainty section. Do not
  read either trend as climate signal without caveats. A possible cause is sensor or orbit drift late in
  the MODIS record, but these scripts did not test that.

## 2. Baselines: the numbers to beat (`02_baseline.py` → `baseline_scores.csv`)

Leave-one-year-out (all 4 regions held out together):

| Model | MAE (d) | RMSE (d) |
|---|---|---|
| **Climatology (region mean of other years)** | **3.22** | **3.95** |
| Region intercept + shared linear trend | 3.30 | 3.99 |
| Region intercept + per-region trend | 3.28 | 3.96 |

By region, the climatology MAE is Raleigh 3.08, Asheville 3.25, Boone 3.16 and Shenandoah 3.40.
**Adding a trend does not help.** The real model has to beat **MAE 3.2 d / RMSE 3.95 d**.

## 3. Anomaly correlations (`03_anomaly_correlations.py`)

These are Pearson r values between within-region anomalies and the peak-date anomaly (n = 100). "Year-mean r" collapses the
data to one point per year (n = 25). That is a conservative check, because the 4 regions share weather and are not independent.

| Feature | Timing | r | p | Year-mean r |
|---|---|---|---|---|
| tmax_Aug | pre-peak | **+0.39** | <0.001 | +0.49 |
| cdd_Sep1_Oct15 | partly | **−0.38** | <0.001 | −0.39 |
| mid_greenup_p50 | pre-peak | **+0.32** | 0.001 | +0.36 |
| precip_AugSep | pre-peak | −0.31 | 0.002 | −0.39 |
| vpd_AugSep | pre-peak | +0.29 | 0.003 | +0.42 |
| srad_Sep | pre-peak | +0.29 | 0.003 | +0.39 |
| cold_nights_Sep1_Oct15 | partly | −0.27 | 0.007 | −0.27 |
| pdsi_AugSep | pre-peak | −0.27 | 0.007 | −0.38 |
| cdd_Sep | pre-peak | −0.26 | 0.008 | −0.25 |
| year | — | −0.24 | 0.016 | −0.27 |
| tmin_SepOct | post | +0.19 | 0.06 | +0.25 |
| srad_SepOct, precip_Oct, tmin_Sep, first_frost_doy | | ≤ 0.14 | n.s. | |

- **What the signs say:** a hot, dry, sunny August–September delays greendown. A wet late summer and early cold
  advance it. A late spring green-up goes with a late autumn. First frost has essentially no relationship with
  the target (r = 0.07), which is expected because the frost usually comes *after* the Oct 5 mid-greendown.
- **All the effects are modest.** The best single feature explains about 15% of the anomaly variance.
- **Helene:** dropping the 2 Helene rows changes every r by ≤ 0.04. The largest change is precip_AugSep
  (−0.31 → −0.27), because Helene's rain gives those two rows the largest precipitation anomalies in the data. Nothing changes rank in a
  meaningful way.
- `first_frost_doy`: **0 of the 100 target rows are censored.** Censoring only affects the 2026 forecast rows
  (Raleigh and Shenandoah).

## 4. Redundancy and shortlists (`04_feature_selection.py`)

**Redundancy cluster** (|r| > 0.8 between anomalies): {cdd_Sep1_Oct15, cdd_Sep, tmin_SepOct,
tmin_Sep}. This is the "September–October coolness" group. Every other feature is a singleton. Other notable pairs
below the threshold are tmax_Aug–vpd_AugSep (0.72) and srad_Sep–srad_SepOct (0.76). With all 15 year-varying features in one model, the VIFs
for the cold cluster are 10–11.5.

**(a) Forecast-safe shortlist** (all known by Oct 1). VIFs are 1.02–1.08. In-sample anomaly R² is 0.41 (a screen, not a validated model).
- `tmax_Aug` (r = +0.39): the strongest pre-peak signal. It stands in for vpd_AugSep.
- `mid_greenup_p50` (r = +0.32): spring timing carries over to autumn. It is nearly uncorrelated with the weather, so it adds unique information.
- `precip_AugSep` (r = −0.31): late-summer moisture. It keeps its effect jointly, while pdsi_AugSep loses its effect once tmax_Aug is in.
- `cdd_Sep` (r = −0.26): the best pre-peak member of the cold cluster, and the hook for a degree-day model.

**(b) Explanatory shortlist.** VIFs are 1.01–1.07. In-sample anomaly R² is 0.48.
- `tmax_Aug`, `mid_greenup_p50`, `precip_AugSep`: same reasons as above.
- `cdd_Sep1_Oct15` (r = −0.38): the best member of the cold cluster. Its window runs past the mean peak date,
  so use it for explanation only.

## 5. Consistency across regions (`05_feature_consistency.py`)

Per-region slopes ± SE (days per unit):

| Feature | Raleigh | Asheville | Boone | Shenandoah | Pooled | Equal-slope p |
|---|---|---|---|---|---|---|
| tmax_Aug (/°C) | +1.22 ± 0.48 | +1.22 ± 0.64 | +1.16 ± 0.54 | +0.94 ± 0.59 | +1.14 | 0.98 |
| cdd_Sep1_Oct15 (/°C·d) | +0.020 ± 0.033 | −0.039 ± 0.015 | −0.038 ± 0.013 | −0.042 ± 0.017 | −0.034 | 0.31 |
| mid_greenup_p50 (/d) | +0.29 ± 0.19 | +0.46 ± 0.15 | +0.37 ± 0.16 | +0.04 ± 0.16 | +0.27 | 0.26 |
| precip_AugSep (/mm) | −0.004 ± 0.008 | −0.018 ± 0.007 | −0.015 ± 0.007 | −0.008 ± 0.006 | −0.011 | 0.44 |
| cdd_Sep (/°C·d) | +0.053 ± 0.066 | −0.056 ± 0.028 | −0.049 ± 0.023 | −0.036 ± 0.033 | −0.041 | 0.45 |

- **A shared slope is reasonable for every feature.** No equal-slope F-test rejects (all p ≥ 0.26).
  The tmax_Aug effect is very consistent, at about +1.1 days per °C.
- **Raleigh is the exception for the cold terms.** Its cdd slopes have the wrong sign, but they are not significant. Raleigh rarely
  gets cold in September, so its cdd range is small. It is also only 11% deciduous forest, with about 2,900 pixels. Consider letting Raleigh's
  cold effect differ, or simply accept the shared slope with ridge shrinkage.
- **Shenandoah shows no mid_greenup effect** (+0.04), while the three NC regions show +0.3 to +0.5.

## 6. Data problems and caveats

- **The two sensors disagree on the trend** (see §1). This is the biggest issue for the "10 years out" question.
- **The between-region ordering does not follow latitude or elevation.** Shenandoah, the northernmost region, is the latest. Raleigh
  (100 m elevation) peaks on the same date as Boone (930 m). Region intercepts are therefore needed. Do not expect the static features to
  explain them, and with only 4 regions they cannot be estimated anyway. This also reinforces that mid-greendown is a *proxy*
  that falls about 1–3 weeks before reported visual peak color.
- **Raleigh is a weak signal region.** It is 11% deciduous, has 2,905 forest pixels, has the lowest MODIS–VIIRS agreement (r = 0.21),
  has the only significant trend, and has the most extreme value (2025).
- **`qa_retention` is ≥ 0.9986 in every row.** The "good-quality" filter removes almost nothing, so it is useless as a weight.
  `peak_spread_days` (5.9–11.6) may serve better as a weight.
- 15% of `peak_doy` values are exact integers, a sign of integer-day pixel dates. This is harmless.
- **`pdsi_AugSep` is usable.** All 100 target rows are populated, so the earlier pull bug does not affect this file.
- **2026 weather runs only to Nov 11** (`weather_daily.csv`). The 2026 `first_frost_doy = 335` "censored" values
  mean only "no frost by Nov 11", not by Nov 30.
- The naive p-values in §3 overstate significance, because the effective n is closer to 25 years than to 100 rows. Use the
  year-mean r or year-clustered SEs when you quote significance.

## Figures (`analysis/figures/`, 200 dpi)

| File | Caption |
|---|---|
| `01_peak_by_year.png` | Peak-color proxy date by year, one panel per region, with the VIIRS overlay, linear trend, and Helene 2024 circled. |
| `01_modis_vs_viirs.png` | MODIS vs VIIRS mid-greendown dates, 2013–2025, against the 1:1 line. VIIRS is consistently about 8 days later. |
| `03_target_corr.png` | Within-region anomaly correlation of each year-varying feature with the peak-date anomaly, sorted by \|r\| (n = 100). |
| `03_target_corr_noHelene.png` | The same, excluding the two Helene 2024 rows (n = 98). |
| `03_feature_corr_matrix.png` | Feature-by-feature correlation of the anomalies. It shows the September–October coolness cluster. |
| `05_feature_consistency.png` | Feature anomaly vs peak-date anomaly by region, with per-region fits and a pooled fit, for the five shortlisted features. |
| `05_slope_forest.png` | Per-region slopes with 95% CIs, showing that the regions agree except for Raleigh's cold terms. |
