# Environmental analysis: how location shapes the timing of peak fall foliage

**Owner:** Environmental (Arin). **Data:** `data/panel_v1` (11 locations × 2001–2025, MODIS mid-greendown `peak_doy`; `pixels.csv` for within-location work).
**Scripts** (run in order from this folder with the `fall-foliage` env; total < 2 min): `01_location_effect_size.py` → `02_spatial_predictors.py` → `03_mu_models.py` → `04_location_sensitivity.py` → `05_within_location_elevation.py` → `06_synthesis.py`. Figures in `figures/` (E01–E24, 200 dpi), tables in `tables/`.
**Role in the model** (`plan.md` §1–2): explain and generalize **μ(location)**, test **interactions** for δ, and refit the **within-location elevation rule**.

---

## Headline results

1. **Location matters as much as weather.** Location explains **45%** of all variation in peak date, year-to-year weather shared across the region explains 36%, and the rest is 19% (two-way ANOVA; location F = 56.7, p ≈ 10⁻⁵⁷; this p treats location-years as independent, so read it as "overwhelming", not as an exact value). With year-pairing, **45 of 55** location pairs differ significantly (Holm-adjusted). Typical peak runs from **Sep 28** (Montpelier, Old Forge) to **Oct 10** (Shenandoah). *(E01–E04)*
2. **Location acts through climate.** The best single descriptor of a site's typical peak is a **thermal date**: the day its climatological mean temperature falls below 15 °C. It gives r = 0.91, R² = 0.82, permutation p = 0.0002, and leave-one-location-out (LOLO) MAE **1.33 d** vs 3.32 d for the grand mean. *Caveat:* the 15 °C threshold and the descriptor were chosen after scanning, so 1.33 d is somewhat optimistic. The genuinely pre-specified model (latitude + elevation) scores **1.50 d**, and choosing the single descriptor inside each fold gives 1.54 d. The result is **not threshold-specific**: r = 0.87–0.91 for every threshold from 8 to 18 °C. Peak comes about **0.33 d later for each day later the site cools**. *(E05–E07, E12)*
3. **Latitude and elevation work together, and agree with Hopkins' bioclimatic law.** Alone, elevation shows nothing (r = 0.21) because the southern sites are also the high ones. Holding latitude fixed: **−1.55 d per °N (p = 0.001) and −1.10 d per +100 m (p = 0.024)**, R² = 0.79. So **1° of latitude ≈ 141 m of elevation**, a *ratio* close to Hopkins' ≈ 122 m (400 ft) per degree. The absolute rates differ: Hopkins' spring rule is ~4 d per degree, ours is 1.55 d per degree for the autumn MODIS proxy. Both effects vanish once the thermal date is in the model, which is *consistent with* (not proof of) mediation through temperature; with n = 11 and r ≈ 0.9 collinearity this cannot be separated statistically. The September-temperature slope × the standard lapse rate predicts −1.29 d per 100 m, close to the −1.10 observed. *(E06, E09, E10, E21, E24)*
4. **Searching for the "best" μ model overfits with n = 11.** Choosing among 136 one- and two-predictor models *inside each fold* gives an honest LOLO MAE of **2.14 d**, worse than either *a priori* model (thermal date 1.33 d; lat + elev 1.50 d). **Recommendation: fix the μ predictors on physical grounds and don't select them by search.** *(E11)*
5. **No evidence that location changes how a site responds to weather** (absence of evidence with ~25 effective years, not proof of identical slopes). All equal-slope F-tests have p ≥ 0.57. None of the 30 interaction tests (6 weather anomalies × 5 location moderators) reaches p < 0.05, even before correction (smallest p = 0.19). The best interaction changes LOYO MAE by only −0.02 d. **Use shared slopes in δ; location belongs entirely in μ.** *(E15–E17)*
6. **Weather anomalies are regional.** Year-to-year anomaly correlation decays with distance, with an e-folding distance of **948 ± 52 km** (Mantel r = −0.84, p = 0.0001): Boone–Asheville r = 0.96, Montpelier–Asheville r ≈ 0.2. Northern, central and southern blocks are visible. *(E13, E14)*
7. **Within a location, elevation matters only on high ground; aspect matters everywhere.**
   - **NC Blue Ridge:** break at **1,300 m** (95% CI 1,250–1,300), flat below (+0.28 d/100 m) and **−2.92 d/100 m above** (CI −3.27 to −2.48). This **replicates the earlier 4-region rule** (1,350 m, −2.93) on new pixels.
   - **North and central sites have almost no forest that high.** Their best group breaks (500 m, 700 m) improve fit by only about 1%, and their upper slopes are about −0.7 d/100 m. Per location, Montpelier has a clear steep segment above **600 m** (−2.5 d/100 m). Other "breaks" are weak or sit at the edge of the search range, and there is no reliable link between break elevation and latitude (r = −0.41, n = 7).
   - **Below the break, within-location slopes are small and change sign between sites** (−0.9 to +1.7 d/100 m).
   - **North-facing slopes peak earlier at all 11 sites**, by 1.4–4.6 d per unit northness, with every 95% CI excluding 0. *(E18–E22)*
8. **What climate doesn't explain:** **Shenandoah peaks +3.7 d later** than its climate predicts (in-sample residual; +4.4 d when Shenandoah is held out of the fit, E12), and **Roanoke 2.5 d earlier**. Every other site is within ±1.3 d. Shenandoah's region is a narrow 1,215 km² ridge strip with 79% deciduous forest; a local forest-type or proxy effect is plausible but untested. *(E12, E23)*

## Recommendations for the combined model

| Use | Recommendation |
|---|---|
| **μ for the 11 sites** | Each site's own historical mean (or the trailing window chosen in the `plan.md` §2 backtest). A site's own mean beats any spatial model for that site. |
| **μ for a new place** (maps, "throughout your area") | `μ ≈ 345.4 − 1.55·lat − 1.10·(elev/100 m)` (DOY; E24) inside the sampled range (35.6–44.3 °N, ~390–960 m median forest elevation), or the thermal-date model if daily temperatures are available: `μ ≈ 0.331 × (DOY when climatological Tmean < 15 °C) + const` (table `E09`). Expect about ±1.5 d error (LOLO). |
| **Interactions in δ** | **None.** No location moderator earns a place under the plan's inclusion rule (§3c). This frees the 2-interaction budget. |
| **Elevation within a location** | Southern Blue Ridge: keep the two-piece rule (flat to 1,300 m, then −2.9 d per 100 m). Northern sites: a steep segment starts lower (Montpelier ~600 m, −2.5 d per 100 m), but only Montpelier shows it clearly, so use it only as a caveat. Elsewhere, don't adjust for elevation below ~1,000 m. |
| **Aspect** | For the editor letter: north-facing slopes color about 2–4 days earlier than south-facing slopes nearby. |

## Hand-back table (DATASET_GUIDE §7)

| feature | role | effect size (± SE) | forecast-safe? | evidence | notes |
|---|---|---|---|---|---|
| thermal date (Tmean < 15 °C, climatology) | explains μ | +0.331 d per d (± 0.051) | yes (static) | r = 0.91, perm p = 0.0002, LOLO MAE 1.33 vs 3.32 | robust over 8–18 °C; new column, see formula below |
| latitude | explains μ | −1.55 d/°N (± 0.29), with elevation | yes | p = 0.001; R² 0.79 with elev; LOLO 1.50 | −1.00 d/°N alone (p = 0.006) |
| median forest elevation | explains μ | −1.10 d/100 m (± 0.40), with latitude | yes | p = 0.024; masked when alone (r = 0.21) | 1° lat ≈ 141 m (Hopkins ≈ 122 m) |
| `norm_tmean_Sep` | explains μ | +1.99 d/°C (± 0.33) | yes | r = 0.89, LOLO 1.47 | collinear with thermal date (r = 0.99) |
| `evt_northern_hardwood` | explains μ (weak) | −11.3 d per unit share (± 2.8) | yes | r = −0.80, LOLO 1.93 | collinear with temperature (r = −0.92); not separable with n = 11 |
| longitude, slope, NLCD shares, area, northness_mean | explains μ | — | yes | LOLO ≥ 2.4 d or worse than grand mean | did not help beyond latitude/temperature |
| weather × {lat, elev, Sep normal, northern hardwood, oak–hickory} | interaction with δ | all \|z\| ≤ 1.3 | — | 30 tests, min p = 0.19; best LOYO change −0.02 d | **rejected**: use shared slopes |
| pixel elevation, South (NC) | within-location | flat < 1,300 m; −2.92 d/100 m above (CI −3.27, −2.48) | yes | year bootstrap, break re-chosen | replicates old rule |
| pixel elevation, North/Central | within-location | ≈ −0.7 d/100 m above 500–700 m; Montpelier −2.5 above 600 m | yes | fit gain ≈ 1% (weak) | caveat only |
| pixel northness | within-location | −1.4 to −4.6 d per unit, all 11 sites | yes | every CI excludes 0 | robust |

**New column formula.** `thermal_date_T` = for each location, average daily `tmean` from `weather_daily.csv` over 2001–2025 by day of year, smooth it with a centered 15-day rolling mean (minimum 8 days), and take the first DOY ≥ 182 where it falls below T °C. Values: `tables/E03_thermal_dates.csv`.

**Data problems found.** (1) In the pixel data, a few retrieval failures give anomalies above 60 d; I dropped them (< 0.01% of pixel-years). (2) Binghamton's within-site elevation slope is strongly *positive* (+1.7 d/100 m) over a narrow 330–530 m range. That probably reflects valley-bottom developed land and floodplains rather than elevation. (3) Shenandoah's +3.7 d residual is the main thing the climate model can't explain. Check it against the news-reported peaks.

**Side finding for the temporal teammate:** `gdd5_MarMay` has a strong pooled slope of **−0.025 d per °C·d (p < 10⁻⁵)**, with equal slopes at every site (p = 0.99): warm springs go with *earlier* fall peaks. `mid_greenup_p50` is +0.40 d/d (p = 0.0005).

## Methods notes

- **Unit of analysis for μ:** the 11 location means, so n = 11. Tests use the t-distribution (df = 9) and 10,000-shuffle permutation tests; slope CIs come from a location bootstrap. LOLO = refit without the site and predict it. WLS (weights 1/SE² of each μ) and excluding Shenandoah leave the conclusions unchanged (table `E12`).
- **Interactions:** pooled within-location anomaly regressions with **year-clustered SE**; moderators are standardized across sites; Holm correction over 30 tests; LOYO with Ridge (α = 1) on standardized anomalies. The base δ model (`tmax_pre60to31`, `precip_pre60`, `mid_greenup_p50`) is provisional; its LOYO MAE is 2.96 d vs 3.34 d for climatology.
- **Pixels:** anomaly = pixel DOY − location-year median. Regressions include location-year fixed effects. CIs come from a block bootstrap over years (1,000 reps; 500 for breaks, with the break re-selected in every replicate). Break grid: 400–1,650 m in 50 m steps, requiring at least 1,000 pixel-years above the break.

## Figure list (for the report)

| # | File | Caption |
|---|---|---|
| E01 | `E01_map_typical_peak.png` | Map of the 11 locations colored by mean peak date. |
| E02 | `E02_peak_by_location.png` | Every year's peak date by location, with mean ± 95% CI. |
| E03 | `E03_variance_decomposition.png` | Location 45%, shared-year weather 36%, residual 19%. |
| E04 | `E04_pairwise_differences.png` | Year-paired differences between all location pairs (Holm stars). |
| E05 | `E05_univariate_predictors.png` | 20 location descriptors ranked by LOLO error. |
| E06 | `E06_mu_vs_descriptors.png` | μ vs latitude, elevation, longitude, Sep normal, northern hardwood, thermal date. |
| E07 | `E07_thermal_threshold_scan.png` | Thermal-date result is stable across 8–18 °C thresholds. |
| E08 | `E08_predictor_correlations.png` | Collinearity among location descriptors. |
| E09 | `E09_coefficient_ladder.png` | Latitude/elevation coefficients across nested models (mediation). |
| E10 | `E10_partial_regression.png` | Added-variable plots for μ ~ latitude + elevation. |
| E11 | `E11_honest_mu_error.png` | Nested-LOLO error: a priori models beat model search. |
| E12 | `E12_thermal_model_fit_residuals.png` | Thermal model fit and held-out residuals (Shenandoah +4.4 d). |
| E13 | `E13_anomaly_coherence_vs_distance.png` | Anomaly correlation decays over ~950 km. |
| E14 | `E14_anomaly_correlation_matrix.png` | Anomaly correlation matrix between locations. |
| E15 | `E15_per_location_slopes.png` | Per-location weather slopes: no location differs. |
| E16 | `E16_interaction_heatmap.png` | 30 interaction tests: none significant. |
| E17 | `E17_loyo_interaction_gain.png` | No interaction improves LOYO error. |
| E18 | `E18_within_location_binned.png` | Pixel anomaly vs elevation, per location. |
| E19 | `E19_within_location_slopes.png` | Within-location elevation and aspect slopes (year bootstrap). |
| E20 | `E20_breakpoint_by_region.png` | Two-piece breakpoint by region group (South: 1,300 m, −2.9 d/100 m). |
| E21 | `E21_elevation_scales.png` | Elevation effect between vs within locations vs lapse-rate prediction. |
| E22 | `E22_breakpoint_by_location.png` | Per-location breakpoints. |
| E23 | `E23_mu_decomposition.png` | Each location's mean split into climate part + local residual. |
| E24 | `E24_lat_elev_lookup.png` | Predicted typical peak over latitude × elevation, with uncertainty. |

Full one-sentence captions for every figure: `figures/captions.tsv`.
