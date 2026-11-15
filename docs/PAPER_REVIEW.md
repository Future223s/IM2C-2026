# PAPER_REVIEW: accuracy and transparency audit and a proposed paper story

Reviewed 2026-11-15, about 1 h before submission. Scope: `docs/plan.md`, `data/panel_v1/*.md`, `analysis/FINDINGS.md`,
`analysis/elevation/ELEVATION_FINDINGS.md`, `analysis/environmental/FINDINGS.md` + `tables/`, `tree_species_analysis/`,
`multi_region_phenology/` reports. **Update (after review): `analysis/delta_model/DELTA_MODEL.md` now exists.** It reports nested-LOYO 2.68 d and backtest 2.83 d as the honest numbers, excludes `mid_greenup_p50` (no 2026 value), adds an all-years-μ sensitivity row (forecasts 0.5–3.1 d later), and flags the opposing spring terms. `analysis/environmental/FINDINGS.md` was also revised: Hopkins is now a ratio-only match, the thermal-date 1.33 d is marked as post-scan (1.50 d is the a-priori figure), the two Shenandoah residuals (+3.7 in-sample, +4.4 held-out) are labeled, and the interaction result is worded as absence of evidence. Its script and tables
(`delta_model.py`, `tables/D1–D10`, `figures/D1–D5`) did exist, and I audited them in §2b. Re-check §2b against the write-up when it lands.
Numbers marked "verified" were re-read from the CSV tables.

---

## 0. Model summary (final model: source text for the paper's abstract, data, model, validation and results sections)

Every number in this section was read from a results table. Sources: `analysis/delta_model/` (script `delta_model.py`, tables `D1–D10`, figures `D1–D5`, write-up `DELTA_MODEL.md`), `analysis/environmental/` (tables `E00–E24`, figures `E01–E24`, `FINDINGS.md`), `data/panel_v1/` (`BUILD_LOG.md`, `DATASET_GUIDE.md`), `tree_species_analysis/` and `multi_region_phenology/`. The whole pipeline reruns in about 3 minutes (fall-foliage conda env).

### 0.1 One-paragraph abstract (draft)

We predict the date of peak autumn foliage at 11 Appalachian locations from Montpelier, VT to Asheville, NC. Each prediction is a location baseline μ (the mean of the previous three seasons) plus a weather-driven departure δ, a weighted sum of five pre-season weather anomalies. The target is the MODIS satellite "mid-greendown" date (2001–2025), the day forest greenness is halfway through its autumn decline. Location explains 45% of the variation in peak date and shared year-to-year weather 36%. Differences between sites track climate: the date a site's mean temperature falls below 15 °C explains 82% of them. Year-to-year departures follow spring warmth (warm spring → earlier), late-September chill (cool → earlier) and August heat (hot → later). None of these responses differs detectably between locations or forest types. On years it had never seen, the model's average error was 2.68 days, against 3.34 days for "predict the usual date" (20% better). Its 80% prediction interval runs from 5.3 days before to 2.9 days after the forecast. For 2026 every site is forecast near normal to slightly early (δ = −0.3 to −3.1 d), from Sep 25 in Vermont to Oct 6 at Shenandoah. These are satellite proxy dates. Peak color as people report it is assumed to come 1–3 weeks later, an offset we have not yet calibrated.

### 0.2 Data

| Item | Detail |
|---|---|
| **Locations (11, north → south)** | Montpelier VT, Old Forge NY, Binghamton NY, Scranton PA, State College PA, Frostburg MD, Valley Head WV, Shenandoah VA, Roanoke VA, Boone NC, Asheville NC. New coordinates checked against OpenStreetMap (all within 0.01°). Shenandoah, Boone and Asheville keep the exact points of the earlier 4-region study. |
| **Region per location** | 40 km radius circle clipped to the EPA Level III ecoregion containing the point. Areas range from 1,064 km² (Scranton) to 4,966 km² (Binghamton). No two regions overlap. |
| **Forest mask** | NLCD 2021 deciduous (41) + mixed (43) forest. A 463 m MODIS pixel counts as forest if ≥ 50% of its 30 m NLCD pixels are forest. Forest pixel counts range from 3,215 (Scranton) to 18,346 (Boone). |
| **Target (`peak_doy`)** | MODIS MCD12Q2 v061 `MidGreendown_1`, QA_Overall ≤ 1. Exact median over forest pixels whose center lies in the region (computed from per-region pixel histograms, not Earth Engine's approximate percentile). 2001–2025: 275 location-years, no gaps. Observed range DOY 262–291 (Sep 19 – Oct 18). |
| **Weather** | gridMET daily tmin, tmax, precipitation, solar radiation and VPD, averaged over each clipped region (not forest-masked), Mar 1 – Nov 30, 2001–2026. 2026 runs to **Nov 11, 2026**. PDSI from GRIDMET/DROUGHT. |
| **Robustness target** | VIIRS VNP22Q2 (2013–2025): +9.8 d later than MODIS on average. Within-location anomaly correlation with MODIS r = 0.67. The VIIRS − MODIS gap grows +7.3 d/decade. |
| **Static descriptors** | SRTM elevation, slope and aspect over forest pixels; NLCD cover shares; LANDFIRE EVT (2016) forest-type groups; FIA species composition (Felix: 11 categories, 30-mile circles); gridMET climate normals; day length (Forsythe 1995). |
| **Validation of the build** | Regression test against the earlier 4-region data: weather identical (max difference about 1e-12). Satellite percentiles within 1.0 d (mean 0.29 d); this difference comes from switching approximate → exact percentiles. All range, shape and missingness checks pass (`data/panel_v1/BUILD_LOG.md`). |

### 0.3 Model equation

```
peak(l, y) = μ(l, y) + δ(l, y)

μ(l, y)  = (1/3) Σ_{j=1..3} peak(l, y − j)                          rolling 3-year mean of the site's own record
δ(l, y)  = Σ_k w_k · [ x_k(l, y) − x̄_k(l) ]                         shared weights, within-site anomalies
x̄_k(l)   = mean of feature k at site l over the training years (2001–2025 for the final fit)
```

- One set of weights w_k applies to all locations (no site-specific slopes).
- A positive δ means a later than usual peak.

### 0.4 The baseline μ

**Why a rolling mean, and why 3 years** (`tables/D5`, figure **D2**). Forward backtest: for every year 2011–2025, use only earlier years to build μ, choose the δ features and fit the weights, then predict that year.

| μ window | μ alone, MAE (d) | μ + δ, MAE (d) |
|---|---|---|
| **3 years** | **3.15** | **2.83** |
| 5 years | 3.23 | 3.07 |
| 10 years | 3.30 | 3.08 |
| All prior years | 3.25 | 3.09 |

- **The 3-year window wins** because the MODIS record has drifted earlier recently (pooled 4-region trend −1.3 d/decade, not significant; clustered p = 0.20), and a short window follows that drift.
- **Cost:** more sampling noise (year-to-year SD ≈ 4 d, so the SE of a 3-year mean ≈ 2.3 d).
- **Caveat:** the drift is not confirmed by VIIRS. **Sensitivity:** an all-years μ moves every 2026 forecast 0.5–3.1 d later (§0.9).

**Explaining μ across space** (environmental analysis, n = 11 sites). Needed for places without a long record and for "why is Vermont earlier than North Carolina?":

| Model of the site mean | Coefficients | R² | Leave-one-location-out MAE |
|---|---|---|---|
| Grand mean (benchmark) | — | — | 3.32 d |
| **Latitude + median forest elevation** (pre-specified) | **−1.55 d per °N (p = 0.001); −1.10 d per +100 m (p = 0.024)** | 0.79 | **1.50 d** |
| Thermal date: first DOY the climatological 15-day-smoothed Tmean < 15 °C | +0.33 d per day later cool-down (p = 0.0001) | 0.82 | 1.33 d (optimistic: threshold chosen by scan; r = 0.87–0.91 for any 8–18 °C threshold) |
| September mean temperature normal | +1.99 d per °C | 0.80 | 1.47 d |
| Northern-hardwood share (EVT) | −11.3 d per unit share | 0.64 | 1.93 d |
| Best model searched inside each fold | — | — | 2.14 d (searching overfits with n = 11) |

- **Geography lookup:** `μ ≈ 345.4 − 1.55·lat − 1.10·(elev/100 m)` (DOY; figure **E24**). Valid inside the sampled range: 35.6–44.3 °N, median forest elevation ~390–960 m.
- **Elevation is masked when tested alone** (r = 0.21), because the southern sites are also the high ones.
- **Hopkins ratio:** 1° latitude ≈ 141 m of elevation, close to Hopkins' bioclimatic ratio (≈ 122 m/°). The rates differ: Hopkins' spring rule is ~4 d/°, ours is 1.55 d/° for this autumn proxy.
- **Mediation:** latitude and elevation lose significance once the thermal date is added. That is consistent with acting through temperature, but n = 11 can't prove it.
- **Largest local residuals:** Shenandoah +3.7 d (in-sample; +4.4 d held out) and Roanoke −2.5 d. Every other site is within ±1.4 d (figure **E23**).
- **Species and μ** (`tables/D9`, Felix's FIA shares, n = 11): later-peaking sites have more oak (r = +0.75), tupelo (+0.76) and yellow-poplar (+0.68). Earlier sites have more beech (r = −0.74), hard maple (−0.72) and evergreen (−0.66). Species composition is collinear with climate, so it can't be separated from temperature with 11 sites.

### 0.5 The departure δ

**Eligible inputs: forecast-safe and available for 2026.** A feature qualifies only if its window ends before the typical peak at **every** location (`feature_timing.csv`) and it has a 2026 value. There were 13 candidates:

- Spring: `tmean_MarApr`, `gdd5_MarMay`
- Summer: `tmax_JunJul`, `precip_JunJul`, `tmax_Aug`
- Windows anchored to the typical peak P: `tmean_pre30`, `tmin_pre30`, `cdd_pre30`, `cold_nights_pre30`, `tmax_pre60to31`, `precip_pre60`, `vpd_pre60`, `srad_pre30`

**Excluded, and why:**
- `mid_greenup_p50` (spring green-up date): predictive (+0.40 d/d, p = 0.0005) but **has no 2026 value**, because MODIS MCD12Q2 2026 is not released.
- All September/October fixed windows: they overlap peak at Montpelier and Old Forge (typical peak Sep 28).
- `first_frost_doy`: usually after peak.
- `dormancy_p50`, `peak_doy_viirs` and other target-adjacent columns.

**Estimator:**
- Ridge regression (α = 1) on standardized within-location anomalies of both features and target.
- Forward selection on leave-one-year-out (LOYO) MAE, all 11 sites held out together.
- A feature enters only if MAE drops by ≥ 0.02 d; at most 5 features (`plan.md` budget: ≤ 5 main effects).

**Selection path** (`tables/D1`; in-sample-selection LOYO MAE, optimistic; honest numbers in §0.6):

| Step | Feature added | LOYO MAE |
|---|---|---|
| 0 | none (climatology) | 3.338 d |
| 1 | Spring warmth GDD5, Mar–May | 2.988 d |
| 2 | Cold degree-days, 30 d before peak | 2.708 d |
| 3 | August max temperature | 2.523 d |
| 4 | Mar–Apr mean temperature | 2.495 d |
| 5 | Jun–Jul max temperature | 2.416 d |

**Stability of selection** (`tables/D3`: how often each feature was chosen when selection was redone in each of the 25 nested folds):
- `gdd5_MarMay` 25/25, `tmax_Aug` 25/25, `tmax_JunJul` 22/25, `cdd_pre30` 17/25
- `precip_JunJul` 13/25, `tmean_MarApr` 11/25, `tmean_pre30` 8/25, `cold_nights_pre30` 2/25, `srad_pre30` 2/25
- The spring-warmth and August-heat signals are rock-solid. The others are interchangeable proxies for "late-season chill" and "early-season warmth".

**Final weights** (fit on all 275 location-years; 95% CIs from a 1,000-replicate bootstrap over years; `tables/D7`, figure **D1**):

| Feature (anomaly vs site mean) | Weight (d per unit) | 95% CI | 1 SD of anomaly | Effect of a 1-SD anomaly | Share of δ variance |
|---|---|---|---|---|---|
| Spring warmth GDD5 Mar 1–May 31 (°C·d above 5 °C) | **−0.0393** | −0.0549 to −0.0254 | 66.3 °C·d | **−2.6 d** | 47% |
| Cold degree-days (base 20 °C), P−30…P−1 | **−0.0546** | −0.0737 to −0.0344 | 32.5 °C·d | **−1.8 d** | 22% |
| August mean daily max temperature (°C) | **+1.326** | +0.731 to +1.859 | 1.23 °C | **+1.6 d** | 19% |
| Mar–Apr mean temperature (°C) | +0.668 | +0.087 to +1.174 | 1.41 °C | +0.9 d | 6% |
| Jun–Jul mean daily max temperature (°C) | −0.832 | −1.409 to +0.025 | 1.09 °C | −0.9 d | 6% |

**Physical interpretation:**
1. **Warm spring → earlier autumn.** Earlier leaf-out and an earlier, longer growing season are followed by earlier senescence. Consistent with the carry-over effect of green-up date (+0.40 d per day of later green-up).
2. **Chill before peak → earlier.** Accumulated cold below 20 °C in the month before the typical peak triggers senescence. This is the mechanism of degree-day senescence models (Delpierre et al. 2009).
3. **Hot August → later.** It extends the active season. The same signal appeared in the 4-region study (+1.1 d/°C) and at every one of the 11 sites (pooled +1.0 d/°C for the peak-anchored version, p = 0.036).

**Collinearity disclosure:**
- The two spring terms have correlated anomalies (r = 0.62) and opposite signs. Holding Mar–Apr fixed, the GDD5 term mainly measures **May** warmth. Read them jointly, not separately.
- **Simpler core:** spring GDD5 + pre-peak cold degree-days + August heat. In-sample-selection LOYO MAE 2.52 vs 2.42 d for the five-term model.

**Tested and rejected for δ:**

| Candidate | Result | Source |
|---|---|---|
| Location-specific slopes | Equal-slope F-tests p ≥ 0.57 for all 6 features tested | env `E15` |
| Weather × {latitude, elevation, Sep normal, northern-hardwood share, oak–hickory share} | 30 tests; none has p < 0.05 (min p = 0.19); best LOYO change −0.02 d | env `E16`, `E17` |
| Spring warmth × FIA species shares (evergreen, oak, hard maple, soft maple, beech, yellow-poplar) | best LOYO change −0.048 d (hard maple); none reaches the 0.05 d bar | `D8` |
| Location-varying temperature sensitivity (PhenoCam) | worsens MAE 4.22 → 4.48 d | Sebastian M4 |
| Linear time trend | did not improve LOYO in the 4-region study (3.30 vs 3.22 d) | `analysis/FINDINGS.md` §2 |
| Precipitation and VPD windows | candidates, not selected | `D1`, `D3` |

### 0.6 Validation (all out-of-sample)

**Protocol 1, nested leave-one-year-out** (2001–2025, n = 275). Hold out a year (all 11 sites). Recompute μ from the other years and redo feature selection and weights on the training years only. Predict.

**Protocol 2, operational forward backtest** (2011–2025, n = 165). Predict each year using only earlier years for μ (3-year), feature selection and weights. This mimics real use.

| Metric | Protocol 1: μ only | Protocol 1: μ + δ | Protocol 2: μ only (3-yr) | Protocol 2: μ + δ |
|---|---|---|---|---|
| MAE | 3.34 d | **2.68 d (−20%)** | 3.15 d | **2.83 d (−10%)** |
| RMSE | 4.13 d | **3.29 d** | — | **3.54 d** |
| Bias (observed − predicted) | 0 | −0.20 d | — | −0.73 d |

- **Correlation of predicted δ with the observed departure** (Protocol 1): r = 0.62. δ captures about 38% of year-to-year variance.
- **Helene 2024** (Asheville, Boone): dropping those rows leaves nested MAE unchanged (2.68 d).

**By location** (Protocol 1, MAE in days; the model beats climatology at all 11 sites):

| Location | μ only | μ + δ |
|---|---|---|
| Montpelier | 2.78 | 2.58 |
| Old Forge | 2.51 | 1.96 |
| Binghamton | 3.27 | 2.46 |
| Scranton | 3.18 | 3.07 |
| State College | 3.39 | 2.58 |
| Frostburg | 4.13 | 3.50 |
| Valley Head | 3.55 | 2.79 |
| Shenandoah | 3.45 | 2.73 |
| Roanoke | 3.90 | 3.22 |
| Boone | 3.27 | 2.32 |
| Asheville | 3.28 | 2.32 |

**By year** (Protocol 2, MAE in days across the 11 sites; 3-yr μ alone vs μ + δ):

| Year | μ only | μ + δ |
|---|---|---|
| 2011 | 3.12 | 4.94 |
| 2012 | 2.30 | 2.83 |
| 2013 | 1.94 | 2.36 |
| 2014 | 1.88 | 1.81 |
| 2015 | 2.94 | 2.74 |
| 2016 | 7.24 | **4.69** |
| 2017 | 2.94 | 2.70 |
| 2018 | 2.70 | 4.46 |
| 2019 | 3.03 | 3.40 |
| 2020 | 2.58 | 1.83 |
| 2021 | 4.39 | **2.67** |
| 2022 | 1.82 | 1.72 |
| 2023 | 1.55 | 1.29 |
| 2024 | 5.21 | **2.91** |
| 2025 | 3.61 | 2.16 |

- **δ earns its keep in anomalous years** (2016, 2021, 2024), cutting errors by 2–3 days.
- **It does slightly worse in some normal years.** 2011, 2018 and 2012 are the misses (2011 also had the least training data).

### 0.7 Uncertainty: the "% certainty" statement

Empirical distribution of the 165 Protocol 2 errors (observed − predicted; `tables/D6`, figure **D3**):

| Quantity | Value |
|---|---|
| 50% interval half-width (median \|error\|) | 2.5 d |
| **80% prediction interval** | **forecast − 5.3 d to forecast + 2.9 d** |
| 90% prediction interval | forecast − 6.3 d to forecast + 5.2 d |
| P(\|error\| ≤ 2 d) | 42% |
| **P(\|error\| ≤ 3 d)** | **61%** |
| P(\|error\| ≤ 5 d) | 82% |
| P(\|error\| ≤ 7 d) | 96% |
| Calibration | the 80% interval from 2011–2018 errors covered **96%** of 2019–2025 cases, so it is conservative |

- **Why the interval is asymmetric:** recent years came earlier than predicted (bias −0.73 d), so the interval reaches further on the early side.
- **Suggested wording:** "We forecast Oct 4 (MODIS proxy) for Asheville. There is an 80% chance the date falls between Sep 29 and Oct 7, and roughly a 60% chance it is within 3 days of Oct 4."

### 0.8 Within-location adjustments (for specific summits and slopes; environmental analysis, 1.08 M pixel-years)

- **Southern Blue Ridge (Boone, Asheville):** flat to **1,300 m** (95% CI 1,250–1,300 m; +0.28 d per 100 m below), then **−2.92 d per 100 m above** (95% CI −3.27 to −2.48; year bootstrap with the break re-selected). This replicates the earlier 4-region rule (1,350 m, −2.93). Example: a 1,500 m ridge near Asheville peaks about 5 d earlier than the regional forecast.
- **North and central sites:** little forest above ~800–1,200 m. Montpelier shows a steep segment above ~600 m (−2.5 d per 100 m). Elsewhere, slopes below the break are small and change sign between sites (−0.9 to +1.7 d per 100 m). **Don't adjust below ~1,000 m.**
- **Aspect:** north-facing slopes peak earlier at all 11 sites, by 1.4–4.6 d per unit northness (all CIs exclude 0).

### 0.9 2026 forecast (MODIS proxy dates; `tables/D10`, figures **D4**, **D5**)

The 2026 weather inputs run to Nov 11, 2026. All selected windows end before each site's typical peak, so these forecasts could have been issued by mid-to-late September.

| Location | μ (2023–25) | δ | **Forecast** | 80% interval | δ contributions (d): spring GDD5 / pre-peak chill / Aug heat / Mar–Apr / Jun–Jul | With all-years μ |
|---|---|---|---|---|---|---|
| Montpelier | Sep 26 | −0.3 | **Sep 25** | Sep 20 – Sep 28 | +0.1 / −1.0 / +0.4 / +0.9 / −0.7 | Sep 28 |
| Old Forge | Sep 27 | −0.5 | **Sep 26** | Sep 21 – Sep 29 | −0.1 / −0.7 / −0.5 / +1.0 / −0.3 | Sep 27 |
| Binghamton | Sep 28 | −1.6 | **Sep 27** | Sep 21 – Sep 30 | −2.0 / +0.2 / −0.6 / +1.4 / −0.6 | Sep 29 |
| Scranton | Oct 4 | −0.8 | **Oct 3** | Sep 27 – Oct 6 | −2.0 / +0.6 / −0.1 / +1.3 / −0.6 | Oct 4 |
| State College | Oct 1 | −1.4 | **Sep 29** | Sep 24 – Oct 2 | −3.2 / +1.1 / −0.7 / +1.8 / −0.3 | Oct 3 |
| Frostburg | Sep 29 | −2.0 | **Sep 27** | Sep 22 – Sep 30 | −5.5 / +1.6 / −0.5 / +2.5 / −0.2 | Sep 30 |
| Valley Head | Sep 29 | −0.4 | **Sep 29** | Sep 24 – Oct 2 | −5.5 / +2.6 / +0.9 / +2.3 / −0.7 | Oct 1 |
| Shenandoah | Oct 9 | −3.1 | **Oct 6** | Sep 30 – Oct 9 | −5.5 / +0.4 / +0.3 / +2.1 / −0.4 | Oct 7 |
| Roanoke | Oct 5 | −1.7 | **Oct 3** | Sep 28 – Oct 6 | −5.2 / +1.0 / +1.1 / +1.7 / −0.4 | Oct 4 |
| Boone | Oct 4 | −1.2 | **Oct 2** | Sep 27 – Oct 5 | −6.3 / +2.1 / +1.2 / +2.1 / −0.4 | Oct 4 |
| Asheville | Oct 6 | −2.1 | **Oct 4** | Sep 29 – Oct 7 | −6.6 / +1.8 / +1.1 / +2.0 / −0.3 | Oct 5 |

**Reading the 2026 forecast:**
- **A very warm 2026 spring,** especially May, pulls every site south of New York earlier by 2–7 d.
- **That is partly cancelled** by a warm Mar–Apr (+1–2.5 d) and by the pre-peak chill and August heat terms.
- **Net δ is small** (−0.3 to −3.1 d) compared with the interval, so 2026 is **near normal to slightly early everywhere**, with the largest early shifts at Shenandoah, Asheville and Frostburg.
- **The 3-year μ** is itself 0.5–3.1 d earlier than the long-run mean, because 2023–2025 were early years.

### 0.10 How the three team analyses enter the model

| Analysis | Role | Key result used |
|---|---|---|
| **Environmental (Arin)** `analysis/environmental/` | explains μ; tests interactions; within-site elevation | Location 45% of variance (ANOVA F = 56.7); climate explains 82% of between-site μ; lat + elev geography model (R² 0.79); no location × weather interaction (30 tests); weather anomalies coherent over ~950 km (Mantel p = 0.0001); NC elevation break 1,300 m, −2.9 d per 100 m; north-facing slopes earlier. |
| **Species (Felix)** `tree_species_analysis/` | explains μ (descriptive only) | Oak/tupelo/yellow-poplar shares go with later sites; beech/hard maple/evergreen with earlier sites. No species × weather interaction improves δ. His ridge model (19 predictors → species %, n = 11) has negative LOOCV R² for 7 of 11 categories, so use observed FIA shares, not predicted ones. |
| **Temporal / PhenoCam (Sebastian)** `multi_region_phenology/` | independent confirmation of the structure | 12 PhenoCam sites, separate data: the location baseline is the main term (M1, MAE 4.95 d); a shared temperature shift adds 14.8% skill (M2, +2.75 d/°C for Sep–Oct; explanatory, not forecast-safe); location-varying slopes hurt (M4, 4.48 d); only timing, not curve shape, is predictable from weather (M5). |

**Points of agreement to emphasize:**
- Location dominates.
- One shared weather response everywhere.
- Warmer late season → later peak.
- About 15–20% skill over climatology, independently, in two datasets with different sensors and sites.

### 0.11 Limitations (state in the paper)

1. **Proxy, not visual peak.** MODIS mid-greendown is the midpoint of the greenness decline. The assumed 1–3-week lag to news-reported peak color is **uncalibrated**. Every date in §0.9 is a proxy date.
2. **2026 is not verifiable yet** against the target (MODIS MCD12Q2 2026 is released around mid-2027). The news-report comparison must be stated as done or pending.
3. **Effective sample size ≈ 25 years.** The 11 sites share weather (anomaly correlation up to 0.96 between neighbors), so the 275 rows are not independent. Weight CIs use a bootstrap over years for this reason.
4. **Modest skill.** δ cuts error by 10–20%. Year-to-year SD is ~4 d, and the 80% interval spans ~8 d.
5. **Sensor dependence and trend.** MODIS and VIIRS disagree on the 2013–2025 trend (gap +7.3 d/decade on these sites). The 3-yr μ follows MODIS drift. MODIS is near end of mission, so an operational system needs a MODIS→VIIRS crosswalk.
6. **n = 11 for the spatial claims.** Climate, latitude, species and elevation are collinear across sites. Causal attribution of μ differences is not possible.
7. **Region-level weather.** gridMET averages are over the whole clipped region, not only forest pixels. Region sizes differ (1,064–4,966 km²).
8. **Weights fit with hindsight.** The final weights use all 2001–2025 data. Protocol 2 gives the honest real-time accuracy.

### 0.12 Ten-year outlook (input for the "10 years out" question)

Error budget at a 10-year horizon:
- **Model error** (Protocol 2): RMSE ≈ 3.5 d.
- **Drift in μ:** the trend CI in the earlier analysis spans −3.2 to +0.7 d/decade, up to ~3 d after 10 years if μ is not updated.
- **Climate moving the inputs outside their training range:** e.g. +1 °C August warming × 1.33 d/°C ≈ +1.3 d.
- **Slow change in forest composition** (species analysis: oak vs maple shares; not quantified).

Combining the first three in quadrature gives ≈ 4.8 d RMSE at 10 years if μ is frozen. With μ re-estimated every year the drift term mostly disappears, leaving ≈ 3.5–3.8 d. **Recommendation:** re-estimate μ annually and refit δ every few years.

### 0.13 Figure and table map for the paper

| Paper section | Figures | Tables |
|---|---|---|
| Study area & data | E01 (map), E02 (peak by location) | `BUILD_LOG.md` checks |
| Why location matters | E03 (variance split), E04 (pairwise), E13 (coherence vs distance) | E01, E02 |
| Explaining μ | E05, E06, E09, E10, E11, E12, E23, E24 | E05, E08–E13 |
| δ model | D1 (weights), D2 (μ window), D5 (2026 breakdown) | D1, D3, D7 |
| No location-specific response | E15, E16, E17 (+ Sebastian M4) | E15–E17, D8 |
| Accuracy & uncertainty | D3 (backtest accuracy) | D2, D4, D5, D6 |
| 2026 forecast | D4 | D10 |
| Within-site elevation | E18, E20, E21 (E19, E22 appendix) | E20–E23 |
| Appendix | E07, E08, E14, species LOOCV table, PhenoCam M1–M6 | `ridge_metrics.csv` |

---

## 1. Executive summary (the honest story)

- **What we model is a satellite proxy, not the date a tourist would call "peak".** The target is MODIS MCD12Q2 mid-greendown (the median over the forest pixels
  in a 40 km circle clipped to its ecoregion, 11 locations, 2001–2025). The claim that it falls "1–3 weeks before reported peak" is **not backed by any table in the repo**.
  It is repeated from earlier notes. The proxy→news offset is the largest unquantified error in the whole pipeline.
- **Where you are matters more than which year it is.** Location explains 45% of the variance in peak date, the shared year explains 36%, and 19% is left over
  (verified, `E01_anova_two_way.csv`). Typical proxy dates run from Sep 28 (Montpelier, Old Forge) to Oct 10 (Shenandoah).
- **Between-site differences follow climate.** The site's climatological "cool-down date" (Tmean < 15 °C) explains μ with r = 0.91 and leave-one-location-out (LOLO) MAE 1.33 d,
  against 3.32 d for the grand mean. Latitude + elevation gives LOLO MAE 1.50 d. These are correlational results from n = 11 sites.
- **Year-to-year departures are small and only modestly predictable.** Year-to-year SD is about 4 d. Forecast-safe late-summer heat gives about +1.0 d/°C
  (`tmax_pre60to31`, pooled p = 0.036, `E15`). The delta model gives a **nested** LOYO MAE of **2.68 d vs 3.34 d** for climatology (−20%).
  An operational (forward-in-time, 2011–2025) backtest gives **2.83 d**, RMSE 3.54, bias −0.73 d, and an 80% interval of −5.3 to +2.9 d (80% of absolute errors ≤ 4.5 d)
  (`analysis/delta_model/tables/D2`, `D6`). The δ signal is mostly **spring warmth → earlier** and **late-summer heat → later**.
  PhenoCam (Sebastian) independently finds a similar skill, −15%, but its predictor is **not forecast-safe** (Sep–Oct temperature).
- **The analyses agree on structure.** The response to weather is the same everywhere: no location-varying slopes (Arin: 30 interaction tests with min p = 0.19;
  Sebastian: M4 with site interactions is *worse* than M2). So location belongs in μ, and δ uses one shared slope.
- **There is no defensible trend claim.** In the 4-region data the pooled MODIS trend is −1.3 d/decade (clustered p = 0.20). VIIRS shows +6.8 d/decade over 2013–2025.
  The sensors diverge by about 7.3–7.7 d/decade. That is a sensor problem, not a climate signal, and it dominates the 10-year question.
- **The species model is not a peak-date model.** It predicts FIA species percentages from the environment, and 7 of 11 categories have negative LOOCV R².
  Its role in the paper is descriptive context only.
- **Still to do (must be stated if not finished):** the calibrated offset, the 2026 forecast scored against 2026 news reports, and the error-vs-horizon budget.

---

## 2. Claim audit

Legend: ✅ supported · ⚠️ needs a caveat or rewording · ❌ wrong or unsupported as written.

| # | Claim (quote/paraphrase) | Source | Verdict | Evidence & fix |
|---|---|---|---|---|
| 1 | MODIS mid-greendown "falls about 1–3 weeks before news-reported peak color" | `DATASET_GUIDE.md` §3, `HOW_TO_USE…` l.17, `plan.md` §5, `analysis/FINDINGS.md` l.129 | ❌ unsupported | No offset table or news data exists (`grep` for news/offset finds only prose; `analysis/elevation/offset_table.csv` is the elevation adjustment). From the definition alone (EVI halfway through its autumn decline), the proxy could plausibly land before *or* near visual peak, so even the sign of the offset needs data. Either quantify it from ≥ 6 location-years of news reports or call it an **assumption** in the Assumptions list. |
| 2 | "Location explains 45%, shared-year weather 36%, residual 19%" | env `FINDINGS.md` #1 | ⚠️ | Verified: 3532 / 2829 / 1494 SS. Caveats: (a) the shares depend on which 11 sites were chosen (a wider latitude range would raise the location share); (b) the "year" term holds *anything* common to a year, including MODIS processing/sensor effects, so don't call it "weather". Rename it "shared year effect". |
| 3 | Location F = 56.7, p ≈ 10⁻⁵⁷; 45/55 pairs differ (Holm) | env #1 | ⚠️ | The numbers are verified. Don't headline p ≈ 10⁻⁵⁷: residuals are spatially correlated (E13: anomaly correlation decays over ~950 km), so the p is overstated. Report the effect size (12-day spread of means) instead. |
| 4 | Thermal date: r = 0.91, R² = 0.82, perm p = 0.0002, LOLO MAE 1.33 vs 3.32 | env #2 | ⚠️ | Verified (`E05`, `E11`). But the 15 °C threshold and the thermal-date descriptor were picked *after* scanning thresholds (E06/E07) and 20 descriptors, so calling it "fixed a priori" in `E11` is not accurate. The honest numbers are **1.54 d** (best single predictor chosen inside each fold) and **1.50 d** (lat + elev, which really was named in `plan.md` §2 beforehand). Report 1.33 d alongside those, with the r = 0.87–0.91 robustness over 8–18 °C. |
| 5 | "Location acts through climate"; lat/elev effects "vanish once thermal date is in the model … consistent with mediation" | env #2–3 | ⚠️ | Correlational, n = 11. In `E09` (lat + elev + thermal), every coefficient has p > 0.2 and VIF is huge (thermal date ≈ f(lat, elev)), so "vanish" also reflects collinearity. Reword: "consistent with temperature as the pathway; n = 11 cannot separate the two". |
| 6 | "1° latitude ≈ 141 m elevation, close to Hopkins' ≈ 122 m (400 ft)" | env #3 | ⚠️ | Only the **ratio** matches. Hopkins' law gives ~4 d per °N. Our −1.55 d/°N is roughly 2.5× weaker. Say both. Lat+elev SEs are ±0.29 and ±0.40, so the ratio's CI is wide (roughly 80–400 m). Give it a range, not one number. |
| 7 | September-temperature slope × lapse rate predicts −1.29 d/100 m vs −1.10 observed | env #3 | ✅ | 1.99 d/°C × 0.65 °C/100 m = 1.29. A nice physical consistency check; keep it. |
| 8 | "Searching for the best μ model overfits" (nested LOLO 2.14 d) | env #4 | ✅ | Verified (`E11`). This is a strength, and the paper should show it (E11). |
| 9 | "Location does not change how a site responds to weather"; "Use shared slopes" | env #5 | ⚠️ | Verified: equal-slope p ≥ 0.57 (`E15`), 30 interactions with min p = 0.19 and all p_holm = 1 (`E16`), best LOYO change +0.005 d (worse, `E17`). But this is **absence of evidence** with about 25 years. Reword: "no detectable difference; interactions would need to exceed ~±0.2 d/°C per SD to be seen". The decision (shared slopes) is still right. |
| 10 | Within-location rule: flat to 1,300–1,350 m, then −2.9 d/100 m | `ELEVATION_FINDINGS.md`, env #7 | ⚠️ | Replicated on new pixels (good). But it rests on 325 pixels above 1,350 m in 2 NC regions, and `frac_forest_above_1350m` ≤ 0.051 at every site. It is a **southern high-peak** rule. Don't apply it north of Virginia. Note that the 4-region single-slope "headline" (−0.18) must not be quoted. |
| 11 | "North-facing slopes peak earlier at all 11 sites, 1.4–4.6 d per unit northness" | env #7 | ⚠️ | Supported for the *proxy*. Possible artifact: north slopes are illuminated differently in the satellite view (BRDF/terrain shading). The letter says "north-facing slopes color 2–4 days earlier". Soften to "satellite data suggest". |
| 12 | Shenandoah peaks +3.7 d later than climate predicts | env #8 | ⚠️ | Inconsistent number: the E12 caption says **+4.4 d** (held-out), and the text says +3.7 (in-sample). Use one number and label which. |
| 13 | μ for a new place: `345.4 − 1.55·lat − 1.10·(elev/100)` ±1.5 d | env recs | ⚠️ | It forecasts the **proxy**, inside 35.6–44.3 °N and 390–960 m (median forest elevation, not point elevation). State both, and that ±1.5 d is LOLO MAE, not a prediction interval. |
| 14 | Provisional δ LOYO MAE 2.96 vs 3.34 climatology | env methods; `E17`, `E18` | ✅/⚠️ | Verified. Features `tmax_pre60to31`, `precip_pre60` and `mid_greenup_p50` are pre-peak. Caveat: `mid_greenup_p50` comes from the **same MODIS product**, which is not released for 2026. Is 2026 green-up available (MCD12Q2 2026 is not)? **If not, the 2026 forecast cannot use it.** Refit without it or source green-up from VIIRS/NRT. Must be resolved in DELTA_MODEL.md. |
| 15 | `gdd5_MarMay` −0.025 d/°C·d, p < 10⁻⁵ | env side finding | ⚠️ | Check whether `pooled_p` in `E15` is year-clustered (the methods say interaction SEs are; not stated for E15). The 4-region `FINDINGS` §3 p-values are explicitly *not* clustered. Quote only clustered or year-mean values. |
| 16 | `tmax_Aug` ≈ +1.1 d/°C "very consistent" | `analysis/FINDINGS.md` §5 | ✅ | 4 regions, equal-slope p = 0.98; the 11-site analogue `tmax_pre60to31` is +1.00 ± 0.48. Report the 11-site number with its SE (CI ≈ 0.06–1.94). It is "consistent" but imprecise. |
| 17 | Correlation table p-values (n = 100) | `analysis/FINDINGS.md` §3 | ⚠️ | Naive p-values; the file itself says so. Don't copy them into the paper. Use the year-mean r column. |
| 18 | "No robust trend" (−1.3 d/decade, clustered p = 0.20, CI −3.2 to +0.7) | `analysis/FINDINGS.md` §1, `plan.md` §6 | ⚠️ | Computed on the **old 4 regions** (incl. Raleigh, now dropped). Not recomputed for the 11-site panel (`locations.csv` has `peak_doy_trend_d_per_decade` per site). Either recompute (1 regression) or say "in the 4-region pilot". |
| 19 | VIIRS–MODIS gap grows "+7.7 d/decade" | `DATASET_GUIDE.md` §3, `plan.md` §6 | ⚠️ | On the 11 sites it is **+7.3** (`BUILD_LOG.md` §6); 7.7 is the 4-region value. Use 7.3 for panel_v1. VIIRS bias is +9.8 d and anomaly r = 0.67 (pooled). |
| 20 | "Between-region ordering does not follow latitude or elevation … static features cannot explain" | `analysis/FINDINGS.md` §6 | ❌ (superseded) | With 11 sites, lat + elev gives R² = 0.79 (`E09`). Paper must use the 11-site result and, if the pilot is mentioned, say that 4 regions could not identify it. |
| 21 | "Helene 2024 is not an outlier" | `analysis/FINDINGS.md` §1 | ✅ | z = −1.2 and −1.1. Keep it as a robustness drop-test. |
| 22 | PhenoCam M2: +2.75 d/°C; MAE 4.95 → 4.22 (14.8% skill) | `staged_curve_models/M2…md`, `staged_model_metrics.csv` | ⚠️ | Verified (n = 158 site-years). The report itself says the Sep–Oct window **overlaps the season**, so this is an *explanatory hindcast*, not a forecast. The paper must not cite 14.8% as forecast skill. +2.75 is a full-data refit (in-sample). |
| 23 | PhenoCam R² 0.49 → 0.63 | M1/M2 reports | ⚠️ | R² pools between-site variance (site means alone give 0.49). Quote MAE and skill, not R². |
| 24 | PhenoCam sites are "Appalachian" | `site_selection_and_readiness.md` | ⚠️ | 12 sites include Michigan (UMBS ×2), Indiana (Morgan–Monroe), NH, MA and NY. Only ~5 are Appalachian (GRSM ×2, Shining Rock, Dolly Sods, Oak Ridge ×2). Call it "eastern deciduous forest PhenoCam network". The target is **RCC90 max** (canopy redness, a single camera footprint), which the report itself says is not validated as peak color. |
| 25 | M4 site interactions don't help (MAE 4.48 > 4.22) | M4 report | ✅ | Independent agreement with Arin #9. Strong cross-check. |
| 26 | Species ridge model (19 predictors, λ = 3, n = 11) | `tree_species_analysis/README.md`, `ridge_metrics.csv` | ❌ if used as evidence | Training R² 0.45–0.97 but LOOCV R² < 0 for **7/11** categories (soft maple −2.14, golden hardwood −4.24, ash −2.39, misc deciduous −0.29, hard maple −0.15, sweetgum −0.03, beech −0.01). Only yellow-poplar (0.86) and evergreen (0.28) predict. 19 predictors > n. Three daylight variables are deterministic functions of latitude (Sep daylight SD = 0.018 h, so ridge standardization amplifies noise). Nov climate is irrelevant to a Sep–Oct event. Uses town-center elevation (ft) and 30-mile FIA circles, not the 40 km ecoregion-clipped forest regions. Fix: present only the **observed FIA composition** per site as context. Report the LOOCV failure honestly, as "species mix is too collinear with climate to separate with 11 sites" (consistent with Arin's `evt_northern_hardwood` r = −0.92 with temperature). |
| 27 | `evt_northern_hardwood` −11.3 d per unit share | env hand-back | ⚠️ | r = −0.92 with temperature. Label it "not separable from climate", which the file already does. Never present it as a species *effect* on timing. |
| 28 | "All Stage 7 checks pass" | `BUILD_LOG.md` | ✅ | Verified. Note the tolerance of 1.5 d on the regression test for `peak_doy`, and the switch to exact percentiles (old values carried ±1 d noise). |
| 29 | Forecast-safe rule (window ends before P = round(mean peak)) | `DATASET_GUIDE.md` §5 | ⚠️ | P is the *mean*. In an early year (peak ~P−8) `*_pre30` windows still overlap the actual peak. Disclose this. Strictly forecast-safe means "issue date fixed before the earliest historical peak" (e.g., Sep 15). `dormancy_p50` and the `*_Sep1_Oct15` / `*_SepOct` / frost features must not appear in any forecast model. |
| 30 | `peak_anom` / `peak_doy_mean` use all years | `HOW_TO_USE…` l.41 | ✅ | Correctly flagged. Verify DELTA_MODEL recomputes μ without the held-out year (E17 climatology 3.338 ≈ BUILD_LOG 3.34 suggests yes). |

## 2b. Delta model audit (write-up now in DELTA_MODEL.md)

| # | Item | Verdict | Evidence & fix |
|---|---|---|---|
| D-a | Candidate features = `pre-peak` at **all** sites, available in 2026, frost excluded (`delta_model.py` l.35–38). MODIS phenology is not a candidate. | ✅ | This resolves claim 14: `mid_greenup_p50` is not used. Say so in the paper. |
| D-b | Forward selection path LOYO MAE 3.34 → **2.42** (`D1`) | ⚠️ | That path is selected on the same LOYO it reports, so it is optimistic. Quote the **nested** figure: **2.68 d / RMSE 3.29** (recomputed from `D2_nested_loyo_predictions.csv`, n = 275). Show `D1` only as the selection path. |
| D-c | Selected features: gdd5_MarMay, cdd_pre30, tmax_Aug, tmean_MarApr, tmax_JunJul (`D7`). | ⚠️ | That is 5 features, inside the plan's budget. But **gdd5_MarMay (−0.039 d/°C·d) and tmean_MarApr (+0.67 d/°C) are both "spring warmth" with opposite signs**. That is a collinearity offset, so don't interpret the two separately (report their combined effect). The `tmax_JunJul` 95% CI reaches +0.02, so it is not clearly nonzero. `D3`: gdd5 and tmax_Aug are picked in 25/25 folds, the others are less stable. Report the selection frequency. |
| D-d | `cdd_pre30` is "pre-peak" | ⚠️ | Pre-peak relative to the *mean* peak P. In early years it overlaps the actual peak, and the 2026 values use weather through P−1 (late Sep/early Oct). So the forecast's issue date is ~P−1 for each site. State the issue date per site. A forecast issued, say, Sep 1 would have to drop `cdd_pre30` (and lose ~0.28 d of MAE, `D1`). |
| D-e | μ window: 3-yr trailing mean wins (model MAE 2.83 vs 3.07–3.09; `D5`) | ⚠️ | The k was chosen on the same 2011–2025 backtest. Also, 3-yr μ beating the long mean is what you'd expect **if the late MODIS record drifts early** (VIIRS disagrees, §1). The choice moves 2026 μ 1–3 d earlier than the long mean (e.g., Montpelier 269.0 vs 271.6). Disclose it, and show the 2026 forecast under "all prior" as a sensitivity row. |
| D-f | 2026 forecast: δ < 0 at all 11 sites (−0.3 to −3.1 d); proxy dates Sep 25 (Montpelier) to Oct 6 (Shenandoah), 80% PI ≈ −5/+3 d (`D10`) | ⚠️ | These are **proxy** dates. MCD12Q2 for 2026 will not be released until mid-2027, so the proxy itself cannot be verified now. Only proxy + offset vs news can be scored. The asymmetric PI comes from the −0.73 d bias, so state that. |
| D-g | Species × weather interactions (`D8`) best −0.05 d | ⚠️ | FIA 30-mile-circle percentages are joined to MODIS ecoregion regions (footprint mismatch), and 11 interactions were tried. A −0.05 d gain is noise-level. Report as "no species interaction earns a place", consistent with E16/E17. |
| D-h | Species vs μ correlation (`D9`, beech r = −0.74, hard maple −0.72) | ⚠️ | Correlational, n = 11, and collinear with latitude/temperature. Northern hardwoods sit in the north. Not a causal species effect. |

---

## 3. Cross-analysis consistency

**Where the three analyses agree (lead with this, it is the robustness story):**
- *Shared response to weather.* Arin: equal-slope p ≥ 0.57, 0/30 interactions (`E15–E17`). Sebastian: M4 (site × temperature) is worse than M2 out of sample.
  Pilot: equal-slope p ≥ 0.26 (`analysis/FINDINGS.md` §5). Three datasets, two sensors (satellite and camera), same conclusion.
- *Warm = late.* All three analyses find warmer late-summer/autumn → later senescence (MODIS `tmax_pre60to31` +1.0 d/°C, `tmean_pre30` +0.9 d/°C;
  PhenoCam Sep–Oct Tmean +2.75 d/°C).
- *Timing shifts, shape doesn't.* Sebastian's M5: only peak timing is predictable from environment, not width or amplitude. This justifies modeling one date per site-year.
- *Climatology is hard to beat.* The gain is −11% (MODIS, forecast-safe) and −15% (PhenoCam, hindcast).
- *Species ≈ climate.* The FIA model and the EVT shares both show that composition is collinear with temperature at n = 11.

**Where the numbers differ, and the honest explanation:**

| Quantity | MODIS panel (Arin) | PhenoCam (Sebastian) | Why they differ |
|---|---|---|---|
| Temperature sensitivity | ~+1.0 d/°C (`tmax_pre60to31`, P−60…P−31) | +2.75 d/°C (Sep–Oct Tmean) | Different windows. Sep–Oct is contemporaneous with the event (partly reverse/common-cause: a warm October *is* a late season), so a larger slope is expected. The targets also differ (greenness halfway-down vs redness max), as do the variables (Tmax vs Tmean), the weather source (gridMET regional mean vs Daymet pixel), the sites, and the n. The 2.75 is in-sample. **Don't average them.** Present them as an upper (hindcast) and lower (forecast-safe) bound. |
| Baseline MAE | 3.34 d (LOYO climatology) | 4.95 d (M1) | PhenoCam is a single camera footprint with daily-max noise and gaps (within-site SD ~6.0 d vs ~4 d for MODIS median over thousands of pixels). Averaging over a region removes local noise. MODIS MAE is lower because of the spatial averaging, not because the model is better. |
| Interannual SD | ~4 d | ~6.0 d | Same reason |
| Skill vs climatology | 11% (E17 provisional) / 20% (delta model, nested; issue date ≈ P−1) | 14.8% (not forecast-safe) | Sebastian's skill would shrink if truncated to a pre-season issue date. |
| Spatial unit | 40 km circle ∩ ecoregion, deciduous+mixed forest | camera ROI (~100 m) | Not commensurable; never pool rows. |
| Species unit (Felix) | EVT share of forest in region | FIA 30-mile (48 km) circle, town elevation | Different footprints; only qualitative comparison. |

Paper sentence: *"Two independent data streams (regional satellite greenness and tower-camera redness), at different sites, agree on the sign of the temperature effect and on the absence of site-specific sensitivity. Their magnitudes differ in the way expected from their windows and footprints."*

---

## 4. Proposed paper structure

1. **Summary sheet (1 p).** The problem, the model (`peak = μ + δ`), the 2026 table (proxy + offset ± interval), the accuracy now and 10 years out, the 3 caveats.
2. **Introduction & problem restatement.** What "peak" means. Visual peak vs satellite proxy (define mid-greendown vs RCC max vs news "peak").
3. **Assumptions** (numbered, each with a justification). See §5 below.
4. **Data.** The 11 regions and how they were built (`DATASET_GUIDE.md` §2–3; `BUILD_LOG.md` checks summary). Fig **E01** (map). A table of sites. MODIS vs VIIRS (`analysis/figures/01_modis_vs_viirs.png`, with 11-site stats from `BUILD_LOG.md` §6).
5. **Model overview.** The μ + δ decomposition, why static factors live in μ (`plan.md` §1). Fig **E03** (variance decomposition) motivates it: 45 / 36 / 19.
6. **μ: why places differ.**
   - 6.1 Site means and differences: **E02**.
   - 6.2 Climate explains μ: **E06**, **E07** (threshold robustness), **E12** (fit + held-out residuals).
   - 6.3 Latitude + elevation and Hopkins: **E10**, **E24** (lookup map for "throughout your area"), with the ratio caveat.
   - 6.4 Honest model selection: **E11** (show that search overfits; a strong transparency point).
   - 6.5 Choosing the μ window (trailing k-year mean backtest, `plan.md` §2), *if DELTA_MODEL produced it*.
   - 6.6 Within-site elevation and aspect: **E20** (South break), **E19**, and the caveat that it holds only for the southern high peaks.
7. **δ: this year's departure.** From `analysis/delta_model/` (write-up pending): the forecast-safe candidate rule; `D1` selection path; **D1_delta_weights.png** (coefficients ± CI); **D2_rolling_window_backtest.png** (μ window); **D3_backtest_accuracy.png** (vs climatology); `D6` uncertainty. Then the shared-slope evidence: **E15, E16, E17** + Sebastian M4. Then the spatial coherence of anomalies: **E13** (why one regional weather signal suffices).
8. **Independent check: PhenoCam.** One subsection with the M1–M6 table from `multi_region_phenology/outputs/reports/staged_curve_models/README.md`, the §3 comparison table, and the explicit "hindcast, not forecast" label.
9. **2026 forecast and news comparison.** **D4_forecast_2026.png**, **D5_delta_breakdown_2026.png** (which weather drove 2026 early). Table per location: μ, δ̂, proxy date, offset, predicted visual peak, 80% PI, reported 2026 peak, error. If the news comparison is incomplete, show what exists and say so.
10. **Sensitivity / robustness.** ±1 °C, ±20% precip, drop-a-feature, drop-a-location, no Helene, VIIRS target (`plan.md` §6).
11. **10-year accuracy.** The error budget (§5 below), the sensor-continuity issue, and the recommendation to re-estimate μ annually.
12. **Strengths & weaknesses.**
13. **Letter to the editor (1 p).**
14. **References + appendix:** species composition (Felix) as descriptive context with its LOOCV table; the full figure list (E04, E05, E08, E09, E14, E18, E21–E23); the PhenoCam site list.

Species material belongs in the **appendix or a "what didn't work" paragraph**, not in the main model.

---

## 5. Transparency checklist (must appear in the paper)

**Assumptions**
- [ ] A1. The MODIS mid-greendown median over deciduous+mixed forest within a 40 km ecoregion-clipped circle is a stable proxy for visual peak, offset by a roughly constant number of days per site.
- [ ] A2. Weather averaged over the whole region (not forest-masked) represents the forest's weather.
- [ ] A3. Location effects are fixed over 2001–2026 (stationary μ), apart from the window choice.
- [ ] A4. The response to weather is the same at every site (tested: §3).
- [ ] A5. 2026 weather inputs (gridMET to Nov 11) are final. gridMET is provisional for recent months.

**Limitations and caveats (the exact statements)**
- [ ] "The proxy→reported-peak offset is assumed (1–3 weeks) / estimated from N reports with SD X d." Pick whichever is true. **If it is not calibrated, say so in the summary, the conclusion and the letter.**
- [ ] "The 2026 comparison with news reports is [complete for k of 11 locations / pending]." Don't imply validation that wasn't done.
- [ ] "n = 11 locations; between-site results (μ models) have 9 df; p-values are indicative."
- [ ] "Locations share weather (anomaly correlation e-folding ~950 km), so the effective sample for weather effects is ~25 years, not 275 rows. Significance uses year-clustered SEs."
- [ ] "Sep–Oct predictors (PhenoCam M2–M6; `*_Sep1_Oct15`) overlap the season: explanatory only."
- [ ] "`mid_greenup_p50` 2026: [source / not available → excluded]."
- [ ] "MODIS (Terra/Aqua) is near end of mission. VIIRS is ~+9.8 d later with a gap growing ~7.3 d/decade, so sensor choice changes the trend sign. No trend claim is made."
- [ ] "The elevation rule is calibrated on ~325 pixels above 1,350 m in western NC."
- [ ] "Species composition could not be separated from climate (LOOCV R² < 0 for 7/11 categories)."
- [ ] "Satellite aspect effect may partly reflect viewing geometry."

**10-year accuracy budget (plan.md §6). Combine the terms in quadrature and show error vs horizon:**
- (i) δ model error: operational backtest RMSE 3.54 d (`D6`), or nested LOYO RMSE 3.29 d.
- (ii) μ drift: trend CI −3.2…+0.7 d/decade → up to ~3 d at 10 y. (This CI comes from the 4-region pilot. Recompute or label it.)
- (iii) out-of-range climate: projected +0.3–0.5 °C/decade × ~1.0 d/°C ≈ 0.3–0.5 d (small; say so).
- (iv) sensor transition: MODIS→VIIRS bias ~10 d with drifting gap. **This is the largest 10-year term unless a crosswalk is built.**
- (v) proxy offset uncertainty: SD of the offset (unknown until calibrated).
- (vi) forest composition change: slow; not quantified (species model cannot inform it).
- Headline phrase: "roughly ±4 d today, ±5–6 d in 10 years if μ is re-estimated annually and a sensor crosswalk exists; worse otherwise". Recompute from actual numbers.

---

## 6. One-page editor letter: content guide

**Key messages (plain language)**
1. Peak color moves down the Appalachians like a wave. It starts in northern Vermont and the Adirondacks, then moves through Pennsylvania and Maryland, and reaches Shenandoah and the North Carolina mountains last. The north-to-south gap is about 1.5–2 weeks.
2. Every place has a "usual" week. In any given year peak lands within about ±4 days of it (2 years in 3).
3. A warm spring tends to bring an earlier autumn. A hot August pushes peak a little later (about 1 day per °C). Early cool nights in September bring it forward. Weather moves the whole region together. In 2026 the model leaned **slightly early** (0–3 days) everywhere.
4. Higher ground colors earlier: above about 4,300 ft in the NC/TN mountains, about 3 days earlier per 330 ft of climb. North-facing slopes may color a few days earlier.
5. 2026 table: location → predicted peak-color week (a range, not a day) → what was reported → hit/miss.
6. Story ideas: the "foliage wave" map (E01/E24); why Shenandoah runs late; satellites disagreeing about whether autumn is getting later.

**Do NOT promise**
- A single exact day. Give a week or a ±range.
- "Peak is getting later/earlier because of climate change." Our sensors disagree.
- That species mix or elevation lets you forecast an individual town without a record (±1.5 d is for the proxy, within the sampled range).
- That the satellite date *is* the peak a reader will see. Say "we convert from satellite timing using past reports".
- Accuracy 10 years out without annual updating and a sensor switchover.

---

## 7. Top 10 fixes before submission (ranked by impact; each < 10 min)

1. **Label the proxy offset honestly.** Wherever "1–3 weeks" appears, either cite the calibration (n reports, mean ± SD) or call it an *assumption*. Put the caveat in the summary sheet and the letter.
2. **State the 2026 news comparison status explicitly** (done for k/11 locations, or pending). The problem asks for it. Don't let a reader assume it was done.
3. **Remove or relabel non-forecast-safe numbers.** PhenoCam 14.8% skill and +2.75 d/°C → "explanatory hindcast". Make sure no `*_Sep1_Oct15`, `*_SepOct`, frost, `dormancy_p50` or `senescence_p50` feature is in the 2026 forecast. Check `mid_greenup_p50` 2026 availability.
4. **Delta model numbers:** quote the nested LOYO (2.68 d) and the operational backtest (2.83 d, ±4.5 d 80% PI), not the selection-path 2.42. Combine the two opposite-signed spring terms. State each site's forecast issue date (`cdd_pre30` needs weather up to ~P−1). Add an "all-prior μ" sensitivity row to the 2026 table (§2b D-e).
5. **Demote the species model** to the appendix with its LOOCV table. Delete any sentence implying it predicts peak timing or that composition was shown to affect timing.
6. **Fix inconsistent numbers:** Shenandoah +3.7 vs +4.4 (pick one, label in-sample/held-out); VIIRS gap 7.7 → 7.3 for 11 sites; drop the 4-region claim that static features can't explain μ.
7. **Hopkins comparison:** add "the ratio matches, but the per-degree effect (−1.55 d/°N) is weaker than Hopkins' ~4 d/°N", and give a range for "141 m".
8. **Soften causal language:** "location acts through climate" → "is consistent with"; "does not change response" → "no detectable difference (power-limited)"; "year = weather" → "shared year effect".
9. **Trend:** either recompute the clustered pooled trend on the 11-site panel (one regression) or label −1.3 d/decade as "4-region pilot". Never quote the naive p = 0.017.
10. **Report honest μ errors:** quote 1.50 d (lat + elev, truly a priori) and 1.54 d (nested) alongside 1.33 d. Note that for the 11 known sites, the site's own mean is used, and the spatial model is for new places only.
