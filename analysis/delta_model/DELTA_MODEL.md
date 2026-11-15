# The combined model: rolling mean μ + weather-driven δ

**Script:** `analysis/delta_model/delta_model.py` (fall-foliage env, ~1 min). **Data:** `data/panel_v1` (11 locations, 2001–2025 MODIS mid-greendown target; 2026 forecast rows). **Inputs combined:** the dataset's forecast-safe weather features, the environmental analysis (`analysis/environmental/FINDINGS.md`), Felix's FIA species composition (`tree_species_analysis/`) and Sebastian's PhenoCam staged models (`multi_region_phenology/`).

## 1. The model

```
predicted peak(l, 2026) = μ(l) + δ(l, 2026)

μ(l)   = mean MODIS peak DOY at location l over the last 3 years (2023–2025)
δ(l,y) = Σ_k w_k · ( x_k(l,y) − x̄_k(l) )      x̄_k(l) = location l's 2001–2025 mean of feature k
```

**δ is the departure from the location's rolling normal, not from last year's date.** It is driven by how unusual this year's weather is *at that location*.

### δ weights (fit on all 275 location-years; `tables/D7_delta_weights.csv`)

| Weather anomaly (forecast-safe) | Weight | Per 1 typical (1-SD) anomaly | 95% CI (year bootstrap) | Share of δ variance |
|---|---|---|---|---|
| Spring warmth, GDD5 Mar 1–May 31 | **−0.039 d per °C·d** | −2.6 d | −0.055 to −0.025 per °C·d | 47% |
| Cold degree-days (base 20 °C), 30 days before typical peak | **−0.055 d per °C·d** | −1.8 d | −0.074 to −0.034 | 22% |
| August mean daily max temperature | **+1.33 d per °C** | +1.6 d | +0.73 to +1.86 | 19% |
| Spring mean temperature, Mar–Apr | **+0.67 d per °C** | +0.9 d | +0.09 to +1.17 | 6% |
| Early-summer max temperature, Jun–Jul | **−0.83 d per °C** | −0.9 d | −1.41 to +0.02 (crosses 0) | 6% |

**What it means physically:**
- **A warm spring brings fall color earlier.** Leaves that start early also finish early.
- **A cool late September** (more cold degree-days before the typical peak) **brings color earlier.** The chill triggers senescence.
- **A hot August delays color.** The forest keeps growing and the season is extended. This is the same +1.1–1.3 d/°C signal seen in the 4-region analysis.

**Transparency notes:**
- **The two spring terms partly offset each other.** Their anomalies correlate at r = 0.62. Holding Mar–Apr fixed, the GDD term effectively measures **May warmth**, which matters most.
- **A simpler 3-feature core** (spring GDD5, pre-peak cold degree-days, August heat) loses only about 0.1 d of accuracy. Use it if a smaller model is preferred.
- **The June–July term's CI includes zero.**
- **Location-specific slopes are not used.** No location × weather interaction improved held-out error. That held in the environmental analysis (30 tests), in Sebastian's M4 (location-varying sensitivity *worsened* MAE, 4.22 → 4.48 d), and for FIA species shares × spring warmth (best gain 0.05 d, below the 0.05 d bar; `tables/D8`).

**How features were chosen:**
- **Eligibility:** only features whose window ends before the typical peak at *every* location (`feature_timing.csv`) and that exist for 2026.
- **Excluded:** `mid_greenup_p50`, which is predictive but has no 2026 value because MODIS hasn't released it, and all September/October fixed windows, which overlap peak in the north.
- **Selection rule:** forward selection on leave-one-year-out MAE, adding a feature only if MAE improves by ≥ 0.02 d, with at most 5 features (`plan.md` budget).

## 2. Is it any good? (all numbers out-of-sample)

| Test | μ only | μ + δ | Improvement |
|---|---|---|---|
| **Nested leave-one-year-out** (feature selection redone without the test year), 2001–2025, n = 275 | 3.34 d | **2.68 d** | **−20%** |
| **Operational backtest** (μ, feature choice and weights use only years before the forecast year), 2011–2025, n = 165, 3-yr μ | 3.15 d | **2.83 d** | −10% |
| Same backtest with an all-years μ | 3.25 d | 3.09 d | −5% |

**Which rolling window?** The backtest compared 3, 5 and 10 years and all prior years (`tables/D5`, figure D2). **3 years won** (2.83 d vs 3.07–3.09 d). A short window follows the recent drift toward earlier MODIS dates. The trade-off is more sampling noise (≈ 4/√3 ≈ 2.3 d), and the MODIS trend itself is uncertain (VIIRS disagrees).

**Comparison with Sebastian's PhenoCam model:** M2 (+2.75 d/°C for Sep–Oct temperature) reaches 14.8% skill over its baseline, but its Sep–Oct window overlaps the peak, so it is explanatory, not a forecast. The forecast-safe analogue in our data, mean temperature in the 30 days before the typical peak, gives about **+0.9 d/°C** (environmental table E15). Both analyses agree on the **direction** (warmer → later), on **shared slopes**, and that **timing, not curve shape, is what weather predicts** (M5).

## 3. "With what certainty?" — prediction intervals

These come from the 165 genuinely out-of-sample backtest errors (`tables/D6_uncertainty.csv`):

| Statement | Value |
|---|---|
| Typical error (MAE) / RMSE | 2.8 d / 3.5 d |
| **80% prediction interval** | forecast **−5.3 d to +2.9 d** |
| 90% interval | −6.3 d to +5.2 d |
| Chance the actual MODIS date is within ±2 / ±3 / ±5 / ±7 days | **42% / 61% / 82% / 96%** |
| Bias | −0.7 d (recent years came earlier than predicted, so the interval is asymmetric) |
| Calibration check | the 80% interval built only from 2011–2018 errors covered **96%** of 2019–2025 cases, so it is conservative |

**Worked example (Asheville 2026):** rolling mean μ = **Oct 6**, δ = **−2.1 d**, so the forecast is **Oct 4**. There is an **80% chance the MODIS proxy date falls between Sep 29 and Oct 7**, and about a **60% chance it is within 3 days of Oct 4**.

## 4. 2026 forecast (MODIS proxy date; `tables/D10_forecast_2026.csv`, figures D4–D5)

| Location | μ (2023–25 mean) | δ | **Forecast** | 80% interval |
|---|---|---|---|---|
| Montpelier | Sep 26 | −0.3 | **Sep 25** | Sep 20 – Sep 28 |
| Old Forge | Sep 27 | −0.5 | **Sep 26** | Sep 21 – Sep 29 |
| Binghamton | Sep 28 | −1.6 | **Sep 27** | Sep 21 – Sep 30 |
| Scranton | Oct 4 | −0.8 | **Oct 3** | Sep 27 – Oct 6 |
| State College | Oct 1 | −1.4 | **Sep 29** | Sep 24 – Oct 2 |
| Frostburg | Sep 29 | −2.0 | **Sep 27** | Sep 22 – Sep 30 |
| Valley Head | Sep 29 | −0.4 | **Sep 29** | Sep 24 – Oct 2 |
| Shenandoah | Oct 9 | −3.1 | **Oct 6** | Sep 30 – Oct 9 |
| Roanoke | Oct 5 | −1.7 | **Oct 3** | Sep 28 – Oct 6 |
| Boone | Oct 4 | −1.2 | **Oct 2** | Sep 27 – Oct 5 |
| Asheville | Oct 6 | −2.1 | **Oct 4** | Sep 29 – Oct 7 |

**Sensitivity to the μ window.** The 3-year μ follows the recent drift toward earlier MODIS dates. Using the all-years (2001–2025) mean instead, with the same δ, every forecast moves **0.5–3.1 days later**:

| Location | 3-yr μ + δ (main) | All-years μ + δ | Difference (d) |
|---|---|---|---|
| Montpelier | **Sep 25** | Sep 28 | -2.6 |
| Old Forge | **Sep 26** | Sep 27 | -1.2 |
| Binghamton | **Sep 27** | Sep 29 | -2.4 |
| Scranton | **Oct 3** | Oct 4 | -0.8 |
| State College | **Sep 29** | Oct 3 | -3.1 |
| Frostburg | **Sep 27** | Sep 30 | -2.5 |
| Valley Head | **Sep 29** | Oct 1 | -2.4 |
| Shenandoah | **Oct 6** | Oct 7 | -1.0 |
| Roanoke | **Oct 3** | Oct 4 | -0.5 |
| Boone | **Oct 2** | Oct 4 | -1.5 |
| Asheville | **Oct 4** | Oct 5 | -1.3 |

The backtest prefers the 3-year window (2.83 vs 3.09 d), but VIIRS does not confirm the MODIS drift (MODIS–VIIRS gap +7.3 d/decade on these 11 sites). **Report both rows; treat the spread as part of the uncertainty.**

**Every location gets a negative δ (an earlier-than-normal 2026).** The main driver is a very warm 2026 spring, especially May (GDD5 term). A cool pre-peak period and warm Mar–Apr pushed partly the other way (figure D5). δ is small (−0.3 to −3.1 d) compared with the ±4–5 d interval, so **2026 is forecast to be near normal to slightly early everywhere**.

**Important for the news comparison:**
- **These are satellite proxy dates.** Mid-greendown is when canopy greenness is halfway through its autumn decline. People usually report "peak color" **about 1–3 weeks later**.
- **That offset has not yet been calibrated** against news reports (`plan.md` §5). Until it is, don't present these as dates for peak color as people see it.

## 5. Where each analysis enters

| Analysis | Contribution to the final model |
|---|---|
| **Environmental (Arin)** | μ is explained by climate (82% of between-site variance via the thermal date; latitude + elevation follow Hopkins' law). No location × weather interactions, so δ uses shared slopes. Within-site elevation rule (NC: flat to 1,300 m, then −2.9 d per 100 m) adjusts the forecast for specific summits. |
| **Species (Felix)** | Composition explains μ: beech, hard maple and evergreen shares go with earlier sites (r ≈ −0.66 to −0.74); oak, yellow-poplar and tupelo with later ones (r ≈ +0.68 to +0.76; `tables/D9`). It is collinear with climate (n = 11) and does not improve δ (`tables/D8`). Felix's ridge model predicts composition from environment; several categories have negative LOOCV R², so use observed FIA shares. |
| **Temporal/PhenoCam (Sebastian)** | An independent dataset (12 PhenoCam sites) confirms the structure: the location baseline is the main term, a shared temperature shift adds about 15% skill, location-varying slopes hurt, and only timing (not curve shape) is predictable. |

## 6. Limits to state in the paper

- The target is a **satellite proxy**. The visual-peak offset is not yet calibrated.
- **There are effectively about 25 independent years.** The 11 locations share each year's weather; anomaly correlation between sites decays over about 950 km.
- **The skill gain is modest (~10–20%) and the interval is wide (80%: −5 to +3 d).** Weather explains part of the year-to-year variation, not all of it.
- **The 3-year μ follows recent drift,** but the MODIS trend is not confirmed by VIIRS. Re-estimate μ every year.
- **The 2026 inputs use weather through Nov 11, 2026.** All selected windows end before each site's typical peak, so the forecast could have been issued by roughly mid-September. The δ weights, however, were fit on data that includes later years.

## Figures

- **D1** `D1_delta_weights.png`: the δ weights, per 1-SD anomaly, with 95% CIs.
- **D2** `D2_rolling_window_backtest.png`: rolling-window backtest; the 3-yr μ + δ is best.
- **D3** `D3_backtest_accuracy.png`: out-of-sample predicted vs observed, and the error distribution with the 80% band.
- **D4** `D4_forecast_2026.png`: the 2026 forecast with 80% intervals and the rolling mean.
- **D5** `D5_delta_breakdown_2026.png`: what drives each site's 2026 δ.
