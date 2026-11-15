# plan.md — From three analyses to one peak-foliage model

**Problem (IMMC "Peak Peeping"):** predict when peak fall foliage occurs **this year (2026)** at locations across the Appalachians, compare the predictions with news reports, test robustness, and say how accurate the model would be **10 years out**. Deliverables are a technical report and a one-page letter to the chief editor.
**Data:** `data/panel_v1/` (see `DATASET_GUIDE.md`), 11 locations × 2001–2026, target = MODIS mid-greendown DOY (`peak_doy`).
**Model form:** `peak(l, y) = μ(l) + δ(l, y)`. Here μ is the location's typical peak and δ is this year's departure from it.

---

## 1. Who owns what, and where each piece enters the model

| Piece | Owner | What it answers | Where it enters |
|---|---|---|---|
| **Mean μ(l)** | shared (Arin builds) | What is "normal" for each location? | The baseline term |
| **Temporal δ** | Temporal teammate | Which year-to-year weather and phenology signals move peak earlier or later? | **Main effects** in δ |
| **Environmental** | Arin | Why do locations differ (elevation, latitude/day length, climate normals, terrain)? | (a) a **spatial model of μ**; (b) **interactions** that change how strongly a location responds to weather; (c) the **within-location elevation adjustment** |
| **Species** | Species teammate | How does the forest mix (maple/birch vs. oak/hickory vs. cove) affect timing? | (a) the spatial model of μ; (b) **interactions** with weather |

**The key structural point:** at a fixed location, elevation and species mix don't change from year to year, so **in a mean + delta model they are completely absorbed by μ(l) and cannot move δ by themselves.** The environmental and species work therefore pays off in only two places:
1. **Explaining and generalizing μ** across space. This answers "why is Montpelier earlier than Roanoke?" and gives a forecast for places without a 25-year record, which matters for "throughout your chosen area".
2. **Interaction terms** that make δ's sensitivity depend on the location. For example, "a hot August delays peak more in oak-dominated forests than in northern hardwoods", or "cold nights matter more at high elevation".

Every teammate hands back the results table in DATASET_GUIDE §7, labeled as main effect, interaction, or explains-μ.

## 2. The mean μ(l)

The original idea was "the mean over the last couple of years". That should be **tested, not assumed**:
- Year-to-year SD is about 4 days, and the earlier analysis found **no robust trend**. A k-year mean has sampling noise of about 4/√k days: about **2.3 d for k = 3** and about **0.8 d for k = 25**. Short windows add noise unless the climate is clearly shifting.
- **Backtest:** for each year y from 2011 to 2025, predict using a trailing mean of the previous k years, with k ∈ {3, 5, 10, all prior}. Compare MAE across k. Pick the k that wins (or a weighted/exponentially-decayed mean if it is close), and report the comparison. It is a nice, defensible figure for the report.
- **Spatial μ model (environmental + species):** `μ(l) ~ lat (or day length) + elev_median_m + frac_forest_above_1350m + evt_northern_hardwood + norm_tmean_Sep`. With n = 11, use **at most 2–3 predictors**. Validate by **leaving one location out**. Use it to (i) explain the differences between locations in the report, and (ii) make maps and forecasts between the 11 sites.
- **Elevation within a location:** reuse the existing two-piece rule (`analysis/elevation/elevation_adjust.py`): flat below about 1,350 m, then about −2.9 d per 100 m. Re-fit the rule on the 11-location `pixels.csv` if it was built; the northern sites will test whether 1,350 m is a southern-Blue-Ridge artifact.

## 3. The delta δ(l, y)

### 3a. Main model (statistical)
Let `x̃` be the within-location anomaly of each feature, relative to that location's training-year mean.

```
δ(l, y) = Σ_k β_k · x̃_k(l, y)                       ← temporal teammate's features (shared slopes)
        + Σ_(k,j) γ_kj · x̃_k(l, y) · s_j(l)         ← interactions from the environmental/species work (few!)
        + ε(l, y)
```
- **Candidate main effects** (forecast-safe only, i.e. `timing == pre-peak` in `feature_timing.csv`). Starting points from the 4-region analysis are `tmax_Aug` (about +1.1 d/°C, very consistent across regions), `mid_greenup_p50`, `precip_AugSep` and `cdd_Sep`. Their relative-window versions (`tmax_pre60to31`, `precip_pre60`, `cdd_pre30`) **should be preferred**, because a fixed "September" window is post-peak in Vermont.
- **Budget:** the effective sample size is about 26 years, so use **≤ 5 main effects + ≤ 2 interactions**. Use ridge regression, or a mixed model with location random slopes (`statsmodels` MixedLM), so that weak location-specific slopes shrink toward the shared one.
- **Standard errors:** clustered by year.

### 3b. Alternative model (process-based, for the report)
A **cold-degree-day threshold model** (Delpierre et al. 2009 style): color reaches its midpoint when cold degree-days, accumulated from a photoperiod-triggered start date, cross a threshold that may depend on elevation or species. Fit 2–3 parameters by minimizing held-out-year error. It gives a mechanistic story and a check on the regression. If both agree, that is a strong robustness argument.

### 3c. Inclusion rule for any feature or interaction
A term stays in the model only if **all** of the following hold:
1. it lowers the **leave-one-year-out MAE** (hold all 11 locations out together);
2. its sign is consistent across most locations (or it is an interaction that explains why not);
3. there is a physical reason for it.

## 4. Validation protocol (everyone uses the same one)

| Test | How | Number to beat |
|---|---|---|
| Leave-one-year-out (LOYO) | Hold out year y for all locations; μ from the other years; fit δ; predict | Climatology LOYO MAE (from `BUILD_LOG.md`; it was 3.2 d / RMSE 3.95 d with 4 regions) |
| Leave-one-location-out | For the spatial μ model only | Grand-mean MAE |
| Train early, test late | Fit 2001–2015, test 2016–2025 | Climatology under the same split. **This is the direct "does it age well?" test.** |
| Per-location breakdown | MAE by location | Flag locations where the model loses to climatology |

Report MAE, RMSE, bias and an 80% prediction interval built from the LOYO residuals.

## 5. 2026 forecast and the news comparison

1. **Forecast the proxy:** `peak_doy_2026(l) = μ(l) + δ̂(l, 2026)`, using 2026 pre-peak features from `panel.csv`.
2. **Bridge proxy → visual peak:** MODIS mid-greendown comes about 1–3 weeks before reported peak color. Collect **reported peak dates** for about 6 locations × the last several years. Possible sources: state foliage reports (Vermont, the NY "I LOVE NY" weekly foliage report, PA DCNR, WV and VA tourism reports), NPS Shenandoah updates, the Appalachian State "fall color" forecasts for Boone, and local news archives. Estimate `offset(l) = reported − peak_doy` (pooled, or north vs. south if they clearly differ), with its spread.
3. **Final prediction** = proxy forecast + offset, with an interval that combines the model and offset uncertainty.
4. **Compare with 2026 reality.** It is mid-November, so peak has already happened at every location. Pull the 2026 news reports per location and score the forecast. This is exactly the comparison the problem statement asks for. Be honest about misses.

## 6. Robustness and the 10-year question

- **Perturbation tests:** shift the 2026 inputs by ±1 °C (temperature) and ±20% (precipitation) and report how far the forecast moves (days). Also drop one feature at a time, drop one location at a time, and exclude Helene 2024.
- **Sensor dependence:** refit with VIIRS as the target (2013–2025). Then state plainly that MODIS and VIIRS disagree on the trend (the gap grows about 7.7 d/decade). MODIS is near end of mission, so a 10-year operational model would need a MODIS→VIIRS crosswalk.
- **10-year accuracy budget.** Combine:
  - (i) model error from LOYO;
  - (ii) drift in μ: the trend CI was −3.2 to +0.7 d/decade, so μ could be off by up to about 3 d after 10 years;
  - (iii) climate inputs moving outside the training range: projected warming × the `tmax_Aug` slope;
  - (iv) slow changes in forest composition (species-teammate input).
  
  Present this as an error that grows with the forecast horizon, and recommend **re-estimating μ every year** (one line of code once new data are in).

## 7. Order of work

| Step | What | Who | Output |
|---|---|---|---|
| 0 | Build the dataset (`BUILD_DATASET.md`), check `BUILD_LOG.md`, send the zip | Arin's Claude → Arin | `panel_v1.zip` |
| 1 | Parallel analyses | all three | the §7 results table from each person |
| 2 | Merge: agree on ≤ 5 main effects + ≤ 2 interactions + the spatial μ predictors | all | frozen feature list |
| 3 | Choose the μ window (§2 backtest), fit δ (§3a) and the CDD model (§3b), run §4 | Arin | scores table and figures |
| 4 | Collect reported peaks → offset → 2026 forecast → compare with 2026 news (§5) | split by location | forecast table and comparison figure |
| 5 | Robustness and 10-year analysis (§6) | Arin + temporal | sensitivity table and error-vs-horizon figure |
| 6 | Write the report (model, assumptions, validation, 2026 results, limits) + the one-page editor letter (plain-language headline per location, how reliable it is, story ideas) | all | final submission |

## 8. Figures the report will need

1. Map of the 11 regions, colored by `peak_doy_mean`.
2. Peak date by year per location, with VIIRS overlaid.
3. μ-window backtest (MAE vs. k).
4. Feature effects (days per unit) with CIs, shared vs. per location.
5. LOYO predicted vs. observed, against climatology.
6. 2026 forecast vs. reported peak, per location.
7. Elevation rule (flat, then about −2.9 d/100 m).
8. Forecast error vs. horizon (10-year question).
