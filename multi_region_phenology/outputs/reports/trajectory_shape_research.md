# Research audit: discover autumn RCC trajectory shapes before explaining them

## Conclusion

The proposed two-stage strategy is scientifically stronger than direct daily RCC-on-weather regression, provided that uncertainty, camera effects, and site structure remain in the model. The current evidence supports **continuous functional modes plus interpretable curve summaries** as the next step. It does not yet justify a mixture-of-experts model with named biological curve classes.

## Evidence from the assembled data

The audit considered 202 platform-years with at least one Sep-Nov RCC90 observation. 177 passed prespecified shape-screening rules for this audit: at least 55 observed dates, no gap longer than 14 days, and a 5th-to-95th percentile RCC90 amplitude of at least 0.025.

A robust six-parameter difference-of-logistics model fitted to daily interpolated, nine-day Savitzky–Golay-smoothed curves achieved R² ≥ 0.70 for 96.6% of eligible curves. This in-sample statistic describes the smoothed representation and is not forecast accuracy. Moreover, 41.2% of fits placed at least one parameter near its optimization bound. The logistic family is therefore a useful summary model, not a universally identified biological mechanism. Derived peak, half-maximum width, and 10–90% durations are safer downstream quantities than raw scale/separation coefficients.

For normalized calendar-aligned curves, the first three functional components explained 71.7% of variation. After aligning curves at their observed RCC peak, the first three components explained 79.0%. The latter analysis isolates rise/decline shape from calendar timing.

Component interpretation against curve summaries:

- `calendar_pc1` was most strongly associated with `peak_doy` (Spearman ρ=0.87); component signs are arbitrary.
- `calendar_pc2` was most strongly associated with `fwhm_days` (Spearman ρ=-0.37); component signs are arbitrary.
- `calendar_pc3` was most strongly associated with `fwhm_days` (Spearman ρ=0.52); component signs are arbitrary.
- `aligned_pc1` was most strongly associated with `fall_90_10_days` (Spearman ρ=0.51); component signs are arbitrary.
- `aligned_pc2` was most strongly associated with `fwhm_days` (Spearman ρ=0.79); component signs are arbitrary.
- `aligned_pc3` was most strongly associated with `peak_intensity_above_baseline` (Spearman ρ=0.30); component signs are arbitrary.

## Are there discrete trajectory types?

K-means solutions with 2–6 clusters were compared on the functional scores. The best silhouette occurred at K=2 with silhouette=0.32, perturbation stability ARI=0.96, and adjusted mutual information with site=0.19.

A silhouette below roughly 0.4 is weak evidence for sharply separated types. Strong site association would additionally suggest that clusters encode camera/species/location identity rather than recurrent weather-driven regimes. Unless both separation and out-of-sample gating performance improve, mixture weights should not be interpreted as biological state probabilities.

## Environmental explanation stage

Environmental summaries used fixed, interpretable windows rather than windows defined by the observed peak: June-August and September-October temperature, precipitation, VPD, and radiation, plus autumn frost days and cooling-degree accumulation. Both predictors and curve targets were demeaned within site before association testing, so the correlations address interannual anomalies rather than permanent differences among locations. Significance was evaluated by permuting environmental anomalies within each site and controlling the false-discovery rate across all tested pairs.

9 associations survived the exploratory 10% false-discovery-rate threshold:

| Curve property | Environmental anomaly | Spearman ρ | q-value | N |
|---|---|---:|---:|---:|
| peak_doy | sep_oct_temperature_mean_c | 0.52 | 0.022 | 158 |
| peak_doy | sep_oct_vpd_kpa | 0.28 | 0.022 | 158 |
| peak_doy | sep_oct_frost_days | -0.23 | 0.022 | 158 |
| peak_doy | sep_oct_cooling_degree_days | -0.39 | 0.022 | 158 |
| calendar_pc1 | sep_oct_temperature_mean_c | 0.58 | 0.022 | 158 |
| calendar_pc1 | sep_oct_vpd_kpa | 0.24 | 0.022 | 158 |
| calendar_pc1 | sep_oct_cooling_degree_days | -0.44 | 0.022 | 158 |
| calendar_pc2 | sep_oct_precipitation_mm | -0.25 | 0.022 | 158 |
| calendar_pc1 | sep_oct_frost_days | -0.28 | 0.040 | 158 |

These rows are not independent mechanisms. September-October mean temperature, cooling accumulation, frost frequency, and VPD are mutually correlated descriptions of the same seasonal weather. Their associations chiefly support a timing signal: warmer autumns tended to have later RCC peaks, whereas greater cooling and more frost tended to accompany earlier peaks. The full September-October windows also overlap the outcome season, so this is explanatory hindcasting, not a real-time forecast experiment. Forecast models must truncate every environmental summary at the issue date.

Leave-one-site-out ridge tests provide a stricter assessment of whether environmental-anomaly relationships transfer to a new location:

| Target | CV R² | CV MAE | Zero-anomaly baseline MAE |
|---|---:|---:|---:|
| peak_doy | 0.28 | 3.87 | 4.58 |
| peak_intensity_above_baseline | -0.04 | 0.02 | 0.02 |
| fwhm_days | -0.08 | 5.63 | 5.34 |
| rise_10_90_days | -0.04 | 6.83 | 6.59 |
| fall_90_10_days | -0.10 | 6.52 | 6.17 |
| shape_asymmetry_days | -0.09 | 9.70 | 9.28 |
| calendar_pc1 | 0.38 | 0.54 | 0.68 |
| calendar_pc2 | -0.01 | 0.52 | 0.53 |
| calendar_pc3 | 0.03 | 0.39 | 0.39 |

These predictive tests are exploratory because Daymet is gridded, sites have unequal record lengths, and quantitative species composition is not yet available.

## Recommended model sequence

1. Fit a hierarchical latent curve to daily RCC90 and GCC90 jointly, with site, platform, ROI, and year effects and observation-level QC. Use a double logistic or monotone-spline rise/fall representation and retain posterior uncertainty.
2. Report derived quantities—peak date, realized amplitude, half-maximum width, rise duration, and fall duration—rather than interpreting all raw nonlinear coefficients.
3. Use calendar-aligned fPCA to discover timing-plus-shape modes, and peak-aligned fPCA to test for residual shape variation.
4. Regress site-demeaned scores/parameters on prespecified environmental anomalies with hierarchical partial pooling. Retain site climatology, elevation, and eventual composition measurements as between-site predictors.
5. Compare against a direct functional weather model under blocked validation. The two-stage method is preferred only if it improves held-out curve reconstruction or scientific interpretability without discarding parameter uncertainty.
6. Attempt a mixture of experts only if stable clusters recur across platforms/sites and a gating model predicts them out of sample better than continuous-score regression.

## Important methodological qualifications

- The two-stage plug-in approach can understate uncertainty if fitted curve parameters are treated as exact. A hierarchical joint model or multiple posterior draws should propagate first-stage uncertainty into the environmental stage.
- RCC90 is a camera color coordinate, not a universal latent foliage-color state. GCC decline, field coloration, leaf fall, and satellite greendown require separate observation functions.
- Amplitude is particularly sensitive to camera exposure, view composition, and species pigments. Cross-site amplitude comparisons require platform effects and preferably raw-image color calibration.
- Functional PCA components are continuous axes with arbitrary sign. Naming them ‘early’, ‘broad’, or ‘sharp’ requires correlation with independent, interpretable metrics.
- Photoperiod is deterministic from latitude and date. Its causal role cannot be separated from seasonal time using a simple within-year correlation; cross-site/species structure or process-based constraints are required.

## Literature basis

- Double-logistic functions are established phenology summaries with interpretable seasonal parameters, but model choice changes transition estimates: Beck et al. (2006), https://doi.org/10.1016/j.rse.2005.10.021; Atkinson et al. (2012), https://doi.org/10.1016/j.rse.2012.04.001.
- PhenoCam transition products use smoothed GCC stages and report amplitude thresholds with confidence intervals rather than treating raw daily extrema as exact: https://daac.ornl.gov/VEGETATION/guides/PhenoCam_V2.html.
- Bayesian hierarchical phenology models demonstrate why sparse observations and parameter uncertainty should be modeled explicitly: Gao, Gray & Reich (2021), https://doi.org/10.1016/j.rse.2021.112484; Babcock, Finley & Looker (2021), https://doi.org/10.1016/j.rse.2021.112471.
- Autumn controls are nonlinear and species/location dependent. Temperature, photoperiod, drought, and prior growing-season conditions can act differently across regimes: Gill et al. (2015), https://doi.org/10.1093/aob/mcv055; Chen et al. (2020), https://doi.org/10.1038/s41558-020-0820-2; Moon et al. (2022), https://doi.org/10.1016/j.agrformet.2022.109026.

## Decision

Proceed with the two-stage functional strategy, but use continuous fPCA scores and uncertainty-aware curve summaries first. Treat mixture-of-experts as a later hypothesis test, not the default architecture.
