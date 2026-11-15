# M5 — Which logistic-curve dimensions are environmentally predictable?

## Motivation

Test whether the selected environmental architecture should modify only timing or should also modify amplitude, width, rise duration, decline duration, and asymmetry.

## Methodology

Apply the best environmental timing architecture from M2–M4 (`M2`) separately to each derived curve property under the same leave-one-year-out protocol. Compare every result with its training-only location-mean baseline.

## Mathematical formula

For each property \(q\),

\[
\theta_{q,s,y}=\overline\theta_{q,s,-y}+f_q(\Delta\mathbf X_{s,y})+\epsilon_{q,s,y}.
\]

An environmental effect is retained only when it reduces held-out MAE by at least 1%. This prespecified practical threshold prevents negligible numerical changes from being interpreted as biological support.

## Results

| Curve property | Environment MAE | Location-baseline MAE | Skill | CV R² | Retain? |
|---|---:|---:|---:|---:|---|
| peak_doy | 4.215 | 4.946 | 14.8% | 0.63 | Yes |
| peak_intensity_above_baseline | 0.018 | 0.017 | -0.9% | 0.29 | No |
| fwhm_days | 5.824 | 5.789 | -0.6% | 0.35 | No |
| rise_10_90_days | 7.139 | 7.177 | 0.5% | 0.04 | No |
| fall_90_10_days | 6.832 | 6.771 | -0.9% | 0.27 | No |
| shape_asymmetry_days | 10.260 | 10.282 | 0.2% | 0.19 | No |

Only timing crosses the practical-support threshold. The purpose of M5 is model selection, not maximizing in-sample fit. Unsupported environmental effects remain random year-to-year deviations around the location baseline.
