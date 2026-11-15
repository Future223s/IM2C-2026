# M2 — Shared temperature-driven timing shift

## Motivation

Test the primary empirical finding: warmer- or colder-than-normal autumns primarily shift the curve in time rather than changing its shape.

## Methodology

Compute September–October mean-temperature anomalies relative to training-only site climatology. Fit one shared coefficient across locations, then shift both logistic transitions by the predicted number of days.

## Mathematical formula

\[
\delta^\tau_{s,y}=\beta_T\Delta T_{s,y},
\]

\[
\widehat C_{s,y}(t)=\overline B_s+\overline A_s\left[L\!\left(\overline k_s^\uparrow[t-\overline\tau_s^\uparrow-\delta^\tau_{s,y}]\right)-L\!\left(\overline k_s^\downarrow[t-\overline\tau_s^\downarrow-\delta^\tau_{s,y}]\right)\right].
\]

## Results

MAE=4.22 days, RMSE=5.27 days, R²=0.63, skill versus M1=14.8%. A positive skill value means temperature improves on the location baseline.

The full-data descriptive coefficient is **+2.75 days per °C**: a 1 °C warmer-than-usual September–October period is associated with a 2.75-day later peak at the same location. The prediction score above, unlike this descriptive refit, is entirely out of sample.

This is an explanatory hindcast because the full September–October temperature window overlaps the foliage season. A real forecast must truncate weather at the issue date or substitute weather forecasts.
