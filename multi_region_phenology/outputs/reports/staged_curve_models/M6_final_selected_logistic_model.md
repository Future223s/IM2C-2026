# M6 — Final selected logistic decomposition

## Motivation

Construct the simplest full-curve model supported by the staged validation results.

## Methodology

The timing architecture with the lowest held-out MAE was `M2`. Use it to shift the location-specific double-logistic baseline. M5 supported no additional curve dimension, so amplitude, width, rise duration, decline duration, and asymmetry retain their location baselines; their annual deviations are represented as uncertainty.

## Mathematical formula

\[
\widehat C_{s,y}(t)=\overline B_s+\overline A_s\left[L\!\left(\overline k_s^\uparrow[t-\overline\tau_s^\uparrow-\widehat\delta^\tau_{s,y}]\right)-L\!\left(\overline k_s^\downarrow[t-\overline\tau_s^\downarrow-\widehat\delta^\tau_{s,y}]\right)\right].
\]

For the selected M2 specification,

\[
\widehat\delta^\tau_{s,y}=2.75\,\Delta T_{s,y}\quad\text{days}.
\]

The selected timing function is the `M2` specification. The two transitions receive the same shift, so the typical width and asymmetry are preserved unless M5 supports relaxing that constraint.

## Results

MAE=4.22 days, RMSE=5.27 days, R²=0.63, skill versus M1=14.8%. This is the selected timing result from the current hindcast experiment. Forecast-safe performance must be re-estimated using environmental inputs available at each issue date.
