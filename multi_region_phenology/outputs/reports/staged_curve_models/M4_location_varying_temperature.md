# M4 — Location-varying temperature sensitivity

## Motivation

Test whether a one-degree temperature anomaly moves the foliage curve by different amounts at different locations.

## Methodology

Extend M3 with temperature-anomaly-by-location interactions. Ridge regularization shrinks unstable site-specific departures toward the shared response, approximating partial pooling. Validation still holds out complete years, so each location remains represented in the training record.

## Mathematical formula

\[
\delta^\tau_{s,y}=(\beta_0+b_s)\Delta T_{s,y}+\boldsymbol{\gamma}^\top\Delta\mathbf{ X}_{s,y,-T},
\]

\[
b_s\approx 0\quad\text{under stronger regularization when the site record is uninformative}.
\]

## Results

MAE=4.48 days, RMSE=5.62 days, R²=0.58, skill versus M1=9.5%. This model is useful only if the location interactions improve held-out-year performance. It is not a cold-start model for a location with no historical curve record.
