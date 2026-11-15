# M3 — Multivariable environmental timing shift

## Motivation

Determine whether precipitation, VPD, radiation, frost, cooling accumulation, and summer history improve timing predictions beyond mean autumn temperature.

## Methodology

Construct training-only within-site anomalies for ten summer and autumn variables. Fit a ridge regression, selecting the penalty from \(0.1,1,10,100\) inside each held-out-year fold using grouped validation by year. Environmental variables still affect only the common timing shift.

## Mathematical formula

\[
\delta^\tau_{s,y}=\boldsymbol{\beta}^\top\Delta\mathbf{ X}_{s,y},
\]

\[
\widehat{\boldsymbol{\beta}}_\lambda=\arg\min_{\boldsymbol{\beta}}\left\{\|\mathbf{ d}-X\boldsymbol{\beta}\|_2^2+\lambda\|\boldsymbol{\beta}\|_2^2\right\}.
\]

## Results

MAE=4.30 days, RMSE=5.37 days, R²=0.61, skill versus M1=13.0%. Compare this directly with M2: added variables are retained only if held-out error declines.
