# M1 — Location-specific baseline curve

## Motivation

Establish the minimum model: each location has a typical logistic trajectory and typical peak timing. No annual environmental correction is used.

## Methodology

For every held-out year, estimate each location's mean peak date using all other years. Predict the held-out peak using that training-only location mean. This is a leave-one-year-out evaluation, so the predicted year never contributes to its own baseline.

## Mathematical formula

\[
\widehat C_{s,y}(t)=\overline C_{s,-y}(t),
\qquad
\widehat\tau_{s,y}=\overline\tau_{s,-y}.
\]

The notation \(-y\) means that year \(y\) is excluded when constructing the location baseline.

## Results

MAE=4.95 days, RMSE=6.15 days, R²=0.49, skill versus M1=0.0%. This is the benchmark that every environmental extension must beat.
