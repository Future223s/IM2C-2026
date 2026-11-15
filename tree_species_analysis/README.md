# IM2C — Tree species composition analysis

This standalone package contains the eleven-location FIA tree-count analysis and the original lambda=3 ridge model. It does not include the separate satellite peak-date model.

## Reproduce

Python 3.10+; install dependencies with `python -m pip install -r requirements.txt`, then run `python src/run_analysis.py` (Windows: `py` instead of `python`). The script validates the raw FIA totals against the environmental workbook and writes category counts, coefficients, predictions and metrics. Original result CSVs are preserved separately from recomputed outputs.

## Contents

- `data/raw_fia/`: eleven original FIA CSVs (30-mile circles), one per location. Asheville's two uploaded copies are identical; copy (1) is retained.
- `data/processed/data.xlsx`: original category percentages workbook.
- `data/processed/data_environment.xlsx`: environmental inputs and Sources worksheet.
- `data/processed/appalachian_model_distributions.csv`: eleven-category percentages, categories by city.
- `data/processed/asheville_model_distribution.csv`: original Asheville example.
- `data/processed/category_counts.csv`: recomputed estimated tree counts and percentages for every city/category.
- `src/model_distribution.py`: original species-group mapping and single-city aggregation.
- `src/run_analysis.py`: reproducible multi-output ridge fit and leave-one-location-out validation.
- `src/predict.py`: predict all eleven percentages from a JSON object of environmental factors.
- `results/tree_ridge_coefficients.csv`, `tree_ridge_predictions.csv`: original fitted coefficients and predictions.
- `results/recomputed_coefficients.csv`, `recomputed_predictions.csv`, `ridge_metrics.csv`: independently reproduced results and evaluation metrics.

## Model

For 19 predictors, impute missing values with training means, standardize using population standard deviations, and fit an unpenalized intercept plus ridge coefficients: B = (Z^T Z + 3 I)^(-1) Z^T (Y - mean(Y)). Outputs are percentages, not fractions. The redundant elevation in meters and constant equinox daylight are excluded. Each validation fold recomputes preprocessing on its ten training locations. CV R² uses the overall observed mean as the reference denominator. The shared penalty preserves the sum of percentages (100), but negative predictions are possible. The optional normalized prediction clips negatives and rescales; reported metrics use the original affine predictions.

Use `python src/predict.py inputs.json` with a JSON object keyed by the exact Predictor names in the coefficients CSV (excluding Intercept). Omitted/null predictors are mean-imputed. Longitude is signed degrees east (negative for these US locations). Temperature is Fahrenheit, elevation feet, precipitation inches, humidity percent, daylight hours.

## Limitations

Only eleven locations; correlated predictors, mixed climate-normal sources/periods, missing humidity/precipitation, town-center rather than forest-wide elevation. See the Sources worksheet. FIA counts are estimates, not canopy area. Broad FIA groups can contain exceptions to the evergreen/deciduous labels. These are exploratory composition regressions, not validated foliage-date predictions. Several categories have negative cross-validation R²; use observed FIA distributions when available. No new data were downloaded for this package.
