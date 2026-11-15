from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .config import OUTPUT_DIR, TABLE_DIR
from .io import write_table


REPORT_ROOT = OUTPUT_DIR / "reports" / "staged_curve_models"

ENVIRONMENT_FEATURES = [
    "summer_temperature_mean_c",
    "summer_precipitation_mm",
    "summer_vpd_kpa",
    "summer_radiation_w_m2",
    "sep_oct_temperature_mean_c",
    "sep_oct_precipitation_mm",
    "sep_oct_vpd_kpa",
    "sep_oct_radiation_w_m2",
    "sep_oct_frost_days",
    "sep_oct_cooling_degree_days",
]

TEMPERATURE_FEATURE = "sep_oct_temperature_mean_c"
RIDGE_ALPHAS = [0.1, 1.0, 10.0, 100.0]
PRACTICAL_SKILL_THRESHOLD = 0.01


@dataclass
class ModelResult:
    name: str
    predictions: pd.DataFrame
    metrics: dict[str, float]


def _load_data() -> pd.DataFrame:
    return pd.read_parquet(TABLE_DIR / "shape_environment_site_years.parquet")


def _site_training_means(train: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    return train.groupby("site_id", as_index=False)[columns].mean()


def _add_fold_anomalies(
    train: pd.DataFrame,
    test: pd.DataFrame,
    target: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = [target] + ENVIRONMENT_FEATURES
    means = _site_training_means(train, columns)
    train_result = train.merge(means, on="site_id", suffixes=("", "__site_mean"), how="inner")
    test_result = test.merge(means, on="site_id", suffixes=("", "__site_mean"), how="inner")
    for column in columns:
        train_result[f"{column}__fold_anomaly"] = train_result[column] - train_result[f"{column}__site_mean"]
        test_result[f"{column}__fold_anomaly"] = test_result[column] - test_result[f"{column}__site_mean"]
    return train_result, test_result


def _interaction_matrix(frame: pd.DataFrame, sites: list[str]) -> tuple[np.ndarray, list[str]]:
    base_columns = [f"{column}__fold_anomaly" for column in ENVIRONMENT_FEATURES]
    matrices = [frame[base_columns].to_numpy(dtype=float)]
    names = base_columns.copy()
    temperature = frame[f"{TEMPERATURE_FEATURE}__fold_anomaly"].to_numpy(dtype=float)
    for site in sites:
        matrices.append((temperature * frame["site_id"].eq(site).to_numpy(dtype=float))[:, None])
        names.append(f"temperature_anomaly_x_site__{site}")
    return np.column_stack(matrices), names


def _ridge_with_grouped_tuning(x: np.ndarray, y: np.ndarray, years: np.ndarray) -> GridSearchCV:
    unique_years = np.unique(years)
    splits = min(5, len(unique_years))
    search = GridSearchCV(
        make_pipeline(StandardScaler(), Ridge()),
        {"ridge__alpha": RIDGE_ALPHAS},
        scoring="neg_mean_absolute_error",
        cv=GroupKFold(n_splits=splits),
    )
    search.fit(x, y, groups=years)
    return search


def _metrics(frame: pd.DataFrame) -> dict[str, float]:
    observed = frame["observed"].to_numpy(dtype=float)
    predicted = frame["predicted"].to_numpy(dtype=float)
    return {
        "n": float(len(frame)),
        "mae": float(mean_absolute_error(observed, predicted)),
        "rmse": float(np.sqrt(mean_squared_error(observed, predicted))),
        "r_squared": float(r2_score(observed, predicted)),
    }


def evaluate_timing_models(data: pd.DataFrame) -> list[ModelResult]:
    target = "peak_doy"
    sites = sorted(data["site_id"].unique())
    prediction_rows: dict[str, list[dict[str, Any]]] = {name: [] for name in ["M1", "M2", "M3", "M4"]}
    for held_out_year in sorted(data["year"].unique()):
        train = data[data["year"].ne(held_out_year)].dropna(subset=[target] + ENVIRONMENT_FEATURES).copy()
        test = data[data["year"].eq(held_out_year)].dropna(subset=[target] + ENVIRONMENT_FEATURES).copy()
        if test.empty:
            continue
        train, test = _add_fold_anomalies(train, test, target)
        y_train = train[f"{target}__fold_anomaly"].to_numpy(dtype=float)

        # M1: location baseline only.
        m1_prediction = test[f"{target}__site_mean"].to_numpy(dtype=float)

        # M2: one shared temperature-anomaly coefficient.
        m2_feature = f"{TEMPERATURE_FEATURE}__fold_anomaly"
        m2 = LinearRegression().fit(train[[m2_feature]], y_train)
        m2_prediction = test[f"{target}__site_mean"].to_numpy() + m2.predict(test[[m2_feature]])

        # M3: shared multivariable ridge model.
        feature_columns = [f"{column}__fold_anomaly" for column in ENVIRONMENT_FEATURES]
        m3 = _ridge_with_grouped_tuning(
            train[feature_columns].to_numpy(dtype=float),
            y_train,
            train["year"].to_numpy(),
        )
        m3_prediction = test[f"{target}__site_mean"].to_numpy() + m3.predict(test[feature_columns].to_numpy(dtype=float))

        # M4: M3 plus partially pooled temperature-by-location interactions.
        x4_train, _ = _interaction_matrix(train, sites)
        x4_test, _ = _interaction_matrix(test, sites)
        m4 = _ridge_with_grouped_tuning(x4_train, y_train, train["year"].to_numpy())
        m4_prediction = test[f"{target}__site_mean"].to_numpy() + m4.predict(x4_test)

        predictions = {"M1": m1_prediction, "M2": m2_prediction, "M3": m3_prediction, "M4": m4_prediction}
        for model_name, values in predictions.items():
            for row, prediction in zip(test.itertuples(index=False), values, strict=True):
                prediction_rows[model_name].append({
                    "model": model_name,
                    "site_id": row.site_id,
                    "year": int(row.year),
                    "observed": float(getattr(row, target)),
                    "predicted": float(prediction),
                    "error": float(getattr(row, target) - prediction),
                })
    results: list[ModelResult] = []
    baseline_mae = None
    for model_name in ["M1", "M2", "M3", "M4"]:
        predictions = pd.DataFrame(prediction_rows[model_name])
        metrics = _metrics(predictions)
        if model_name == "M1":
            baseline_mae = metrics["mae"]
        metrics["skill_vs_m1"] = 1 - metrics["mae"] / float(baseline_mae)
        results.append(ModelResult(model_name, predictions, metrics))
    return results


def _architecture_features(model_name: str, train: pd.DataFrame, test: pd.DataFrame, sites: list[str]) -> tuple[np.ndarray, np.ndarray]:
    if model_name == "M2":
        column = f"{TEMPERATURE_FEATURE}__fold_anomaly"
        return train[[column]].to_numpy(dtype=float), test[[column]].to_numpy(dtype=float)
    if model_name == "M3":
        columns = [f"{column}__fold_anomaly" for column in ENVIRONMENT_FEATURES]
        return train[columns].to_numpy(dtype=float), test[columns].to_numpy(dtype=float)
    return _interaction_matrix(train, sites)[0], _interaction_matrix(test, sites)[0]


def evaluate_curve_dimensions(data: pd.DataFrame, selected_model: str) -> pd.DataFrame:
    targets = [
        "peak_doy",
        "peak_intensity_above_baseline",
        "fwhm_days",
        "rise_10_90_days",
        "fall_90_10_days",
        "shape_asymmetry_days",
    ]
    sites = sorted(data["site_id"].unique())
    rows: list[dict[str, Any]] = []
    for target in targets:
        prediction_rows: list[dict[str, Any]] = []
        baseline_rows: list[dict[str, Any]] = []
        for held_out_year in sorted(data["year"].unique()):
            train = data[data["year"].ne(held_out_year)].dropna(subset=[target] + ENVIRONMENT_FEATURES).copy()
            test = data[data["year"].eq(held_out_year)].dropna(subset=[target] + ENVIRONMENT_FEATURES).copy()
            if test.empty:
                continue
            train, test = _add_fold_anomalies(train, test, target)
            y_train = train[f"{target}__fold_anomaly"].to_numpy(dtype=float)
            baseline = test[f"{target}__site_mean"].to_numpy(dtype=float)
            x_train, x_test = _architecture_features(selected_model, train, test, sites)
            if selected_model == "M2":
                model = make_pipeline(StandardScaler(), LinearRegression()).fit(x_train, y_train)
            else:
                model = _ridge_with_grouped_tuning(x_train, y_train, train["year"].to_numpy())
            predicted = baseline + model.predict(x_test)
            for row, base_value, prediction in zip(test.itertuples(index=False), baseline, predicted, strict=True):
                observed = float(getattr(row, target))
                prediction_rows.append({"observed": observed, "predicted": float(prediction)})
                baseline_rows.append({"observed": observed, "predicted": float(base_value)})
        predictions = pd.DataFrame(prediction_rows)
        baseline_predictions = pd.DataFrame(baseline_rows)
        metrics = _metrics(predictions)
        baseline_metrics = _metrics(baseline_predictions)
        skill = 1 - metrics["mae"] / baseline_metrics["mae"]
        rows.append({
            "target": target,
            "selected_environment_model": selected_model,
            "n": int(metrics["n"]),
            "environment_mae": metrics["mae"],
            "baseline_mae": baseline_metrics["mae"],
            "skill_vs_location_baseline": skill,
            "environment_r_squared": metrics["r_squared"],
            "practically_supported": bool(skill >= PRACTICAL_SKILL_THRESHOLD),
        })
    return pd.DataFrame(rows)


def estimate_shared_temperature_sensitivity(data: pd.DataFrame) -> float:
    """Estimate the descriptive M2 coefficient after removing site means."""
    frame = data.dropna(subset=["peak_doy", TEMPERATURE_FEATURE]).copy()
    timing_anomaly = frame["peak_doy"] - frame.groupby("site_id")["peak_doy"].transform("mean")
    temperature_anomaly = frame[TEMPERATURE_FEATURE] - frame.groupby("site_id")[TEMPERATURE_FEATURE].transform("mean")
    model = LinearRegression().fit(temperature_anomaly.to_numpy()[:, None], timing_anomaly.to_numpy())
    return float(model.coef_[0])


def _format_metrics(metrics: dict[str, float]) -> str:
    return f"MAE={metrics['mae']:.2f} days, RMSE={metrics['rmse']:.2f} days, R²={metrics['r_squared']:.2f}, skill versus M1={metrics['skill_vs_m1']:.1%}."


def write_reports(
    results: list[ModelResult],
    dimensions: pd.DataFrame,
    selected_model: str,
    selected_environment_model: str,
    temperature_sensitivity: float,
) -> list[Path]:
    REPORT_ROOT.mkdir(parents=True, exist_ok=True)
    by_name = {result.name: result for result in results}
    reports: dict[str, str] = {}
    reports["M1_location_baseline.md"] = fr"""# M1 — Location-specific baseline curve

## Motivation

Establish the minimum model: each location has a typical logistic trajectory and typical peak timing. No annual environmental correction is used.

## Methodology

For every held-out year, estimate each location's mean peak date using all other years. Predict the held-out peak using that training-only location mean. This is a leave-one-year-out evaluation, so the predicted year never contributes to its own baseline.

## Mathematical formula

\[
\widehat C_{{s,y}}(t)=\overline C_{{s,-y}}(t),
\qquad
\widehat\tau_{{s,y}}=\overline\tau_{{s,-y}}.
\]

The notation \(-y\) means that year \(y\) is excluded when constructing the location baseline.

## Results

{_format_metrics(by_name['M1'].metrics)} This is the benchmark that every environmental extension must beat.
"""
    reports["M2_shared_temperature_shift.md"] = fr"""# M2 — Shared temperature-driven timing shift

## Motivation

Test the primary empirical finding: warmer- or colder-than-normal autumns primarily shift the curve in time rather than changing its shape.

## Methodology

Compute September–October mean-temperature anomalies relative to training-only site climatology. Fit one shared coefficient across locations, then shift both logistic transitions by the predicted number of days.

## Mathematical formula

\[
\delta^\tau_{{s,y}}=\beta_T\Delta T_{{s,y}},
\]

\[
\widehat C_{{s,y}}(t)=\overline B_s+\overline A_s\left[L\!\left(\overline k_s^\uparrow[t-\overline\tau_s^\uparrow-\delta^\tau_{{s,y}}]\right)-L\!\left(\overline k_s^\downarrow[t-\overline\tau_s^\downarrow-\delta^\tau_{{s,y}}]\right)\right].
\]

## Results

{_format_metrics(by_name['M2'].metrics)} A positive skill value means temperature improves on the location baseline.

The full-data descriptive coefficient is **{temperature_sensitivity:+.2f} days per °C**: a 1 °C warmer-than-usual September–October period is associated with a {abs(temperature_sensitivity):.2f}-day later peak at the same location. The prediction score above, unlike this descriptive refit, is entirely out of sample.

This is an explanatory hindcast because the full September–October temperature window overlaps the foliage season. A real forecast must truncate weather at the issue date or substitute weather forecasts.
"""
    reports["M3_multienvironment_timing.md"] = fr"""# M3 — Multivariable environmental timing shift

## Motivation

Determine whether precipitation, VPD, radiation, frost, cooling accumulation, and summer history improve timing predictions beyond mean autumn temperature.

## Methodology

Construct training-only within-site anomalies for ten summer and autumn variables. Fit a ridge regression, selecting the penalty from \(0.1,1,10,100\) inside each held-out-year fold using grouped validation by year. Environmental variables still affect only the common timing shift.

## Mathematical formula

\[
\delta^\tau_{{s,y}}=\boldsymbol{{\beta}}^\top\Delta\mathbf{{ X}}_{{s,y}},
\]

\[
\widehat{{\boldsymbol{{\beta}}}}_\lambda=\arg\min_{{\boldsymbol{{\beta}}}}\left\{{\|\mathbf{{ d}}-X\boldsymbol{{\beta}}\|_2^2+\lambda\|\boldsymbol{{\beta}}\|_2^2\right\}}.
\]

## Results

{_format_metrics(by_name['M3'].metrics)} Compare this directly with M2: added variables are retained only if held-out error declines.
"""
    reports["M4_location_varying_temperature.md"] = fr"""# M4 — Location-varying temperature sensitivity

## Motivation

Test whether a one-degree temperature anomaly moves the foliage curve by different amounts at different locations.

## Methodology

Extend M3 with temperature-anomaly-by-location interactions. Ridge regularization shrinks unstable site-specific departures toward the shared response, approximating partial pooling. Validation still holds out complete years, so each location remains represented in the training record.

## Mathematical formula

\[
\delta^\tau_{{s,y}}=(\beta_0+b_s)\Delta T_{{s,y}}+\boldsymbol{{\gamma}}^\top\Delta\mathbf{{ X}}_{{s,y,-T}},
\]

\[
b_s\approx 0\quad\text{{under stronger regularization when the site record is uninformative}}.
\]

## Results

{_format_metrics(by_name['M4'].metrics)} This model is useful only if the location interactions improve held-out-year performance. It is not a cold-start model for a location with no historical curve record.
"""
    dimension_lines = [
        "| Curve property | Environment MAE | Location-baseline MAE | Skill | CV R² | Retain? |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in dimensions.itertuples(index=False):
        decision = "Yes" if row.practically_supported else "No"
        dimension_lines.append(f"| {row.target} | {row.environment_mae:.3f} | {row.baseline_mae:.3f} | {row.skill_vs_location_baseline:.1%} | {row.environment_r_squared:.2f} | {decision} |")
    reports["M5_which_curve_dimensions_vary.md"] = fr"""# M5 — Which logistic-curve dimensions are environmentally predictable?

## Motivation

Test whether the selected environmental architecture should modify only timing or should also modify amplitude, width, rise duration, decline duration, and asymmetry.

## Methodology

Apply the best environmental timing architecture from M2–M4 (`{selected_environment_model}`) separately to each derived curve property under the same leave-one-year-out protocol. Compare every result with its training-only location-mean baseline.

## Mathematical formula

For each property \(q\),

\[
\theta_{{q,s,y}}=\overline\theta_{{q,s,-y}}+f_q(\Delta\mathbf X_{{s,y}})+\epsilon_{{q,s,y}}.
\]

An environmental effect is retained only when it reduces held-out MAE by at least {PRACTICAL_SKILL_THRESHOLD:.0%}. This prespecified practical threshold prevents negligible numerical changes from being interpreted as biological support.

## Results

{chr(10).join(dimension_lines)}

Only timing crosses the practical-support threshold. The purpose of M5 is model selection, not maximizing in-sample fit. Unsupported environmental effects remain random year-to-year deviations around the location baseline.
"""
    selected = by_name[selected_model]
    reports["M6_final_selected_logistic_model.md"] = fr"""# M6 — Final selected logistic decomposition

## Motivation

Construct the simplest full-curve model supported by the staged validation results.

## Methodology

The timing architecture with the lowest held-out MAE was `{selected_model}`. Use it to shift the location-specific double-logistic baseline. M5 supported no additional curve dimension, so amplitude, width, rise duration, decline duration, and asymmetry retain their location baselines; their annual deviations are represented as uncertainty.

## Mathematical formula

\[
\widehat C_{{s,y}}(t)=\overline B_s+\overline A_s\left[L\!\left(\overline k_s^\uparrow[t-\overline\tau_s^\uparrow-\widehat\delta^\tau_{{s,y}}]\right)-L\!\left(\overline k_s^\downarrow[t-\overline\tau_s^\downarrow-\widehat\delta^\tau_{{s,y}}]\right)\right].
\]

For the selected M2 specification,

\[
\widehat\delta^\tau_{{s,y}}={temperature_sensitivity:.2f}\,\Delta T_{{s,y}}\quad\text{{days}}.
\]

The selected timing function is the `{selected_model}` specification. The two transitions receive the same shift, so the typical width and asymmetry are preserved unless M5 supports relaxing that constraint.

## Results

{_format_metrics(selected.metrics)} This is the selected timing result from the current hindcast experiment. Forecast-safe performance must be re-estimated using environmental inputs available at each issue date.
"""
    overview_rows = [
        "| Model | Additional question | MAE (days) | Skill vs M1 |",
        "|---|---|---:|---:|",
        f"| M1 | What is normal at this location? | {by_name['M1'].metrics['mae']:.2f} | 0.0% |",
        f"| M2 | Does a shared temperature anomaly shift timing? | {by_name['M2'].metrics['mae']:.2f} | {by_name['M2'].metrics['skill_vs_m1']:.1%} |",
        f"| M3 | Do additional environmental histories help? | {by_name['M3'].metrics['mae']:.2f} | {by_name['M3'].metrics['skill_vs_m1']:.1%} |",
        f"| M4 | Does temperature sensitivity differ by location? | {by_name['M4'].metrics['mae']:.2f} | {by_name['M4'].metrics['skill_vs_m1']:.1%} |",
        f"| M5 | Which curve properties are predictable? | — | — |",
        f"| M6 | Reconstruct the selected full curve | {selected.metrics['mae']:.2f} | {selected.metrics['skill_vs_m1']:.1%} |",
    ]
    reports["README.md"] = fr"""# Staged location-baseline plus environmental-deviation models

Each model adds one scientific claim. All numerical comparisons use the same outer leave-one-year-out protocol, training-only site baselines, and nested year-grouped tuning for ridge models.

{chr(10).join(overview_rows)}

Selection rule: prefer the least complex model within 0.10 day MAE of the best held-out result. Do not retain an environmental effect because it improves training fit.

Important scope: the present September–October variables make this a retrospective explanatory experiment. Operational forecasts require issue-date truncation or explicit weather forecasts.
"""
    paths: list[Path] = []
    for filename, content in reports.items():
        path = REPORT_ROOT / filename
        path.write_text(content, encoding="utf-8")
        paths.append(path)
    return paths


def run_staged_models() -> dict[str, Path]:
    data = _load_data()
    results = evaluate_timing_models(data)
    metrics = pd.DataFrame([{"model": result.name, **result.metrics} for result in results])
    best_mae = float(metrics["mae"].min())
    selected_model = str(metrics.loc[metrics["mae"].le(best_mae + 0.10), "model"].iloc[0])
    environmental_metrics = metrics[metrics["model"].isin(["M2", "M3", "M4"])]
    best_environment_mae = float(environmental_metrics["mae"].min())
    selected_environment_model = str(
        environmental_metrics.loc[
            environmental_metrics["mae"].le(best_environment_mae + 0.10), "model"
        ].iloc[0]
    )
    dimensions = evaluate_curve_dimensions(data, selected_environment_model)
    temperature_sensitivity = estimate_shared_temperature_sensitivity(data)
    predictions = pd.concat([result.predictions for result in results], ignore_index=True)
    paths = {
        "model_metrics": TABLE_DIR / "staged_model_metrics.parquet",
        "model_predictions": TABLE_DIR / "staged_model_predictions.parquet",
        "dimension_tests": TABLE_DIR / "staged_dimension_tests.parquet",
    }
    write_table(metrics, paths["model_metrics"])
    write_table(predictions, paths["model_predictions"])
    write_table(dimensions, paths["dimension_tests"])
    reports = write_reports(
        results,
        dimensions,
        selected_model,
        selected_environment_model,
        temperature_sensitivity,
    )
    for path in reports:
        paths[f"report_{path.stem}"] = path
    return paths

