from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares
from scipy.signal import savgol_filter
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score, r2_score, silhouette_score
from sklearn.model_selection import LeaveOneGroupOut, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge

from .config import FIGURE_DIR, REPORT_DIR, TABLE_DIR
from .io import write_table


AUTUMN_GRID = np.arange(244, 335, dtype=float)  # Sep 1-Nov 30 in a non-leap calendar
RELATIVE_GRID = np.arange(-20, 21, dtype=float)


@dataclass
class CurveRecord:
    site_id: str
    platform_id: str
    year: int
    observed_days: int
    longest_gap_days: int
    raw_amplitude: float
    raw_peak_doy: float
    daily_curve: np.ndarray


def _read(name: str) -> pd.DataFrame:
    return pd.read_parquet(TABLE_DIR / f"{name}.parquet")


def _longest_gap(doys: np.ndarray) -> int:
    observed = set(int(value) for value in doys)
    current = longest = 0
    for day in AUTUMN_GRID.astype(int):
        if day in observed:
            current = 0
        else:
            current += 1
            longest = max(longest, current)
    return longest


def prepare_curves() -> tuple[list[CurveRecord], pd.DataFrame]:
    observations = _read("phenology_observations")
    data = observations[
        observations["observation_type"].eq("rcc_90")
        & observations["quality_usable"]
        & observations["value"].notna()
        & observations["timestamp"].dt.month.isin([9, 10, 11])
    ].copy()
    data["year"] = data["timestamp"].dt.year
    # Overlapping ROI sequences on a platform-date remain in the source table but
    # are collapsed to a median here so a configuration transition is not double-weighted.
    data = data.groupby(["site_id", "platform_id", "year", "timestamp"], as_index=False)["value"].median()
    curves: list[CurveRecord] = []
    audit_rows: list[dict[str, Any]] = []
    for (site_id, platform_id, year), group in data.groupby(["site_id", "platform_id", "year"]):
        group = group.sort_values("timestamp")
        doys = group["timestamp"].dt.day_of_year.to_numpy(dtype=float)
        values = group["value"].to_numpy(dtype=float)
        observed_days = len(np.unique(doys))
        longest_gap = _longest_gap(doys)
        amplitude = float(np.quantile(values, 0.95) - np.quantile(values, 0.05))
        eligible = observed_days >= 55 and longest_gap <= 14 and amplitude >= 0.025
        reason = "eligible"
        if observed_days < 55:
            reason = "fewer_than_55_days"
        elif longest_gap > 14:
            reason = "gap_over_14_days"
        elif amplitude < 0.025:
            reason = "rcc_amplitude_below_0.025"
        audit_rows.append({
            "site_id": site_id,
            "platform_id": platform_id,
            "year": int(year),
            "observed_days": observed_days,
            "longest_gap_days": longest_gap,
            "raw_rcc90_5_95_amplitude": amplitude,
            "shape_analysis_eligible": eligible,
            "exclusion_reason": reason,
        })
        if not eligible:
            continue
        interpolator = PchipInterpolator(doys, values, extrapolate=False)
        daily = interpolator(AUTUMN_GRID)
        if np.isnan(daily).any():
            # Only edge gaps are extrapolated, using the nearest observed value.
            daily[AUTUMN_GRID < doys.min()] = values[0]
            daily[AUTUMN_GRID > doys.max()] = values[-1]
        daily = savgol_filter(daily, window_length=9, polyorder=2, mode="interp")
        raw_peak = float(AUTUMN_GRID[int(np.nanargmax(daily))])
        curves.append(CurveRecord(site_id, platform_id, int(year), observed_days, longest_gap, amplitude, raw_peak, daily))
    return curves, pd.DataFrame(audit_rows)


def _logistic(x: np.ndarray) -> np.ndarray:
    return 1 / (1 + np.exp(np.clip(-x, -50, 50)))


def double_logistic(parameters: np.ndarray, day: np.ndarray) -> np.ndarray:
    baseline, scale, onset, separation, rise_rate, fall_rate = parameters
    decline = onset + separation
    return baseline + scale * (_logistic(rise_rate * (day - onset)) - _logistic(fall_rate * (day - decline)))


def _crossing(day: np.ndarray, curve: np.ndarray, threshold: float, rising: bool, peak_index: int) -> float:
    if rising:
        x, y = day[: peak_index + 1], curve[: peak_index + 1]
        candidates = np.where((y[:-1] <= threshold) & (y[1:] >= threshold))[0]
        index = candidates[-1] if len(candidates) else None
    else:
        x, y = day[peak_index:], curve[peak_index:]
        candidates = np.where((y[:-1] >= threshold) & (y[1:] <= threshold))[0]
        index = candidates[0] if len(candidates) else None
    if index is None:
        return np.nan
    x0, x1 = x[index], x[index + 1]
    y0, y1 = y[index], y[index + 1]
    if y1 == y0:
        return float((x0 + x1) / 2)
    return float(x0 + (threshold - y0) * (x1 - x0) / (y1 - y0))


def fit_curves(curves: list[CurveRecord]) -> tuple[pd.DataFrame, np.ndarray]:
    rows: list[dict[str, Any]] = []
    fitted_calendar: list[np.ndarray] = []
    fine_grid = np.linspace(235, 350, 1151)
    for record in curves:
        values = record.daily_curve
        baseline0 = float(np.quantile(values, 0.10))
        amplitude0 = max(float(np.quantile(values, 0.90) - baseline0), 0.02)
        peak0 = record.raw_peak_doy
        initial = np.array([baseline0, min(amplitude0 * 1.4, 0.4), max(245, peak0 - 12), 24, 0.15, 0.15])
        lower = np.array([max(0, baseline0 - 0.15), 0.005, 240, 3, 0.015, 0.015])
        upper = np.array([min(0.7, baseline0 + 0.15), 0.8, 325, 80, 1.5, 1.5])
        fit = least_squares(
            lambda parameters: double_logistic(parameters, AUTUMN_GRID) - values,
            initial,
            bounds=(lower, upper),
            loss="soft_l1",
            f_scale=0.01,
            max_nfev=5000,
        )
        prediction = double_logistic(fit.x, AUTUMN_GRID)
        residual = prediction - values
        rmse = float(np.sqrt(np.mean(residual**2)))
        denominator = float(np.sum((values - values.mean()) ** 2))
        r_squared = float(1 - np.sum(residual**2) / denominator) if denominator > 0 else np.nan
        fine_curve = double_logistic(fit.x, fine_grid)
        peak_index = int(np.argmax(fine_curve))
        peak_day = float(fine_grid[peak_index])
        baseline = float(fit.x[0])
        realized_amplitude = float(fine_curve[peak_index] - baseline)
        thresholds = {fraction: baseline + fraction * realized_amplitude for fraction in (0.1, 0.5, 0.9)}
        rise10 = _crossing(fine_grid, fine_curve, thresholds[0.1], True, peak_index)
        rise50 = _crossing(fine_grid, fine_curve, thresholds[0.5], True, peak_index)
        rise90 = _crossing(fine_grid, fine_curve, thresholds[0.9], True, peak_index)
        fall90 = _crossing(fine_grid, fine_curve, thresholds[0.9], False, peak_index)
        fall50 = _crossing(fine_grid, fine_curve, thresholds[0.5], False, peak_index)
        fall10 = _crossing(fine_grid, fine_curve, thresholds[0.1], False, peak_index)
        near_bound = bool(np.any(np.isclose(fit.x, lower, rtol=0, atol=0.02)) or np.any(np.isclose(fit.x, upper, rtol=0, atol=0.02)))
        rows.append({
            "site_id": record.site_id,
            "platform_id": record.platform_id,
            "year": record.year,
            "observed_days": record.observed_days,
            "longest_gap_days": record.longest_gap_days,
            "raw_amplitude": record.raw_amplitude,
            "raw_peak_doy": record.raw_peak_doy,
            "fit_baseline": fit.x[0],
            "fit_scale_parameter": fit.x[1],
            "fit_onset_parameter": fit.x[2],
            "fit_logistic_separation": fit.x[3],
            "fit_rise_rate_parameter": fit.x[4],
            "fit_fall_rate_parameter": fit.x[5],
            "fit_rmse": rmse,
            "fit_r_squared": r_squared,
            "fit_parameter_near_bound": near_bound,
            "peak_doy": peak_day,
            "peak_intensity_above_baseline": realized_amplitude,
            "rise_10_doy": rise10,
            "rise_50_doy": rise50,
            "rise_90_doy": rise90,
            "fall_90_doy": fall90,
            "fall_50_doy": fall50,
            "fall_10_doy": fall10,
            "fwhm_days": fall50 - rise50,
            "rise_10_90_days": rise90 - rise10,
            "fall_90_10_days": fall10 - fall90,
            "shape_asymmetry_days": (fall10 - fall90) - (rise90 - rise10),
            "fit_suitable_for_derived_metrics": bool(r_squared >= 0.5 and realized_amplitude >= 0.02),
        })
        fitted_calendar.append(prediction)
    return pd.DataFrame(rows), np.vstack(fitted_calendar)


def functional_pca(curves: list[CurveRecord], parameters: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    eligible = parameters["fit_suitable_for_derived_metrics"].to_numpy()
    selected_records = [record for record, keep in zip(curves, eligible, strict=True) if keep]
    selected_parameters = parameters.loc[eligible].reset_index(drop=True)
    calendar_curves = np.vstack([record.daily_curve for record in selected_records])
    lower = np.quantile(calendar_curves, 0.05, axis=1, keepdims=True)
    upper = np.quantile(calendar_curves, 0.95, axis=1, keepdims=True)
    calendar_normalized = np.clip((calendar_curves - lower) / np.maximum(upper - lower, 1e-6), -0.25, 1.25)

    aligned_curves: list[np.ndarray] = []
    aligned_indices: list[int] = []
    for index, (record, row) in enumerate(zip(selected_records, selected_parameters.itertuples(index=False), strict=True)):
        peak = record.raw_peak_doy
        target = peak + RELATIVE_GRID
        if target.min() < AUTUMN_GRID.min() or target.max() > AUTUMN_GRID.max():
            continue
        curve = PchipInterpolator(AUTUMN_GRID, calendar_normalized[index])(target)
        aligned_curves.append(curve)
        aligned_indices.append(index)
    aligned_matrix = np.vstack(aligned_curves)

    calendar_pca = PCA(n_components=min(10, len(calendar_normalized), calendar_normalized.shape[1])).fit(calendar_normalized)
    aligned_pca = PCA(n_components=min(10, len(aligned_matrix), aligned_matrix.shape[1])).fit(aligned_matrix)
    calendar_scores = calendar_pca.transform(calendar_normalized)
    aligned_scores = aligned_pca.transform(aligned_matrix)

    score_rows = selected_parameters[["site_id", "platform_id", "year"]].copy()
    for index in range(min(6, calendar_scores.shape[1])):
        score_rows[f"calendar_pc{index + 1}"] = calendar_scores[:, index]
    aligned_frame = selected_parameters.loc[aligned_indices, ["site_id", "platform_id", "year"]].copy()
    for index in range(min(6, aligned_scores.shape[1])):
        aligned_frame[f"aligned_pc{index + 1}"] = aligned_scores[:, index]
    score_rows = score_rows.merge(aligned_frame, on=["site_id", "platform_id", "year"], how="left")

    loading_rows: list[dict[str, Any]] = []
    for representation, grid, pca in [
        ("calendar_normalized", AUTUMN_GRID, calendar_pca),
        ("peak_aligned_normalized", RELATIVE_GRID, aligned_pca),
    ]:
        for component in range(min(6, pca.components_.shape[0])):
            for day, loading in zip(grid, pca.components_[component], strict=True):
                loading_rows.append({
                    "representation": representation,
                    "component": component + 1,
                    "day": day,
                    "loading": loading,
                    "explained_variance_ratio": pca.explained_variance_ratio_[component],
                })
    metadata = {
        "calendar_mean": calendar_pca.mean_,
        "calendar_components": calendar_pca.components_,
        "calendar_explained": calendar_pca.explained_variance_ratio_,
        "calendar_scores": calendar_scores,
        "calendar_normalized": calendar_normalized,
        "aligned_mean": aligned_pca.mean_,
        "aligned_components": aligned_pca.components_,
        "aligned_explained": aligned_pca.explained_variance_ratio_,
        "aligned_scores": aligned_scores,
        "aligned_matrix": aligned_matrix,
        "selected_parameters": selected_parameters,
    }
    return score_rows, pd.DataFrame(loading_rows), metadata


def cluster_diagnostics(metadata: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = metadata["calendar_scores"]
    cumulative = np.cumsum(metadata["calendar_explained"])
    component_count = min(max(int(np.searchsorted(cumulative, 0.90) + 1), 2), 8)
    x = scores[:, :component_count]
    sites = metadata["selected_parameters"]["site_id"].to_numpy()
    diagnostics: list[dict[str, Any]] = []
    assignments = metadata["selected_parameters"][["site_id", "platform_id", "year"]].copy()
    rng = np.random.default_rng(20261115)
    best_silhouette = -np.inf
    best_labels: np.ndarray | None = None
    best_k = 0
    for clusters in range(2, 7):
        model = KMeans(n_clusters=clusters, n_init=50, random_state=42).fit(x)
        labels = model.labels_
        silhouette = float(silhouette_score(x, labels))
        stability: list[float] = []
        for _ in range(50):
            perturbation = rng.normal(0, 0.05, size=x.shape) * np.maximum(x.std(axis=0), 1e-6)
            perturbed_labels = KMeans(n_clusters=clusters, n_init=20, random_state=int(rng.integers(1_000_000))).fit_predict(x + perturbation)
            stability.append(adjusted_rand_score(labels, perturbed_labels))
        site_association = float(adjusted_mutual_info_score(sites, labels))
        diagnostics.append({
            "cluster_count": clusters,
            "pca_component_count": component_count,
            "cumulative_variance_used": float(cumulative[component_count - 1]),
            "silhouette_score": silhouette,
            "mean_perturbation_adjusted_rand": float(np.mean(stability)),
            "adjusted_mutual_information_with_site": site_association,
        })
        assignments[f"cluster_k{clusters}"] = labels + 1
        if silhouette > best_silhouette:
            best_silhouette = silhouette
            best_labels = labels + 1
            best_k = clusters
    assignments["best_silhouette_cluster"] = best_labels
    assignments["best_cluster_count"] = best_k
    return pd.DataFrame(diagnostics), assignments


def environmental_summaries() -> pd.DataFrame:
    environment = _read("environment_daily")
    environment["year"] = environment["date"].dt.year
    rows: list[dict[str, Any]] = []
    for (site_id, year), group in environment.groupby(["site_id", "year"]):
        summer = group[group["date"].dt.month.isin([6, 7, 8])]
        autumn = group[group["date"].dt.month.isin([9, 10])]
        rows.append({
            "site_id": site_id,
            "year": int(year),
            "summer_temperature_mean_c": summer["temperature_mean_c"].mean(),
            "summer_precipitation_mm": summer["precipitation_mm_day"].sum(),
            "summer_vpd_kpa": summer["vpd_kpa"].mean(),
            "summer_radiation_w_m2": summer["shortwave_radiation_w_m2"].mean(),
            "sep_oct_temperature_mean_c": autumn["temperature_mean_c"].mean(),
            "sep_oct_precipitation_mm": autumn["precipitation_mm_day"].sum(),
            "sep_oct_vpd_kpa": autumn["vpd_kpa"].mean(),
            "sep_oct_radiation_w_m2": autumn["shortwave_radiation_w_m2"].mean(),
            "sep_oct_frost_days": autumn["frost_day"].sum(),
            "sep_oct_cooling_degree_days": autumn["cooling_degree_days_base15_c"].sum(),
        })
    return pd.DataFrame(rows)


def _benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    values = p_values.to_numpy(dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0, 1)
    return pd.Series(result, index=p_values.index)


def environment_associations(parameters: pd.DataFrame, scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    shape = parameters[parameters["fit_suitable_for_derived_metrics"]].merge(scores, on=["site_id", "platform_id", "year"], how="left")
    targets = [
        "peak_doy", "peak_intensity_above_baseline", "fwhm_days",
        "rise_10_90_days", "fall_90_10_days", "shape_asymmetry_days",
        "calendar_pc1", "calendar_pc2", "calendar_pc3",
    ]
    # Avoid pseudo-replication: environmental summaries occur once per site-year.
    site_year = shape.groupby(["site_id", "year"], as_index=False)[targets].median()
    environment = environmental_summaries()
    merged = site_year.merge(environment, on=["site_id", "year"], how="inner")
    predictors = [column for column in environment.columns if column not in {"site_id", "year"}]
    for column in predictors + targets:
        merged[f"{column}_anomaly"] = merged[column] - merged.groupby("site_id")[column].transform("mean")

    rows: list[dict[str, Any]] = []
    rng = np.random.default_rng(20261115)
    for target in targets:
        for predictor in predictors:
            subset = merged[[f"{target}_anomaly", f"{predictor}_anomaly", "site_id"]].dropna()
            rho, p_value = spearmanr(subset.iloc[:, 0], subset.iloc[:, 1])
            x = subset.iloc[:, 1].to_numpy()
            y = subset.iloc[:, 0].to_numpy()
            site_values = subset["site_id"].to_numpy()
            site_indices = [np.where(site_values == site)[0] for site in np.unique(site_values)]
            exceedances = 0
            permutation_count = 499
            for _ in range(permutation_count):
                permuted = x.copy()
                for indices in site_indices:
                    permuted[indices] = rng.permutation(permuted[indices])
                permutation_rho = spearmanr(y, permuted).statistic
                exceedances += int(abs(permutation_rho) >= abs(rho))
            permutation_p = (exceedances + 1) / (permutation_count + 1)
            rows.append({
                "target": target,
                "predictor": predictor,
                "n_site_years": len(subset),
                "site_count": subset["site_id"].nunique(),
                "within_site_spearman_rho": rho,
                "nominal_p_value": p_value,
                "within_site_permutation_p_value": permutation_p,
            })
    associations = pd.DataFrame(rows)
    associations["fdr_q_value"] = _benjamini_hochberg(associations["within_site_permutation_p_value"])
    associations["exploratory_fdr_0_10"] = associations["fdr_q_value"].le(0.10)

    x_columns = [f"{column}_anomaly" for column in predictors]
    predictive_rows: list[dict[str, Any]] = []
    logo = LeaveOneGroupOut()
    for target in targets:
        subset = merged[["site_id", f"{target}_anomaly"] + x_columns].dropna()
        x = subset[x_columns].to_numpy()
        y = subset[f"{target}_anomaly"].to_numpy()
        groups = subset["site_id"].to_numpy()
        if len(np.unique(groups)) < 3:
            continue
        model = make_pipeline(StandardScaler(), Ridge(alpha=10.0))
        prediction = cross_val_predict(model, x, y, groups=groups, cv=logo)
        predictive_rows.append({
            "target": target,
            "n_site_years": len(subset),
            "site_count": len(np.unique(groups)),
            "leave_one_site_out_r_squared": r2_score(y, prediction),
            "mae": float(np.mean(np.abs(y - prediction))),
            "baseline_zero_anomaly_mae": float(np.mean(np.abs(y))),
            "ridge_alpha": 10.0,
        })
    return associations, pd.DataFrame(predictive_rows), merged


def component_interpretation(parameters: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    merged = parameters.merge(scores, on=["site_id", "platform_id", "year"], how="inner")
    metrics = ["peak_doy", "peak_intensity_above_baseline", "fwhm_days", "rise_10_90_days", "fall_90_10_days", "shape_asymmetry_days"]
    pcs = [column for column in scores if column.startswith("calendar_pc") or column.startswith("aligned_pc")]
    rows: list[dict[str, Any]] = []
    for pc in pcs:
        for metric in metrics:
            subset = merged[[pc, metric]].dropna()
            rho, p_value = spearmanr(subset[pc], subset[metric])
            rows.append({"component": pc, "metric": metric, "spearman_rho": rho, "p_value": p_value, "n": len(subset)})
    return pd.DataFrame(rows)


def create_figures(
    curves: list[CurveRecord],
    parameters: pd.DataFrame,
    metadata: dict[str, Any],
    diagnostics: pd.DataFrame,
    associations: pd.DataFrame,
) -> list[Path]:
    output = FIGURE_DIR / "shape_research"
    output.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    paths: list[Path] = []

    examples = parameters.sort_values("fit_r_squared").iloc[[0, 1, -4, -3, -2, -1]]
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=True)
    curve_lookup = {(record.site_id, record.platform_id, record.year): record for record in curves}
    for ax, row in zip(axes.ravel(), examples.itertuples(index=False), strict=True):
        record = curve_lookup[(row.site_id, row.platform_id, row.year)]
        fit_parameters = np.array([row.fit_baseline, row.fit_scale_parameter, row.fit_onset_parameter, row.fit_logistic_separation, row.fit_rise_rate_parameter, row.fit_fall_rate_parameter])
        ax.plot(AUTUMN_GRID, record.daily_curve, color="#555555", label="smoothed observation")
        ax.plot(AUTUMN_GRID, double_logistic(fit_parameters, AUTUMN_GRID), color="#d95f02", label="double logistic")
        ax.set_title(f"{row.site_id}, {row.year}\n$R^2$={row.fit_r_squared:.2f}")
        ax.set_xlabel("Day of year")
        ax.set_ylabel("RCC90")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Double-logistic fit examples: weakest and strongest fits")
    fig.tight_layout()
    path = output / "double_logistic_fit_examples.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    explained = metadata["calendar_explained"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes[0, 0].bar(np.arange(1, min(8, len(explained)) + 1), explained[:8] * 100, color="#4575b4")
    axes[0, 0].set(title="Calendar-aligned fPCA", xlabel="Component", ylabel="Variance explained (%)")
    mean = metadata["calendar_mean"]
    for component in range(3):
        scale = math.sqrt(metadata["calendar_explained"][component] * len(AUTUMN_GRID))
        axes[0, 1].plot(AUTUMN_GRID, mean + scale * metadata["calendar_components"][component], label=f"PC{component + 1} +")
        axes[0, 1].plot(AUTUMN_GRID, mean - scale * metadata["calendar_components"][component], linestyle="--", alpha=0.7)
    axes[0, 1].set(title="Calendar modes (± direction)", xlabel="Day of year", ylabel="Normalized RCC90")
    axes[0, 1].legend(fontsize=8)
    aligned_explained = metadata["aligned_explained"]
    axes[1, 0].bar(np.arange(1, min(8, len(aligned_explained)) + 1), aligned_explained[:8] * 100, color="#66a61e")
    axes[1, 0].set(title="Peak-aligned fPCA", xlabel="Component", ylabel="Variance explained (%)")
    aligned_mean = metadata["aligned_mean"]
    for component in range(3):
        scale = math.sqrt(aligned_explained[component] * len(RELATIVE_GRID))
        axes[1, 1].plot(RELATIVE_GRID, aligned_mean + scale * metadata["aligned_components"][component], label=f"PC{component + 1} +")
        axes[1, 1].plot(RELATIVE_GRID, aligned_mean - scale * metadata["aligned_components"][component], linestyle="--", alpha=0.7)
    axes[1, 1].set(title="Peak-aligned shape modes", xlabel="Days relative to observed RCC peak", ylabel="Normalized RCC90")
    axes[1, 1].legend(fontsize=8)
    fig.tight_layout()
    path = output / "functional_pca_modes.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    axes[0].plot(diagnostics["cluster_count"], diagnostics["silhouette_score"], marker="o")
    axes[0].axhline(0.4, color="grey", linestyle="--", linewidth=1)
    axes[0].set(title="Cluster separation", xlabel="Number of clusters", ylabel="Silhouette")
    axes[1].plot(diagnostics["cluster_count"], diagnostics["mean_perturbation_adjusted_rand"], marker="o", color="#1b9e77")
    axes[1].set(title="Perturbation stability", xlabel="Number of clusters", ylabel="Adjusted Rand")
    axes[2].plot(diagnostics["cluster_count"], diagnostics["adjusted_mutual_information_with_site"], marker="o", color="#d95f02")
    axes[2].set(title="Association with site", xlabel="Number of clusters", ylabel="Adjusted mutual information")
    fig.tight_layout()
    path = output / "cluster_diagnostics.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    heat = associations.pivot(index="target", columns="predictor", values="within_site_spearman_rho")
    fig, ax = plt.subplots(figsize=(13, 7))
    image = ax.imshow(heat.to_numpy(), cmap="coolwarm", vmin=-0.5, vmax=0.5, aspect="auto")
    ax.set_xticks(np.arange(len(heat.columns)), heat.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(heat.index)), heat.index)
    fig.colorbar(image, ax=ax, label="Within-site Spearman rho")
    ax.set_title("Exploratory environmental associations with curve properties")
    fig.tight_layout()
    path = output / "environment_shape_associations.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)
    return paths


def write_report(
    audit: pd.DataFrame,
    parameters: pd.DataFrame,
    metadata: dict[str, Any],
    diagnostics: pd.DataFrame,
    associations: pd.DataFrame,
    predictive: pd.DataFrame,
    interpretation: pd.DataFrame,
) -> Path:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    eligible = int(audit["shape_analysis_eligible"].sum())
    suitable = parameters["fit_suitable_for_derived_metrics"]
    good_70 = float(parameters["fit_r_squared"].ge(0.70).mean())
    bound_fraction = float(parameters["fit_parameter_near_bound"].mean())
    cumulative = np.cumsum(metadata["calendar_explained"])
    cumulative_aligned = np.cumsum(metadata["aligned_explained"])
    best = diagnostics.loc[diagnostics["silhouette_score"].idxmax()]
    significant = associations[associations["exploratory_fdr_0_10"]].sort_values("fdr_q_value")
    pc_labels: list[str] = []
    for pc in ["calendar_pc1", "calendar_pc2", "calendar_pc3", "aligned_pc1", "aligned_pc2", "aligned_pc3"]:
        subset = interpretation[interpretation["component"].eq(pc)].copy()
        if subset.empty:
            continue
        strongest = subset.loc[subset["spearman_rho"].abs().idxmax()]
        pc_labels.append(f"- `{pc}` was most strongly associated with `{strongest.metric}` (Spearman ρ={strongest.spearman_rho:.2f}); component signs are arbitrary.")
    lines = [
        "# Research audit: discover autumn RCC trajectory shapes before explaining them",
        "",
        "## Conclusion",
        "",
        "The proposed two-stage strategy is scientifically stronger than direct daily RCC-on-weather regression, provided that uncertainty, camera effects, and site structure remain in the model. The current evidence supports **continuous functional modes plus interpretable curve summaries** as the next step. It does not yet justify a mixture-of-experts model with named biological curve classes.",
        "",
        "## Evidence from the assembled data",
        "",
        f"The audit considered {len(audit)} platform-years with at least one Sep-Nov RCC90 observation. {eligible} passed prespecified shape-screening rules for this audit: at least 55 observed dates, no gap longer than 14 days, and a 5th-to-95th percentile RCC90 amplitude of at least 0.025.",
        "",
        f"A robust six-parameter difference-of-logistics model fitted to daily interpolated, nine-day Savitzky–Golay-smoothed curves achieved R² ≥ 0.70 for {good_70:.1%} of eligible curves. This in-sample statistic describes the smoothed representation and is not forecast accuracy. Moreover, {bound_fraction:.1%} of fits placed at least one parameter near its optimization bound. The logistic family is therefore a useful summary model, not a universally identified biological mechanism. Derived peak, half-maximum width, and 10–90% durations are safer downstream quantities than raw scale/separation coefficients.",
        "",
        f"For normalized calendar-aligned curves, the first three functional components explained {cumulative[2]:.1%} of variation. After aligning curves at their observed RCC peak, the first three components explained {cumulative_aligned[2]:.1%}. The latter analysis isolates rise/decline shape from calendar timing.",
        "",
        "Component interpretation against curve summaries:",
        "",
        *pc_labels,
        "",
        "## Are there discrete trajectory types?",
        "",
        f"K-means solutions with 2–6 clusters were compared on the functional scores. The best silhouette occurred at K={int(best.cluster_count)} with silhouette={best.silhouette_score:.2f}, perturbation stability ARI={best.mean_perturbation_adjusted_rand:.2f}, and adjusted mutual information with site={best.adjusted_mutual_information_with_site:.2f}.",
        "",
        "A silhouette below roughly 0.4 is weak evidence for sharply separated types. Strong site association would additionally suggest that clusters encode camera/species/location identity rather than recurrent weather-driven regimes. Unless both separation and out-of-sample gating performance improve, mixture weights should not be interpreted as biological state probabilities.",
        "",
        "## Environmental explanation stage",
        "",
        "Environmental summaries used fixed, interpretable windows rather than windows defined by the observed peak: June-August and September-October temperature, precipitation, VPD, and radiation, plus autumn frost days and cooling-degree accumulation. Both predictors and curve targets were demeaned within site before association testing, so the correlations address interannual anomalies rather than permanent differences among locations. Significance was evaluated by permuting environmental anomalies within each site and controlling the false-discovery rate across all tested pairs.",
        "",
    ]
    if significant.empty:
        lines.append("No curve-property/environment association survived a 10% false-discovery-rate threshold. That result favors shrinkage and out-of-sample evaluation over selecting variables from nominal correlations.")
    else:
        lines += [
            f"{len(significant)} associations survived the exploratory 10% false-discovery-rate threshold:",
            "",
            "| Curve property | Environmental anomaly | Spearman ρ | q-value | N |",
            "|---|---|---:|---:|---:|",
        ]
        for row in significant.itertuples(index=False):
            lines.append(f"| {row.target} | {row.predictor} | {row.within_site_spearman_rho:.2f} | {row.fdr_q_value:.3f} | {row.n_site_years} |")
    lines += [
        "",
        "These rows are not independent mechanisms. September-October mean temperature, cooling accumulation, frost frequency, and VPD are mutually correlated descriptions of the same seasonal weather. Their associations chiefly support a timing signal: warmer autumns tended to have later RCC peaks, whereas greater cooling and more frost tended to accompany earlier peaks. The full September-October windows also overlap the outcome season, so this is explanatory hindcasting, not a real-time forecast experiment. Forecast models must truncate every environmental summary at the issue date.",
        "",
        "Leave-one-site-out ridge tests provide a stricter assessment of whether environmental-anomaly relationships transfer to a new location:",
        "",
        "| Target | CV R² | CV MAE | Zero-anomaly baseline MAE |",
        "|---|---:|---:|---:|",
    ]
    for row in predictive.itertuples(index=False):
        lines.append(f"| {row.target} | {row.leave_one_site_out_r_squared:.2f} | {row.mae:.2f} | {row.baseline_zero_anomaly_mae:.2f} |")
    lines += [
        "",
        "These predictive tests are exploratory because Daymet is gridded, sites have unequal record lengths, and quantitative species composition is not yet available.",
        "",
        "## Recommended model sequence",
        "",
        "1. Fit a hierarchical latent curve to daily RCC90 and GCC90 jointly, with site, platform, ROI, and year effects and observation-level QC. Use a double logistic or monotone-spline rise/fall representation and retain posterior uncertainty.",
        "2. Report derived quantities—peak date, realized amplitude, half-maximum width, rise duration, and fall duration—rather than interpreting all raw nonlinear coefficients.",
        "3. Use calendar-aligned fPCA to discover timing-plus-shape modes, and peak-aligned fPCA to test for residual shape variation.",
        "4. Regress site-demeaned scores/parameters on prespecified environmental anomalies with hierarchical partial pooling. Retain site climatology, elevation, and eventual composition measurements as between-site predictors.",
        "5. Compare against a direct functional weather model under blocked validation. The two-stage method is preferred only if it improves held-out curve reconstruction or scientific interpretability without discarding parameter uncertainty.",
        "6. Attempt a mixture of experts only if stable clusters recur across platforms/sites and a gating model predicts them out of sample better than continuous-score regression.",
        "",
        "## Important methodological qualifications",
        "",
        "- The two-stage plug-in approach can understate uncertainty if fitted curve parameters are treated as exact. A hierarchical joint model or multiple posterior draws should propagate first-stage uncertainty into the environmental stage.",
        "- RCC90 is a camera color coordinate, not a universal latent foliage-color state. GCC decline, field coloration, leaf fall, and satellite greendown require separate observation functions.",
        "- Amplitude is particularly sensitive to camera exposure, view composition, and species pigments. Cross-site amplitude comparisons require platform effects and preferably raw-image color calibration.",
        "- Functional PCA components are continuous axes with arbitrary sign. Naming them ‘early’, ‘broad’, or ‘sharp’ requires correlation with independent, interpretable metrics.",
        "- Photoperiod is deterministic from latitude and date. Its causal role cannot be separated from seasonal time using a simple within-year correlation; cross-site/species structure or process-based constraints are required.",
        "",
        "## Literature basis",
        "",
        "- Double-logistic functions are established phenology summaries with interpretable seasonal parameters, but model choice changes transition estimates: Beck et al. (2006), https://doi.org/10.1016/j.rse.2005.10.021; Atkinson et al. (2012), https://doi.org/10.1016/j.rse.2012.04.001.",
        "- PhenoCam transition products use smoothed GCC stages and report amplitude thresholds with confidence intervals rather than treating raw daily extrema as exact: https://daac.ornl.gov/VEGETATION/guides/PhenoCam_V2.html.",
        "- Bayesian hierarchical phenology models demonstrate why sparse observations and parameter uncertainty should be modeled explicitly: Gao, Gray & Reich (2021), https://doi.org/10.1016/j.rse.2021.112484; Babcock, Finley & Looker (2021), https://doi.org/10.1016/j.rse.2021.112471.",
        "- Autumn controls are nonlinear and species/location dependent. Temperature, photoperiod, drought, and prior growing-season conditions can act differently across regimes: Gill et al. (2015), https://doi.org/10.1093/aob/mcv055; Chen et al. (2020), https://doi.org/10.1038/s41558-020-0820-2; Moon et al. (2022), https://doi.org/10.1016/j.agrformet.2022.109026.",
        "",
        "## Decision",
        "",
        "Proceed with the two-stage functional strategy, but use continuous fPCA scores and uncertainty-aware curve summaries first. Treat mixture-of-experts as a later hypothesis test, not the default architecture.",
    ]
    path = REPORT_DIR / "trajectory_shape_research.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run_shape_research() -> dict[str, Path]:
    curves, audit = prepare_curves()
    parameters, _ = fit_curves(curves)
    scores, loadings, metadata = functional_pca(curves, parameters)
    diagnostics, assignments = cluster_diagnostics(metadata)
    associations, predictive, merged = environment_associations(parameters, scores)
    interpretation = component_interpretation(parameters, scores)
    outputs = {
        "shape_curve_audit": audit,
        "curve_shape_parameters": parameters,
        "functional_pca_scores": scores,
        "functional_pca_loadings": loadings,
        "cluster_diagnostics": diagnostics,
        "cluster_assignments": assignments,
        "environment_shape_associations": associations,
        "environment_shape_predictive_cv": predictive,
        "shape_environment_site_years": merged,
        "component_interpretation": interpretation,
    }
    paths: dict[str, Path] = {}
    for name, frame in outputs.items():
        path = TABLE_DIR / f"{name}.parquet"
        write_table(frame, path)
        paths[name] = path
    figures = create_figures(curves, parameters, metadata, diagnostics, associations)
    paths["report"] = write_report(audit, parameters, metadata, diagnostics, associations, predictive, interpretation)
    paths.update({f"figure_{index + 1}": path for index, path in enumerate(figures)})
    return paths

