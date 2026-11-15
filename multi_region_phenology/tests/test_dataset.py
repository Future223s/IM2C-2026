from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "outputs" / "tables"


def read(name: str) -> pd.DataFrame:
    path = TABLES / f"{name}.parquet"
    assert path.exists(), f"Run run_pipeline.py first; missing {path}"
    return pd.read_parquet(path)


def test_required_relational_tables_exist() -> None:
    required = {
        "sites",
        "site_years",
        "environment_daily",
        "phenology_observations",
        "phenology_events",
        "species_composition",
        "imagery_metadata",
        "data_provenance",
    }
    assert required <= {path.stem for path in TABLES.glob("*.parquet")}


def test_primary_keys_and_foreign_keys() -> None:
    sites = read("sites")
    platforms = read("platforms")
    rois = read("camera_rois")
    observations = read("phenology_observations")
    assert sites["site_id"].is_unique
    assert platforms["platform_id"].is_unique
    assert rois["roi_id"].is_unique
    assert observations["observation_id"].is_unique
    site_ids = set(sites["site_id"])
    assert set(platforms["site_id"]) <= site_ids
    assert set(rois["site_id"]) <= site_ids
    assert set(observations["site_id"]) <= site_ids


def test_asynchronous_records_are_not_forced_into_daily_wide_rows() -> None:
    observations = read("phenology_observations")
    assert {"timestamp", "observation_type", "platform_id", "spatial_support_id"} <= set(observations)
    assert set(observations["observation_type"].unique()) == {"gcc_90", "rcc_90"}
    assert observations["timestamp"].notna().all()


def test_coverage_is_record_derived_and_bounded() -> None:
    coverage = read("coverage_audit")
    assert coverage["phenocam_usable_days_any_platform"].between(0, 91).all()
    assert coverage["phenocam_usable_fraction"].between(0, 1).all()
    expected = (
        coverage["phenocam_usable_days_any_platform"].ge(45)
        & coverage["phenocam_usable_fraction"].ge(0.60)
        & coverage["phenocam_longest_calendar_gap_days"].le(14)
    )
    assert expected.equals(coverage["phenocam_autumn_sufficient"])


def test_daymet_has_365_day_calendar_and_unique_site_dates() -> None:
    environment = read("environment_daily")
    assert not environment.duplicated(["site_id", "date"]).any()
    counts = environment.groupby(["site_id", environment["date"].dt.year]).size()
    assert counts.eq(365).all()


def test_trailing_temperature_window_has_no_future_leakage() -> None:
    environment = read("environment_daily")
    site = environment[environment["site_id"].eq("harvard_forest")].sort_values("date").reset_index(drop=True)
    row_index = 50
    expected = site.loc[row_index - 7:row_index - 1, "temperature_mean_c"].mean()
    actual = site.loc[row_index, "temperature_mean_c_roll7_lag1"]
    assert np.isclose(actual, expected)


def test_events_only_use_trajectory_ready_site_years() -> None:
    events = read("phenology_events")
    coverage = read("coverage_audit")
    ready = set(map(tuple, coverage.loc[coverage["phenocam_autumn_sufficient"], ["site_id", "year"]].to_numpy()))
    event_keys = set(map(tuple, events[["site_id", "year"]].to_numpy()))
    assert event_keys <= ready
    assert events["biological_caution"].str.contains("not asserted").all()


def test_high_and_low_elevation_smoky_sites_are_not_merged() -> None:
    sites = read("sites").set_index("site_id")
    high = sites.loc["grsm_purchase_knob"]
    low = sites.loc["grsm_neon"]
    assert high["elevation_m"] - low["elevation_m"] > 900
    assert (high["latitude"], high["longitude"]) != (low["latitude"], low["longitude"])


def test_species_metadata_is_not_mislabeled_as_fraction() -> None:
    species = read("species_composition")
    assert species["composition_measure"].eq("listed_dominant_species_presence").all()
    assert species["composition_unit"].eq("binary_metadata_indicator_not_abundance").all()


def test_pending_sources_have_no_fabricated_completeness() -> None:
    inventory = read("dataset_inventory")
    pending = inventory[~inventory["availability"].eq("downloaded")]
    assert pending["autumn_completeness"].isna().all()


def test_shape_research_preserves_platform_year_identity() -> None:
    parameters = read("curve_shape_parameters")
    assert not parameters.duplicated(["site_id", "platform_id", "year"]).any()
    assert parameters["fit_r_squared"].between(-1, 1).all()
    assert parameters["peak_doy"].between(244, 335).all()


def test_shape_clusters_are_diagnostics_not_species_weights() -> None:
    diagnostics = read("cluster_diagnostics")
    assignments = read("cluster_assignments")
    assert diagnostics["cluster_count"].tolist() == [2, 3, 4, 5, 6]
    assert diagnostics["silhouette_score"].between(-1, 1).all()
    assert "best_silhouette_cluster" in assignments


def test_environment_shape_analysis_avoids_platform_pseudoreplication() -> None:
    site_years = read("shape_environment_site_years")
    associations = read("environment_shape_associations")
    assert not site_years.duplicated(["site_id", "year"]).any()
    assert associations["within_site_permutation_p_value"].between(0, 1).all()
    assert associations["fdr_q_value"].between(0, 1).all()


def test_staged_model_predictions_are_unique_and_complete() -> None:
    predictions = read("staged_model_predictions")
    assert set(predictions["model"]) == {"M1", "M2", "M3", "M4"}
    assert not predictions.duplicated(["model", "site_id", "year"]).any()
    assert predictions.groupby("model").size().eq(158).all()


def test_m1_uses_training_only_location_mean() -> None:
    source = read("shape_environment_site_years")
    predictions = read("staged_model_predictions")
    row = predictions[predictions["model"].eq("M1")].iloc[0]
    expected = source.loc[
        source["site_id"].eq(row["site_id"]) & source["year"].ne(row["year"]),
        "peak_doy",
    ].mean()
    assert np.isclose(row["predicted"], expected)


def test_staged_environmental_model_improves_timing_only() -> None:
    metrics = read("staged_model_metrics").set_index("model")
    dimensions = read("staged_dimension_tests").set_index("target")
    assert metrics.loc["M2", "mae"] < metrics.loc["M1", "mae"]
    supported = set(dimensions.index[dimensions["practically_supported"]])
    assert supported == {"peak_doy"}


def test_all_staged_model_reports_exist() -> None:
    report_dir = ROOT / "outputs" / "reports" / "staged_curve_models"
    expected = {
        "M1_location_baseline.md",
        "M2_shared_temperature_shift.md",
        "M3_multienvironment_timing.md",
        "M4_location_varying_temperature.md",
        "M5_which_curve_dimensions_vary.md",
        "M6_final_selected_logistic_model.md",
        "README.md",
    }
    assert expected <= {path.name for path in report_dir.glob("*.md")}

