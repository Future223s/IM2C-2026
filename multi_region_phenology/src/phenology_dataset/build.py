from __future__ import annotations

from datetime import date
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import CACHE_DIR, OUTPUT_DIR, SHARED_PHENOCAM_CACHE, TABLE_DIR, load_config
from .io import CachedSession, sha256, write_table
from .sources import (
    daymet_point,
    derive_environment,
    load_phenocam_catalog,
    phenocam_roi_summary,
    raw_image_url,
)


def _platform_map(config: dict[str, Any]) -> dict[str, str]:
    return {
        platform: site["site_id"]
        for site in config["sites"]
        for platform in site["platforms"]
    }


def _quality(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    for column in ["image_count", "gcc_90", "rcc_90", "snow_flag", "outlierflag_gcc_90"]:
        if column not in result:
            result[column] = np.nan
    result["processed_available"] = (
        pd.to_numeric(result["image_count"], errors="coerce").fillna(0).gt(0)
        & pd.to_numeric(result["gcc_90"], errors="coerce").notna()
        & pd.to_numeric(result["rcc_90"], errors="coerce").notna()
    )
    snow = pd.to_numeric(result["snow_flag"], errors="coerce")
    outlier = pd.to_numeric(result["outlierflag_gcc_90"], errors="coerce")
    result["quality_usable"] = result["processed_available"] & ~snow.isin([1, 5]) & ~outlier.eq(1)
    result["quality_flag"] = "not_flagged"
    result.loc[~result["processed_available"], "quality_flag"] = "processed_value_missing"
    result.loc[snow.eq(1), "quality_flag"] = "snow_or_obscured_manual_flag_1"
    result.loc[snow.eq(5), "quality_flag"] = "snow_on_trees_manual_flag_5"
    result.loc[outlier.eq(1), "quality_flag"] = "gcc90_automated_outlier"
    return result


def build_sites(config: dict[str, Any], cameras: list[dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    camera_by_id = {camera["Sitename"]: camera for camera in cameras}
    site_rows: list[dict[str, Any]] = []
    platform_rows: list[dict[str, Any]] = []
    for site in config["sites"]:
        canonical = camera_by_id[site["canonical_platform"]]
        metadata = canonical.get("sitemetadata") or {}
        site_rows.append({
            "site_id": site["site_id"],
            "site_name": site["name"],
            "research_region": site["research_region"],
            "state": site["state"],
            "candidate_status": site["status"],
            "latitude": float(canonical["Lat"]),
            "longitude": float(canonical["Lon"]),
            "elevation_m": pd.to_numeric(canonical.get("Elev"), errors="coerce"),
            "primary_vegetation_code": metadata.get("primary_veg_type"),
            "secondary_vegetation_code": metadata.get("secondary_veg_type"),
            "forest_description": metadata.get("site_description"),
            "dominant_species_verbatim": metadata.get("dominant_species"),
            "network_ids_json": json.dumps(site.get("network_ids", {}), sort_keys=True),
            "canonical_environment_support": "nearest Daymet 1 km grid cell to canonical platform coordinate",
            "coordinate_source": "PhenoCam camera API",
        })
        for platform_id in site["platforms"]:
            camera = camera_by_id[platform_id]
            camera_metadata = camera.get("sitemetadata") or {}
            platform_rows.append({
                "platform_id": platform_id,
                "site_id": site["site_id"],
                "platform_type": "fixed_repeat_photography_camera",
                "latitude": float(camera["Lat"]),
                "longitude": float(camera["Lon"]),
                "elevation_m": pd.to_numeric(camera.get("Elev"), errors="coerce"),
                "first_observation_date": camera.get("date_first"),
                "last_observation_date": camera.get("date_last"),
                "active_at_catalog_download": bool(camera.get("active")),
                "camera_description": camera_metadata.get("camera_description"),
                "camera_orientation": camera_metadata.get("camera_orientation"),
                "footprint_description": camera_metadata.get("site_description"),
                "source_dataset_id": "phenocam_camera_api",
            })
    return pd.DataFrame(site_rows), pd.DataFrame(platform_rows)


def build_species(sites: pd.DataFrame) -> pd.DataFrame:
    replacements = {"Liriodendron tilipifera": "Liriodendron tulipifera"}
    rows: list[dict[str, Any]] = []
    for site in sites.itertuples(index=False):
        names = [value.strip() for value in str(site.dominant_species_verbatim or "").split(",") if value.strip()]
        for raw_name in names:
            rows.append({
                "species_composition_id": f"{site.site_id}__{len(rows) + 1:04d}",
                "site_id": site.site_id,
                "species_scientific_name": replacements.get(raw_name, raw_name),
                "species_name_verbatim": raw_name,
                "composition_measure": "listed_dominant_species_presence",
                "composition_value": 1.0,
                "composition_unit": "binary_metadata_indicator_not_abundance",
                "measurement_date": pd.NaT,
                "sampling_footprint": "camera site metadata; exact surveyed footprint not supplied",
                "uncertainty": "No basal area, tree-count fraction, or camera-visible fraction is implied.",
                "source_dataset_id": "phenocam_camera_api",
            })
    return pd.DataFrame(rows)


def build_phenocam_tables(
    config: dict[str, Any],
    rois: list[dict[str, Any]],
    client: CachedSession,
    force: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    platform_to_site = _platform_map(config)
    selected = [
        roi for roi in rois
        if roi["site"] in platform_to_site
        and roi["roitype"] == "DB"
        and int(roi["sequence_number"]) >= 1000
    ]
    observations: list[pd.DataFrame] = []
    imagery: list[pd.DataFrame] = []
    roi_rows: list[dict[str, Any]] = []
    provenance: list[dict[str, Any]] = []
    for roi in sorted(selected, key=lambda row: (row["site"], int(row["sequence_number"]))):
        frame, path = phenocam_roi_summary(client, roi, force)
        frame = _quality(frame)
        frame = frame[frame["date"].dt.year <= int(config["analysis_end_year"])].copy()
        site_id = platform_to_site[roi["site"]]
        frame["site_id"] = site_id
        roi_rows.append({
            "roi_id": roi["roi_name"],
            "platform_id": roi["site"],
            "site_id": site_id,
            "roi_type": roi["roitype"],
            "roi_sequence": int(roi["sequence_number"]),
            "description": roi.get("description"),
            "owner_name": roi.get("owner_name"),
            "first_date_catalog": roi.get("first_date"),
            "last_date_catalog": roi.get("last_date"),
            "first_date_downloaded": frame["date"].min() if not frame.empty else pd.NaT,
            "last_date_downloaded": frame["date"].max() if not frame.empty else pd.NaT,
            "roi_page_url": roi.get("roi_page"),
            "one_day_summary_url": roi.get("one_day_summary"),
            "three_day_summary_url": roi.get("three_day_summary"),
            "spatial_support": "camera region of interest; polygon/mask downloadable from ROI page",
        })
        image_frame = pd.DataFrame({
            "image_id": [f"{roi['roi_name']}__{value:%Y%m%d}" for value in frame["date"]],
            "site_id": site_id,
            "platform_id": roi["site"],
            "timestamp_local_standard_date": frame["date"],
            "timestamp_precision": "daily_summary_with_midday_reference_image",
            "timezone_convention": "camera local standard time; no daylight saving shift",
            "roi_id": roi["roi_name"],
            "roi_sequence": int(roi["sequence_number"]),
            "midday_filename": frame.get("midday_filename"),
            "image_count_in_daily_summary": frame.get("image_count"),
            "snow_flag": frame.get("snow_flag"),
            "outlierflag_gcc_90": frame.get("outlierflag_gcc_90"),
            "quality_usable": frame["quality_usable"],
            "quality_flag": frame["quality_flag"],
        })
        image_frame["raw_image_reference"] = [
            raw_image_url(roi["site"], day, filename)
            for day, filename in zip(frame["date"], frame.get("midday_filename"), strict=False)
        ]
        imagery.append(image_frame)
        for variable, definition in [
            ("gcc_90", "90th percentile daily green chromatic coordinate"),
            ("rcc_90", "90th percentile daily red chromatic coordinate"),
        ]:
            obs = pd.DataFrame({
                "site_id": site_id,
                "platform_id": roi["site"],
                "timestamp": frame["date"],
                "timestamp_precision": "day",
                "timezone_convention": "camera local standard date",
                "observation_type": variable,
                "observation_definition": definition,
                "value": pd.to_numeric(frame[variable], errors="coerce"),
                "unit": "unitless_chromatic_coordinate",
                "species_scientific_name": pd.NA,
                "individual_id": pd.NA,
                "spatial_support_type": "camera_roi",
                "spatial_support_id": roi["roi_name"],
                "observation_uncertainty": "No per-observation statistical uncertainty supplied; use QC flags and cross-image dispersion.",
                "quality_usable": frame["quality_usable"],
                "quality_flag": frame["quality_flag"],
                "source_dataset_id": "phenocam_1day_roi_summary",
                "source_record_reference": roi["one_day_summary"],
            })
            obs["observation_id"] = [f"{roi['roi_name']}__{day:%Y%m%d}__{variable}" for day in frame["date"]]
            observations.append(obs)
        provenance.append({
            "dataset_id": f"phenocam_1day__{roi['roi_name']}",
            "dataset_name": "PhenoCam one-day ROI summary",
            "source_url": roi["one_day_summary"],
            "version": "live archive snapshot",
            "download_date": date.today().isoformat(),
            "local_cache_path": str(path),
            "sha256": sha256(path),
            "processing_steps": "Parsed official CSV; retained raw daily values; added non-destructive QC flags; no smoothing or temporal interpolation.",
            "license_or_policy": "PhenoCam Fair Use Data policy; verify attribution requirements before publication",
            "availability_status": "downloaded",
        })
    return (
        pd.concat(observations, ignore_index=True),
        pd.concat(imagery, ignore_index=True),
        pd.DataFrame(roi_rows),
        provenance,
    )


def build_environment(
    config: dict[str, Any],
    sites: pd.DataFrame,
    client: CachedSession,
    force: bool,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    frames: list[pd.DataFrame] = []
    provenance: list[dict[str, Any]] = []
    for site in sites.itertuples(index=False):
        frame, metadata, path, url = daymet_point(
            client,
            site.site_id,
            float(site.latitude),
            float(site.longitude),
            int(config["environment_start_year"]),
            int(config["analysis_end_year"]),
            force,
        )
        frames.append(frame)
        provenance.append({
            "dataset_id": f"daymet_v4r1__{site.site_id}",
            "dataset_name": "Daymet daily surface weather, nearest 1 km pixel",
            "source_url": url,
            "version": metadata.get("Daymet Software Version", "4.0; service snapshot"),
            "download_date": date.today().isoformat(),
            "local_cache_path": str(path),
            "sha256": sha256(path),
            "processing_steps": "Parsed official point extraction; derived Tmean, VPD proxy, RH proxy, frost, GDD, cooling, and trailing lagged windows. No future dates in rolling windows.",
            "license_or_policy": "NASA Earthdata/ORNL DAAC data policy; cite DOI 10.3334/ORNLDAAC/2129",
            "availability_status": "downloaded",
        })
    return derive_environment(pd.concat(frames, ignore_index=True)), provenance


def _gap_stats(dates: pd.Series, start: pd.Timestamp, end: pd.Timestamp) -> tuple[int, int]:
    observed = set(pd.to_datetime(dates).dt.normalize())
    longest = current = 0
    gaps = 0
    for day in pd.date_range(start, end, freq="D"):
        if day in observed:
            if current:
                gaps += int(current > 3)
            current = 0
        else:
            current += 1
            longest = max(longest, current)
    if current:
        gaps += int(current > 3)
    return longest, gaps


def build_coverage(
    config: dict[str, Any],
    sites: pd.DataFrame,
    observations: pd.DataFrame,
    environment: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric = observations[observations["observation_type"].eq("rcc_90")].copy()
    metric["year"] = metric["timestamp"].dt.year
    rows: list[dict[str, Any]] = []
    for site in sites.itertuples(index=False):
        for year in range(int(config["environment_start_year"]), int(config["analysis_end_year"]) + 1):
            start = pd.Timestamp(f"{year}-{config['autumn_start']}")
            end = pd.Timestamp(f"{year}-{config['autumn_end']}")
            group = metric[
                metric["site_id"].eq(site.site_id)
                & metric["timestamp"].between(start, end)
            ]
            usable = group[group["quality_usable"] & group["value"].notna()]
            # Multiple co-located platforms/ROIs are evidence, not extra calendar days.
            unique_days = usable["timestamp"].drop_duplicates()
            longest, gaps = _gap_stats(unique_days, start, end)
            env = environment[
                environment["site_id"].eq(site.site_id)
                & environment["date"].between(start, end)
            ]
            rows.append({
                "site_id": site.site_id,
                "year": year,
                "season": "autumn_sep01_nov30",
                "expected_calendar_days": 91,
                "phenocam_observed_days_any_platform": int(group["timestamp"].nunique()),
                "phenocam_usable_days_any_platform": int(unique_days.nunique()),
                "phenocam_usable_fraction": float(unique_days.nunique() / 91),
                "phenocam_observation_records": int(len(group)),
                "phenocam_platform_count": int(group["platform_id"].nunique()),
                "phenocam_longest_calendar_gap_days": longest,
                "phenocam_gap_runs_over_3_days": gaps,
                "daymet_days": int(env["date"].nunique()),
                "daymet_fraction": float(env["date"].nunique() / 91),
                "phenocam_autumn_sufficient": bool(
                    unique_days.nunique() >= int(config["minimum_usable_autumn_days"])
                    and unique_days.nunique() / 91 >= float(config["minimum_autumn_fraction"])
                    and longest <= int(config["maximum_long_gap_days"])
                ),
            })
    coverage = pd.DataFrame(rows)
    site_years = coverage.copy()
    site_years["site_year_id"] = site_years["site_id"] + "__" + site_years["year"].astype(str)
    site_years["available_modalities"] = np.where(
        site_years["phenocam_observed_days_any_platform"].gt(0),
        "phenocam_color_indices;daymet_environment",
        "daymet_environment",
    )
    site_years["data_completeness_class"] = np.select(
        [
            site_years["phenocam_autumn_sufficient"],
            site_years["phenocam_usable_days_any_platform"].gt(0),
        ],
        ["camera_trajectory_ready", "camera_sparse_evidence"],
        default="environment_only",
    )
    return coverage, site_years


def build_events(observations: pd.DataFrame, coverage: pd.DataFrame) -> pd.DataFrame:
    eligible = coverage[coverage["phenocam_autumn_sufficient"]][["site_id", "year"]]
    rcc = observations[
        observations["observation_type"].eq("rcc_90")
        & observations["quality_usable"]
        & observations["value"].notna()
    ].copy()
    rcc["year"] = rcc["timestamp"].dt.year
    rcc = rcc.merge(eligible, on=["site_id", "year"], how="inner")
    rcc = rcc[rcc["timestamp"].dt.month.isin([9, 10, 11])]
    rows: list[dict[str, Any]] = []
    for (site_id, platform_id, year), group in rcc.groupby(["site_id", "platform_id", "year"]):
        group = group.sort_values("timestamp")
        peak = group.loc[group["value"].idxmax()]
        dates = group["timestamp"].drop_duplicates().sort_values().reset_index(drop=True)
        pos = int((dates == peak["timestamp"]).to_numpy().argmax())
        before = dates.iloc[pos - 1] if pos > 0 else peak["timestamp"]
        after = dates.iloc[pos + 1] if pos + 1 < len(dates) else peak["timestamp"]
        rows.append({
            "phenology_event_id": f"{site_id}__{platform_id}__{year}__rcc90_observed_max",
            "site_id": site_id,
            "platform_id": platform_id,
            "year": int(year),
            "event_definition": "rcc90_observed_autumn_maximum",
            "event_date": peak["timestamp"],
            "event_value": float(peak["value"]),
            "event_value_unit": "unitless_chromatic_coordinate",
            "estimation_methodology": "Date of maximum usable observed daily RCC90 during Sep-Nov; no interpolation or smoothing.",
            "uncertainty_interval_start": before,
            "uncertainty_interval_end": after,
            "uncertainty_interpretation": "Adjacent observed dates bracket sampling resolution only; does not include camera/biological measurement error.",
            "source_dataset_id": "phenocam_1day_roi_summary",
            "biological_caution": "RCC maximum is a camera color metric and is not asserted to equal landscape peak visible color.",
        })
    return pd.DataFrame(rows)


def verified_inventory(sites: pd.DataFrame, rois: pd.DataFrame, coverage: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for site in sites.itertuples(index=False):
        site_rois = rois[rois["site_id"].eq(site.site_id)]
        site_cov = coverage[coverage["site_id"].eq(site.site_id)]
        observed = site_cov[site_cov["phenocam_observed_days_any_platform"].gt(0)]
        rows.append({
            "site_id": site.site_id,
            "dataset_id": "phenocam_1day_roi_summary",
            "dataset_name": "PhenoCam one-day DB ROI summaries",
            "source_url": "https://phenocam.nau.edu/webcam/network/table/",
            "temporal_coverage": (
                f"{site_rois['first_date_downloaded'].min().date()} to {site_rois['last_date_downloaded'].max().date()}"
                if not site_rois.empty else None
            ),
            "frequency": "daily summary from subdaily images",
            "spatial_resolution": "camera ROI",
            "availability": "downloaded",
            "autumn_completeness": f"{int(site_cov['phenocam_autumn_sufficient'].sum())} trajectory-ready site-years; calculated from records",
            "missingness": f"{int(len(observed) - site_cov['phenocam_autumn_sufficient'].sum())} observed but sparse autumn site-years",
            "quality": "illumination, snow, view/ROI changes, and camera footprint; record-level flags retained",
            "relevance": "canopy color trajectory proxy (GCC/RCC), not species-specific",
        })
        rows.append({
            "site_id": site.site_id,
            "dataset_id": "daymet_v4r1_single_pixel",
            "dataset_name": "Daymet daily surface weather V4 R1",
            "source_url": "https://daymet.ornl.gov/web_services",
            "temporal_coverage": "2000-2025 downloaded",
            "frequency": "daily",
            "spatial_resolution": "nearest 1 km grid cell",
            "availability": "downloaded",
            "autumn_completeness": f"{int((site_cov['daymet_fraction'] == 1).sum())}/{len(site_cov)} autumns complete",
            "missingness": "No missing Sep-Nov dates in downloaded period expected; checked from records",
            "quality": "modeled/interpolated grid, not a tower measurement; secondary variables are derived",
            "relevance": "environmental state and lagged histories",
        })
        networks = json.loads(site.network_ids_json)
        if "neon" in networks:
            rows.append({
                "site_id": site.site_id,
                "dataset_id": "NEON_DP1.10055.001",
                "dataset_name": "NEON plant phenology observations",
                "source_url": "https://data.neonscience.org/data-products/DP1.10055.001/RELEASE-2026",
                "temporal_coverage": "site-specific records not yet downloaded; product catalog spans Aug 2013-Dec 2024",
                "frequency": "irregular individual-plant survey bouts",
                "spatial_resolution": "tagged individual or patch",
                "availability": "verified accessible; download requires NEON API token",
                "autumn_completeness": None,
                "missingness": "not calculated; no records fabricated from catalog range",
                "quality": "status/intensity protocol; discrete visits bracket transition dates",
                "relevance": "species- and individual-specific coloration/senescence evidence",
            })
        if "ameriflux" in networks:
            rows.append({
                "site_id": site.site_id,
                "dataset_id": f"ameriflux_{networks['ameriflux']}",
                "dataset_name": "AmeriFlux BASE/FLUXNET candidate",
                "source_url": "https://ameriflux.lbl.gov/sites/site-search/",
                "temporal_coverage": "not downloaded or record-audited in this release",
                "frequency": "half-hourly to daily product dependent",
                "spatial_resolution": "tower footprint",
                "availability": "site identifier verified; product/license must be checked per site",
                "autumn_completeness": None,
                "missingness": "not calculated",
                "quality": "tower footprint and gap-filling vary by product",
                "relevance": "measured meteorology, fluxes, productivity context",
            })
        for dataset_id, name, url, resolution, relevance in [
            ("MCD12Q2_061", "MODIS Land Cover Dynamics V6.1", "https://lpdaac.usgs.gov/products/mcd12q2v061/", "500 m annual phenometrics", "satellite greenness senescence/dormancy; not visible color peak"),
            ("landsat_c2_l2", "Landsat Collection 2 Level-2", "https://www.usgs.gov/landsat-missions/landsat-collection-2", "30 m irregular clear-sky scenes", "surface reflectance and vegetation-index trajectories"),
            ("sentinel2_l2a", "Sentinel-2 Level-2A", "https://dataspace.copernicus.eu/explore-data/data-collections/sentinel-data/sentinel-2", "10-20 m irregular clear-sky scenes", "surface reflectance, red-edge and vegetation-index trajectories"),
        ]:
            rows.append({
                "site_id": site.site_id,
                "dataset_id": dataset_id,
                "dataset_name": name,
                "source_url": url,
                "temporal_coverage": "global product accessible; site-pixel records not extracted in this release",
                "frequency": "annual" if dataset_id.startswith("MCD") else "irregular revisit subject to cloud",
                "spatial_resolution": resolution,
                "availability": "verified accessible; pending site-footprint extraction",
                "autumn_completeness": None,
                "missingness": "not calculated; catalog availability is not treated as clear-sky coverage",
                "quality": "cloud/snow, mixed pixels, view geometry, and land-cover heterogeneity require QC",
                "relevance": relevance,
            })
    return pd.DataFrame(rows)


def build_provenance_catalog(downloaded: list[dict[str, Any]]) -> pd.DataFrame:
    camera_catalog = CACHE_DIR / "phenocam" / "api" / "cameras.json"
    roi_catalog = CACHE_DIR / "phenocam" / "api" / "roilists.json"
    if not camera_catalog.exists():
        camera_catalog = SHARED_PHENOCAM_CACHE / "api" / "cameras.json"
    if not roi_catalog.exists():
        roi_catalog = SHARED_PHENOCAM_CACHE / "api" / "roilists.json"
    catalog = [
        {
            "dataset_id": "phenocam_camera_api",
            "dataset_name": "PhenoCam camera and ROI catalogs",
            "source_url": "https://phenocam.nau.edu/api/",
            "version": "live API snapshot",
            "download_date": date.today().isoformat(),
            "local_cache_path": str(camera_catalog),
            "sha256": sha256(camera_catalog),
            "processing_steps": "Selected configured camera platforms.",
            "license_or_policy": "PhenoCam Fair Use Data policy",
            "availability_status": "downloaded_or_reused_verified_cache",
        },
        {
            "dataset_id": "phenocam_roi_api",
            "dataset_name": "PhenoCam ROI catalog",
            "source_url": "https://phenocam.nau.edu/api/roilists/?limit=3000",
            "version": "live API snapshot",
            "download_date": date.today().isoformat(),
            "local_cache_path": str(roi_catalog),
            "sha256": sha256(roi_catalog),
            "processing_steps": "Selected deciduous-broadleaf ROI sequences >=1000; retained sequence boundaries.",
            "license_or_policy": "PhenoCam Fair Use Data policy",
            "availability_status": "downloaded_or_reused_verified_cache",
        },
        {
            "dataset_id": "NEON_DP1.10055.001",
            "dataset_name": "NEON plant phenology observations RELEASE-2026",
            "source_url": "https://data.neonscience.org/data-products/DP1.10055.001/RELEASE-2026",
            "version": "RELEASE-2026",
            "download_date": None,
            "local_cache_path": None,
            "sha256": None,
            "processing_steps": "Not downloaded; use /data/query with API token and preserve phe_statusintensity individual records.",
            "license_or_policy": "CC BY 4.0",
            "availability_status": "verified_not_downloaded_requires_token",
        },
        {
            "dataset_id": "MCD12Q2_061",
            "dataset_name": "MODIS Land Cover Dynamics V6.1",
            "source_url": "https://lpdaac.usgs.gov/products/mcd12q2v061/",
            "version": "061",
            "download_date": None,
            "local_cache_path": None,
            "sha256": None,
            "processing_steps": "Not extracted; future pipeline must retain QA and pixel support.",
            "license_or_policy": "NASA Earth science data policy",
            "availability_status": "verified_not_downloaded",
        },
    ]
    return pd.DataFrame(downloaded + catalog).drop_duplicates("dataset_id", keep="first")


def run(force: bool = False) -> dict[str, Path]:
    config = load_config()
    client = CachedSession()
    cameras, roi_catalog = load_phenocam_catalog(client, force=force)
    sites, platforms = build_sites(config, cameras)
    species = build_species(sites)
    observations, imagery, rois, phenocam_provenance = build_phenocam_tables(config, roi_catalog, client, force)
    environment, daymet_provenance = build_environment(config, sites, client, force)
    coverage, site_years = build_coverage(config, sites, observations, environment)
    events = build_events(observations, coverage)
    inventory = verified_inventory(sites, rois, coverage)
    provenance = build_provenance_catalog(phenocam_provenance + daymet_provenance)
    tables = {
        "sites": sites,
        "platforms": platforms,
        "camera_rois": rois,
        "site_years": site_years,
        "environment_daily": environment,
        "phenology_observations": observations,
        "phenology_events": events,
        "species_composition": species,
        "imagery_metadata": imagery,
        "data_provenance": provenance,
        "coverage_audit": coverage,
        "dataset_inventory": inventory,
    }
    paths: dict[str, Path] = {}
    for name, frame in tables.items():
        path = TABLE_DIR / f"{name}.parquet"
        write_table(frame, path)
        paths[name] = path
    manifest = {
        "release": "0.1.0",
        "built_at": pd.Timestamp.now(tz="UTC").isoformat(),
        "analysis_end_year": config["analysis_end_year"],
        "table_rows": {name: int(len(frame)) for name, frame in tables.items()},
        "tables": {name: str(path) for name, path in paths.items()},
        "notes": [
            "PhenoCam and Daymet records are downloaded; other modalities remain explicit verified/pending inventory entries.",
            "Observations retain platform, ROI, timestamp, footprint, quality, and provenance identifiers.",
            "No modality is discarded for lack of same-day observations in another modality.",
        ],
    }
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = OUTPUT_DIR / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    paths["manifest"] = manifest_path
    return paths

