from __future__ import annotations

from io import BytesIO, StringIO
import json
from pathlib import Path
import re
from typing import Any

import numpy as np
import pandas as pd

from .config import CACHE_DIR, SHARED_PHENOCAM_CACHE
from .io import CachedSession


PHENOCAM_API = "https://phenocam.nau.edu/api"
PHENOCAM_ARCHIVE = "https://phenocam.nau.edu/data/archive"
DAYMET_API = "https://daymet.ornl.gov/single-pixel/api/data"


def _metadata_path(name: str, force: bool = False) -> Path:
    shared = SHARED_PHENOCAM_CACHE / "api" / name
    return shared if shared.exists() and not force else CACHE_DIR / "phenocam" / "api" / name


def load_phenocam_catalog(client: CachedSession, force: bool = False) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    camera_path = _metadata_path("cameras.json", force)
    roi_path = _metadata_path("roilists.json", force)
    cameras = client.get_json(f"{PHENOCAM_API}/cameras/?limit=2000", camera_path, force)
    rois = client.get_json(f"{PHENOCAM_API}/roilists/?limit=3000", roi_path, force)
    return cameras["results"], rois["results"]


def phenocam_roi_summary(client: CachedSession, roi: dict[str, Any], force: bool = False) -> tuple[pd.DataFrame, Path]:
    name = roi["roi_name"]
    shared = SHARED_PHENOCAM_CACHE / "roi_summaries" / f"{name}_1day.csv"
    path = shared if shared.exists() and not force else CACHE_DIR / "phenocam" / "roi_summaries" / f"{name}_1day.csv"
    content = client.get_bytes(roi["one_day_summary"], path, force)
    frame = pd.read_csv(BytesIO(content), comment="#", na_values=["NA"])
    frame["date"] = pd.to_datetime(frame["date"])
    frame["roi_id"] = name
    frame["roi_sequence"] = int(roi["sequence_number"])
    frame["platform_id"] = roi["site"]
    return frame, path


def _parse_daymet_csv(content: bytes) -> tuple[pd.DataFrame, dict[str, str]]:
    text = content.decode("utf-8")
    lines = text.splitlines()
    header_index = next(i for i, line in enumerate(lines) if line.startswith("year,yday"))
    metadata: dict[str, str] = {}
    for line in lines[:header_index]:
        if ":" in line:
            key, value = line.split(":", 1)
            metadata[key.strip()] = value.strip()
    frame = pd.read_csv(StringIO("\n".join(lines[header_index:])))
    frame.columns = [re.sub(r"\s*\([^)]*\)", "", column).strip() for column in frame.columns]
    return frame, metadata


def daymet_point(
    client: CachedSession,
    site_id: str,
    latitude: float,
    longitude: float,
    start_year: int,
    end_year: int,
    force: bool = False,
) -> tuple[pd.DataFrame, dict[str, str], Path, str]:
    years = ",".join(str(year) for year in range(start_year, end_year + 1))
    url = f"{DAYMET_API}?lat={latitude}&lon={longitude}&vars=dayl,prcp,srad,swe,tmax,tmin,vp&years={years}"
    path = CACHE_DIR / "daymet" / f"{site_id}_{start_year}_{end_year}.csv"
    content = client.get_bytes(url, path, force)
    frame, metadata = _parse_daymet_csv(content)
    frame["site_id"] = site_id
    # Daymet uses a 365-day calendar: leap-year day 365 corresponds to Dec 30.
    frame["date"] = pd.to_datetime(frame["year"].astype(str), format="%Y") + pd.to_timedelta(frame["yday"] - 1, unit="D")
    return frame, metadata, path, url


def derive_environment(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy().sort_values(["site_id", "date"])
    result = result.rename(columns={
        "dayl": "daylength_s",
        "prcp": "precipitation_mm_day",
        "srad": "shortwave_radiation_w_m2",
        "swe": "snow_water_equivalent_kg_m2",
        "tmax": "temperature_max_c",
        "tmin": "temperature_min_c",
        "vp": "vapor_pressure_pa",
    })
    result["temperature_mean_c"] = (result["temperature_min_c"] + result["temperature_max_c"]) / 2
    saturation_kpa = 0.6108 * np.exp(17.27 * result["temperature_mean_c"] / (result["temperature_mean_c"] + 237.3))
    result["vpd_kpa"] = (saturation_kpa - result["vapor_pressure_pa"] / 1000).clip(lower=0)
    result["relative_humidity_proxy_pct"] = (100 * (result["vapor_pressure_pa"] / 1000) / saturation_kpa).clip(0, 100)
    result["frost_day"] = result["temperature_min_c"].lt(0)
    result["growing_degree_days_base5_c"] = (result["temperature_mean_c"] - 5).clip(lower=0)
    result["cooling_degree_days_base15_c"] = (15 - result["temperature_mean_c"]).clip(lower=0)
    for _, indices in result.groupby("site_id", sort=False).groups.items():
        idx = list(indices)
        site = result.loc[idx].copy()
        for window in (3, 7, 14, 21, 30):
            result.loc[idx, f"temperature_mean_c_roll{window}_lag1"] = site["temperature_mean_c"].shift(1).rolling(window, min_periods=max(2, window // 2)).mean().to_numpy()
            result.loc[idx, f"cooling_degree_days_base15_c_sum{window}_lag1"] = site["cooling_degree_days_base15_c"].shift(1).rolling(window, min_periods=max(2, window // 2)).sum().to_numpy()
        for window in (7, 30):
            result.loc[idx, f"precipitation_mm_sum{window}_lag1"] = site["precipitation_mm_day"].shift(1).rolling(window, min_periods=max(2, window // 2)).sum().to_numpy()
            result.loc[idx, f"vpd_kpa_roll{window}_lag1"] = site["vpd_kpa"].shift(1).rolling(window, min_periods=max(2, window // 2)).mean().to_numpy()
    result["source_dataset_id"] = "daymet_v4r1_single_pixel"
    result["measurement_basis"] = "modeled_gridded_from_station_observations"
    result["spatial_support"] = "nearest_1km_daymet_grid_cell"
    result["quality_flag"] = np.where(result.isna().any(axis=1), "contains_missing", "not_flagged")
    return result


def raw_image_url(platform_id: str, date: pd.Timestamp, filename: str | float) -> str | None:
    if pd.isna(filename):
        return None
    return f"{PHENOCAM_ARCHIVE}/{platform_id}/{date.year}/{date.month:02d}/{filename}"

