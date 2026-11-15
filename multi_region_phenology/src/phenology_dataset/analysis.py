from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config import FIGURE_DIR, REPORT_DIR, TABLE_DIR
from .io import write_table


def _read(name: str) -> pd.DataFrame:
    return pd.read_parquet(TABLE_DIR / f"{name}.parquet")


def site_ranking() -> pd.DataFrame:
    sites = _read("sites")
    platforms = _read("platforms")
    coverage = _read("coverage_audit")
    summary = coverage.groupby("site_id", as_index=False).agg(
        first_camera_year=("year", lambda x: int(coverage.loc[x.index][coverage.loc[x.index, "phenocam_observed_days_any_platform"].gt(0)]["year"].min())),
        last_camera_year=("year", lambda x: int(coverage.loc[x.index][coverage.loc[x.index, "phenocam_observed_days_any_platform"].gt(0)]["year"].max())),
        observed_camera_years=("phenocam_observed_days_any_platform", lambda x: int((x > 0).sum())),
        trajectory_ready_years=("phenocam_autumn_sufficient", "sum"),
        median_usable_fraction=("phenocam_usable_fraction", lambda x: float(x[x > 0].median())),
        median_longest_gap_days=("phenocam_longest_calendar_gap_days", lambda x: float(x[x < 91].median())),
    )
    active = platforms.groupby("site_id")["active_at_catalog_download"].any().rename("has_active_camera")
    result = sites.merge(summary, on="site_id").merge(active, on="site_id")
    result["verified_neon_field_candidate"] = result["network_ids_json"].str.contains('"neon"')
    max_years = max(float(result["trajectory_ready_years"].max()), 1)
    result["priority_score"] = (
        45 * result["trajectory_ready_years"] / max_years
        + 25 * result["median_usable_fraction"].fillna(0)
        + 15 * result["has_active_camera"].astype(int)
        + 15 * result["verified_neon_field_candidate"].astype(int)
    )
    result = result.sort_values(["priority_score", "trajectory_ready_years"], ascending=False).reset_index(drop=True)
    result.insert(0, "rank", np.arange(1, len(result) + 1))
    return result


def variability_summary() -> tuple[pd.DataFrame, pd.DataFrame]:
    events = _read("phenology_events")
    events["event_doy"] = events["event_date"].dt.day_of_year
    site_year = events.groupby(["site_id", "year"], as_index=False).agg(
        event_doy=("event_doy", "median"),
        platform_count=("platform_id", "nunique"),
        platform_min_doy=("event_doy", "min"),
        platform_max_doy=("event_doy", "max"),
    )
    site_year["platform_spread_days"] = site_year["platform_max_doy"] - site_year["platform_min_doy"]
    by_site = site_year.groupby("site_id", as_index=False).agg(
        site_year_count=("year", "nunique"),
        mean_event_doy=("event_doy", "mean"),
        median_event_doy=("event_doy", "median"),
        within_site_sd_days=("event_doy", "std"),
        earliest_event_doy=("event_doy", "min"),
        latest_event_doy=("event_doy", "max"),
        median_platform_spread_days=("platform_spread_days", "median"),
        max_platform_spread_days=("platform_spread_days", "max"),
    )
    overall = pd.DataFrame([{
        "metric": "rcc90_observed_autumn_maximum",
        "site_count": int(by_site["site_id"].nunique()),
        "site_year_count": int(len(site_year)),
        "median_within_site_sd_days": float(by_site["within_site_sd_days"].median()),
        "between_site_sd_of_site_means_days": float(by_site["mean_event_doy"].std()),
        "between_site_range_of_site_means_days": float(by_site["mean_event_doy"].max() - by_site["mean_event_doy"].min()),
        "median_same_site_year_platform_spread_days": float(site_year.loc[site_year["platform_count"] > 1, "platform_spread_days"].median()),
        "caution": "RCC90 maxima are operational camera metrics; platform spread is measurement/footprint disagreement, not pure biological uncertainty.",
    }])
    return by_site, overall


def create_figures() -> list[Path]:
    plt.style.use("seaborn-v0_8-whitegrid")
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    coverage = _read("coverage_audit")
    sites = _read("sites")
    events = _read("phenology_events")
    observations = _read("phenology_observations")

    matrix = coverage.pivot(index="site_id", columns="year", values="phenocam_usable_fraction")
    fig, ax = plt.subplots(figsize=(15, 7))
    image = ax.imshow(matrix.to_numpy(), aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_yticks(np.arange(len(matrix.index)), labels=matrix.index)
    tick_step = max(1, len(matrix.columns) // 13)
    tick_positions = np.arange(0, len(matrix.columns), tick_step)
    ax.set_xticks(tick_positions, labels=matrix.columns[tick_positions], rotation=45, ha="right")
    fig.colorbar(image, ax=ax, label="usable autumn-day fraction")
    ax.set(title="Record-derived PhenoCam autumn coverage (Sep 1–Nov 30)", xlabel="Year", ylabel="Ecological site")
    fig.tight_layout()
    path = FIGURE_DIR / "phenocam_autumn_coverage.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    events["event_doy"] = events["event_date"].dt.day_of_year
    order = events.groupby("site_id")["event_doy"].median().sort_values().index
    fig, ax = plt.subplots(figsize=(12, 7))
    values = [events.loc[events["site_id"].eq(site), "event_doy"].to_numpy() for site in order]
    ax.boxplot(values, vert=False, tick_labels=order, patch_artist=True, boxprops={"facecolor": "#c76b37", "alpha": 0.75})
    rng = np.random.default_rng(42)
    for position, site_values in enumerate(values, start=1):
        ax.scatter(site_values, position + rng.normal(0, 0.035, len(site_values)), color="black", s=10, alpha=0.45)
    ax.set(title="Observed autumn RCC90 maxima by platform-year", xlabel="Day of year", ylabel="Ecological site")
    ax.text(0.01, -0.12, "Camera color metric only; not asserted to equal landscape peak visible color.", transform=ax.transAxes, fontsize=9)
    fig.tight_layout()
    path = FIGURE_DIR / "rcc90_maxima_variation.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    selected_sites = ["harvard_forest", "bartlett_forest", "morgan_monroe", "grsm_neon", "umbs_flux"]
    seasonal = observations[
        observations["site_id"].isin(selected_sites)
        & observations["observation_type"].eq("rcc_90")
        & observations["quality_usable"]
        & observations["timestamp"].dt.year.isin([2010, 2015, 2020, 2025])
        & observations["timestamp"].dt.month.isin([9, 10, 11])
    ].copy()
    seasonal["day_of_year"] = seasonal["timestamp"].dt.day_of_year
    seasonal = seasonal.groupby(["site_id", seasonal["timestamp"].dt.year.rename("year"), "day_of_year"], as_index=False)["value"].median()
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
    axes = axes.ravel()
    for ax, site_id in zip(axes, selected_sites, strict=False):
        subset = seasonal[seasonal["site_id"].eq(site_id)]
        for year, group in subset.groupby("year"):
            ax.plot(group["day_of_year"], group["value"], label=str(year), linewidth=1.4)
        ax.set_title(site_id)
        ax.set_xlabel("Day of year")
        ax.set_ylabel("RCC90")
    axes[-1].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower right", title="Year")
    fig.suptitle("Selected observed autumn RCC90 trajectories")
    fig.tight_layout()
    path = FIGURE_DIR / "seasonal_rcc90_curves.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    return paths


def write_report(ranking: pd.DataFrame, by_site: pd.DataFrame, overall: pd.DataFrame) -> Path:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    sites = _read("sites")
    inventory = _read("dataset_inventory")
    manifest = json.loads((TABLE_DIR.parent / "manifest.json").read_text(encoding="utf-8"))
    merged = ranking.merge(by_site, on="site_id", how="left")
    lines = [
        "# Multi-region autumn phenology dataset — release 0.1",
        "",
        "## What is actually assembled",
        "",
        f"This release contains {manifest['table_rows']['phenology_observations']:,} PhenoCam GCC90/RCC90 observations from {manifest['table_rows']['camera_rois']} deciduous ROIs on {manifest['table_rows']['platforms']} camera platforms, plus {manifest['table_rows']['environment_daily']:,} Daymet daily environmental rows for {manifest['table_rows']['sites']} ecological sites. Coverage ends with the last complete calendar year, 2025.",
        "",
        "NEON individual-plant phenology, AmeriFlux, Landsat, Sentinel-2, and MODIS records are **not** represented as downloaded observations. They are retained in the inventory as verified-accessible or site-identifier candidates with null completeness until their records are ingested and audited.",
        "",
        "## Recommended core and complementary set",
        "",
        "The score ranks the current release only: 45% trajectory-ready camera years, 25% median camera completeness, 15% an active camera, and 15% a verified NEON site identifier. It is a transparent triage score, not a scientific weighting.",
        "",
        "| Rank | Site | Elev. m | Ready camera years | Observed years | Median usable | Active | NEON candidate | Priority score |",
        "|---:|---|---:|---:|---:|---:|:---:|:---:|---:|",
    ]
    for row in merged.itertuples(index=False):
        lines.append(f"| {row.rank} | {row.site_name} | {row.elevation_m:.0f} | {int(row.trajectory_ready_years)} | {int(row.observed_camera_years)} | {row.median_usable_fraction:.1%} | {'yes' if row.has_active_camera else 'no'} | {'yes' if row.verified_neon_field_candidate else 'no'} | {row.priority_score:.1f} |")
    lines += [
        "",
        "A defensible first hierarchical core is Harvard, Bartlett, Morgan–Monroe, UMBS, and the NEON Great Smoky/Oak Ridge sites. Purchase Knob, Shining Rock, and Dolly Sods add high-elevation historical contrast but lack a current camera and currently have no ingested independent field modality. Cary is a strong camera-only secondary record. UMBS FASET is best treated as an experimental complement rather than pooled with the nearby UMBS flux tower.",
        "",
        "## Verified coordinates, forest metadata, and limitations",
        "",
        "Coordinates/elevations and dominant-species lists below come from the downloaded PhenoCam camera catalog. Species are presence-style metadata, not composition fractions.",
        "",
        "| Site | Latitude | Longitude | Elev. m | Vegetation | Dominant species (verbatim) |",
        "|---|---:|---:|---:|---|---|",
    ]
    for row in sites.itertuples(index=False):
        primary = "" if pd.isna(row.primary_vegetation_code) else str(row.primary_vegetation_code)
        secondary = "" if pd.isna(row.secondary_vegetation_code) else str(row.secondary_vegetation_code)
        species = "" if pd.isna(row.dominant_species_verbatim) else str(row.dominant_species_verbatim)
        lines.append(f"| {row.site_name} | {row.latitude:.5f} | {row.longitude:.5f} | {row.elevation_m:.0f} | {primary}/{secondary} | {species} |")
    metric = overall.iloc[0]
    lines += [
        "",
        "## Preliminary modeling-readiness findings",
        "",
        f"Across site-year median RCC90 maxima, the median within-site interannual SD is {metric['median_within_site_sd_days']:.1f} days. The SD of site mean dates is {metric['between_site_sd_of_site_means_days']:.1f} days and their range is {metric['between_site_range_of_site_means_days']:.1f} days. For site-years with multiple camera platforms, the median platform spread is {metric['median_same_site_year_platform_spread_days']:.1f} days.",
        "",
        "This supports a hierarchical model: there is material temporal variation within sites and material spatial structure between sites. But the platform spread also shows that some apparent variability is measurement-footprint disagreement. A latent-state model should therefore include modality/platform observation functions rather than treating each date as interchangeable truth.",
        "",
        "The environmental table includes raw Daymet variables plus trailing, lagged cooling summaries at 3, 7, 14, 21, and 30 days. Here, ‘several scales’ means several look-back windows over the same daily temperature process. Every rolling feature ends one day before its row date, preventing future leakage for date-indexed forecasts.",
        "",
        "## Critical limitations",
        "",
        "- RCC90 maximum is an operational red-chromatic-coordinate event, not a validated universal peak-color date.",
        "- Daymet is modeled/interpolated at 1 km; it is not a direct canopy or tower measurement. Its vapor pressure and radiation are secondary estimates.",
        "- Camera ROIs, satellite pixels, tagged trees, and tower footprints are not coextensive. The schema preserves their identifiers so future joins can express proximity rather than false co-location.",
        "- Camera view/ROI changes are retained as separate `platform_id` and `roi_id` values. Adjacent platforms are not averaged in the observation table.",
        "- Species metadata say which taxa are described as dominant but do not quantify basal area, stem fraction, or camera-visible canopy fraction.",
        "- Satellite clear-sky completeness and NEON individual-tree completeness remain uncalculated until records are downloaded.",
        "",
        "## Next ingestion priority",
        "",
        "1. Download NEON DP1.10055.001 for HARV, BART, GRSM, and ORNL with an API token; retain `individualID`, taxon, phenophase status/intensity, and visit dates.",
        "2. Extract Landsat/Sentinel-2 observations over explicit camera ROI/tower buffers; retain cloud masks, pixel IDs, acquisition times, and footprint geometry.",
        "3. Extract MCD12Q2 V6.1 pixel-level senescence/dormancy metrics and QA, labeling them as greenness phenometrics rather than color peaks.",
        "4. Audit site-specific AmeriFlux licenses and daily meteorological coverage before replacing or supplementing Daymet.",
        "5. Add quantitative vegetation inventories (NEON/FIA/LTER) with composition measure and footprint explicit.",
        "",
        "## Inventory status counts",
        "",
    ]
    for status, count in inventory["availability"].value_counts().items():
        lines.append(f"- {status}: {count}")
    path = REPORT_DIR / "site_selection_and_readiness.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run_analysis() -> dict[str, Path]:
    ranking = site_ranking()
    by_site, overall = variability_summary()
    paths = {
        "site_ranking": TABLE_DIR / "site_selection_ranking.parquet",
        "variability_by_site": TABLE_DIR / "variability_by_site.parquet",
        "variability_overall": TABLE_DIR / "variability_overall.parquet",
    }
    write_table(ranking, paths["site_ranking"])
    write_table(by_site, paths["variability_by_site"])
    write_table(overall, paths["variability_overall"])
    figure_paths = create_figures()
    paths["report"] = write_report(ranking, by_site, overall)
    paths.update({f"figure_{index + 1}": path for index, path in enumerate(figure_paths)})
    return paths

