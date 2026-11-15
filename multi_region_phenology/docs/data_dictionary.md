# Data dictionary

Release 0.1 stores every table as both Parquet and CSV under `outputs/tables/`. Parquet is canonical because it preserves dates, booleans, and nulls. IDs are stable within this release.

## Relationship map

```text
sites ──< platforms ──< camera_rois
  │           │              │
  │           └──────< imagery_metadata
  │                          │
  ├────< site_years          └────< phenology_observations
  ├────< environment_daily                │
  ├────< species_composition              └────> phenology_events (derived)
  └────< dataset_inventory

data_provenance identifies source snapshots used by all tables.
```

No table requires observations from different modalities to occur on the same date.

## `sites`

One ecological site or intentionally distinct experimental/elevation unit.

| Field | Definition |
|---|---|
| `site_id` | Primary key. |
| `latitude`, `longitude`, `elevation_m` | Canonical location from the configured PhenoCam platform. |
| `research_region`, `state` | Geographic grouping. |
| `primary_vegetation_code`, `secondary_vegetation_code` | PhenoCam vegetation codes; `DB` means deciduous broadleaf. |
| `dominant_species_verbatim` | Catalog text, not a quantitative inventory. |
| `network_ids_json` | Known network identifiers; presence does not imply that records were downloaded. |
| `canonical_environment_support` | Spatial interpretation of the Daymet join. |

## `platforms`

One physical camera platform. A site can have multiple platforms, and nearby experimental or elevation units can remain separate sites.

`platform_id` is the primary key; `site_id` is a foreign key. Coordinates, elevation, activity, observation bounds, orientation, and footprint description are platform-specific.

## `camera_rois`

One documented deciduous-broadleaf region of interest and sequence. `roi_id` is the primary key. ROI sequence changes are never erased. Catalog and actually downloaded date bounds are separate.

## `site_years`

One site × calendar year from 2000–2025. `site_year_id` is the primary key. `data_completeness_class` is:

- `camera_trajectory_ready`: at least 45 unique usable autumn days, at least 60% of Sep 1–Nov 30, and no gap longer than 14 days;
- `camera_sparse_evidence`: at least one camera observation but thresholds not met;
- `environment_only`: no downloaded camera observation in autumn.

Multiple platforms can fill different dates, but their records remain separate in the observation table.

## `environment_daily`

One site × Daymet date. The composite key is (`site_id`, `date`). Daymet supplies a nearest-pixel modeled 1 km grid estimate rather than a point measurement.

| Field | Unit / meaning |
|---|---|
| `temperature_min_c`, `temperature_max_c`, `temperature_mean_c` | °C; mean is arithmetic mean of daily min/max. |
| `precipitation_mm_day` | mm day⁻¹. |
| `shortwave_radiation_w_m2` | W m⁻², daylight-period average in Daymet. |
| `vapor_pressure_pa` | Pa. |
| `snow_water_equivalent_kg_m2` | kg m⁻². |
| `daylength_s` | seconds. |
| `vpd_kpa` | kPa; derived as saturation vapor pressure at daily mean temperature minus Daymet vapor pressure, clipped at zero. |
| `relative_humidity_proxy_pct` | Derived daily proxy, not a directly measured RH series. |
| `frost_day` | `temperature_min_c < 0`. |
| `growing_degree_days_base5_c` | `max(Tmean − 5, 0)`. |
| `cooling_degree_days_base15_c` | `max(15 − Tmean, 0)`. |
| `temperature_mean_c_roll{3,7,14,21,30}_lag1` | Mean over the named trailing window ending one day before `date`. |
| `cooling_degree_days_base15_c_sum{3,7,14,21,30}_lag1` | Cooling accumulated over the named trailing window ending one day before `date`. |
| `precipitation_mm_sum{7,30}_lag1`, `vpd_kpa_roll{7,30}_lag1` | Trailing histories, also ending one day before `date`. |

“Several cooling scales” therefore means multiple trailing windows over the same temperature process, not different instruments or incompatible variables.

Daymet uses a 365-day calendar. In leap years its day 365 maps to December 30; no artificial December 31 value is inserted.

## `phenology_observations`

One measured variable at one original daily timestamp, platform, and ROI. `observation_id` is the primary key. Current types are `gcc_90` and `rcc_90`.

`quality_usable` is a non-destructive analysis flag. The row remains present when false. Snow codes 1 and 5 and the official GCC90 automated outlier flag are recorded in `quality_flag`. `species_scientific_name` and `individual_id` are null for camera ROIs; future NEON records will populate them.

## `phenology_events`

Derived operational events. Release 0.1 includes `rcc90_observed_autumn_maximum`: the date of the maximum usable observed RCC90 between September and November, without smoothing or interpolation. The adjacent observation dates form a sampling-resolution interval. This event is **not** asserted to be landscape peak visible color.

## `species_composition`

One species listed in camera-site metadata. `composition_measure = listed_dominant_species_presence` and value 1 means presence on the dominant-species list only. It is not tree-count fraction, basal-area fraction, or camera-visible canopy fraction.

## `imagery_metadata`

One ROI-day summary with a referenced midday filename and raw archive URL, image count, ROI, timestamps, and quality flags. Raw JPEGs are referenced but not bulk-downloaded.

## `coverage_audit`

One site × year, calculated from downloaded records. Camera completeness uses unique calendar dates so overlapping platforms do not inflate the fraction. `phenocam_observation_records` still records repeated evidence. Daymet completeness is calculated independently.

## `dataset_inventory`

One site × candidate dataset. `availability` distinguishes `downloaded` from verified-but-pending sources. Null completeness means “not calculated,” never zero and never inferred from a catalog range.

## `data_provenance`

Source URL, version/snapshot, download date, cache path, SHA-256 digest, processing steps, policy, and access status. Token-requiring or pending sources have no local digest.

