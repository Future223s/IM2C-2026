# Multi-region autumn phenology dataset

This project builds an incremental relational dataset for asynchronous autumn phenology observations. Release 0.1 contains actual PhenoCam daily GCC90/RCC90 observations and actual Daymet daily environmental estimates through 2025. NEON, AmeriFlux, MODIS, Landsat, and Sentinel-2 are represented only in the verified inventory until records are downloaded and audited.

## Current release

- 12 ecological sites and 16 camera platforms
- 30 deciduous camera ROIs
- 146,692 long-format camera observations
- 113,880 site-days of Daymet environment
- 312 audited site-years
- 179 platform-year operational RCC90 maxima
- Parquet and CSV for every table

The schema preserves site, platform, ROI, timestamp, footprint, QC, and provenance. It does not require cross-modal same-day overlap.

## Reproduce

From this directory, using Python 3.11+:

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt
PYTHONPATH=src .venv/bin/python run_pipeline.py
PYTHONPATH=src .venv/bin/python run_analysis.py
PYTHONPATH=src .venv/bin/python run_shape_research.py
PYTHONPATH=src .venv/bin/python run_staged_models.py
PYTHONPATH=src .venv/bin/pytest -q
```

Use `--force` with `run_pipeline.py` only when a fresh source snapshot is wanted. Otherwise the pipeline uses its cache, including the existing verified PhenoCam catalog cache in the parent workspace when available.

Launch the dashboard with:

```bash
.venv/bin/streamlit run dashboard.py
```

## Deliverables

- [Site selection and readiness report](outputs/reports/site_selection_and_readiness.md)
- [Trajectory-shape research audit](outputs/reports/trajectory_shape_research.md)
- [Staged environmental model ladder](outputs/reports/staged_curve_models/README.md)
- [Dataset inventory](outputs/tables/dataset_inventory.csv)
- [Data dictionary](docs/data_dictionary.md)
- [Exploratory notebook](notebooks/exploratory_analysis.ipynb)
- [Coverage dashboard](dashboard.py)
- [Coverage figure](outputs/figures/phenocam_autumn_coverage.png)
- [RCC variation figure](outputs/figures/rcc90_maxima_variation.png)
- [Seasonal curves](outputs/figures/seasonal_rcc90_curves.png)
- [Functional trajectory modes](outputs/figures/shape_research/functional_pca_modes.png)
- [Curve-cluster diagnostics](outputs/figures/shape_research/cluster_diagnostics.png)
- [Build manifest](outputs/manifest.json)

## Source documentation

- PhenoCam network and archive: <https://phenocam.nau.edu/webcam/network/table/>
- PhenoCam access guidance and Fair Use reminder: <https://phenocam.nau.edu/education/PhenoCam_Access_Guide.pdf>
- Daymet single-pixel service: <https://daymet.ornl.gov/web_services>
- Daymet V4 documentation and method: <https://daac.ornl.gov/DAYMET/guides/Daymet_Daily_V4.html>
- NEON plant phenology DP1.10055.001: <https://data.neonscience.org/data-products/DP1.10055.001/RELEASE-2026>
- NEON data-query API: <https://data.neonscience.org/data-api/endpoints/data-query/>
- MODIS MCD12Q2 V6.1 user guide: <https://lpdaac.usgs.gov/documents/1417/MCD12Q2_User_Guide_V61.pdf>
- USA-NPN data services: <https://data.usanpn.org/>
- AmeriFlux data policies: <https://ameriflux.lbl.gov/resource/ameriflux-data-use-policies/>
- Harvard Forest data archive: <https://harvardforest.fas.harvard.edu/data-archives/data-archive/>

## NEON next step

`scripts/download_neon_phenology.py` provides the token-authenticated official retrieval path. It preserves original files; a subsequent normalizer should map `phe_statusintensity`, `phe_perindividual`, and taxonomy records into the long observation table without averaging individuals.

```bash
export NEON_API_TOKEN='...'
.venv/bin/python scripts/download_neon_phenology.py --sites HARV,BART,GRSM,ORNL
```

Do not commit the token. The release deliberately does not claim NEON records were constructed without this step.

