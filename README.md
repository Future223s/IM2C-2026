# IM2C 2026 — Peak Peeping: predicting peak fall foliage in the Appalachians

We predict the date of peak autumn foliage at 11 Appalachian locations, from Montpelier, VT to Asheville, NC, and forecast 2026.

**Model:** `peak = μ + δ`
- **μ:** the location's mean peak over the previous 3 years.
- **δ:** a weather-driven departure, a weighted sum of five pre-season weather anomalies:
  - spring warmth → earlier;
  - late-September chill → earlier;
  - August heat → later;
  - plus two smaller spring and early-summer terms.

**Target:** the MODIS satellite "mid-greendown" date, 2001–2025. This is a proxy: peak color as people report it is assumed to come about 1–3 weeks later, an offset that is not calibrated.

**Key results:**
- Location explains 45% of the variation in peak date, and its differences track climate: latitude plus elevation, or the date a site cools below 15 °C.
- On years it never saw, the model averages **2.68 days** of error, against 3.34 days for "predict the usual date" (20% better).
- The 80% prediction interval runs from 5.3 days before to 2.9 days after the forecast.
- **2026:** near normal to slightly early everywhere, from Sep 25 (Montpelier) to Oct 6 (Shenandoah).

## Try the interactive app (Peak Foliage Explorer)

You need only Python 3; nothing to install, and it works offline.

```bash
git clone -b panel-v1-dataset https://github.com/Future223s/IM2C-2026.git
cd IM2C-2026
python3 app/serve.py          # Windows: py app\serve.py
```

Your browser opens **http://localhost:8000**; if it doesn't, paste that address in. Press Ctrl+C in the terminal to stop the app.
- If port 8000 is busy: `python3 app/serve.py --port 8080`
- Open the page through the server, not by double-clicking `app/index.html`.

In the app:
1. Pick a location (dropdown or map) and a year (2001–2026).
2. See the rolling mean μ, the weather adjustment δ, and the prediction with its 80% interval.
3. Compare with the actual satellite date (2001–2025), see which weather factors drove δ, and view monthly weather vs normal.

Every past year shown is a held-out test: μ uses only earlier years, and the δ weights are refit without that year. Details: [`app/README.md`](app/README.md).

## Repository map

| Path | What it is |
|---|---|
| [`docs/PAPER_REVIEW.md`](docs/PAPER_REVIEW.md) | **Start here for the paper.** §0 is the full model summary; also a claim audit, paper outline and caveats |
| [`analysis/delta_model/`](analysis/delta_model/DELTA_MODEL.md) | Final μ + δ model, validation, uncertainty and 2026 forecast (figures D1–D5) |
| [`analysis/environmental/`](analysis/environmental/FINDINGS.md) | How location, elevation, latitude and climate affect peak timing (figures E01–E24) |
| [`tree_species_analysis/`](tree_species_analysis/README.md) | FIA tree-species composition vs environment |
| [`multi_region_phenology/`](multi_region_phenology/README.md) | PhenoCam staged curve models M1–M6 (independent dataset) |
| [`data/panel_v1/`](data/panel_v1/HOW_TO_USE_THIS_DATASET.md) | The 11-location dataset (2001–2026) plus usage guide and build log |
| `data/build_panel.py` | Rebuilds the dataset from Google Earth Engine |
| `docs/plan.md` | Team plan |
| `analysis/*.py`, `analysis/elevation/` | Earlier 4-region pilot analysis |

## Reproduce the analyses

Use the `fall-foliage` conda environment (Python 3.12 with pandas, numpy, statsmodels, scikit-learn and matplotlib). Rebuilding the dataset also needs `earthengine-api` and Earth Engine access.

```bash
cd analysis/environmental && for s in 01 02 03 04 05 06; do python ${s}_*.py; done   # ~2 min
cd ../delta_model && python delta_model.py                                         # ~1 min
python app/build_data.py                                                           # refresh app data
```

## Main caveats

- Dates are a **satellite proxy**. The offset to visually reported peak color is not calibrated, and the 2026 MODIS product isn't released yet.
- The 275 location-years carry only about 25 independent years of weather. Spatial claims rest on 11 sites.
- MODIS and VIIRS disagree on the trend (a gap of about +7 days per decade), so we make no trend claim. Re-estimate μ every year.
