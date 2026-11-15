# Peak Foliage Explorer: interactive demo of our model

An interactive page to explore our peak-foliage model at the 11 Appalachian study locations.

## Run it (any computer with Python 3; no installs needed)

```bash
python3 app/serve.py          # opens http://localhost:8000 in your browser
python3 app/serve.py --port 8080   # if port 8000 is busy
```

Windows: `py app\serve.py`. Stop the server with Ctrl+C.

## What you can do

- **Pick a location** from the dropdown or by clicking the map (north → south, Vermont to North Carolina).
- **Pick a year** (2001–2026) with the slider or the quick buttons.
- **See the prediction:** the rolling mean μ (average of the 3 previous years), the weather adjustment δ, the forecast μ + δ, and its 80% prediction interval.
- **Compare with what actually happened:** the observed MODIS peak date (and VIIRS from 2013), the error, whether the result fell inside the 80% interval, and whether δ beat the rolling mean alone.
- **Every year at a glance:** a timeline of predicted vs. actual for the selected location. Click any year to jump to it.
- **Why δ is what it is:** each weather factor's value that year, its normal, the anomaly, its weight, and its contribution in days.
- **That season's weather:** monthly mean temperature and precipitation (Mar–Nov) against the 2001–2025 normal.

## How the numbers are made (honest hindcasts)

- **Past years:** μ uses only the 3 years *before* the selected year, and the δ weights are refit *without* the selected year. Every past prediction is therefore genuinely out-of-sample. Across 2002–2025 this gives MAE 3.0 d, vs 3.7 d for the rolling mean alone. The paper's headline numbers come from `analysis/delta_model` (nested leave-one-year-out 2.68 d; forward backtest 2.83 d).
- **2026:** uses the published model (weights fit on 2001–2025; μ = 2023–2025 mean). This is the forecast in `analysis/delta_model/DELTA_MODEL.md`.
- **Interval:** the 80% range of 165 out-of-sample backtest errors (−5.3 to +2.9 d).

## Caveats

- **Dates are the MODIS mid-greendown satellite proxy,** the day forest greenness is halfway through its autumn decline over a 40 km region. People usually report peak color about 1–3 weeks later; that offset is an assumption and has not been calibrated.
- **2026 has no observed date yet,** because the 2026 MODIS product has not been released.

## Files

- `index.html`: the app (plain HTML/JS/SVG, no external libraries; works offline).
- `data.json`: precomputed predictions and weather.
- `build_data.py`: rebuilds `data.json` from `data/panel_v1` and `analysis/delta_model/tables`. Run it with the fall-foliage env after the model changes.
- `serve.py`: a tiny local web server using only the Python standard library.
