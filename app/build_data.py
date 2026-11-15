"""Precompute everything the interactive app needs -> app/data.json (run once with the fall-foliage env).

For every location-year the prediction is made the way it would have been made in practice:
  mu    = mean MODIS peak of the 3 previous years (needs >= 1 prior year)
  delta = fixed 5-feature model (analysis/delta_model) refit WITHOUT that year (leave-one-year-out weights)
2026 uses weights fit on all 2001-2025 (the published forecast).
Interval = forecast + [q10, q90] of the out-of-sample backtest errors (analysis/delta_model/tables/D6).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
P = ROOT / "data" / "panel_v1"
D = ROOT / "analysis" / "delta_model" / "tables"

panel = pd.read_csv(P / "panel.csv")
locs = pd.read_csv(P / "locations.csv")
wx = pd.read_csv(P / "weather_daily.csv", parse_dates=["date"])
W = pd.read_csv(D / "D7_delta_weights.csv")
U = pd.read_csv(D / "D6_uncertainty.csv", index_col=0).iloc[:, 0]
FEATS = list(W.feature)
ORDER = list(locs.location)
train = panel[panel.is_forecast_row == 0]
ALPHA = 1.0


def fit(tr):
    xm = tr.groupby("location")[FEATS].mean()
    ym = tr.groupby("location").peak_doy.mean()
    X = tr[FEATS].values - xm.loc[tr.location].values
    sd = X.std(0)
    Z = X / sd
    y = tr.peak_doy.values - ym.loc[tr.location].values
    b = np.linalg.solve(Z.T @ Z + ALPHA * np.eye(len(FEATS)), Z.T @ y)
    return b / sd, xm


rows = []
for y in range(2001, 2027):
    tr = train[train.year != y]
    w, xm_tr = fit(tr)
    xm_all = train.groupby("location")[FEATS].mean()      # normals shown in the app
    for loc in ORDER:
        r = panel[(panel.location == loc) & (panel.year == y)].iloc[0]
        prior = train[(train.location == loc) & (train.year < y) & (train.year >= y - 3)]
        mu = prior.peak_doy.mean() if len(prior) else None
        anom = {f: float(r[f] - xm_tr.loc[loc, f]) for f in FEATS}
        contrib = {f: float(w[i] * anom[f]) for i, f in enumerate(FEATS)}
        delta = float(sum(contrib.values()))
        rows.append(dict(location=loc, year=y, mu=None if mu is None else float(mu), n_prior=len(prior),
                         delta=delta, contrib=contrib, weights={f: float(w[i]) for i, f in enumerate(FEATS)},
                         values={f: float(r[f]) for f in FEATS},
                         normals={f: float(xm_all.loc[loc, f]) for f in FEATS},
                         actual=None if pd.isna(r.peak_doy) else float(r.peak_doy),
                         viirs=None if pd.isna(r.peak_doy_viirs) else float(r.peak_doy_viirs),
                         helene=int(r.helene_2024)))

# monthly weather (Mar-Nov) per location-year + 2001-2025 climatology
wx["year"], wx["month"] = wx.date.dt.year, wx.date.dt.month
mon = wx.groupby(["location", "year", "month"]).agg(tmean=("tmean", "mean"), tmin=("tmin", "mean"),
                                                     tmax=("tmax", "mean"), pr=("pr", "sum"),
                                                     n=("tmean", "size")).reset_index()
clim = mon[mon.year <= 2025].groupby(["location", "month"])[["tmean", "tmin", "tmax", "pr"]].mean().reset_index()
monthly = {}
for (loc, y), g in mon.groupby(["location", "year"]):
    monthly.setdefault(loc, {})[str(y)] = [dict(m=int(r.month), tmean=round(r.tmean, 2), tmin=round(r.tmin, 2),
                                                tmax=round(r.tmax, 2), pr=round(r.pr, 1), n=int(r.n))
                                           for r in g.itertuples()]
climatology = {loc: [dict(m=int(r.month), tmean=round(r.tmean, 2), tmin=round(r.tmin, 2), tmax=round(r.tmax, 2),
                          pr=round(r.pr, 1)) for r in g.itertuples()] for loc, g in clim.groupby("location")}

out = dict(
    locations=[dict(name=r.location, state=r.state, lat=r.lat, lon=r.lon, ecoregion=r.ecoregion,
                    elev=round(r.elev_median_m), mean_peak=r.peak_doy_mean, sd=r.peak_doy_sd,
                    hardwood=r.frac_hardwood, nh=r.evt_northern_hardwood, oak=r.evt_oak_hickory)
               for r in locs.itertuples()],
    features=[dict(key=r.feature, label=r.label, unit=r.unit, sd=r.anomaly_SD, weight=r.weight_d_per_unit)
              for r in W.itertuples()],
    interval=dict(q10=U["q10"], q90=U["q90"], q05=U["q05"], q95=U["q95"], mae=U["MAE"],
                  p3=U["P_within_3d"], p5=U["P_within_5d"]),
    rows=rows, monthly=monthly, climatology=climatology,
    last_weather_date=str(wx.date.max().date()),
)
(Path(__file__).resolve().parent / "data.json").write_text(json.dumps(out, separators=(",", ":")))
ev = pd.DataFrame([r for r in rows if r["actual"] is not None and r["mu"] is not None])
e = ev.actual - ev.mu - ev.delta
print(f"wrote data.json: {len(rows)} location-years; app-consistent MAE (2002-2025) = {e.abs().mean():.2f} d, "
      f"mu-only {(ev.actual - ev.mu).abs().mean():.2f} d")
