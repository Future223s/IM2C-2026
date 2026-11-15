"""Task 2: leave-one-YEAR-out baselines (climatology and linear trend)."""
import numpy as np
import pandas as pd

from common import REGIONS, TARGET, load, ROOT

df = load()
years = sorted(df.year.unique())
preds = []
for y in years:
    train, test = df[df.year != y], df[df.year == y]
    for _, row in test.iterrows():
        tr = train[train.region == row.region]
        clim = tr[TARGET].mean()
        # region intercept + shared slope, fit on training years only
        X = pd.get_dummies(train.region).astype(float)
        X["year_c"] = train.year - 2013
        beta, *_ = np.linalg.lstsq(X.values, train[TARGET].values, rcond=None)
        b = dict(zip(X.columns, beta))
        trend_shared = b[row.region] + b["year_c"] * (row.year - 2013)
        # per-region slope
        p = np.polyfit(tr.year, tr[TARGET], 1)
        trend_region = np.polyval(p, row.year)
        preds.append({"region": row.region, "year": y, "obs": row[TARGET],
                      "helene_2024": row.helene_2024, "climatology": clim,
                      "trend_shared_slope": trend_shared,
                      "trend_per_region_slope": trend_region})
P = pd.DataFrame(preds)
models = ["climatology", "trend_shared_slope", "trend_per_region_slope"]


def score(sub, label):
    out = []
    for m in models:
        e = sub[m] - sub.obs
        out.append({"subset": label, "model": m, "n": len(sub),
                    "MAE": e.abs().mean(), "RMSE": np.sqrt((e ** 2).mean()),
                    "bias": e.mean()})
    return out


rows = []
for r in REGIONS:
    rows += score(P[P.region == r], r)
rows += score(P, "Overall")
rows += score(P[P.helene_2024 == 0], "Overall, no Helene")
rows += score(P[P.year >= 2016], "Overall, 2016-2025 only")
res = pd.DataFrame(rows)
print(res.round(2).to_string(index=False))
res.round(3).to_csv(ROOT / "baseline_scores.csv", index=False)
P.round(3).to_csv(ROOT / "baseline_predictions.csv", index=False)
print("saved analysis/baseline_scores.csv, analysis/baseline_predictions.csv")
