"""Delta model: peak(l, y) = mu(l, y) + delta(l, y), with calibrated uncertainty and a 2026 forecast.

mu    = trailing k-year mean of the location's own peak (k chosen by an operational backtest)
delta = sum_k beta_k * (x_k(l, y) - mean_l x_k)  shared slopes (no location interactions; see
        analysis/environmental FINDINGS #5 and multi_region_phenology M4), forecast-safe features only.

Validation:
  A. nested leave-one-year-out (feature selection repeated inside every fold) -> honest delta skill
  B. operational backtest 2011-2025: everything (mu, feature selection, weights) uses only years < y
  C. residual quantiles from B -> prediction intervals; coverage checked
Also tests FIA species-composition interactions (tree_species_analysis) and compares with PhenoCam M2.
Run with the fall-foliage env from this folder.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "environmental"))
from common import (INK, INK2, MUTED, GRID, BLUE, ORANGE, AQUA, RED, SURFACE, plt, doy_label)  # noqa

ROOT = Path(__file__).resolve().parent
PANEL = ROOT.parents[1] / "data" / "panel_v1"
FIG, TAB = ROOT / "figures", ROOT / "tables"
FIG.mkdir(exist_ok=True); TAB.mkdir(exist_ok=True)

panel = pd.read_csv(PANEL / "panel.csv")
locs = pd.read_csv(PANEL / "locations.csv")
timing = pd.read_csv(PANEL / "feature_timing.csv")
ORDER = list(locs.location)
train = panel[panel.is_forecast_row == 0].copy()
f26 = panel[panel.year == 2026].copy()

# ---------------------------------------------------------------- candidate features (forecast-safe + available 2026)
safe = timing.groupby("feature").timing.apply(lambda s: (s == "pre-peak").all())
CANDS = [f for f in safe[safe].index
         if f in panel.columns and f26[f].notna().all() and f != "first_frost_doy"]
LABEL = {"tmean_MarApr": "Spring temp, Mar–Apr", "gdd5_MarMay": "Spring warmth GDD5, Mar–May",
         "tmax_JunJul": "Early-summer heat, Jun–Jul", "precip_JunJul": "Early-summer rain, Jun–Jul",
         "tmax_Aug": "August max temp", "tmean_pre30": "Mean temp, 30 d before typical peak",
         "tmin_pre30": "Min temp, 30 d before typical peak", "cdd_pre30": "Cold degree-days, 30 d before",
         "cold_nights_pre30": "Cold nights (<5 °C), 30 d before", "tmax_pre60to31": "Late-summer heat, 60–31 d before",
         "precip_pre60": "Rain, 60 d before typical peak", "vpd_pre60": "Dryness (VPD), 60 d before",
         "srad_pre30": "Sunshine, 30 d before typical peak"}
UNIT = {"tmean_MarApr": "°C", "gdd5_MarMay": "°C·d", "tmax_JunJul": "°C", "precip_JunJul": "mm", "tmax_Aug": "°C",
        "tmean_pre30": "°C", "tmin_pre30": "°C", "cdd_pre30": "°C·d", "cold_nights_pre30": "nights",
        "tmax_pre60to31": "°C", "precip_pre60": "mm", "vpd_pre60": "kPa", "srad_pre30": "W/m²"}
print("candidates:", CANDS)
ALPHA = 1.0


def fit_delta(tr, feats):
    """Ridge on standardized within-location anomalies (training rows only). Returns predictor function + coefs."""
    xm = tr.groupby("location")[feats].mean()
    ym = tr.groupby("location").peak_doy.mean()
    X = tr[feats].values - xm.loc[tr.location].values
    sd = X.std(0)
    Z = X / sd
    y = tr.peak_doy.values - ym.loc[tr.location].values
    b = np.linalg.solve(Z.T @ Z + ALPHA * np.eye(len(feats)), Z.T @ y)

    def predict(part):
        Xp = part[feats].values - xm.loc[part.location].values
        return (Xp / sd) @ b
    return predict, b / sd, b          # raw-unit weights (d per unit), standardized weights (d per SD)


def loyo_mae(df, feats):
    errs = []
    for y in df.year.unique():
        tr, te = df[df.year != y], df[df.year == y]
        mu = tr.groupby("location").peak_doy.mean().loc[te.location].values
        d = fit_delta(tr, feats)[0](te) if feats else 0
        errs.append(te.peak_doy.values - mu - d)
    return np.mean(np.abs(np.concatenate(errs)))


def forward_select(df, cands, max_k=5, min_gain=0.02):
    chosen, best = [], loyo_mae(df, [])
    path = [("(climatology)", best)]
    while len(chosen) < max_k:
        trial = {f: loyo_mae(df, chosen + [f]) for f in cands if f not in chosen}
        f, m = min(trial.items(), key=lambda kv: kv[1])
        if best - m < min_gain:
            break
        chosen.append(f); best = m; path.append((f, m))
    return chosen, path


# ---------------------------------------------------------------- 1. selection on all years (for the final weights)
FEATS, PATH = forward_select(train, CANDS)
print("selected:", FEATS)
pd.DataFrame(PATH, columns=["added", "LOYO_MAE"]).to_csv(TAB / "D1_forward_selection_path.csv", index=False)

# ---------------------------------------------------------------- A. nested LOYO (honest)
rows = []
for y in sorted(train.year.unique()):
    tr, te = train[train.year != y], train[train.year == y]
    feats_y, _ = forward_select(tr, CANDS)
    mu = tr.groupby("location").peak_doy.mean().loc[te.location].values
    d = fit_delta(tr, feats_y)[0](te) if feats_y else 0
    rows.append(te[["location", "year", "peak_doy"]].assign(mu=mu, delta=d, pred=mu + d,
                                                             feats="+".join(feats_y)))
NL = pd.concat(rows, ignore_index=True)
NL.to_csv(TAB / "D2_nested_loyo_predictions.csv", index=False)
nl_mae_model = (NL.peak_doy - NL.pred).abs().mean()
nl_mae_clim = (NL.peak_doy - NL.mu).abs().mean()
print(f"nested LOYO: model {nl_mae_model:.3f} vs climatology {nl_mae_clim:.3f}")
sel_freq = pd.Series("+".join(NL.drop_duplicates("year").feats).split("+")).value_counts()
sel_freq.to_csv(TAB / "D3_selection_frequency.csv")

# ---------------------------------------------------------------- B. operational backtest (only the past is used)
KS = {"3 yr": 3, "5 yr": 5, "10 yr": 10, "All prior": 99}
rows = []
for y in range(2011, 2026):
    past, te = train[train.year < y], train[train.year == y]
    feats_y, _ = forward_select(past, CANDS)
    d = fit_delta(past, feats_y)[0](te) if feats_y else np.zeros(len(te))
    for name, k in KS.items():
        mu = past[past.year >= y - k].groupby("location").peak_doy.mean().loc[te.location].values
        rows.append(te[["location", "year", "peak_doy"]].assign(window=name, mu=mu, delta=d, pred=mu + d))
BT = pd.concat(rows, ignore_index=True)
BT.to_csv(TAB / "D4_operational_backtest.csv", index=False)
W = (BT.assign(e_clim=(BT.peak_doy - BT.mu).abs(), e_model=(BT.peak_doy - BT.pred).abs())
       .groupby("window")[["e_clim", "e_model"]].mean().loc[list(KS)])
W.to_csv(TAB / "D5_mu_window_backtest.csv")
print(W.round(3))
BEST_K = W.e_model.idxmin()
KBEST = KS[BEST_K]
print("best rolling window:", BEST_K)

# ---------------------------------------------------------------- C. uncertainty from operational residuals
res = BT[BT.window == BEST_K].assign(err=lambda d: d.peak_doy - d.pred)
q = res.err.quantile([0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95])
absq = res.err.abs().quantile([0.5, 0.8, 0.9])
p_within = {k: (res.err.abs() <= k).mean() for k in (2, 3, 5, 7)}
# coverage check of the 80% interval with a split: quantiles from 2011-2018, coverage on 2019-2025
early, late = res[res.year <= 2018], res[res.year > 2018]
lo80, hi80 = early.err.quantile([0.10, 0.90])
cov_late = ((late.err >= lo80) & (late.err <= hi80)).mean()
UNC = pd.Series({"n_backtest": len(res), "MAE": res.err.abs().mean(), "RMSE": np.sqrt((res.err ** 2).mean()),
                 "bias": res.err.mean(), "half_width_50": absq[0.5], "half_width_80": absq[0.8],
                 "half_width_90": absq[0.9], **{f"P_within_{k}d": v for k, v in p_within.items()},
                 "q10": q[0.10], "q90": q[0.90], "q05": q[0.05], "q95": q[0.95],
                 "coverage_80_on_2019_2025_using_2011_2018_quantiles": cov_late})
UNC.to_csv(TAB / "D6_uncertainty.csv")
print(UNC.round(3))

# ---------------------------------------------------------------- final weights on all 2001-2025
predict, w_raw, w_std = fit_delta(train, FEATS)
SDX = (train[FEATS] - train.groupby("location")[FEATS].transform("mean")).std()
WT = pd.DataFrame({"feature": FEATS, "label": [LABEL.get(f, f) for f in FEATS],
                   "unit": [UNIT.get(f, "") for f in FEATS], "weight_d_per_unit": w_raw,
                   "weight_d_per_SD": w_std, "anomaly_SD": SDX.values})
# bootstrap over years for weight CIs
rng = np.random.default_rng(3)
yrs = train.year.unique()
boots = []
for _ in range(1000):
    pick = rng.choice(yrs, len(yrs), replace=True)
    bt = pd.concat([train[train.year == y].assign(year=1000 + i) for i, y in enumerate(pick)])
    boots.append(fit_delta(bt, FEATS)[1])
boots = np.array(boots)
WT["lo95"], WT["hi95"] = np.percentile(boots, 2.5, 0), np.percentile(boots, 97.5, 0)
WT["share_of_delta_variance"] = (w_std ** 2) / (w_std ** 2).sum()
WT.to_csv(TAB / "D7_delta_weights.csv", index=False)
print(WT.round(4).to_string())

# ---------------------------------------------------------------- species (FIA) interaction check
fia_path = ROOT.parents[1] / "tree_species_analysis" / "data" / "processed" / "appalachian_model_distributions.csv"
SPEC = None
if fia_path.exists():
    fia = pd.read_csv(fia_path).set_index("Category").T
    fia.index = [c.rsplit(" ", 1)[0] for c in fia.index]          # "Montpelier VT" -> "Montpelier"
    fia = fia.loc[ORDER]
    rows = []
    base = loyo_mae(train, FEATS)
    top = FEATS[0]
    for cat in ["Evergreen", "Oak", "Hard maple", "Soft maple", "Beech", "Yellow-poplar"]:
        s = (fia[cat] - fia[cat].mean()) / fia[cat].std()
        t = train.assign(**{f"ix_{top}": (train[top] - train.groupby("location")[top].transform("mean"))
                            * train.location.map(s)})
        rows.append(dict(interaction=f"{top} × FIA {cat} %", LOYO_MAE=loyo_mae(t, FEATS + [f"ix_{top}"]),
                         base_MAE=base))
    SPEC = pd.DataFrame(rows).assign(change=lambda d: d.LOYO_MAE - d.base_MAE)
    SPEC.to_csv(TAB / "D8_species_interactions.csv", index=False)
    print(SPEC.round(3).to_string())
    # composition vs mu (cross-sectional, n = 11)
    mu_all = train.groupby("location").peak_doy.mean().loc[ORDER]
    cs = pd.DataFrame({c: [np.corrcoef(fia[c], mu_all)[0, 1]] for c in fia.columns}).T.rename(columns={0: "r_with_mu"})
    cs.sort_values("r_with_mu").to_csv(TAB / "D9_species_vs_mu.csv")

# ---------------------------------------------------------------- 2026 forecast
mu26 = train[train.year >= 2026 - KBEST].groupby("location").peak_doy.mean().loc[f26.location].values
d26 = predict(f26)
F = f26[["location"] + FEATS].copy()
F["mu_rolling"] = mu26
F["delta"] = d26
F["pred_doy"] = mu26 + d26
F["pred_date"] = [doy_label(v, 2026) for v in F.pred_doy]
F["lo80"] = [doy_label(v + q[0.10], 2026) for v in F.pred_doy]
F["hi80"] = [doy_label(v + q[0.90], 2026) for v in F.pred_doy]
F["mu_date"] = [doy_label(v, 2026) for v in mu26]
for f in FEATS:   # contribution of each feature to delta (days)
    F[f"contrib_{f}"] = (f26[f].values - train.groupby("location")[f].mean().loc[f26.location].values) * \
        WT.set_index("feature").loc[f, "weight_d_per_unit"]
F.to_csv(TAB / "D10_forecast_2026.csv", index=False)
print(F[["location", "mu_date", "delta", "pred_date", "lo80", "hi80"]].round(2).to_string())

# ---------------------------------------------------------------- figures
# D1 weights
fig, ax = plt.subplots(figsize=(9, 0.9 + 0.65 * len(FEATS)))
yy = np.arange(len(WT))[::-1]
ax.errorbar(WT.weight_d_per_SD, yy, xerr=[WT.weight_d_per_SD - WT.lo95 * WT.anomaly_SD,
                                          WT.hi95 * WT.anomaly_SD - WT.weight_d_per_SD],
            fmt="o", color=BLUE, capsize=3, ms=8)
for yi, r in zip(yy, WT.itertuples()):
    ax.annotate(f"{r.weight_d_per_unit:+.3g} d per {r.unit}  (1 SD = {r.anomaly_SD:.3g} {r.unit})",
                (r.weight_d_per_SD, yi), xytext=(0, 10), textcoords="offset points", ha="center", fontsize=8,
                color=INK)
ax.axvline(0, color=INK2, lw=1)
ax.set_yticks(yy, WT.label)
ax.set_ylim(-0.6, len(WT) - 0.3)
ax.set_xlabel("Days of peak shift for a 1-SD anomaly (95% CI, year bootstrap)")
ax.grid(axis="y", visible=False)
ax.set_title("The δ function: weights of each forecast-safe weather anomaly")
fig.savefig(FIG / "D1_delta_weights.png"); plt.close(fig)

# D2 rolling-window backtest
fig, ax = plt.subplots(figsize=(8, 3.8))
x = np.arange(len(W))
ax.bar(x - 0.18, W.e_clim, 0.34, color="#d3d2cc", label="μ only (rolling mean)")
ax.bar(x + 0.18, W.e_model, 0.34, color=BLUE, label="μ + δ")
for xi, (a, b) in enumerate(zip(W.e_clim, W.e_model)):
    ax.text(xi - 0.18, a + 0.04, f"{a:.2f}", ha="center", fontsize=8, color=INK)
    ax.text(xi + 0.18, b + 0.04, f"{b:.2f}", ha="center", fontsize=8, color=INK)
ax.set_xticks(x, W.index)
ax.set_ylabel("MAE, 2011–2025 (days)")
ax.set_xlabel("Rolling-mean window for μ")
ax.legend(fontsize=8.5, loc="upper right")
ax.grid(axis="x", visible=False)
ax.set_title("Operational backtest: train on past years only, predict the next year")
fig.savefig(FIG / "D2_rolling_window_backtest.png"); plt.close(fig)

# D3 predicted vs observed (operational backtest)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw=dict(width_ratios=[1, 1.1]))
a1.plot([258, 296], [258, 296], color=MUTED, ls="--", lw=1)
a1.scatter(res.pred, res.peak_doy, s=16, color=BLUE, alpha=0.6, lw=0)
tk = [262, 270, 278, 286, 294]
a1.set_xticks(tk, [doy_label(v) for v in tk]); a1.set_yticks(tk, [doy_label(v) for v in tk])
a1.set_xlabel("Predicted (μ + δ)"); a1.set_ylabel("Observed MODIS peak")
a1.set_title(f"Backtest 2011–2025, n = {len(res)}  (MAE {res.err.abs().mean():.2f} d)", fontsize=10)
bins = np.arange(-14, 15, 1)
a2.hist(res.err, bins=bins, color=BLUE, alpha=0.8, edgecolor=SURFACE)
for qq, lab in [(q[0.10], "10%"), (q[0.90], "90%")]:
    a2.axvline(qq, color=ORANGE, lw=1.6, ls="--")
    a2.text(qq, a2.get_ylim()[1] * 0.92, f" {lab}: {qq:+.1f} d", color=INK2, fontsize=8)
a2.set_xlabel("Error: observed − predicted (days)")
a2.set_ylabel("Location-years")
a2.set_title(f"80% of errors fall in [{q[0.10]:+.1f}, {q[0.90]:+.1f}] d;  P(|error| ≤ 3 d) = "
             f"{p_within[3]:.0%}", fontsize=10)
fig.suptitle("How accurate is μ + δ when forecasting a year it has never seen?", x=0.01, ha="left",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "D3_backtest_accuracy.png"); plt.close(fig)

# D4 2026 forecast
fig, ax = plt.subplots(figsize=(10, 5.2))
yy = np.arange(len(F))[::-1]
ax.errorbar(F.pred_doy, yy, xerr=[[-q[0.10]] * len(F), [q[0.90]] * len(F)], fmt="o", color=BLUE, ms=8,
            capsize=3, lw=2, label="μ + δ forecast, 80% interval")
ax.scatter(F.mu_rolling, yy, marker="|", s=220, color=INK, lw=2, label=f"Rolling mean μ ({BEST_K})", zorder=4)
for yi, r in zip(yy, F.itertuples()):
    ax.text(297.5, yi, f"{r.pred_date}  (δ {r.delta:+.1f} d)", va="center", fontsize=8.5, color=INK)
ax.set_yticks(yy, F.location)
tk = [266, 270, 274, 278, 282, 286, 290]
ax.set_xticks(tk, [doy_label(v, 2026) for v in tk])
ax.set_xlim(263, 307)
ax.set_xlabel("2026 date (MODIS mid-greendown proxy; reported peak color is typically 1–3 weeks later)")
ax.legend(loc="upper center", bbox_to_anchor=(0.4, -0.13), ncol=2, fontsize=8.5)
ax.grid(axis="y", visible=False)
ax.set_title("2026 forecast: rolling mean + weather-driven δ")
fig.savefig(FIG / "D4_forecast_2026.png"); plt.close(fig)

# D5 contribution breakdown 2026
C = F.set_index("location")[[f"contrib_{f}" for f in FEATS]].loc[ORDER]
fig, ax = plt.subplots(figsize=(10, 5.2))
cols = [BLUE, ORANGE, AQUA, "#eda100", "#e87ba4"]
left_pos = np.zeros(len(C)); left_neg = np.zeros(len(C))
yy = np.arange(len(C))[::-1]
for k, f in enumerate(FEATS):
    v = C[f"contrib_{f}"].values
    lp = np.where(v >= 0, left_pos, left_neg)
    ax.barh(yy, v, left=lp, color=cols[k % len(cols)], height=0.6, label=LABEL.get(f, f),
            edgecolor=SURFACE, linewidth=1.5)
    left_pos += np.where(v >= 0, v, 0); left_neg += np.where(v < 0, v, 0)
ax.scatter(C.sum(axis=1), yy, color=INK, zorder=4, s=30, label="δ (total)")
ax.axvline(0, color=INK2, lw=1)
ax.set_yticks(yy, C.index)
ax.set_xlabel("Contribution to 2026 δ (days; + = later than normal)")
ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3)
ax.grid(axis="y", visible=False)
ax.set_title("Why 2026 is early or late at each location: δ broken down by weather driver")
fig.savefig(FIG / "D5_delta_breakdown_2026.png"); plt.close(fig)
print("figures written")
