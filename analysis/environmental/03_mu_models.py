"""Step 3: Multivariable models of the location mean (mu), n = 11.

(a) geography model mu ~ lat + elev (+ lon): partial effects, Hopkins-style equivalence
(b) mediation: do lat / elev act through temperature?  mu ~ lat + elev + thermal
(c) every 1- and 2-predictor model, scored by LOLO; honest error via NESTED LOLO (model chosen inside each fold)
(d) influence diagnostics (Cook's D, leave-Shenandoah-out), WLS robustness
"""
import itertools

import numpy as np
import pandas as pd
import statsmodels.api as sm

from common import *  # noqa

panel, locs, train, ORDER = load()
P = pd.read_csv(TAB / "E04_location_predictors.csv").set_index("location").loc[ORDER]
y = P["peak_doy_mean"].values
n = len(y)
sd_mu = train.groupby("location")["peak_doy"].std().loc[ORDER].values / 5   # SE of each mu (sqrt 25)


def ols(cols, data=P, w=None):
    X = sm.add_constant(data[cols].astype(float))
    return (sm.WLS(y, X, weights=w) if w is not None else sm.OLS(y, X)).fit()


def lolo(cols):
    err = []
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        X = sm.add_constant(P[cols].astype(float).values, has_constant="add")
        b = np.linalg.lstsq(X[m], y[m], rcond=None)[0]
        err.append(y[i] - X[i] @ b)
    return np.array(err)


# ---------------------------------------------------------------- (a)+(b) nested model ladder
ladder = {
    "Latitude": ["lat"],
    "Latitude + elevation": ["lat", "elev_median_m"],
    "Lat + elev + longitude": ["lat", "elev_median_m", "lon"],
    "Thermal date (15 °C)": ["thermal_date_15C"],
    "Lat + elev + thermal date": ["lat", "elev_median_m", "thermal_date_15C"],
    "Sep temp normal": ["norm_tmean_Sep"],
}
rows, coefs = [], []
for name, cols in ladder.items():
    f = ols(cols)
    e = lolo(cols)
    k = len(cols) + 1
    aicc = f.aic + 2 * k * (k + 1) / (n - k - 1)
    rows.append(dict(model=name, k=len(cols), R2=f.rsquared, adjR2=f.rsquared_adj, AICc=aicc,
                     LOLO_MAE=np.abs(e).mean(), LOLO_RMSE=np.sqrt((e ** 2).mean())))
    for c in cols:
        coefs.append(dict(model=name, term=c, coef=f.params[c], se=f.bse[c], p=f.pvalues[c],
                          lo=f.conf_int().loc[c, 0], hi=f.conf_int().loc[c, 1]))
L = pd.DataFrame(rows)
CO = pd.DataFrame(coefs)
L.to_csv(TAB / "E08_mu_model_ladder.csv", index=False)
CO.to_csv(TAB / "E09_mu_model_coefficients.csv", index=False)
print(L.round(3).to_string())
print(CO.round(4).to_string())

geo = ols(["lat", "elev_median_m"])
b_lat, b_el = geo.params["lat"], geo.params["elev_median_m"]
equiv_m_per_deg = b_lat / b_el if b_el != 0 else np.nan     # metres of elevation equal to 1° latitude
print(f"\nGeography model: {b_lat:+.2f} d/°lat, {b_el * 100:+.2f} d/100 m -> 1° lat ≈ {abs(equiv_m_per_deg):.0f} m")

# ---------------------------------------------------------------- Fig E09: coefficients across the ladder
show = [("lat", "Latitude (d per °N)", 1), ("elev_median_m", "Elevation (d per 100 m)", 100),
        ("thermal_date_15C", "Thermal date (d per d)", 1)]
fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
models = list(ladder)
for ax, (term, lab, sc) in zip(axes, show):
    sub = CO[CO.term == term].set_index("model").reindex(models).dropna()
    yy = np.arange(len(sub))[::-1]
    ax.errorbar(sub.coef * sc, yy, xerr=[(sub.coef - sub.lo) * sc, (sub.hi - sub.coef) * sc], fmt="o",
                color=BLUE, ecolor=BLUE, capsize=3, ms=7)
    for yi, (c, p) in zip(yy, zip(sub.coef * sc, sub.p)):
        ax.annotate(f"{c:+.2f} (p={p:.3f})", (c, yi), xytext=(0, 8), textcoords="offset points",
                    ha="center", fontsize=7.8, color=INK)
    ax.axvline(0, color=INK2, lw=1)
    ax.set_yticks(yy, sub.index)
    ax.set_ylim(-0.6, len(sub) - 0.3)
    ax.set_title(lab, fontsize=10)
    ax.grid(axis="y", visible=False)
fig.suptitle("How each effect changes as other variables enter the model (95% CI, n = 11)", x=0.01, ha="left",
             fontsize=12, fontweight="bold")
fig.tight_layout()
save(fig, "E09_coefficient_ladder",
     "Coefficients of latitude, elevation and the 15 °C thermal date in nested models of the location mean. "
     "Elevation only matters once latitude is held fixed; both lose their effect once the thermal date enters "
     "(mediation through temperature).")

# ---------------------------------------------------------------- (b) partial-regression (added-variable) plots
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
for ax, (term, other, lab, sc) in zip(axes, [("elev_median_m", "lat", "Elevation (m), latitude removed", 1),
                                             ("lat", "elev_median_m", "Latitude (°), elevation removed", 1)]):
    ry = sm.OLS(y, sm.add_constant(P[[other]])).fit().resid
    rx = sm.OLS(P[term], sm.add_constant(P[[other]])).fit().resid
    b = np.polyfit(rx, ry, 1)
    xs = np.linspace(rx.min(), rx.max(), 20)
    ax.plot(xs, np.polyval(b, xs), color=BLUE)
    ax.scatter(rx, ry, s=55, color=INK, edgecolor="white", lw=1.5, zorder=3)
    off = {"elev_median_m": {"Boone": (-4, -11, "right"), "Frostburg": (-4, 3, "right")},
           "lat": {"Frostburg": (4, -11, "left"), "Binghamton": (-4, 3, "right")}}[term]
    for xi, yi, nm in zip(rx, ry, ORDER):
        dx_, dy_, ha_ = off.get(nm, (4, 3, "left"))
        ax.annotate(nm, (xi, yi), xytext=(dx_, dy_), textcoords="offset points", fontsize=7.5, color=INK2,
                    ha=ha_)
    ax.axhline(0, color=GRID)
    ax.set_xlabel(lab)
    ax.set_ylabel("Mean peak, other variable removed (d)")
    unit = "100 m" if term == "elev_median_m" else "°"
    mult = 100 if term == "elev_median_m" else 1
    ax.set_title(f"Partial effect: {b[0] * mult:+.2f} d per {unit}", fontsize=10)
fig.suptitle("Added-variable plots for the geography model  μ ~ latitude + elevation", x=0.01, ha="left",
             fontsize=12, fontweight="bold")
fig.tight_layout()
save(fig, "E10_partial_regression",
     f"Added-variable plots. Holding latitude fixed, higher sites peak later? The partial slopes are "
     f"{b_lat:+.2f} d per degree latitude and {b_el * 100:+.2f} d per 100 m of forest elevation.")

# ---------------------------------------------------------------- (c) all 1- and 2-predictor models
cands = ["lat", "lon", "elev_median_m", "slope_mean_deg", "norm_tmean_Sep", "norm_tmin_Oct",
         "norm_precip_AugSep", "norm_gdd5_MarMay", "frac_hardwood", "frac_evergreen", "frac_developed",
         "evt_northern_hardwood", "evt_oak_hickory", "evt_cove_mesic", "thermal_date_12C", "thermal_date_15C"]
models = [[c] for c in cands] + [list(p) for p in itertools.combinations(cands, 2)]
res = []
for cols in models:
    e = lolo(cols)
    f = ols(cols)
    k = len(cols) + 1
    res.append(dict(model=" + ".join(cols), k=len(cols), LOLO_MAE=np.abs(e).mean(), R2=f.rsquared,
                    AICc=f.aic + 2 * k * (k + 1) / (n - k - 1), max_vif=(
                        1 / (1 - np.corrcoef(P[cols[0]], P[cols[1]])[0, 1] ** 2) if len(cols) == 2 else 1)))
AM = pd.DataFrame(res).sort_values("LOLO_MAE").reset_index(drop=True)
AM.to_csv(TAB / "E10_all_subsets_mu.csv", index=False)
print("\nTop 12 models by LOLO MAE:\n", AM.head(12).round(3).to_string())


# honest error: nested LOLO (choose the best model on the 10 training sites, by inner LOLO)
def nested(model_set):
    errs, chosen = [], []
    for i in range(n):
        keep = np.array([j for j in range(n) if j != i])
        best, best_mae = None, np.inf
        for cols in model_set:
            X = sm.add_constant(P[cols].astype(float).values, has_constant="add")
            inner = []
            for j in keep:
                m = np.isin(np.arange(n), keep) & (np.arange(n) != j)
                b = np.linalg.lstsq(X[m], y[m], rcond=None)[0]
                inner.append(abs(y[j] - X[j] @ b))
            if np.mean(inner) < best_mae:
                best, best_mae = cols, np.mean(inner)
        X = sm.add_constant(P[best].astype(float).values, has_constant="add")
        b = np.linalg.lstsq(X[keep], y[keep], rcond=None)[0]
        errs.append(y[i] - X[i] @ b)
        chosen.append(" + ".join(best))
    return np.array(errs), chosen


grand = np.array([y[i] - np.delete(y, i).mean() for i in range(n)])
e1, ch1 = nested([[c] for c in cands])
e2, ch2 = nested(models)
fixed_geo = lolo(["lat", "elev_median_m"])
fixed_th = lolo(["thermal_date_15C"])
H = pd.DataFrame([
    dict(approach="Grand mean (no predictor)", MAE=np.abs(grand).mean()),
    dict(approach="Latitude + elevation (fixed a priori)", MAE=np.abs(fixed_geo).mean()),
    dict(approach="Thermal date 15 °C (fixed a priori)", MAE=np.abs(fixed_th).mean()),
    dict(approach="Best single predictor, chosen in each fold", MAE=np.abs(e1).mean()),
    dict(approach="Best 1- or 2-predictor model, chosen in each fold", MAE=np.abs(e2).mean()),
    dict(approach="Best model in-sample LOLO (optimistic)", MAE=AM.LOLO_MAE.iloc[0]),
])
H.to_csv(TAB / "E11_honest_mu_error.csv", index=False)
print("\nHonest errors:\n", H.round(3).to_string())
print("models chosen in nested folds (2-pred search):", pd.Series(ch2).value_counts().to_dict())

fig, ax = plt.subplots(figsize=(9, 3.6))
yy = np.arange(len(H))[::-1]
cols_ = ["#d3d2cc", BLUE, BLUE, BLUE, BLUE, "#9ec5f4"]
ax.barh(yy, H.MAE, color=cols_, height=0.6)
for yi, v in zip(yy, H.MAE):
    ax.text(v + 0.04, yi, f"{v:.2f} d", va="center", fontsize=8.5, color=INK)
ax.set_yticks(yy, H.approach)
ax.set_xlabel("Leave-one-location-out MAE of the location mean (days)")
ax.grid(axis="y", visible=False)
ax.set_title("Honest out-of-sample error for predicting a NEW location's typical peak")
ax.text(0, -0.3, "Light bar: picking the best of 136 models on all 11 sites, then scoring it, is optimistic. "
                 "'Chosen in each fold' repeats the search without the held-out site.",
        transform=ax.transAxes, fontsize=8, color=INK2)
save(fig, "E11_honest_mu_error",
     "Leave-one-location-out error for predicting the typical peak of a location not used in fitting. Nested "
     "selection (search repeated inside every fold) gives the honest number for a model chosen from data.")

# ---------------------------------------------------------------- (d) influence + robustness
final_cols = ["thermal_date_15C"]
f = ols(final_cols)
infl = f.get_influence()
cook = infl.cooks_distance[0]
noS = P.index != "Shenandoah"
f_noS = sm.OLS(y[noS], sm.add_constant(P.loc[noS, final_cols])).fit()
f_geo_noS = sm.OLS(y[noS], sm.add_constant(P.loc[noS, ["lat", "elev_median_m"]])).fit()
wls = ols(final_cols, w=1 / sd_mu ** 2)
rob = pd.DataFrame([
    dict(fit="Thermal date, OLS, all 11", slope=f.params.iloc[1], se=f.bse.iloc[1], R2=f.rsquared),
    dict(fit="Thermal date, WLS (1/SE² of each μ)", slope=wls.params.iloc[1], se=wls.bse.iloc[1], R2=wls.rsquared),
    dict(fit="Thermal date, OLS, without Shenandoah", slope=f_noS.params.iloc[1], se=f_noS.bse.iloc[1],
         R2=f_noS.rsquared),
    dict(fit="Lat + elev, without Shenandoah: lat", slope=f_geo_noS.params["lat"], se=f_geo_noS.bse["lat"],
         R2=f_geo_noS.rsquared),
    dict(fit="Lat + elev, without Shenandoah: elev/100 m", slope=f_geo_noS.params["elev_median_m"] * 100,
         se=f_geo_noS.bse["elev_median_m"] * 100, R2=f_geo_noS.rsquared),
])
rob.to_csv(TAB / "E12_mu_robustness.csv", index=False)
print("\n", rob.round(3).to_string())

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw=dict(width_ratios=[1.1, 1]))
pred = f.fittedvalues
lo_ = lolo(final_cols)
a1.plot([268, 286], [268, 286], color=MUTED, ls="--", lw=1)
a1.scatter(pred, y, s=60, color=INK, edgecolor="white", lw=1.5, zorder=3)
off12 = {"Old Forge": (5, -9, "left"), "Montpelier": (-5, 2, "right"), "Valley Head": (-5, 2, "right"),
         "Frostburg": (5, 2, "left"), "Scranton": (-5, 2, "right")}
for xi, yi, nm in zip(pred, y, ORDER):
    dx_, dy_, ha_ = off12.get(nm, (5, -3, "left"))
    a1.annotate(nm, (xi, yi), xytext=(dx_, dy_), textcoords="offset points", fontsize=7.5, color=INK2, ha=ha_)
tk = [270, 274, 278, 282, 286]
a1.set_xticks(tk, [doy_label(v) for v in tk]); a1.set_yticks(tk, [doy_label(v) for v in tk])
a1.set_xlabel("Predicted from thermal date (15 °C)")
a1.set_ylabel("Observed mean peak")
a1.set_title(f"Observed vs. predicted (R² = {f.rsquared:.2f}, LOLO MAE {np.abs(lo_).mean():.2f} d)", fontsize=10)
yy = np.arange(n)[::-1]
a2.barh(yy, lo_, color=[RED if v > 0 else BLUE for v in lo_], height=0.6)
for yi, v, c in zip(yy, lo_, cook):
    a2.text(max(v, 0) + 0.15, yi, f"{v:+.1f} d  (Cook's D {c:.2f})", va="center", ha="left",
            fontsize=7.8, color=INK)
a2.axvline(0, color=INK2, lw=1)
a2.set_yticks(yy, ORDER)
a2.set_xlim(-4.5, 9)
a2.set_xlabel("Held-out residual: observed − predicted (days)")
a2.set_title("Leave-one-location-out residuals (red = later than climate predicts)", fontsize=10)
a2.grid(axis="y", visible=False)
fig.suptitle("Climate model of the location mean: fit and what it leaves unexplained", x=0.01, ha="left",
             fontsize=12, fontweight="bold")
fig.tight_layout()
save(fig, "E12_thermal_model_fit_residuals",
     "Left: observed location mean vs. the value predicted from the date the site's climatological mean "
     "temperature falls below 15 °C. Right: held-out residuals with Cook's distance; Shenandoah is ~5 d later "
     "than its climate implies.")
pd.Series(dict(b_lat=b_lat, b_elev_per100=b_el * 100, equiv_m_per_deg=equiv_m_per_deg,
               geo_R2=geo.rsquared, geo_p_lat=geo.pvalues["lat"], geo_p_elev=geo.pvalues["elev_median_m"])
          ).to_csv(TAB / "E13_geography_model_summary.csv")
