"""Step 6: Synthesis figures.
E23: each location's mean = grand mean + climate part (thermal-date model) + local residual
E24: predicted typical peak over latitude x elevation (geography model mu ~ lat + elev), with prediction SE
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm

from common import *  # noqa

panel, locs, train, ORDER = load()
P = pd.read_csv(TAB / "E04_location_predictors.csv").set_index("location").loc[ORDER]
y = P.peak_doy_mean.values
gm = y.mean()

# ---------------------------------------------------------------- E23 decomposition
th = sm.OLS(y, sm.add_constant(P[["thermal_date_15C"]])).fit()
climate = th.fittedvalues.values - gm
local = y - th.fittedvalues.values
D = pd.DataFrame(dict(location=ORDER, mean=y, climate=climate, local=local))
D.to_csv(TAB / "E24_mu_decomposition.csv", index=False)
fig, ax = plt.subplots(figsize=(10, 5))
yy = np.arange(len(ORDER))[::-1]
h = 0.36
ax.barh(yy + h / 2, climate, height=h, color=BLUE, label="Climate part (thermal-date model)")
ax.barh(yy - h / 2, local, height=h, color=ORANGE, label="Local residual (not explained by climate)")
for yi, c, l, m in zip(yy, climate, local, y):
    ax.text(6.6, yi, f"{doy_label(m)}  ({m - gm:+.1f} d)", va="center", fontsize=8.5, color=INK)
ax.axvline(0, color=INK2, lw=1)
ax.set_yticks(yy, ORDER)
ax.set_xlim(-8, 9.5)
ax.text(6.6, len(ORDER) - 0.35, "Mean peak (vs. all-site mean)", fontsize=8, color=INK2, fontweight="bold")
ax.set_xlabel(f"Days relative to the 11-location mean ({doy_label(gm)})")
ax.legend(loc="lower left", fontsize=8.5)
ax.grid(axis="y", visible=False)
ax.set_title("Why each location peaks when it does: climate vs. local factors")
save(fig, "E23_mu_decomposition",
     "Each location's mean peak relative to the 11-site mean, split into the part predicted by its climate "
     "(date the climatological mean temperature falls below 15 °C) and the local residual. Climate explains "
     "the north–south pattern; Shenandoah (+4 d) and Roanoke (−3 d) carry the largest local effects.")

# ---------------------------------------------------------------- E24 lat x elev lookup
geo = sm.OLS(y, sm.add_constant(P[["lat", "elev_median_m"]])).fit()
lat = np.linspace(35, 45, 101)
el = np.linspace(200, 1400, 121)
LA, EL = np.meshgrid(lat, el)
Xg = sm.add_constant(np.column_stack([LA.ravel(), EL.ravel()]), has_constant="add")
pr = geo.get_prediction(Xg).summary_frame(alpha=0.05)
Z = pr["mean"].values.reshape(LA.shape)
SE = pr["mean_se"].values.reshape(LA.shape)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw=dict(width_ratios=[1.25, 1]))
levels = np.arange(262, 292, 2)
cf = a1.contourf(LA, EL, Z, levels=levels, cmap=SEQ_CMAP, extend="both")
cs = a1.contour(LA, EL, Z, levels=levels[::2], colors="white", linewidths=0.8)
a1.clabel(cs, fmt=lambda v: doy_label(v), fontsize=7.5)
a1.scatter(P.lat, P.elev_median_m, s=50, color=INK, edgecolor="white", lw=1.5, zorder=3)
for nm, la, e in zip(ORDER, P.lat, P.elev_median_m):
    dy_ = -11 if nm in ("State College", "Binghamton") else 4
    a1.annotate(nm, (la, e), xytext=(4, dy_), textcoords="offset points", fontsize=7.5, color=INK)
cb = fig.colorbar(cf, ax=a1, shrink=0.85)
cb.set_ticks(levels[::3]); cb.set_ticklabels([doy_label(v) for v in levels[::3]])
cb.set_label("Predicted typical peak")
a1.set_xlabel("Latitude (°N)")
a1.set_ylabel("Median forest elevation (m)")
a1.set_title(f"μ = {geo.params['const']:.0f} {geo.params['lat']:+.2f}·lat {geo.params['elev_median_m'] * 100:+.2f}·elev/100 m"
             f"   (R² {geo.rsquared:.2f})", fontsize=9.5)
a1.grid(False)
cf2 = a2.contourf(LA, EL, 1.96 * SE, levels=np.arange(0, 9, 1), cmap=SEQ_CMAP, extend="max")
a2.scatter(P.lat, P.elev_median_m, s=30, color=INK, edgecolor="white", lw=1, zorder=3)
cb2 = fig.colorbar(cf2, ax=a2, shrink=0.85)
cb2.set_label("± days (95% CI of the mean)")
a2.set_xlabel("Latitude (°N)")
a2.set_title("Uncertainty: grows away from the 11 sites", fontsize=9.5)
a2.grid(False)
fig.suptitle("Typical peak for any Appalachian place from latitude and elevation (geography model)", x=0.01,
             ha="left", fontsize=12, fontweight="bold")
fig.tight_layout()
save(fig, "E24_lat_elev_lookup",
     "Predicted typical MODIS peak date over latitude and median forest elevation from μ ~ lat + elev (n = 11), "
     "with the 95% confidence half-width. Use inside the sampled range (dots); add the within-location elevation "
     "rule for individual summits.")
print(D.round(2).to_string())
print(geo.summary().tables[1])
