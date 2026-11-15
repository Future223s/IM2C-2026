# How to use this dataset — instructions for you and your AI agent

**Give this file, `DATASET_GUIDE.md`, and the unzipped `panel_v1/` folder to your AI agent at the start of the session.**
This file says *what to do and how*. `DATASET_GUIDE.md` is the column dictionary (what every column means). If the two ever
disagree about a column, `DATASET_GUIDE.md` wins; if they disagree about method, this file wins.

> **Agent: read this whole file before writing any code.** Then read `DATASET_GUIDE.md` §2–§5.
> Follow the shared rules in §3 exactly, even where a different method seems better. Our three analyses
> only combine into one model if everyone uses the same anomalies, the same validation and the same hand-back table.
> If you think a rule is wrong, finish the analysis the agreed way and add your concern in the notes section of your hand-back.

---

## 1. The project in five sentences

1. IMMC "Peak Peeping": predict the **date of peak fall foliage in 2026** at 11 Appalachian locations (Vermont → North Carolina), compare with news reports, and say how accurate the model would be **10 years out**.
2. The target is `peak_doy`: the MODIS satellite "mid-greendown" day of year (DOY 274 = Oct 1), a **proxy** that falls about 1–3 weeks before the peak color people report. Model it as it is; the offset to news-reported peak is handled separately at the end.
3. The agreed model is **`peak(location, year) = μ(location) + δ(location, year)`**: μ is the location's typical peak, δ is this year's departure from it.
4. Year-to-year SD is only about **4 days**, and the honest benchmark is **climatology** (predict μ) with leave-one-year-out **MAE 3.34 d / RMSE 4.13 d**. Every claim of "this helps" means "this beats 3.34 d".
5. Static things (elevation, species mix, terrain) **never change within a location**, so they cannot move δ on their own. They either explain why μ differs between locations, or they change how strongly a location responds to weather (an interaction).

## 2. What is in the folder

| File | Grain | Use it for |
|---|---|---|
| `panel.csv` | 11 locations × 2001–2026 = 286 rows | **Main table.** Target, weather features, satellite phenology. 2026 rows have no target (`is_forecast_row = 1`). |
| `locations.csv` | 11 rows | Static features: geography, elevation, terrain, land cover, species groups (`evt_*`), day length, climate normals, typical peak. |
| `weather_daily.csv` | location × day, Mar 1–Nov 30 | Build your own weather windows. Ends 2026-11-11. |
| `pdsi.csv` | location × 5-day step | Drought index, if you need custom windows. |
| `feature_timing.csv` | feature × location | Is each feature known **before** that location's typical peak? (`pre-peak` / `partly` / `post`) |
| `evt_class_map.csv` | LANDFIRE class | Which vegetation class went into which `evt_*` group (species teammate: audit or regroup). |
| `DATASET_GUIDE.md` | — | Column dictionary, units, caveats. |
| `BUILD_LOG.md` | — | How the data were built and checked (regression test, ranges, baseline per location). |
| `pixels.csv` *(separate zip: `panel_v1_pixels.zip`)* | forest pixel × year | Pixel-level peak date + elevation. Environmental work only. |

Keys: `location` (string, exact spelling as in `locations.csv`) and `year` (int). Always order tables north → south using
`locations.csv` row order.

**Facts an agent tends to get wrong:**
- 275 rows have a target (2001–2025). Effective sample size for year-to-year questions is **about 25 years**, not 275 rows: all 11 locations share each year's weather.
- `peak_doy_mean` and `peak_anom` in `panel.csv` use **all** 2001–2025 years. Fine for exploring; **leaks** when scoring. When you score a model, recompute the mean without the held-out year (§4 code does this).
- Weather is a mean over the whole region (not forest-only). Temperatures °C, precipitation mm totals, srad W/m², VPD kPa.
- `first_frost_doy = 335` with `frost_censored = 1` means "no frost seen" (in 2026: no frost by Nov 11).

## 3. Shared rules (everyone, every analysis)

1. **Work in within-location anomalies.** For any year-varying column `x`: `x̃ = x − mean of x for that location over the training years`. Do the same to the target. Relate anomalies to anomalies.
2. **Never use these as predictors:** `dormancy_p50`, `peak_doy_p25`, `peak_doy_p75`, `peak_spread_days`, `peak_doy_viirs`, `mid_greendown_count`, `peak_date`, `peak_anom`, `peak_doy_mean` (as a δ predictor). They are measured at or after the target, or are the target.
3. **Forecasting claims use only forecast-safe features**: `timing == "pre-peak"` in `feature_timing.csv` for that location. September features are only `partly` known for Montpelier and Old Forge (their peak ≈ Sep 28). Prefer the relative-window versions (`*_pre30`, `tmax_pre60to31`, `*_pre60`), which are pre-peak everywhere by construction. Explanatory claims ("why was 2019 late?") may use any feature; say which kind of claim you are making.
4. **Significance uses year-clustered standard errors** (or year-averaged data). Never report naive OLS p-values on 275 rows.
5. **Validate by leaving one year out** (hold out all 11 locations of that year together) and compare against climatology on the same folds. In-sample R² is a screen, not evidence.
6. **Budget:** the final δ model gets **≤ 5 main effects + ≤ 2 interactions** in total across the team. Be selective; rank your candidates.
7. **A term earns its place only if** it lowers LOYO MAE, its sign is consistent across most locations (or an interaction explains why not), and there is a physical reason for it.
8. **Report Helene sensitivity** (drop `helene_2024 == 1` rows) for any headline number. In the old data it changed nothing; confirm it here.
9. Don't modify the CSVs. If you create a new column, write down its exact formula (it goes in your hand-back).

## 4. Copy-paste code (tested against this data)

```python
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

P = "panel_v1/"                      # path to the unzipped folder
panel = pd.read_csv(P + "panel.csv")
locs = pd.read_csv(P + "locations.csv")
timing = pd.read_csv(P + "feature_timing.csv")
ORDER = list(locs["location"])       # north -> south
train = panel[panel["is_forecast_row"] == 0].copy()

def add_anoms(df, cols):
    """Within-location anomalies over the rows in df (use training rows only when scoring)."""
    df = df.copy()
    for c in cols:
        df[c + "_a"] = df[c] - df.groupby("location")[c].transform("mean")
    return df

def forecast_safe(feature):
    """True if the feature is pre-peak at every location."""
    t = timing[timing["feature"] == feature]["timing"]
    return len(t) > 0 and (t == "pre-peak").all()
```

**Effect of one feature on the anomaly, year-clustered SE:**
```python
d = add_anoms(train, ["peak_doy", "tmax_Aug"])
fit = smf.ols("peak_doy_a ~ tmax_Aug_a", data=d).fit(
    cov_type="cluster", cov_kwds={"groups": d["year"]})
print(fit.params["tmax_Aug_a"], fit.bse["tmax_Aug_a"], fit.pvalues["tmax_Aug_a"])
# days of peak shift per +1 unit of the feature, ± SE, p
```

**Per-location slopes (is the sign consistent?):**
```python
rows = []
for loc, g in d.groupby("location"):
    f = smf.ols("peak_doy_a ~ tmax_Aug_a", data=g).fit()
    rows.append(dict(location=loc, slope=f.params.iloc[1], se=f.bse.iloc[1]))
slopes = pd.DataFrame(rows).set_index("location").loc[ORDER]
```

**Leave-one-year-out score vs climatology (the test every result must pass):**
```python
from sklearn.linear_model import Ridge

def loyo(df, features, alpha=1.0, interactions=()):
    """Hold out each year (all locations). mu and anomalies use training years only.
    interactions: list of (weather_feature, static_column) pairs -> x̃ * centered static.
    Returns per-row predictions with columns pred_clim and pred_model."""
    out = []
    for y in sorted(df["year"].unique()):
        tr, te = df[df["year"] != y], df[df["year"] == y]
        mu = tr.groupby("location")["peak_doy"].mean()
        def design(part):
            if not features:
                return np.zeros((len(part), 0))
            xm = tr.groupby("location")[features].mean()
            X = part[features].values - xm.loc[part["location"]].values
            cols = [X]
            for w, s in interactions:
                st = locs.set_index("location")[s]
                s_c = (st - st.mean()).loc[part["location"]].values
                xw = part[w].values - tr.groupby("location")[w].mean().loc[part["location"]].values
                cols.append((xw * s_c)[:, None])
            return np.hstack(cols)
        y_tr = tr["peak_doy"].values - mu.loc[tr["location"]].values
        Xtr, Xte = design(tr), design(te)
        delta = Ridge(alpha=alpha).fit(Xtr, y_tr).predict(Xte) if Xtr.shape[1] else 0.0
        base = mu.loc[te["location"]].values
        out.append(te[["location", "year", "peak_doy"]].assign(pred_clim=base, pred_model=base + delta))
    return pd.concat(out, ignore_index=True)

def score(res):
    e_c, e_m = res["peak_doy"] - res["pred_clim"], res["peak_doy"] - res["pred_model"]
    return pd.Series({"MAE_clim": e_c.abs().mean(), "MAE_model": e_m.abs().mean(),
                      "RMSE_clim": np.sqrt((e_c**2).mean()), "RMSE_model": np.sqrt((e_m**2).mean()),
                      "bias_model": e_m.mean()})

res = loyo(train, ["tmax_pre60to31", "precip_pre60"])
print(score(res))                                            # overall
print(res.groupby("location").apply(score, include_groups=False).loc[ORDER].round(2))   # per location
```
`MAE_clim` reproduces the 3.34 d baseline. Scale matters for Ridge: with very different feature units, standardize the
anomalies first or use `alpha` small (e.g. 0.1) and check results are stable.

**Train early, test late** (the "does it age well?" test): fit on 2001–2015 and score on 2016–2025, with climatology
computed from 2001–2015 only. Adapt `loyo` by replacing the loop with a single split.

## 5. Your assignment

Find your role below and do those steps in order. Each ends with the hand-back table in §6.

### 5a. Temporal teammate — what moves a given year earlier or later (the δ main effects)

Goal: a ranked shortlist of **≤ 5 forecast-safe weather/phenology features** for δ, with shared slopes.

1. Screen every year-varying column in `panel.csv` §4b–4c (skip the forbidden list in §3.2). For each: anomaly correlation with `peak_doy` anomaly (pooled), year-mean correlation (collapse to 25 points), year-clustered slope ± SE.
2. Compare each **fixed** window with its **relative** counterpart (`tmax_Aug` vs `tmax_pre60to31`, `cdd_Sep` vs `cdd_pre30`, `precip_AugSep` vs `precip_pre60`, `srad_Sep` vs `srad_pre30`, `tmin_Sep` vs `tmin_pre30`). The plan prefers relative windows because they are pre-peak everywhere; keep a fixed window only if it is clearly better in LOYO.
3. Check redundancy (|r| > 0.8 between anomalies → keep one per cluster) and per-location slope consistency (equal-slope F-test or a mixed model with random slopes, `smf.mixedlm`).
4. Forward-select with `loyo()`: add the feature that most lowers LOYO MAE; stop when nothing lowers it by ≥ 0.05 d.
5. Test the spring features (`tmean_MarApr`, `gdd5_MarMay`, `mid_greenup_p50`) — a late spring went with a late fall before (r ≈ +0.32).
6. Starting points from the old 4-region analysis: `tmax_Aug` ≈ **+1.1 d/°C** (very consistent), `mid_greenup_p50` +0.27 d/d, `precip_AugSep` −0.011 d/mm, `cdd_Sep` −0.04 d/°C·d. Report whether they hold across 11 locations, especially the 7 new northern ones.
7. Optionally, fit the process alternative (plan §3b): color midpoint when cold degree-days accumulated from a photoperiod-triggered start cross a threshold. Use `weather_daily.csv` and `doy_daylen_below_11h` from `locations.csv`; fit 2–3 parameters by minimizing LOYO error.

### 5b. Environmental teammate (Arin) — why locations differ (μ) and where sensitivity changes

Goal: (a) a spatial model of μ with ≤ 2–3 predictors; (b) ≤ 2 candidate interactions; (c) the within-location elevation rule refit.

1. **μ across space, n = 11.** Regress `peak_doy_mean` on one or two of: `lat` or `daylen_Oct1_h`, `elev_median_m`, `norm_tmean_Sep`, `norm_tmin_Oct`, `evt_northern_hardwood`. Validate by **leaving one location out** (score against the grand mean). Never more than 3 predictors.
   *Heads-up:* latitude alone is not the whole story (Frostburg and Valley Head at Oct 2 are earlier than Scranton/State College at Oct 4; Shenandoah, Oct 10, is latest). `frac_forest_above_1350m` is ≈ 0 except Boone (0.03) and Asheville (0.05) — useless as a cross-location predictor.
2. **Interactions:** does the response to a weather anomaly depend on a static feature? Fit `peak_doy_a ~ x_a + x_a:s_c` with year-clustered SE, where `s_c` is the static feature centered across locations (e.g. does the `tmax_pre60to31` effect grow with `elev_median_m` or with `norm_tmean_Sep`?). Keep only interactions that lower LOYO MAE via `loyo(..., interactions=[(x, s)])`.
3. **Elevation within a location** (`pixels.csv`): refit the old two-piece rule (flat below ~1,350 m, about −2.9 d per 100 m above; see `analysis/elevation/` in the repo). Compute each pixel's anomaly vs its location-year median, regress on elevation, block-bootstrap over years. Key question: is 1,350 m a southern-Blue-Ridge artifact? Only Boone and Asheville have forest above 1,350 m, so test whether the northern sites show a break at a *lower* elevation.
4. `northness_mean` and `slope_mean_deg` are location averages near zero/similar; northness is more useful at pixel level (it had −4.3 d per unit before).

### 5c. Species teammate — how forest composition changes timing

Goal: show whether species mix explains μ differences and/or changes weather sensitivity.

1. **Audit the groups first.** Open `evt_class_map.csv` (113 classes present, 60 forest). The six `evt_*` columns are shares of all LANDFIRE tree pixels and sum to 1. Grouping: explicit northern-hardwood names first, otherwise the first-named genus (oak-pine → oak). Ruderal/plantation/floodplain/swamp/urban forest → `evt_other`. **State College's `evt_other` (0.37) is mostly "Ruderal Forest – Northern and Central Hardwood and Conifer"** (young regrowth). If you regroup, publish your new mapping and the new column formulas.
2. LANDFIRE EVT here is the **2016** vegetation map (the only version in Earth Engine). NLCD `frac_*` columns are 2021 and are shares of *total region area*, not of forest. Don't mix the denominators.
3. **μ:** with n = 11, test one species variable at a time against `peak_doy_mean` (leave-one-location-out), and check whether it adds anything beyond latitude/normal temperature (species mix is strongly tied to latitude here: northern hardwood 0.74 in VT/NY vs 0.04 in NC).
4. **Interactions:** the key question for δ is whether species changes sensitivity: e.g. does `tmax_pre60to31` delay peak more where `evt_oak_hickory` is high? Same method as §5b step 2.
5. If you bring in FIA / USFS forest-type data, compute it over the **same regions** (40 km circle clipped to the EPA L3 ecoregion; ask Arin for the geometries) and document it exactly.

## 6. What to hand back (all roles)

A short Markdown file with this table, one row per feature or term you tested — including the ones that failed:

| feature | role (δ main effect / interaction with ___ / explains μ) | effect size (days per unit, ± SE) | forecast-safe? | evidence (r, p with year-clustered SE, LOYO MAE change vs climatology 3.34) | notes |
|---|---|---|---|---|---|

Then:
- **Ranked recommendation:** your top candidates for the final model (respect the ≤ 5 + 2 budget).
- **Did not help:** features you tried that failed, and how badly.
- **Data problems** you found (with location/year).
- **New columns** with the exact formula.
- **Helene check:** did dropping the 2024 Asheville/Boone rows change any headline number?

## 7. Things the agent should not do

- Don't predict 2026 with features that are `partly` or `post` for that location.
- Don't use VIIRS as a target or predictor. MODIS and VIIRS disagree on the trend (VIIRS runs ≈ +10 d later here, and the gap grows ≈ +7 d/decade). It is a robustness check only.
- Don't claim a trend ("peak is getting later/earlier") from one sensor. Per-location trends are in `locations.csv` (`peak_doy_trend_d_per_decade`, all between −2.4 and −0.1 d/decade, none tested for significance here).
- Don't treat 275 rows as 275 independent observations.
- Don't add a time trend term unless it lowers LOYO MAE (it didn't in the old data).
- Don't overwrite or "fix" the data. Report problems to Arin.
