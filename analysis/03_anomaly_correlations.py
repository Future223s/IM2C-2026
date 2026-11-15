"""Task 3: correlations of within-region anomalies with the target anomaly."""
import numpy as np
import pandas as pd
from scipy import stats

from common import FEATURES, TARGET, TIMING, load, save, plt, ROOT

FEATS = FEATURES + ["year"]
full = load()
print(f"first_frost censored among target rows: {int(full.frost_censored.sum())} of {len(full)}")
print(f"pdsi_AugSep non-null among target rows: {full.pdsi_AugSep.notna().sum()}")


def target_corr(df):
    A = anomalies_with_year(df)
    rows = []
    # year-mean anomalies: one point per year (n=25). Regions share weather, so this
    # is a conservative check on the effective sample size.
    Y = A.groupby("year")[FEATS[:-1] + [TARGET]].mean().reset_index()
    for f in FEATS:
        ok = A[[f, TARGET]].dropna()
        r, p = stats.pearsonr(ok[f], ok[TARGET])
        rs, _ = stats.spearmanr(ok[f], ok[TARGET])
        yv = Y["year"] if f == "year" else Y[f]
        ry, py = stats.pearsonr(yv, Y[TARGET])
        rows.append({"feature": f, "timing": TIMING[f], "n": len(ok), "r": r, "p": p,
                     "spearman": rs, "r_yearmean": ry, "p_yearmean": py})
    return pd.DataFrame(rows).sort_values("r", key=np.abs, ascending=False)


def anomalies_with_year(df):
    A = df[["region", "year", "helene_2024"]].copy()
    for c in FEATURES + [TARGET]:
        A[c] = df[c] - df.groupby("region")[c].transform("mean")
    A["year"] = df.year  # year needs no de-meaning (same years in every region)
    return A


res_all = target_corr(full)
res_noh = target_corr(full[full.helene_2024 == 0])
print("\nAll rows (n=100)\n", res_all.round(3).to_string(index=False))
print("\nWithout Helene 2024 rows (n=98)\n", res_noh.round(3).to_string(index=False))
cmp = res_all[["feature", "r"]].merge(res_noh[["feature", "r"]], on="feature",
                                      suffixes=("_all", "_noHelene"))
cmp["delta"] = cmp.r_noHelene - cmp.r_all
print("\nChange in r when Helene rows dropped\n", cmp.round(3).to_string(index=False))
res_all.round(4).to_csv(ROOT / "anomaly_corr_target.csv", index=False)
res_noh.round(4).to_csv(ROOT / "anomaly_corr_target_noHelene.csv", index=False)


def target_heatmap(res, title, fname):
    fig, ax = plt.subplots(figsize=(4.2, 0.32 * len(res) + 1))
    vals = res[["r"]].values
    ax.imshow(vals, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    labels = [f"{f}  [{t}]" for f, t in zip(res.feature, res.timing)]
    ax.set_yticks(range(len(res)), labels)
    ax.set_xticks([0], ["r with peak-date anomaly"])
    for i, (r, p) in enumerate(zip(res.r, res.p)):
        star = "**" if p < 0.01 else ("*" if p < 0.05 else "")
        ax.text(0, i, f"{r:+.2f}{star}", ha="center", va="center", fontsize=8,
                color="white" if abs(r) > 0.55 else "#1a1a19")
    ax.grid(False)
    ax.set_title(title, fontsize=9.5)
    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-1, 1))
    fig.colorbar(sm, ax=ax, fraction=0.08, pad=0.04, label="Pearson r")
    fig.text(0.0, -0.01, "Within-region anomalies; * p<0.05, ** p<0.01 (naive n)",
             fontsize=7, color="#5f5f5a")
    save(fig, fname)


target_heatmap(res_all, "Feature vs peak-date anomaly (n=100)", "03_target_corr.png")
target_heatmap(res_noh, "Feature vs peak-date anomaly, no Helene (n=98)",
               "03_target_corr_noHelene.png")

# ---- feature x feature matrix -------------------------------------------
A = anomalies_with_year(full)
order = list(res_all.feature)
C = A[order].corr()
C.round(3).to_csv(ROOT / "anomaly_corr_features.csv")
fig, ax = plt.subplots(figsize=(8.5, 7.5))
ax.imshow(C.values, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(len(order)), order, rotation=60, ha="right")
ax.set_yticks(range(len(order)), order)
for i in range(len(order)):
    for j in range(len(order)):
        v = C.values[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6.5,
                color="white" if abs(v) > 0.55 else "#1a1a19")
ax.grid(False)
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-1, 1))
fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.03, label="Pearson r")
ax.set_title("Feature-by-feature correlation of within-region anomalies (n=100)")
save(fig, "03_feature_corr_matrix.png")
