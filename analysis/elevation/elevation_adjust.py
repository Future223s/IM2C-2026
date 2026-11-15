"""Turn a regional peak-color forecast into an elevation-specific one.

Reference elevation = the region's `elev_median_m` from static_features.csv (same as the
main model). Slopes come from results.json (02_elevation_analysis.py).

Two rules:
- adjust_for_elevation: the linear rule from the handoff (one slope).
- adjust_for_elevation_two_piece: RECOMMENDED. The data are flat below ~1,350 m and steep
  above it, so a single slope understates high-elevation shifts and overstates low ones.

Run this file to print/save the offset table (offset_table.csv).
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
STATIC = HERE.parents[1] / "data" / "data" / "static_features.csv"


def adjust_for_elevation(region_forecast_doy, target_elev_m, region_median_elev_m,
                         slope_days_per_100m):
    return region_forecast_doy + slope_days_per_100m * (target_elev_m - region_median_elev_m) / 100


def _two_piece(elev_m, hinge_m, below, above):
    return (below * min(elev_m, hinge_m) + above * max(elev_m - hinge_m, 0)) / 100


def adjust_for_elevation_two_piece(region_forecast_doy, target_elev_m, region_median_elev_m,
                                   hinge_m, slope_below, slope_above):
    """Slopes in days per 100 m below / above hinge_m."""
    return region_forecast_doy + (_two_piece(target_elev_m, hinge_m, slope_below, slope_above)
                                  - _two_piece(region_median_elev_m, hinge_m, slope_below, slope_above))


def load_params():
    res = json.load(open(HERE / "results.json"))
    tp = res["two_piece_best"]
    return {"linear": res["headline_mountain"]["elev100"]["est"],
            "hinge_m": tp["hinge_m"], "below": tp["below100"]["est"], "above": tp["above100"]["est"],
            "region_specific": {r: {"below": res[f"two_piece_{r}"]["below"]["est"],
                                    "above": res[f"two_piece_{r}"]["above"]["est"]}
                                for r in ["Asheville", "Boone"]},
            "max_sampled_elev_m": 1819}


def median_elev():
    s = pd.read_csv(STATIC).set_index("name").elev_median_m
    return s.to_dict()


GRID = [500, 800, 1000, 1200, 1500, 1800]
PLACES = {"Asheville": [(650, "Asheville city"), (2037, "Mount Mitchell")],
          "Boone": [(1800, "Grandfather Mountain")]}


def offset_table():
    p, med = load_params(), median_elev()
    rows = []
    for r in ["Asheville", "Boone"]:
        rs = p["region_specific"][r]
        targets = sorted({**{e: "" for e in GRID}, **dict(PLACES[r])}.items())
        for e, label in targets:
            rows.append({
                "region": r, "region_median_elev_m": round(med[r]), "target_elev_m": e, "place": label,
                "offset_two_piece_d": adjust_for_elevation_two_piece(0, e, med[r], p["hinge_m"], p["below"], p["above"]),
                "offset_linear_d": adjust_for_elevation(0, e, med[r], p["linear"]),
                "offset_region_specific_d": adjust_for_elevation_two_piece(0, e, med[r], p["hinge_m"], rs["below"], rs["above"]),
                "extrapolated": e > p["max_sampled_elev_m"],
            })
    return pd.DataFrame(rows)


if __name__ == "__main__":
    t = offset_table()
    t.to_csv(HERE / "offset_table.csv", index=False)
    print(t.round(1).to_string(index=False))
