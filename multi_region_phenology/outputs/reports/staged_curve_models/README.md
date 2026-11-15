# Staged location-baseline plus environmental-deviation models

Each model adds one scientific claim. All numerical comparisons use the same outer leave-one-year-out protocol, training-only site baselines, and nested year-grouped tuning for ridge models.

| Model | Additional question | MAE (days) | Skill vs M1 |
|---|---|---:|---:|
| M1 | What is normal at this location? | 4.95 | 0.0% |
| M2 | Does a shared temperature anomaly shift timing? | 4.22 | 14.8% |
| M3 | Do additional environmental histories help? | 4.30 | 13.0% |
| M4 | Does temperature sensitivity differ by location? | 4.48 | 9.5% |
| M5 | Which curve properties are predictable? | — | — |
| M6 | Reconstruct the selected full curve | 4.22 | 14.8% |

Selection rule: prefer the least complex model within 0.10 day MAE of the best held-out result. Do not retain an environmental effect because it improves training fit.

Important scope: the present September–October variables make this a retrospective explanatory experiment. Operational forecasts require issue-date truncation or explicit weather forecasts.
