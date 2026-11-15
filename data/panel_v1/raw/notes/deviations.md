### Percentiles: exact instead of approximate (decided by Arin, 2026-11-15)

- First build: the regression test failed on `peak_doy` / `peak_doy_viirs` (up to 1.2 d) and Boone `elev_median_m` (0.001 m);
  weather, frost and NLCD columns matched to ~1e-12.
- Cause (verified): `ee.Reducer.percentile` is a bucketed approximation whose result for one region depends on the *other*
  regions in the same `reduceRegions` request. The ported code with the notebook's original 5 regions reproduced the old
  2001 medians exactly (Boone 275.354009, Asheville 278.252745, Shenandoah 283.411840); with the 11 regions they became
  275.312 / 277.553 / 284.000, and Boone alone 275.450. The old `model_df.csv` values therefore carry ~±1 d of this noise.
- Fix: every phenology and elevation percentile is now computed **exactly** from an unweighted per-region pixel histogram
  (dates and SRTM are integers), using numpy's 'linear' percentile definition. A location's value no longer depends on the
  other locations. Raw histograms: `raw/phenology_histograms.csv`, `raw/elevation_histograms.csv`.
- Side effect: "unweighted" counts a pixel only if its **center** is inside the clipped region. The old weighted reducers also
  counted edge pixels that partly overlap the region. `forest_px_count` and `*_count` are therefore 2–5% lower than in
  `model_df.csv` (Boone −359, Asheville −299, Shenandoah −298; the same every year).
- Regression test now: exact (< 1e-6) for weather, frost and `frac_hardwood`; tolerance 1.5 d for `peak_doy` /
  `peak_doy_viirs` and 5 m for `elev_median_m`.
- `raw/phenology_approx_percentile.csv` and `raw/static_approx_percentile.csv` keep the first (approximate) build for reference.
