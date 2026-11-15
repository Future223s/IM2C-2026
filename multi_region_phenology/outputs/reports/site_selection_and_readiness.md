# Multi-region autumn phenology dataset — release 0.1

## What is actually assembled

This release contains 146,692 PhenoCam GCC90/RCC90 observations from 30 deciduous ROIs on 16 camera platforms, plus 113,880 Daymet daily environmental rows for 12 ecological sites. Coverage ends with the last complete calendar year, 2025.

NEON individual-plant phenology, AmeriFlux, Landsat, Sentinel-2, and MODIS records are **not** represented as downloaded observations. They are retained in the inventory as verified-accessible or site-identifier candidates with null completeness until their records are ingested and audited.

## Recommended core and complementary set

The score ranks the current release only: 45% trajectory-ready camera years, 25% median camera completeness, 15% an active camera, and 15% a verified NEON site identifier. It is a transparent triage score, not a scientific weighting.

| Rank | Site | Elev. m | Ready camera years | Observed years | Median usable | Active | NEON candidate | Priority score |
|---:|---|---:|---:|---:|---:|:---:|:---:|---:|
| 1 | Bartlett Experimental Forest | 268 | 18 | 18 | 100.0% | yes | yes | 100.0 |
| 2 | Harvard Forest EMS Tower | 340 | 17 | 18 | 99.5% | yes | yes | 97.4 |
| 3 | Morgan–Monroe State Forest | 275 | 17 | 18 | 97.8% | yes | no | 82.0 |
| 4 | Cary Institute of Ecosystem Studies | 127 | 16 | 18 | 98.9% | yes | no | 79.7 |
| 5 | University of Michigan Biological Station | 230 | 15 | 15 | 100.0% | yes | no | 77.5 |
| 6 | Great Smoky Mountains — NEON tower | 589 | 9 | 9 | 100.0% | yes | yes | 77.5 |
| 7 | Oak Ridge — NEON tower | 351 | 9 | 9 | 96.7% | yes | yes | 76.7 |
| 8 | UMBS FASET experimental tower | 240 | 14 | 18 | 95.6% | yes | no | 73.9 |
| 9 | Great Smoky Mountains — Purchase Knob | 1550 | 16 | 19 | 94.5% | no | no | 63.6 |
| 10 | Shining Rock Wilderness | 1500 | 13 | 19 | 84.6% | no | no | 53.7 |
| 11 | Dolly Sods / Canaan Valley | 1133 | 8 | 11 | 92.3% | no | no | 43.1 |
| 12 | Oak Ridge — Chestnut Ridge | 371 | 6 | 8 | 100.0% | no | no | 40.0 |

A defensible first hierarchical core is Harvard, Bartlett, Morgan–Monroe, UMBS, and the NEON Great Smoky/Oak Ridge sites. Purchase Knob, Shining Rock, and Dolly Sods add high-elevation historical contrast but lack a current camera and currently have no ingested independent field modality. Cary is a strong camera-only secondary record. UMBS FASET is best treated as an experimental complement rather than pooled with the nearby UMBS flux tower.

## Verified coordinates, forest metadata, and limitations

Coordinates/elevations and dominant-species lists below come from the downloaded PhenoCam camera catalog. Species are presence-style metadata, not composition fractions.

| Site | Latitude | Longitude | Elev. m | Vegetation | Dominant species (verbatim) |
|---|---:|---:|---:|---|---|
| Harvard Forest EMS Tower | 42.53780 | -72.17150 | 340 | DB/EN | Quercus rubra, Acer rubrum, Pinus strobus |
| Bartlett Experimental Forest | 44.06460 | -71.28810 | 268 | DB/EN | Acer rubrum, Betula papyrifera, Fagus grandifolia |
| Great Smoky Mountains — Purchase Knob | 35.58600 | -83.07330 | 1550 | DB/GR | Betula alleghaniensis, Fagus grandifolia, Acer rubrum, Liriodendron tulipifera |
| Great Smoky Mountains — NEON tower | 35.68896 | -83.50195 | 589 | DB/ | Liriodendron tulipifera, Acer rubrum, Acer pensylvanicum |
| Morgan–Monroe State Forest | 39.32310 | -86.41310 | 275 | DB/ | Acer saccharum, Quercus rubra, Quercus alba, Liriodendron tilipifera, Sassafras albidum, Carya glabra, Fagus grandifolia |
| University of Michigan Biological Station | 45.55984 | -84.71382 | 230 | DB/EN | Populus grandidentata, Acer rubrum, Quercus rubra, Populus tremuloides, Betula papyrifera, Pinus strobus, Fagus grandifolia |
| UMBS FASET experimental tower | 45.56250 | -84.69760 | 240 | DB/EN | Acer rubrum, Quercus rubra, Pinus strobus, Acer saccharum, Fagus grandifolia (all Populus and Betula killed by stem girdling in 2008) |
| Cary Institute of Ecosystem Studies | 41.78390 | -73.73410 | 127 | DB/ | Quercus rubra |
| Oak Ridge — Chestnut Ridge | 35.93110 | -84.33230 | 371 | DB/ | Quercus alba, Quercus prinus, Carya ovata, Acer rubrum, Liriodendron tulipifera |
| Oak Ridge — NEON tower | 35.96413 | -84.28259 | 351 | DB/EN | Quercus montana, Liriodendron tulipifera, Cornus florida |
| Dolly Sods / Canaan Valley | 39.09953 | -79.42704 | 1133 | DB/ | Acer saccharum, Acer rubrum, Fagus grandifolia |
| Shining Rock Wilderness | 35.39016 | -82.77497 | 1500 | DB/ | Betula alleghaniensis, Fagus grandifolia, Acer rubrum, Liriodendron tulipifera |

## Preliminary modeling-readiness findings

Across site-year median RCC90 maxima, the median within-site interannual SD is 6.0 days. The SD of site mean dates is 7.6 days and their range is 24.4 days. For site-years with multiple camera platforms, the median platform spread is 1.0 days.

This supports a hierarchical model: there is material temporal variation within sites and material spatial structure between sites. But the platform spread also shows that some apparent variability is measurement-footprint disagreement. A latent-state model should therefore include modality/platform observation functions rather than treating each date as interchangeable truth.

The environmental table includes raw Daymet variables plus trailing, lagged cooling summaries at 3, 7, 14, 21, and 30 days. Here, ‘several scales’ means several look-back windows over the same daily temperature process. Every rolling feature ends one day before its row date, preventing future leakage for date-indexed forecasts.

## Critical limitations

- RCC90 maximum is an operational red-chromatic-coordinate event, not a validated universal peak-color date.
- Daymet is modeled/interpolated at 1 km; it is not a direct canopy or tower measurement. Its vapor pressure and radiation are secondary estimates.
- Camera ROIs, satellite pixels, tagged trees, and tower footprints are not coextensive. The schema preserves their identifiers so future joins can express proximity rather than false co-location.
- Camera view/ROI changes are retained as separate `platform_id` and `roi_id` values. Adjacent platforms are not averaged in the observation table.
- Species metadata say which taxa are described as dominant but do not quantify basal area, stem fraction, or camera-visible canopy fraction.
- Satellite clear-sky completeness and NEON individual-tree completeness remain uncalculated until records are downloaded.

## Next ingestion priority

1. Download NEON DP1.10055.001 for HARV, BART, GRSM, and ORNL with an API token; retain `individualID`, taxon, phenophase status/intensity, and visit dates.
2. Extract Landsat/Sentinel-2 observations over explicit camera ROI/tower buffers; retain cloud masks, pixel IDs, acquisition times, and footprint geometry.
3. Extract MCD12Q2 V6.1 pixel-level senescence/dormancy metrics and QA, labeling them as greenness phenometrics rather than color peaks.
4. Audit site-specific AmeriFlux licenses and daily meteorological coverage before replacing or supplementing Daymet.
5. Add quantitative vegetation inventories (NEON/FIA/LTER) with composition measure and footprint explicit.

## Inventory status counts

- verified accessible; pending site-footprint extraction: 36
- downloaded: 24
- site identifier verified; product/license must be checked per site: 8
- verified accessible; download requires NEON API token: 4
