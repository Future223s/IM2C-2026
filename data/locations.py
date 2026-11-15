"""Location config for the 11-location Appalachian peak-foliage panel (v1).

Every build / analysis step imports from here. Coordinates of the 8 new towns were checked
against OpenStreetMap Nominatim on 2026-11-15 (all within 0.01 deg, no corrections needed).
Boone, Asheville and Shenandoah keep the exact points used in data/notebooks/DataPreprocessing.ipynb.
"""

LOCATIONS = {  # name: (lon, lat, state)
    "Montpelier":    (-72.576, 44.260, "VT"),
    "Old Forge":     (-74.972, 43.710, "NY"),
    "Binghamton":    (-75.918, 42.099, "NY"),
    "Scranton":      (-75.662, 41.409, "PA"),
    "State College": (-77.860, 40.793, "PA"),
    "Frostburg":     (-78.928, 39.658, "MD"),
    "Valley Head":   (-80.037, 38.546, "WV"),
    "Shenandoah":    (-78.40685152551191, 38.499092552797826, "VA"),  # keep exact old point
    "Roanoke":       (-79.941, 37.271, "VA"),
    "Boone":         (-81.68516978460896, 36.21834246106182, "NC"),   # keep exact old point
    "Asheville":     (-82.55831963941482, 35.633546529345516, "NC"),  # keep exact old point
}
ORDER = list(LOCATIONS)     # north -> south; use this order in every table and figure
BUFFER_M = 40_000
FOREST_THRESHOLD = 0.5
PHENO_YEARS = range(2001, 2026)
WX_YEARS = range(2001, 2027)
WX_START, WX_END = "03-01", "12-01"   # was 07-01; extended so spring/summer features exist
HELENE = {("Asheville", 2024), ("Boone", 2024)}

PROJECT = "fall-foliage-project"
CARRIED_OVER = ["Shenandoah", "Boone", "Asheville"]   # regression-test locations
