"""Supplementary: a summer day in each district, from MODIS Terra + Aqua (1 km).

No satellite watches a district for 24 hours. MODIS on Terra (~10:30 / ~22:30)
and Aqua (~13:30 / ~01:30) gives four surface-temperature snapshots a day.
For summer 2024 (Jun-Sep) we average each pass per district and for lowland
forest cells, using the per-pixel view time to place each pass on the clock.

Kept on its own 1 km grid: it is never combined with the 30 m Landsat statistics.
Writes data/processed/diurnal.json.
"""
import json

import numpy as np
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import from_origin
from shapely.geometry import mapping

import common as C

KM = 1000.0
t = from_origin(np.floor(C.TRANSFORM.c / KM) * KM, np.ceil(C.TRANSFORM.f / KM) * KM, KM, KM)
w, h = int(np.ceil(C.WIDTH * C.RES / KM)) + 1, int(np.ceil(C.HEIGHT * C.RES / KM)) + 1
cat = C.catalog()
items = list(cat.search(collections=["modis-11A2-061"], bbox=C.BBOX, datetime="2024-06-01/2024-09-30").items())

passes = {"terra_day": ("MOD11A2", "LST_Day_1km", "Day_view_time"),
          "terra_night": ("MOD11A2", "LST_Night_1km", "Night_view_time"),
          "aqua_day": ("MYD11A2", "LST_Day_1km", "Day_view_time"),
          "aqua_night": ("MYD11A2", "LST_Night_1km", "Night_view_time")}
grids = {}
for key, (prod, lst_asset, vt_asset) in passes.items():
    lsts, vts = [], []
    for it in [i for i in items if i.id.startswith(prod)]:
        a, _ = C.read_to_grid(it.assets[lst_asset].href, Resampling.nearest, t, w, h, "float32")
        v, _ = C.read_to_grid(it.assets[vt_asset].href, Resampling.nearest, t, w, h, "float32")
        ok = a > 0
        lsts.append(np.where(ok, a * 0.02 - 273.15, np.nan))
        vts.append(np.where(ok & (v < 255), v * 0.1, np.nan))   # local solar time, hours
    with np.errstate(all="ignore"):
        grids[key] = (np.nanmean(np.stack(lsts), 0), np.nanmean(np.stack(vts), 0))
    print(key, len(lsts), "composites")

# districts and lowland forest on the 1 km grid
d_fc = C.load_geojson("hk_districts.geojson")["features"]
dist = rasterize([(mapping(C.geom_to_utm(f["geometry"])), i) for i, f in enumerate(d_fc, 1)],
                 out_shape=(h, w), transform=t, fill=0, dtype="int16")
S = C.load("static")
# share of each 1 km cell that is >= 80 % tree below 200 m (from the 30 m layers)
rows, cols = np.indices(C.SHAPE)
xs = C.TRANSFORM.c + (cols + 0.5) * C.RES
ys = C.TRANSFORM.f - (rows + 0.5) * C.RES
ci = ((xs - t.c) // KM).astype(int); ri = ((t.f - ys) // KM).astype(int)
inside = (ci >= 0) & (ci < w) & (ri >= 0) & (ri < h)
flat = ri[inside] * w + ci[inside]
rural_px = np.bincount(flat, weights=S["rural_ref"][inside].astype(float), minlength=w * h)
all_px = np.bincount(flat, minlength=w * h)
forest = (rural_px / np.maximum(all_px, 1)).reshape(h, w) >= 0.6

out = {"passes": {}, "districts": {}, "forest": {}}
for key, (lst, vt) in grids.items():
    land = np.isfinite(lst) & (dist > 0)
    out["passes"][key] = round(float(np.nanmean(vt[land])), 2)
    out["forest"][key] = round(float(np.nanmean(lst[forest & np.isfinite(lst)])), 2)
    for i, f in enumerate(d_fc, 1):
        m = (dist == i) & np.isfinite(lst)
        name = f["properties"]["name"]
        out["districts"].setdefault(name, {})[key] = round(float(np.nanmean(lst[m])), 2) if m.sum() else None
out["forest_cells"] = int(forest.sum())
json.dump(out, open(C.PROC / "diurnal.json", "w"), indent=1)
print(json.dumps(out["passes"]), "forest cells", out["forest_cells"])
print("forest", out["forest"])
for n in ("Yau Tsim Mong", "Yuen Long", "Sai Kung"):
    print(n, out["districts"][n])
