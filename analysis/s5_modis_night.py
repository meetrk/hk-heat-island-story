"""Supplementary, VISUAL ONLY: MODIS MOD11A2 night-time LST (1 km), summer 2024.

MODIS is 1 km, Landsat 30 m. Following the course feedback, the two are never
combined in any statistic: this layer only feeds a side-by-side picture of the
night-time heat pattern. It is kept on its own 1 km grid on purpose.
"""
import numpy as np
from rasterio.enums import Resampling
from rasterio.transform import from_origin

import common as C

KM = 1000.0
t = from_origin(np.floor(C.TRANSFORM.c / KM) * KM, np.ceil(C.TRANSFORM.f / KM) * KM, KM, KM)
w, h = int(np.ceil(C.WIDTH * C.RES / KM)) + 1, int(np.ceil(C.HEIGHT * C.RES / KM)) + 1
cat = C.catalog()
items = list(cat.search(collections=["modis-11A2-061"], bbox=C.BBOX,
                        datetime="2024-06-01/2024-09-30").items())
items = [i for i in items if i.id.startswith("MOD11A2")]   # Terra only (~10:30 / ~22:30)
night, day = [], []
for it in items:
    for key, store in (("LST_Night_1km", night), ("LST_Day_1km", day)):
        a, _ = C.read_to_grid(it.assets[key].href, Resampling.nearest, t, w, h, "float32")
        a[a == 0] = np.nan
        store.append(a * 0.02 - 273.15)
with np.errstate(all="ignore"):
    n = np.nanmean(np.stack(night), axis=0)
    d = np.nanmean(np.stack(day), axis=0)
np.savez_compressed(C.PROC / "modis_2024.npz", night=n, day=d,
                    transform=np.array(t)[:6], n_items=len(items))
print("MODIS composites", len(items), "items; night range", np.nanpercentile(n, [2, 98]))
