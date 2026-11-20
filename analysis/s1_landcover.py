"""Step 1 - static layers on the 30 m grid.

* ESA WorldCover 2021 (10 m) -> share of each 30 m pixel that is built-up,
  tree, water ... (block-mean of 3x3 sub-pixels; never used at 10 m)
* Copernicus DEM 30 m -> elevation (rural reference must sit low, otherwise
  "rural is cooler" partly means "mountains are cooler")
* Impact Observatory annual land cover 2017 & 2023 (10 m -> 30 m shares),
  an independent check on the Landsat change map in step 4
* Hong Kong boundary, 18 districts, focus zones
"""
import numpy as np
from rasterio.enums import Resampling
from rasterio.transform import from_origin

import common as C

cat = C.catalog()
fine_t = from_origin(C.TRANSFORM.c, C.TRANSFORM.f, 10, 10)
FW, FH = C.WIDTH * 3, C.HEIGHT * 3


def block_share(fine, classes):
    m = np.isin(fine, classes).astype("float32")
    return m.reshape(C.HEIGHT, 3, C.WIDTH, 3).mean(axis=(1, 3))


def mosaic_fine(items, asset):
    out = np.zeros((FH, FW), "uint8")
    for it in items:
        arr, _ = C.read_to_grid(it.assets[asset].href, Resampling.nearest, fine_t, FW, FH, "uint8")
        out = np.where(out == 0, arr, out)
    return out

# --- WorldCover 2021 (v200)
wc_items = [i for i in cat.search(collections=["esa-worldcover"], bbox=C.BBOX).items()
            if "2021" in i.id]
wc = mosaic_fine(wc_items, "map")
shares = {
    "built": block_share(wc, [50]), "tree": block_share(wc, [10]),
    "grass": block_share(wc, [20, 30]), "crop": block_share(wc, [40]),
    "water": block_share(wc, [80]), "wetland": block_share(wc, [90, 95]),
    "bare": block_share(wc, [60]),
}
print("worldcover tiles", len(wc_items))

# --- Impact Observatory annual LULC (class 7 = built area)
io = {}
for year in (2017, 2023):
    items = list(cat.search(collections=["io-lulc-annual-v02"], bbox=C.BBOX,
                            datetime=f"{year}-01-01/{year}-12-31").items())
    io[year] = block_share(mosaic_fine(items, "data"), [7])
    print("io-lulc", year, len(items))

# --- DEM
dem = np.full(C.SHAPE, np.nan, "float32")
for it in cat.search(collections=["cop-dem-glo-30"], bbox=C.BBOX).items():
    arr, nd = C.read_to_grid(it.assets["data"].href, Resampling.bilinear)
    arr[arr == (nd if nd is not None else -32767)] = np.nan
    arr[arr < -100] = np.nan
    dem = np.where(np.isnan(dem), arr, dem)

# --- Boundaries & zones
hk_geom = C.load_geojson("hk_boundary.geojson")["features"][0]["geometry"]
hk = C.rasterize_geoms([hk_geom]).astype(bool)
d_fc = C.load_geojson("hk_districts.geojson")["features"]
district = C.rasterize_geoms([f["geometry"] for f in d_fc], list(range(1, len(d_fc) + 1)))
nm = C.rasterize_geoms([C.NORTHERN_METROPOLIS]).astype(bool)

land = hk & (shares["water"] < 0.5)
urban_ref = land & (shares["built"] >= 0.8)
rural_ref = land & (shares["tree"] >= 0.8) & (dem < 200)

C.save("static", hk=hk, land=land, district=district, nm=nm & land, dem=dem,
       urban_ref=urban_ref, rural_ref=rural_ref,
       io_built_2017=io[2017], io_built_2023=io[2023], **{f"share_{k}": v for k, v in shares.items()})
print("grid", C.SHAPE, "land px", land.sum(), "urban ref", urban_ref.sum(), "rural ref", rural_ref.sum(),
      "NM px", (nm & land).sum())
