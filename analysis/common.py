"""Shared configuration, the common 30 m grid, and raster helpers.

Every quantitative layer in this project lives on ONE grid: 30 m pixels in
UTM zone 50N (EPSG:32650), covering Hong Kong. This is the fix the course
feedback asked for - LST, NDVI, NDBI and the land-cover shares are compared
pixel by pixel only because they share the same resolution and alignment.
"""
from __future__ import annotations

import json
import os
import pathlib

import numpy as np
import planetary_computer as pc
import pystac_client
import rasterio
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from shapely.geometry import shape, mapping
from shapely.ops import transform as shp_transform

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
ASSETS = ROOT / "presentation" / "assets"
for p in (RAW, PROC, ASSETS):
    p.mkdir(parents=True, exist_ok=True)

os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "4")
os.environ.setdefault("GDAL_HTTP_RETRY_DELAY", "2")
os.environ.setdefault("GDAL_HTTP_TIMEOUT", "40")
os.environ.setdefault("GDAL_HTTP_CONNECTTIMEOUT", "15")
os.environ.setdefault("GDAL_HTTP_MULTIPLEX", "YES")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif,.TIF,.tiff")

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
BBOX = [113.82, 22.14, 114.45, 22.57]          # lon/lat, Hong Kong SAR
CRS = "EPSG:32650"
RES = 30.0

# Three summer epochs, two summers each (Hong Kong's wet season is cloudy)
EPOCHS = {"2014-15": (2014, 2015), "2019-20": (2019, 2020), "2024-25": (2024, 2025)}
SUMMER_MONTHS = (5, 10)                        # May-October

# Thermal constants for band 10 (from the Landsat MTL files)
K_CONST = {"landsat-8": (774.8853, 1321.0789), "landsat-9": (799.0284, 1329.2405)}

# Rough outline of the Northern Metropolis development zone
# (Yuen Long, San Tin, Kwu Tung, Fanling North). Clipped to HK land later.
NORTHERN_METROPOLIS = {"type": "Polygon", "coordinates": [[
    [113.99, 22.425], [114.165, 22.425], [114.165, 22.545], [113.99, 22.545], [113.99, 22.425]]]}


def _grid():
    t = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    xs, ys = t.transform([BBOX[0], BBOX[2], BBOX[0], BBOX[2]], [BBOX[1], BBOX[1], BBOX[3], BBOX[3]])
    x0, x1 = np.floor(min(xs) / RES) * RES, np.ceil(max(xs) / RES) * RES
    y0, y1 = np.floor(min(ys) / RES) * RES, np.ceil(max(ys) / RES) * RES
    w, h = int((x1 - x0) / RES), int((y1 - y0) / RES)
    return from_origin(x0, y1, RES, RES), w, h


TRANSFORM, WIDTH, HEIGHT = _grid()
SHAPE = (HEIGHT, WIDTH)
TO_UTM = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)


def catalog():
    return pystac_client.Client.open(STAC_URL)


def read_to_grid(href: str, resampling=Resampling.nearest, transform=None, width=None, height=None,
                 dtype="float32", band=1):
    """Read any COG, warped onto the common grid (or a finer aligned one)."""
    transform = transform or TRANSFORM
    width, height = width or WIDTH, height or HEIGHT
    with rasterio.open(pc.sign(href)) as src:
        with WarpedVRT(src, crs=CRS, transform=transform, width=width, height=height,
                       resampling=resampling) as vrt:
            arr = vrt.read(band, out_dtype=dtype)
            nodata = src.nodata
    return arr, nodata


def geom_to_utm(geom):
    return shp_transform(lambda x, y, z=None: TO_UTM.transform(x, y), shape(geom))


def rasterize_geoms(geoms, values=None, fill=0, dtype="int16"):
    items = [(mapping(geom_to_utm(g)), v) for g, v in zip(geoms, values or [1] * len(geoms))]
    return rasterize(items, out_shape=SHAPE, transform=TRANSFORM, fill=fill, dtype=dtype)


def load_geojson(name):
    return json.load(open(RAW / name))


def save(name, **arrays):
    np.savez_compressed(PROC / f"{name}.npz", **arrays)


def load(name):
    return dict(np.load(PROC / f"{name}.npz"))


def write_json(name, obj):
    json.dump(obj, open(PROC / name, "w"), indent=1, default=float)
