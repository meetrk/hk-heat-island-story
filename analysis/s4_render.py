"""Step 4 - render the story's visuals.

PNG maps (transparent sea, subtle hillshade texture) into presentation/assets/
and presentation/assets/data.js (results + overlays in image-pixel space).
"""
import json

import numpy as np
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
from shapely.geometry import mapping, shape

import common as C

S = C.load("static"); D = C.load("derived")
R = json.load(open(C.PROC / "results.json"))
EP = {k: C.load(f"composites_{k}") for k in C.EPOCHS}
LAND = S["land"]
OUT_W = 1800                                   # web width; keeps aspect of the grid
scale = OUT_W / C.WIDTH
OUT_H = int(round(C.HEIGHT * scale))

HEAT = colormaps["inferno"]
GREEN = LinearSegmentedColormap.from_list("g", ["#0d1410", "#1d3b2a", "#2f6b45", "#57a66b", "#b6e3a1"])
CONCRETE = LinearSegmentedColormap.from_list("c", ["#0f1720", "#33404d", "#7d8a96", "#d9dee3"])
DIVERGE = LinearSegmentedColormap.from_list("d", ["#3a7bd5", "#9cc3f0", "#20242b", "#f4a261", "#e63946"])

# hillshade texture from the DEM
dem = np.nan_to_num(S["dem"], nan=0.0)
gy, gx = np.gradient(dem, C.RES)
slope = np.pi / 2 - np.arctan(np.hypot(gx, gy))
aspect = np.arctan2(-gx, gy)
az, alt = np.radians(315), np.radians(40)
shade = np.sin(alt) * np.sin(slope) + np.cos(alt) * np.cos(slope) * np.cos(az - aspect)
shade = np.clip((shade - shade.min()) / (np.ptp(shade) + 1e-9), 0, 1)


def to_png(name, rgb, alpha):
    rgba = np.dstack([np.clip(rgb, 0, 1), alpha.astype(float)])
    img = Image.fromarray((rgba * 255).astype("uint8"), "RGBA").resize((OUT_W, OUT_H), Image.LANCZOS)
    img.save(C.ASSETS / f"{name}.png", optimize=True)


def colour(arr, cmap, vmin, vmax, mask, texture=0.18):
    t = np.clip((arr - vmin) / (vmax - vmin), 0, 1)
    rgb = cmap(np.nan_to_num(t))[..., :3]
    rgb = rgb * (1 - texture + texture * shade[..., None] * 1.6)
    return rgb, mask & np.isfinite(arr)


lst24 = EP["2024-25"]["lst"]
lo, hi = np.nanpercentile(lst24[D["valid_2024-25"]], [2, 98])
lo, hi = float(np.floor(lo)), float(np.ceil(hi))
for k in C.EPOCHS:
    to_png(f"lst_{k}", *colour(EP[k]["lst"], HEAT, lo, hi, LAND))
to_png("ndvi_2024", *colour(EP["2024-25"]["ndvi"], GREEN, 0.0, 0.85, LAND))
to_png("ndbi_2024", *colour(EP["2024-25"]["ndbi"], CONCRETE, -0.35, 0.15, LAND))
to_png("lst_tutorial_2024", *colour(D["lst_tut_2024"], HEAT, lo, hi, LAND))

# hot spots: glowing cores only (z > 1), rest a dim ghost of the land
z = D["z_2024-25"]
ghost = np.dstack([shade * 0.18 + 0.05] * 3)
hot = np.clip((z - 1) / 2.0, 0, 1)
rgb = np.where((z > 1)[..., None], HEAT(0.55 + 0.45 * np.nan_to_num(hot))[..., :3], ghost)
to_png("hotspots_2024", rgb, LAND)

# change map, 2014-15 -> 2024-25
ch = D["change"]
base = np.dstack([shade * 0.16 + 0.05] * 3)
pal = {1: (0.25, 0.27, 0.30), 2: (1.0, 0.24, 0.45), 3: (0.30, 0.85, 0.70)}
rgb = base.copy()
for c, col in pal.items():
    rgb[ch == c] = col
to_png("change", rgb, LAND)
to_png("dlst", *colour(D["d_lst"], DIVERGE, -4, 4, LAND))
# Northern Metropolis crop of the change map
nm_rows, nm_cols = np.where(S["nm"])
r0, r1, c0, c1 = nm_rows.min(), nm_rows.max(), nm_cols.min(), nm_cols.max()
crop = Image.fromarray((np.dstack([np.clip(rgb, 0, 1), LAND]) * 255).astype("uint8"), "RGBA") \
    .crop((c0, r0, c1, r1))
crop = crop.resize((1400, int(1400 * (r1 - r0) / (c1 - c0))), Image.NEAREST)
crop.save(C.ASSETS / "change_nm.png", optimize=True)
lst_rgb, _ = colour(EP["2014-15"]["lst"], HEAT, lo, hi, LAND)
for k in ("2014-15", "2024-25"):
    rgbk, _ = colour(EP[k]["lst"], HEAT, lo, hi, LAND)
    im = Image.fromarray((np.dstack([np.clip(rgbk, 0, 1), LAND]) * 255).astype("uint8"), "RGBA") \
        .crop((c0, r0, c1, r1)).resize(crop.size, Image.LANCZOS)
    im.save(C.ASSETS / f"nm_lst_{k}.png", optimize=True)

# MODIS night/day (visual only, own 1 km grid, nearest-neighbour blocks on purpose)
M = np.load(C.PROC / "modis_2024.npz")
from rasterio.features import rasterize as _rz
from affine import Affine
_mt = Affine(*M["transform"])
_hk = _rz([(mapping(C.geom_to_utm(C.load_geojson("hk_boundary.geojson")["features"][0]["geometry"])), 1)],
          out_shape=M["night"].shape, transform=_mt, fill=0, all_touched=True).astype(bool)
for key in ("night", "day"):
    a = np.where(_hk, M[key], np.nan)
    v = np.nanpercentile(a, [3, 97])
    t = np.clip((a - v[0]) / (v[1] - v[0]), 0, 1)
    rgba = HEAT(np.nan_to_num(t)); rgba[..., 3] = np.isfinite(a)
    Image.fromarray((rgba * 255).astype("uint8"), "RGBA").resize(
        (a.shape[1] * 24, a.shape[0] * 24), Image.NEAREST).save(C.ASSETS / f"modis_{key}.png", optimize=True)
R["modis_range"] = {k: [round(float(x), 1) for x in np.nanpercentile(M[k], [3, 97])] for k in ("night", "day")}


# ---- overlays in image-pixel space
def to_px(geom):
    g = C.geom_to_utm(geom).simplify(60)
    inv = ~C.TRANSFORM
    def conv(x, y, z=None):
        c, r = inv * (x, y)
        return c * scale, r * scale
    from shapely.ops import transform
    return transform(conv, g)


def path_d(g):
    polys = [g] if g.geom_type == "Polygon" else list(g.geoms)
    d = ""
    for p in polys:
        for ring in [p.exterior, *p.interiors]:
            pts = list(ring.coords)
            d += "M" + "L".join(f"{x:.1f},{y:.1f}" for x, y in pts) + "Z"
    return d

overlays = {"districts": [], "parks": [], "w": OUT_W, "h": OUT_H}
for f in C.load_geojson("hk_districts.geojson")["features"]:
    g = to_px(f["geometry"])
    c = g.representative_point()
    overlays["districts"].append({"name": f["properties"]["name"], "d": path_d(g), "cx": c.x, "cy": c.y})
hkg = to_px(C.load_geojson("hk_boundary.geojson")["features"][0]["geometry"])
overlays["outline"] = path_d(hkg)
for f in C.load_geojson("hk_parks.geojson")["features"]:
    c = to_px(f["geometry"]).centroid
    overlays["parks"].append({"name": f["properties"]["name"], "x": c.x, "y": c.y})
nmg = to_px(C.NORTHERN_METROPOLIS)
overlays["nm"] = path_d(nmg)
overlays["nm_box"] = [c0 * scale, r0 * scale, c1 * scale, r1 * scale]
# ---- downsampled grids for hover readouts (120 m cells), one binary file
F = 4
gh, gw = C.HEIGHT // F, C.WIDTH // F
def block(a, how="mean"):
    a = a[:gh * F, :gw * F].reshape(gh, F, gw, F)
    with np.errstate(all="ignore"):
        if how == "mean":
            return np.nanmean(a, axis=(1, 3))
        flat = a.transpose(0, 2, 1, 3).reshape(gh, gw, F * F)
        return np.apply_along_axis(lambda v: np.bincount(v, minlength=4).argmax(), 2, flat)
def i16(a):
    return np.where(np.isfinite(a), np.round(a * 100), -32768).astype("<i2")
landf = block(LAND.astype("float32")) >= 0.5
lst14 = np.where(landf, block(EP["2014-15"]["lst"]), np.nan)
lst24 = np.where(landf, block(EP["2024-25"]["lst"]), np.nan)
dl = np.where(landf, block(D["d_lst"]), np.nan)
ndvi = np.where(landf, block(EP["2024-25"]["ndvi"]), np.nan)
ndbi = np.where(landf, block(EP["2024-25"]["ndbi"]), np.nan)
dist = block(S["district"].astype("int64"), "mode").astype("uint8")
chg = block(D["change"].astype("int64"), "mode").astype("uint8")
u8 = lambda a: np.where(np.isfinite(a), np.clip(np.round((a + 1) * 100), 0, 254), 255).astype("uint8")
with open(C.ASSETS / "grid.bin", "wb") as fh:
    for a in (i16(lst14), i16(lst24), i16(dl)):
        fh.write(a.tobytes())
    for a in (u8(ndvi), u8(ndbi), dist, chg):
        fh.write(a.tobytes())
R["grid"] = {"w": int(gw), "h": int(gh), "cell_m": int(C.RES * F),
             "full_w": int(C.WIDTH), "full_h": int(C.HEIGHT),
             "nm_crop": [int(c0), int(r0), int(c1), int(r1)],   # col0,row0,col1,row1 on the 30 m grid
             "districts": [f["properties"]["name"] for f in C.load_geojson("hk_districts.geojson")["features"]]}
R["modis"] = {k: [[None if not np.isfinite(v) else round(float(v), 1) for v in row]
                  for row in np.where(_hk, M[k], np.nan)] for k in ("night", "day")}
_dp = C.PROC / "diurnal.json"
if _dp.exists():
    R["diurnal"] = json.load(open(_dp))
R["lst_scale"] = [lo, hi]
R["overlays"] = overlays
open(C.ASSETS / "data.js", "w").write("window.UHI = " + json.dumps(R, default=float) + ";\n")
print("rendered", OUT_W, OUT_H, "LST scale", lo, hi)
