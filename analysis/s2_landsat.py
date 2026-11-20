"""Step 2 - Landsat 8/9 Collection 2 Level-2, summers 2014-2025.

For every summer scene (May-October) we read, on the common 30 m grid:
  red, nir08, swir16  -> NDVI, NDBI
  lwir11 (ST_B10)     -> LST, USGS surface-temperature product
  trad (ST_TRAD)      -> at-sensor thermal radiance -> brightness temperature
  qa_pixel            -> cloud / shadow mask

Outputs
  composites_<epoch>.npz : May-Oct median LST, NDVI, NDBI, BT for each epoch
  scene_stats.json       : per-scene urban & rural mean LST (for the trend)

The tutorial method (NDVI -> fractional vegetation -> emissivity -> LST) is
applied to the composites in step 3, as a check on the USGS product.
"""
import concurrent.futures as cf
import datetime as dt
import sys

import numpy as np
from rasterio.enums import Resampling

import common as C

S = C.load("static")
LAND, URB, RUR = S["land"], S["urban_ref"], S["rural_ref"]
cat = C.catalog()


def search(y0, y1):
    items = cat.search(collections=["landsat-c2-l2"], bbox=C.BBOX,
                       datetime=f"{y0}-05-01/{y1}-10-31",
                       query={"platform": {"in": ["landsat-8", "landsat-9"]},
                              "eo:cloud_cover": {"lt": 70}}).item_collection()
    return [i for i in items
            if C.SUMMER_MONTHS[0] <= i.datetime.month <= C.SUMMER_MONTHS[1]]


def read_scene(item, full):
    qa, _ = C.read_to_grid(item.assets["qa_pixel"].href, Resampling.nearest, dtype="uint16")
    clear = ((qa & 0b11011) == 0) & (qa != 0)          # fill, dilated cloud, cloud, shadow
    frac = (clear & LAND).sum() / LAND.sum()
    if frac < 0.10:
        return None
    st, _ = C.read_to_grid(item.assets["lwir11"].href, Resampling.nearest)
    lst = st * 0.00341802 + 149.0 - 273.15
    ok = clear & (st > 0) & (lst > 10) & (lst < 70)
    lst[~ok] = np.nan
    out = {"id": item.id, "date": item.datetime.date().isoformat(),
           "platform": item.properties["platform"], "clear_frac": float(frac),
           "urban_clear": float(np.isfinite(lst[URB]).mean()),
           "rural_clear": float(np.isfinite(lst[RUR]).mean()),
           "urban_lst": float(np.nanmean(lst[URB])) if np.isfinite(lst[URB]).any() else None,
           "rural_lst": float(np.nanmean(lst[RUR])) if np.isfinite(lst[RUR]).any() else None}
    if not full:
        return out, None
    bands = {}
    for b in ("red", "nir08", "swir16"):
        a, _ = C.read_to_grid(item.assets[b].href, Resampling.nearest)
        bands[b] = a * 0.0000275 - 0.2
    trad, _ = C.read_to_grid(item.assets["trad"].href, Resampling.nearest)
    L = trad * 0.001                                     # W / (m2 sr um)
    k1, k2 = C.K_CONST[item.properties["platform"]]
    with np.errstate(divide="ignore", invalid="ignore"):
        bt = k2 / np.log(k1 / L + 1.0) - 273.15          # brightness temperature, deg C
        ndvi = (bands["nir08"] - bands["red"]) / (bands["nir08"] + bands["red"])
        ndbi = (bands["swir16"] - bands["nir08"]) / (bands["swir16"] + bands["nir08"])
    stack = np.stack([lst, ndvi, ndbi, bt]).astype("float16")
    stack[:, ~ok] = np.nan
    return out, stack


def run(years, full, workers=8):
    stats, stacks = [], []
    items = search(*years)
    print(years, "items", len(items), flush=True)
    with cf.ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(read_scene, it, full): it for it in items}
        for n, f in enumerate(cf.as_completed(futs), 1):
            try:
                r = f.result()
            except Exception as e:                        # network hiccup - skip scene
                print("  skip", futs[f].id, type(e).__name__, flush=True)
                continue
            if r is None:
                continue
            stats.append(r[0])
            if r[1] is not None:
                stacks.append(r[1])
            if n % 10 == 0:
                print(f"  {n}/{len(items)} kept {len(stats)}", flush=True)
    return stats, stacks


if __name__ == "__main__":
    import json
    stats_path = C.PROC / "scene_stats.json"
    all_stats = json.load(open(stats_path)) if stats_path.exists() else []
    done_ids = {s["id"] for s in all_stats}
    epoch_years = {y for ys in C.EPOCHS.values() for y in ys}
    which = sys.argv[1] if len(sys.argv) > 1 else "all"

    if which in ("all", "epochs"):
        for name, (y0, y1) in C.EPOCHS.items():
            if (C.PROC / f"composites_{name}.npz").exists():
                continue
            stats, stacks = run((y0, y1), full=True)
            all_stats += [s for s in stats if s["id"] not in done_ids]
            done_ids |= {s["id"] for s in stats}
            arr = np.stack(stacks).astype("float32")      # scenes x 4 x H x W
            with np.errstate(all="ignore"):
                med = np.nanmedian(arr, axis=0)
                count = np.isfinite(arr[:, 0]).sum(axis=0).astype("int16")
            del arr, stacks
            C.save(f"composites_{name}", lst=med[0], ndvi=med[1], ndbi=med[2], bt=med[3], count=count,
                   n_scenes=len(stats))
            C.write_json("scene_stats.json", sorted(all_stats, key=lambda s: s["date"]))
            print(name, "scenes used", len(stats), flush=True)

    if which in ("all", "trend"):
        tpath = C.PROC / "scene_stats_trend.json"
        all_stats = json.load(open(tpath)) if tpath.exists() else []
        done_ids = {s["id"] for s in all_stats}
        for y in range(2014, 2026):
            if y in epoch_years:
                continue
            stats, _ = run((y, y), full=False)
            all_stats += [s for s in stats if s["id"] not in done_ids]
            done_ids |= {s["id"] for s in stats}
            C.write_json("scene_stats_trend.json", sorted(all_stats, key=lambda s: s["date"]))
    print("done", len(all_stats), "scene records")
