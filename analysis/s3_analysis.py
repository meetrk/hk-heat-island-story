"""Step 3 - the analysis behind the three research questions.

Everything here works on the common 30 m grid (see common.py).

  Q1  UHI intensity, hot-spot extent, UTFVI per epoch + 2014-2025 trend
  Q2  LST vs NDVI / NDBI: stratified sample, correlation, multiple regression,
      park cooling (park interior vs 300 m street ring) + cooling-distance curve
  Q3  Northern Metropolis: pixels newly developed between 2014-15 and 2024-25
      vs pixels that stayed unchanged; difference in their LST change

Also re-derives LST with the tutorial method (NDVI -> fractional vegetation ->
emissivity -> LST from brightness temperature) and checks it against the USGS
ST_B10 product.

Writes data/processed/results.json and data/processed/derived.npz.
"""
import json

import numpy as np
from scipy import ndimage, stats

import common as C

rng = np.random.default_rng(9083)
S = C.load("static")
LAND, URB, RUR, NM, DIST = S["land"], S["urban_ref"], S["rural_ref"], S["nm"], S["district"]
BUILT, TREE = S["share_built"], S["share_tree"]
PX_KM2 = (C.RES / 1000) ** 2
DNAMES = [f["properties"]["name"] for f in C.load_geojson("hk_districts.geojson")["features"]]
EP = {k: C.load(f"composites_{k}") for k in C.EPOCHS}
R = {"grid": {"crs": C.CRS, "res_m": C.RES, "shape": list(C.SHAPE),
              "land_km2": float(LAND.sum() * PX_KM2)}}


def r2(x):
    return None if x is None or not np.isfinite(x) else round(float(x), 3)


def valid_mask(e):
    return LAND & np.isfinite(e["lst"]) & np.isfinite(e["ndvi"]) & (e["count"] >= 3)


# ---------------------------------------------------------------- tutorial LST
def tutorial_lst(e, v):
    """Avdan & Jovanovska (2016), as used in the Earth Engine LST tutorials."""
    nmin, nmax = np.nanpercentile(e["ndvi"][v], [1, 99])
    fv = np.clip((e["ndvi"] - nmin) / (nmax - nmin), 0, 1) ** 2
    em = 0.004 * fv + 0.986
    btk = e["bt"] + 273.15
    lst = btk / (1 + (10.895e-6 * btk / 1.4388e-2) * np.log(em)) - 273.15
    return lst, float(nmin), float(nmax)


derived = {}
R["epochs"] = {}
for name, e in EP.items():
    v = valid_mask(e)
    lst = e["lst"]
    mu, sd = float(np.nanmean(lst[v])), float(np.nanstd(lst[v]))
    z = (lst - mu) / sd
    utfvi = ((lst + 273.15) - (mu + 273.15)) / (lst + 273.15)
    edges = [-np.inf, 0, 0.005, 0.010, 0.015, 0.020, np.inf]
    labels = ["excellent", "good", "normal", "bad", "worse", "worst"]
    ut_share = {lab: float(((utfvi[v] >= lo) & (utfvi[v] < hi)).mean())
                for lab, lo, hi in zip(labels, edges[:-1], edges[1:])}
    lst_t, nmin, nmax = tutorial_lst(e, v)
    both = v & np.isfinite(lst_t)
    urban = float(np.nanmean(lst[URB & v])); rural = float(np.nanmean(lst[RUR & v]))
    core = np.isin(DIST, [DNAMES.index("Yau Tsim Mong") + 1, DNAMES.index("Sham Shui Po") + 1]) & v
    districts = []
    for i, dn in enumerate(DNAMES, 1):
        m = (DIST == i) & v
        if m.sum() > 50:
            mb = m & (BUILT >= 0.5)                      # where people live and work
            districts.append({"name": dn, "lst": r2(np.nanmean(lst[m])), "ndvi": r2(np.nanmean(e["ndvi"][m])),
                              "built": r2(BUILT[m].mean()), "land_km2": r2(((DIST == i) & LAND).sum() * PX_KM2),
                              "built_km2": r2(mb.sum() * PX_KM2),
                              "lst_built": r2(np.nanmean(lst[mb])) if mb.sum() > 50 else None,
                              "hotspot_km2": r2(((z > 1) & m).sum() * PX_KM2)})
    R["epochs"][name] = {
        "n_scenes": int(e["n_scenes"]), "valid_km2": r2(v.sum() * PX_KM2),
        "lst_mean": r2(mu), "lst_sd": r2(sd),
        "urban_lst": r2(urban), "rural_lst": r2(rural), "uhi_intensity": r2(urban - rural),
        "core_lst": r2(np.nanmean(lst[core])), "core_minus_rural": r2(np.nanmean(lst[core]) - rural),
        "hotspot_km2": r2(((z > 1) & v).sum() * PX_KM2),
        "severe_hotspot_km2": r2(((z > 2) & v).sum() * PX_KM2),
        "utfvi_share": {k: r2(x) for k, x in ut_share.items()},
        "tutorial_check": {"ndvi_min": r2(nmin), "ndvi_max": r2(nmax),
                           "r_vs_usgs": r2(np.corrcoef(lst[both], lst_t[both])[0, 1]),
                           "mean_diff": r2(np.nanmean(lst_t[both] - lst[both]))},
        "districts": sorted(districts, key=lambda d: -d["lst"]),
    }
    derived[f"z_{name}"] = np.where(v, z, np.nan).astype("float32")
    derived[f"valid_{name}"] = v

# ------------------------------------------------------------------- Q1 trend
def scene_series():
    recs = json.load(open(C.PROC / "scene_stats.json"))
    p = C.PROC / "scene_stats_trend.json"
    if p.exists():
        recs += json.load(open(p))
    seen, out = set(), []
    for s in recs:
        if s["id"] in seen or s["urban_lst"] is None or s["rural_lst"] is None:
            continue
        seen.add(s["id"])
        if s["urban_clear"] >= 0.6 and s["rural_clear"] >= 0.6:
            out.append(s)
    return sorted(out, key=lambda s: s["date"])

scenes = scene_series()
years = {}
for s in scenes:
    years.setdefault(int(s["date"][:4]), []).append(s)
yearly = []
for y in sorted(years):
    u = np.array([s["urban_lst"] for s in years[y]]); r = np.array([s["rural_lst"] for s in years[y]])
    yearly.append({"year": y, "n": len(u), "uhi": r2(np.median(u - r)),
                   "urban": r2(np.median(u)), "rural": r2(np.median(r))})
yy = np.array([d["year"] for d in yearly]); uh = np.array([d["uhi"] for d in yearly])
fit = stats.linregress(yy, uh) if len(yy) >= 4 else None
R["q1_trend"] = {
    "scenes_used": len(scenes), "yearly": yearly,
    "slope_per_decade": r2(fit.slope * 10) if fit else None,
    "p_value": r2(fit.pvalue) if fit else None,
    "scenes": [{"date": s["date"], "uhi": r2(s["urban_lst"] - s["rural_lst"])} for s in scenes],
}

# -------------------------------------------------------------- Q2 relationship
e = EP["2024-25"]; v = derived["valid_2024-25"]
bins = np.digitize(BUILT, [0.2, 0.4, 0.6, 0.8])
idx = []
for b in range(5):
    cand = np.flatnonzero((v & (bins == b)).ravel())
    idx.append(rng.choice(cand, size=min(400, cand.size), replace=False))
idx = np.concatenate(idx)
L_ = e["lst"].ravel()[idx]; N_ = e["ndvi"].ravel()[idx]; B_ = e["ndbi"].ravel()[idx]
ok = np.isfinite(L_) & np.isfinite(N_) & np.isfinite(B_)
L_, N_, B_ = L_[ok], N_[ok], B_[ok]
X = np.column_stack([np.ones_like(N_), N_, B_])
coef, *_ = np.linalg.lstsq(X, L_, rcond=None)
pred = X @ coef
r2_ols = 1 - ((L_ - pred) ** 2).sum() / ((L_ - L_.mean()) ** 2).sum()
zs = lambda a: (a - a.mean()) / a.std()
std_coef, *_ = np.linalg.lstsq(np.column_stack([zs(N_), zs(B_)]), zs(L_), rcond=None)
allv = v & np.isfinite(e["ndbi"])
ndvi_bins = np.arange(-0.1, 0.95, 0.05)
curve = []
for lo in ndvi_bins:
    m = allv & (e["ndvi"] >= lo) & (e["ndvi"] < lo + 0.05)
    if m.sum() > 200:
        curve.append({"ndvi": r2(lo + 0.025), "lst": r2(np.nanmean(e["lst"][m])), "n": int(m.sum())})
R["q2"] = {
    "sample_n": int(L_.size),
    "r_ndvi": r2(stats.pearsonr(N_, L_)[0]), "r_ndbi": r2(stats.pearsonr(B_, L_)[0]),
    "r_ndvi_all_pixels": r2(np.corrcoef(e["ndvi"][allv], e["lst"][allv])[0, 1]),
    "r_ndbi_all_pixels": r2(np.corrcoef(e["ndbi"][allv], e["lst"][allv])[0, 1]),
    "ols": {"intercept": r2(coef[0]), "b_ndvi": r2(coef[1]), "b_ndbi": r2(coef[2]), "r2": r2(r2_ols),
            "std_b_ndvi": r2(std_coef[0]), "std_b_ndbi": r2(std_coef[1])},
    "ndvi_curve": curve,
    "sample": [[r2(a), r2(b), r2(c)] for a, b, c in zip(N_, B_, L_)],
}
# also check the same regression for 2014-15 (robustness)
e0 = EP["2014-15"]; v0 = derived["valid_2014-15"] & np.isfinite(e0["ndbi"])
R["q2"]["r_ndvi_all_pixels_2014"] = r2(np.corrcoef(e0["ndvi"][v0], e0["lst"][v0])[0, 1])

# park cooling
parks = C.load_geojson("hk_parks.geojson")["features"]
park_rows, profiles = [], []
dbins = np.arange(0, 630, 30)
for f in parks:
    pm = C.rasterize_geoms([f["geometry"]]).astype(bool) & v
    if pm.sum() < 15:
        continue
    dist = ndimage.distance_transform_edt(~pm) * C.RES
    ring = (dist > 0) & (dist <= 300) & v & (S["share_water"] < 0.2) & (BUILT >= 0.3)
    if ring.sum() < 30:
        continue
    pin = float(np.nanmean(e["lst"][pm])); pring = float(np.nanmean(e["lst"][ring]))
    park_rows.append({"name": f["properties"]["name"], "area_ha": r2(pm.sum() * 0.09),
                      "park_lst": r2(pin), "ring_lst": r2(pring), "cooling": r2(pring - pin),
                      "ndvi": r2(np.nanmean(e["ndvi"][pm]))})
    prof = []
    for lo in dbins:
        m = (dist > lo) & (dist <= lo + 30) & v & (S["share_water"] < 0.2) if lo > 0 else pm
        prof.append(np.nanmean(e["lst"][m]) - pin if m.sum() > 3 else np.nan)
    profiles.append(prof)
prof = np.nanmean(np.array(profiles), axis=0)
R["q2"]["parks"] = sorted(park_rows, key=lambda p: -p["cooling"])
R["q2"]["park_mean_cooling"] = r2(np.mean([p["cooling"] for p in park_rows]))
R["q2"]["cooling_profile"] = [{"d": int(d), "dT": r2(t)} for d, t in zip(dbins, prof)]

# ------------------------------------------------------------- Q3 development
a, b = EP["2014-15"], EP["2024-25"]
both = derived["valid_2014-15"] & derived["valid_2024-25"] & np.isfinite(a["ndbi"]) & np.isfinite(b["ndbi"])
d_lst, d_ndvi, d_ndbi = b["lst"] - a["lst"], b["ndvi"] - a["ndvi"], b["ndbi"] - a["ndbi"]
developed = both & (d_ndbi > 0.08) & (d_ndvi < -0.10)
greened = both & (d_ndvi > 0.10) & (d_ndbi < -0.05)
stable = both & (np.abs(d_ndbi) < 0.03) & (np.abs(d_ndvi) < 0.05)
io_new = (S["io_built_2023"] - S["io_built_2017"]) > 0.3


def zone_stats(zone, label):
    dev, stb, grn = developed & zone, stable & zone, greened & zone
    out = {"zone": label, "developed_km2": r2(dev.sum() * PX_KM2), "stable_km2": r2(stb.sum() * PX_KM2),
           "greened_km2": r2(grn.sum() * PX_KM2),
           "dlst_developed": r2(np.nanmean(d_lst[dev])), "dlst_stable": r2(np.nanmean(d_lst[stb])),
           "dlst_greened": r2(np.nanmean(d_lst[grn])) if grn.sum() > 30 else None}
    out["development_warming"] = r2(out["dlst_developed"] - out["dlst_stable"])
    if out["dlst_greened"] is not None:
        out["greening_cooling"] = r2(out["dlst_greened"] - out["dlst_stable"])
    out["io_agreement"] = r2((io_new[dev] | (S["io_built_2023"][dev] > 0.5)).mean()) if dev.sum() else None
    # spatially-blocked bootstrap (1.5 km blocks) for an honest interval
    blk = (np.arange(C.HEIGHT)[:, None] // 50) * 1000 + (np.arange(C.WIDTH)[None, :] // 50)
    db, sb = blk[dev], blk[stb]
    ub = np.unique(np.concatenate([db, sb]))
    dvals, svals = d_lst[dev], d_lst[stb]
    boots = []
    for _ in range(300):
        pick = rng.choice(ub, size=ub.size, replace=True)
        w = np.bincount(np.searchsorted(ub, pick), minlength=ub.size)
        wd, ws = w[np.searchsorted(ub, db)], w[np.searchsorted(ub, sb)]
        if wd.sum() and ws.sum():
            boots.append((dvals * wd).sum() / wd.sum() - (svals * ws).sum() / ws.sum())
    out["ci95"] = [r2(np.percentile(boots, 2.5)), r2(np.percentile(boots, 97.5))]
    return out

R["q3"] = {"northern_metropolis": zone_stats(NM, "Northern Metropolis"),
           "hong_kong": zone_stats(LAND, "All Hong Kong")}
dist_change = []
for i, dn in enumerate(DNAMES, 1):
    m = both & (DIST == i)
    if m.sum() > 200:
        dist_change.append({"name": dn, "dlst": r2(np.nanmean(d_lst[m])),
                            "developed_km2": r2((developed & m).sum() * PX_KM2)})
R["q3"]["districts"] = sorted(dist_change, key=lambda d: -d["developed_km2"])

# Q1 (where): per-district warming 2014-15 -> 2024-25. 2024-25 summers were hotter
# everywhere, so the headline is warming relative to the Hong Kong-wide mean change.
hk_mean = float(np.nanmean(d_lst[both]))
hk_built = both & (BUILT >= 0.5)
hk_mean_built = float(np.nanmean(d_lst[hk_built]))
warm = []
for i, dn in enumerate(DNAMES, 1):
    m = both & (DIST == i); mb = m & (BUILT >= 0.5)
    if m.sum() > 200:
        warm.append({"name": dn, "dlst": r2(np.nanmean(d_lst[m])),
                     "rel": r2(np.nanmean(d_lst[m]) - hk_mean),
                     "dlst_built": r2(np.nanmean(d_lst[mb])) if mb.sum() > 50 else None,
                     "rel_built": r2(np.nanmean(d_lst[mb]) - hk_mean_built) if mb.sum() > 50 else None,
                     "developed_km2": r2((developed & m).sum() * PX_KM2),
                     "greened_km2": r2((greened & m).sum() * PX_KM2),
                     "land_km2": r2(((DIST == i) & LAND).sum() * PX_KM2)})
R["q1_districts"] = {"hk_mean_dlst": r2(hk_mean), "hk_mean_dlst_built": r2(hk_mean_built),
                     "districts": sorted(warm, key=lambda d: -d["rel"])}

change = np.zeros(C.SHAPE, "uint8")
change[stable] = 1; change[developed] = 2; change[greened] = 3
derived.update(d_lst=np.where(both, d_lst, np.nan).astype("float32"), change=change,
               lst_tut_2024=tutorial_lst(b, derived["valid_2024-25"])[0].astype("float32"))
C.save("derived", **derived)
C.write_json("results.json", R)
print(json.dumps({k: R[k] for k in ("q1_trend",)}, default=float)[:600])
for k, ep in R["epochs"].items():
    print(k, {x: ep[x] for x in ("n_scenes", "uhi_intensity", "core_minus_rural", "hotspot_km2", "tutorial_check")})
print("Q2", {k: R["q2"][k] for k in ("r_ndvi", "r_ndbi", "r_ndvi_all_pixels", "ols", "park_mean_cooling")})
print("Q3", R["q3"]["northern_metropolis"])
