"""Fetch Hong Kong territory + 18 district polygons from OpenStreetMap (Nominatim).
Run once; writes data/raw/hk_boundary.geojson and data/raw/hk_districts.geojson."""
import json, time, urllib.request, urllib.parse, pathlib
RAW = pathlib.Path(__file__).resolve().parents[1] / "data" / "raw"
UA = {"User-Agent": "uhi-coursework-script/1.0"}
DISTRICTS = ["Central and Western", "Wan Chai", "Eastern", "Southern", "Yau Tsim Mong",
             "Sham Shui Po", "Kowloon City", "Wong Tai Sin", "Kwun Tong", "Kwai Tsing",
             "Tsuen Wan", "Tuen Mun", "Yuen Long", "North", "Tai Po", "Sha Tin", "Sai Kung", "Islands"]

def nominatim(q):
    url = "https://nominatim.openstreetmap.org/search?" + urllib.parse.urlencode(
        {"q": q, "format": "geojson", "polygon_geojson": 1, "polygon_threshold": 0.0005, "limit": 5})
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA)) as r:
        fc = json.load(r)
    time.sleep(1.2)
    for f in fc["features"]:
        if f["geometry"]["type"] in ("Polygon", "MultiPolygon") and f["properties"].get("category") == "boundary":
            return f
    raise RuntimeError("no polygon for " + q)

hk = nominatim("Hong Kong")
json.dump({"type": "FeatureCollection", "features": [hk]}, open(RAW / "hk_boundary.geojson", "w"))
feats = []
for d in DISTRICTS:
    f = nominatim(f"{d} District, Hong Kong")
    f["properties"] = {"name": d}
    feats.append(f); print("ok", d)
json.dump({"type": "FeatureCollection", "features": feats}, open(RAW / "hk_districts.geojson", "w"))
