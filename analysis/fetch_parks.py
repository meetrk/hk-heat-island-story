"""Fetch outlines of large urban parks (OSM / Nominatim) for the park-cooling test."""
import json, time, urllib.request, urllib.parse
import common as C
UA = {"User-Agent": "uhi-coursework-script/1.0"}
PARKS = ["Kowloon Park", "Victoria Park", "Hong Kong Park", "Kowloon Walled City Park",
         "Lai Chi Kok Park", "Morse Park, Wong Tai Sin", "Sha Tin Park", "Tuen Mun Park", "Yuen Long Park",
         "Kwun Tong Promenade", "Hong Kong Zoological and Botanical Gardens, Central", "Shek Kip Mei Park",
         "Tsing Yi Park", "Ma On Shan Park", "Tai Po Waterfront Park", "Sham Shui Po Park"]
feats = []
for p in PARKS:
    url = "https://nominatim.openstreetmap.org/search?" + urllib.parse.urlencode(
        {"q": p + ", Hong Kong", "format": "geojson", "polygon_geojson": 1, "limit": 5})
    fc = json.load(urllib.request.urlopen(urllib.request.Request(url, headers=UA)))
    time.sleep(1.2)
    poly = [f for f in fc["features"] if f["geometry"]["type"] in ("Polygon", "MultiPolygon")]
    if poly:
        poly[0]["properties"] = {"name": p.split(",")[0]}
        feats.append(poly[0]); print("ok", p)
    else:
        print("--", p)
json.dump({"type": "FeatureCollection", "features": feats}, open(C.RAW / "hk_parks.geojson", "w"))
