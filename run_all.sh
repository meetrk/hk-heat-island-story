#!/usr/bin/env bash
# Full pipeline. Boundaries/parks are fetched once; delete data/processed to recompute.
set -euo pipefail
cd "$(dirname "$0")/analysis"
PY=../.venv/bin/python
[ -f ../data/raw/hk_districts.geojson ] || $PY fetch_boundaries.py
[ -f ../data/raw/hk_parks.geojson ] || $PY fetch_parks.py
$PY s1_landcover.py
$PY s2_landsat.py epochs
$PY s2_landsat.py trend
$PY s5_modis_night.py
$PY s6_diurnal.py
$PY s3_analysis.py
$PY s4_render.py
echo "Open ../presentation/index.html"
