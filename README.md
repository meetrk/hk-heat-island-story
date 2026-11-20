# Hotter Streets, Greener Relief?

**Hong Kong's urban heat island, 2014–2025, mapped from Landsat 8/9.**
Group project for HKU CCST9083 *Earth as Seen by Satellite*.

**▶ Live story deck: https://hk-heat-island-story.vercel.app** (arrow keys to move, `F` for fullscreen, hover any map or chart)

## Research questions

1. How have the intensity and spatial extent of Hong Kong's daytime surface urban heat island changed between 2014 and 2025?
2. How strongly is land surface temperature (LST) related to greenness (NDVI) and built-up intensity (NDBI), and how much cooler are parks than the streets around them?
3. Has land conversion in the Northern Metropolis produced measurable surface warming compared with land that stayed unchanged?

## Findings

| | Result |
|---|---|
| Heat island (built-up minus lowland forest) | 6.2 °C (2014–15) → 6.5 °C (2019–20) → 6.8 °C (2024–25). The 12-year scene-level trend is +0.56 °C/decade but **not significant** (p = 0.37), so the gap held rather than clearly grew |
| Hot-spot area (LST > city mean + 1 SD) | 173 → 181 km² |
| LST vs greenness | r = −0.56 in a 2,000-pixel stratified sample (−0.73 over all land pixels) |
| LST vs concrete | r = +0.59. In a joint regression (R² = 0.35), concrete is the stronger driver |
| Park cooling | parks are 1.65 °C cooler than the streets within 300 m, and the effect fades within ~100 m |
| Districts | on built-up ground only, Yuen Long (43.1 °C), Islands (airport) and Kwai Tsing (port) run hottest. Yuen Long also warmed most: +0.9 °C beyond the Hong Kong-wide change |
| Northern Metropolis | 10.6 km² turned from green to built. That land warmed **+2.1 °C more** than unchanged land (95 % CI 1.6–2.8, spatial block bootstrap) |

## Method

Every quantitative layer sits on **one 30 m grid** (UTM zone 50N), so LST, NDVI and NDBI are compared pixel by pixel:

- **Landsat 8/9 Collection 2 Level-2**, May–October median composites for 2014–15, 2019–20 and 2024–25 (50, 34 and 81 scenes), plus per-scene statistics for the 2014–2025 trend. LST comes from `ST_B10`. NDVI = (NIR − Red)/(NIR + Red), NDBI = (SWIR1 − NIR)/(SWIR1 + NIR).
- **Tutorial LST check:** NDVI → fractional vegetation → emissivity → LST from band-10 brightness temperature (Avdan & Jovanovska, 2016). It reproduces the USGS pattern (r = 0.94) but reads ~11 °C lower, because it starts from at-sensor radiance with no atmospheric correction.
- **Land cover:** ESA WorldCover and Impact Observatory (10 m) are *averaged up* to 30 m shares, never used at 10 m. Urban reference = ≥ 80 % built-up. Rural reference = ≥ 80 % tree cover below 200 m elevation (Copernicus DEM), so "rural is cooler" does not just mean "mountains are cooler".
- **MODIS (1 km)** is used only for the day/night pictures and the 24-hour slide (four passes a day: Terra ≈10:00/21:30, Aqua ≈14:00/02:30). It never enters the 30 m statistics.

**Limitations:** satellite LST is surface, not air, temperature. Landsat observes only at ~10:30. The thermal band is natively 100 m. The Northern Metropolis outline is an approximation, not the official planning boundary.

## Repository layout

```
analysis/        Python pipeline (data from Microsoft Planetary Computer, no login needed)
  common.py        config + the common 30 m grid
  fetch_*.py       one-off OpenStreetMap downloads (boundary, districts, parks)
  s1_landcover.py  land-cover shares, DEM, zones
  s2_landsat.py    Landsat composites (epochs) and per-scene trend stats
  s3_analysis.py   statistics for all three questions -> data/processed/results.json
  s4_render.py     maps + data for the deck
  s5_modis_night.py, s6_diurnal.py   MODIS day/night and 24-hour cycle (visual only)
gee/uhi_hongkong.js  the same method as a Google Earth Engine Code Editor script
data/raw/        OSM boundaries and parks
data/processed/  results.json, diurnal.json, scene stats (large rasters are not committed)
presentation/    the story deck (static HTML + D3), deployed on Vercel
run_all.sh       runs the whole pipeline
```

## Run it

```bash
python3 -m venv .venv
.venv/bin/pip install pystac-client planetary-computer rasterio pyproj shapely numpy scipy matplotlib pillow
./run_all.sh            # ~30-60 min, mostly downloading Landsat scenes
cd presentation && python3 -m http.server 8765   # then open http://localhost:8765
```

## Team

Feng Yujie · Meet Kachhadiya · Tony Tse · Song Yixuan · Liu Yangfan · Wong Ching Him

## Credits and licences

- Code: MIT (see `LICENSE`).
- Data: USGS Landsat Collection 2, ESA WorldCover (CC BY 4.0), Impact Observatory LULC (CC BY 4.0), Copernicus DEM, NASA MODIS, all via Microsoft Planetary Computer. Boundaries and parks © OpenStreetMap contributors (ODbL).
- Fonts: Fraunces, Inter, JetBrains Mono (SIL Open Font License), via Fontsource.
- References: Avdan & Jovanovska (2016) *J. Sensors*; Nichol et al. (2009) *Atmos. Res.*; Oke (1982) *QJRMS*; Voogt & Oke (2003) *RSE*; Weng, Lu & Schubring (2004) *RSE*; Zha, Gao & Ni (2003) *IJRS*.
