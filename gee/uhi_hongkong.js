// =====================================================================
// Hotter Streets, Greener Relief? - Hong Kong surface UHI 2014-2025
// Google Earth Engine version of the analysis (paste into code.earthengine.google.com)
//
// Same method as the Python pipeline in ../analysis, which produced the
// numbers in the presentation (it ran on the identical Landsat Collection 2
// data via Microsoft Planetary Computer).
// Every quantitative layer sits on ONE 30 m grid (course feedback).
// =====================================================================

var hk   = ee.Geometry.Rectangle([113.82, 22.14, 114.45, 22.57]);
var grid = ee.Projection('EPSG:32650').atScale(30);           // common 30 m grid
var land = ee.Image('ESA/WorldCover/v200/2021').select('Map').neq(80);

// ---------- 1. Landsat 8/9: LST (USGS), NDVI, NDBI, brightness temperature
function prep(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(parseInt('11011', 2)).eq(0);      // fill, dilated cloud, cloud, shadow
  var lst = img.select('ST_B10').multiply(0.00341802).add(149.0).subtract(273.15).rename('LST');
  var sr  = img.select(['SR_B4', 'SR_B5', 'SR_B6']).multiply(0.0000275).add(-0.2);
  var ndvi = sr.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
  var ndbi = sr.normalizedDifference(['SR_B6', 'SR_B5']).rename('NDBI');
  // at-sensor radiance (ST_TRAD) -> brightness temperature with band-10 constants
  var isL9 = ee.String(img.get('SPACECRAFT_ID')).equals('LANDSAT_9');
  var k1 = ee.Number(ee.Algorithms.If(isL9, 799.0284, 774.8853));
  var k2 = ee.Number(ee.Algorithms.If(isL9, 1329.2405, 1321.0789));
  var L  = img.select('ST_TRAD').multiply(0.001);
  var bt = ee.Image(k2).divide(ee.Image(k1).divide(L).add(1).log()).rename('BT');  // Kelvin
  return lst.addBands([ndvi, ndbi, bt]).updateMask(clear).updateMask(lst.gt(10).and(lst.lt(70)));
}
function summer(y0, y1) {
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(hk).filterDate(y0 + '-05-01', y1 + '-10-31')
    .filter(ee.Filter.calendarRange(5, 10, 'month'))
    .filter(ee.Filter.lt('CLOUD_COVER', 70))
    .map(prep).median().reproject(grid).clip(hk).updateMask(land);
}
var e14 = summer(2014, 2015), e19 = summer(2019, 2020), e24 = summer(2024, 2025);

// ---------- 2. Tutorial LST (NDVI -> FV -> emissivity -> LST), as a check
function tutorialLST(img) {
  var p = img.select('NDVI').reduceRegion({reducer: ee.Reducer.percentile([1, 99]),
            geometry: hk, scale: 30, maxPixels: 1e10});
  var nmin = ee.Number(p.get('NDVI_p1')), nmax = ee.Number(p.get('NDVI_p99'));
  var fv = img.select('NDVI').subtract(nmin).divide(nmax.subtract(nmin)).clamp(0, 1).pow(2);
  var em = fv.multiply(0.004).add(0.986);
  var bt = img.select('BT');
  return bt.divide(bt.multiply(10.895e-6 / 1.4388e-2).multiply(em.log()).add(1))
           .subtract(273.15).rename('LST_tutorial');
}
var lstTut24 = tutorialLST(e24);

// ---------- 3. Land-cover shares, aggregated 10 m -> 30 m (never used at 10 m)
var wc = ee.Image('ESA/WorldCover/v200/2021').select('Map');
function share(cls) {
  return wc.eq(cls).reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 64}).reproject(grid);
}
var dem   = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM').mosaic().reproject(grid);
var urban = share(50).gte(0.8);                       // >= 80 % built-up
var rural = share(10).gte(0.8).and(dem.lt(200));      // >= 80 % trees, below 200 m

// ---------- Q1: UHI intensity per epoch
function uhi(img, label) {
  var u = img.select('LST').updateMask(urban).reduceRegion(ee.Reducer.mean(), hk, 30, null, null, false, 1e10);
  var r = img.select('LST').updateMask(rural).reduceRegion(ee.Reducer.mean(), hk, 30, null, null, false, 1e10);
  print(label + ' UHI intensity (deg C)', ee.Number(u.get('LST')).subtract(r.get('LST')));
}
uhi(e14, '2014-15'); uhi(e19, '2019-20'); uhi(e24, '2024-25');

// Hot spots: LST more than 1 SD above the land mean
var st = e24.select('LST').reduceRegion(ee.Reducer.mean().combine(ee.Reducer.stdDev(), null, true), hk, 30, null, null, false, 1e10);
var z24 = e24.select('LST').subtract(ee.Number(st.get('LST_mean'))).divide(ee.Number(st.get('LST_stdDev')));
// Urban Thermal Field Variance Index (tutorial): (Ts - Tmean) / Ts, in Kelvin
var utfvi = e24.select('LST').add(273.15).subtract(ee.Number(st.get('LST_mean')).add(273.15))
              .divide(e24.select('LST').add(273.15)).rename('UTFVI');

// ---------- Q2: LST vs NDVI / NDBI, stratified sample on the 30 m grid
var strata = share(50).multiply(5).floor().clamp(0, 4).int().rename('stratum');
var sample = e24.select(['LST', 'NDVI', 'NDBI']).addBands(strata).stratifiedSample({
  numPoints: 400, classBand: 'stratum', region: hk, scale: 30, projection: grid, seed: 9083, geometries: false});
print('LST vs NDVI', ui.Chart.feature.byFeature(sample, 'NDVI', 'LST').setChartType('ScatterChart'));
var reg = sample.map(function(f) { return f.set('one', 1); })
  .reduceColumns(ee.Reducer.linearRegression(3, 1), ['one', 'NDVI', 'NDBI', 'LST']);
print('OLS LST = a + b*NDVI + c*NDBI  [a, b, c]', reg.get('coefficients'));

// ---------- Q3: Northern Metropolis change, 2014-15 -> 2024-25
var nm = ee.Geometry.Rectangle([113.99, 22.425, 114.165, 22.545]);
var dNDBI = e24.select('NDBI').subtract(e14.select('NDBI'));
var dNDVI = e24.select('NDVI').subtract(e14.select('NDVI'));
var dLST  = e24.select('LST').subtract(e14.select('LST'));
var developed = dNDBI.gt(0.08).and(dNDVI.lt(-0.10));
var stable    = dNDBI.abs().lt(0.03).and(dNDVI.abs().lt(0.05));
var dDev = dLST.updateMask(developed).reduceRegion(ee.Reducer.mean(), nm, 30, null, null, false, 1e10);
var dStb = dLST.updateMask(stable).reduceRegion(ee.Reducer.mean(), nm, 30, null, null, false, 1e10);
print('Development warming, Northern Metropolis (deg C)', ee.Number(dDev.get('LST')).subtract(dStb.get('LST')));

// ---------- Supplementary, VISUAL ONLY: MODIS night LST at 1 km
var night = ee.ImageCollection('MODIS/061/MOD11A2').filterDate('2024-06-01', '2024-09-30')
  .select('LST_Night_1km').mean().multiply(0.02).subtract(273.15).clip(hk);

// ---------- Maps
var heat = ['#08306b', '#2171b5', '#6baed6', '#fee08b', '#fdae61', '#f46d43', '#a50026'];
Map.centerObject(hk, 11);
Map.addLayer(e24.select('LST'), {min: 25, max: 45, palette: heat}, 'LST 2024-25');
Map.addLayer(lstTut24, {min: 25, max: 45, palette: heat}, 'LST 2024-25 (tutorial method)', false);
Map.addLayer(e24.select('NDVI'), {min: 0, max: 0.8, palette: ['white', 'darkgreen']}, 'NDVI', false);
Map.addLayer(e24.select('NDBI'), {min: -0.3, max: 0.2, palette: ['darkgreen', 'white', 'black']}, 'NDBI', false);
Map.addLayer(z24.gt(1).selfMask(), {palette: ['red']}, 'Hot spots (z > 1)', false);
Map.addLayer(utfvi, {min: -0.02, max: 0.03, palette: heat}, 'UTFVI', false);
Map.addLayer(dLST, {min: -4, max: 4, palette: ['#2166ac', 'white', '#b2182b']}, 'LST change', false);
Map.addLayer(developed.selfMask(), {palette: ['magenta']}, 'Newly developed', false);
Map.addLayer(night, {min: 24, max: 30, palette: heat}, 'MODIS night LST 1 km (visual only)', false);
