// Lake Tanganyika (Bujumbura / Rusizi delta) water study: export script for the Earth Engine Code Editor.
// Paste into https://code.earthengine.google.com , press Run, then press "Run" on the 3 tasks in the Tasks tab.
// No API keys are needed: you are signed in with your own Google account.

// ---------------- settings you may change ----------------
var AOI = ee.Geometry.Rectangle([29.20, -3.55, 29.55, -3.15]);   // lon/lat box around Bujumbura and the Rusizi plain
var START_YEAR = 2017, END_YEAR = 2025;
var S1_THRESH_DB = -16;     // Sentinel-1 VV below this = water (radar-dark)
var MNDWI_THRESH = 0.0;     // Sentinel-2 MNDWI above this = water
var SCALE = 30;             // metres, keeps the files small
var FOLDER = 'TanganyikaStudy';   // Google Drive folder that will be created

// ---------------- helpers ----------------
var slope = ee.Terrain.slope(ee.Image('USGS/SRTMGL1_003'));
var flat = slope.lt(5);                       // removes radar-shadow false water on hillsides
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(AOI)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .select('VV');
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(AOI);
var jrc = ee.ImageCollection('JRC/GSW1_4/MonthlyHistory').filterBounds(AOI);
var pxkm2 = ee.Image.pixelArea().divide(1e6);

function areaOf(maskImg) {   // km2 of non-zero pixels in maskImg inside AOI (null if image fully masked)
  var r = pxkm2.updateMask(maskImg).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: AOI, scale: SCALE, maxPixels: 1e10, tileScale: 4});
  return r.get('area');
}

var nMonths = (END_YEAR - START_YEAR + 1) * 12;
var months = ee.List.sequence(0, nMonths - 1).map(function (i) {
  var d0 = ee.Date.fromYMD(START_YEAR, 1, 1).advance(ee.Number(i), 'month');
  var d1 = d0.advance(1, 'month');
  // Sentinel-1: monthly median backscatter
  var s1m = s1.filterDate(d0, d1);
  var n1 = s1m.size();
  var s1img = ee.Image(ee.Algorithms.If(n1.gt(0), s1m.median(),
      ee.Image.constant(0).rename('VV').updateMask(ee.Image.constant(0))));
  var s1water = s1img.lt(S1_THRESH_DB).and(flat).selfMask();
  // Sentinel-2: monthly median MNDWI on clear pixels (SCL 4 vegetation, 5 bare, 6 water)
  var s2m = s2.filterDate(d0, d1).map(function (im) {
    var scl = im.select('SCL');
    var clear = scl.eq(4).or(scl.eq(5)).or(scl.eq(6));
    return im.normalizedDifference(['B3', 'B11']).rename('MNDWI').updateMask(clear);
  });
  var n2 = s2m.size();
  var s2img = ee.Image(ee.Algorithms.If(n2.gt(0), s2m.median(),
      ee.Image.constant(0).rename('MNDWI').updateMask(ee.Image.constant(0))));
  var s2water = s2img.gt(MNDWI_THRESH).selfMask();
  var s2valid = s2img.mask().selfMask();
  // JRC monthly water history (version 1.4 ends in 2021): 2 = water, 1 = land, 0 = no data
  var jm = jrc.filterDate(d0, d1);
  var nj = jm.size();
  var jimg = ee.Image(ee.Algorithms.If(nj.gt(0), jm.first().select('water'),
      ee.Image.constant(0).rename('water')));
  return ee.Feature(null, {
    month: d0.format('YYYY-MM'),
    s1_scenes: n1, s1_water_km2: areaOf(s1water),
    s2_scenes: n2, s2_water_km2: areaOf(s2water), s2_valid_km2: areaOf(s2valid),
    jrc_water_km2: areaOf(jimg.eq(2).selfMask()), jrc_valid_km2: areaOf(jimg.gte(1).selfMask())
  });
});

Export.table.toDrive({
  collection: ee.FeatureCollection(months),
  description: 'tanganyika_monthly_water_timeseries',
  folder: FOLDER, fileFormat: 'CSV'
});

// ---------------- yearly Sentinel-1 water frequency (percent of scenes classed as water) ----------------
var yearImages = ee.List.sequence(START_YEAR, END_YEAR).map(function (y) {
  y = ee.Number(y);
  var col = s1.filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y, 12, 31).advance(1, 'day'))
    .map(function (im) { return im.lt(S1_THRESH_DB).rename('w'); });
  return col.mean().multiply(100).updateMask(flat).unmask(255).toUint8()
    .rename(ee.String('y').cat(y.format('%d')));
});
var yearly = ee.Image(yearImages.iterate(function (cur, acc) {
  return ee.Image(acc).addBands(ee.Image(cur));
}, ee.Image([]))).clip(AOI);

Export.image.toDrive({
  image: yearly, description: 'tanganyika_s1_water_frequency_by_year',
  folder: FOLDER, region: AOI, scale: SCALE, crs: 'EPSG:32735', maxPixels: 1e10, fileFormat: 'GeoTIFF'
});   // value 255 = no data / hillside excluded

// ---------------- reference layers ----------------
var ref = ee.Image.cat([
  ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0).rename('jrc_occurrence'),
  ee.ImageCollection('ESA/WorldCover/v200').first().select('Map').rename('landcover'),
  ee.Image('USGS/SRTMGL1_003').rename('elevation_m'),
  slope.rename('slope_deg')
]).toInt16().clip(AOI);

Export.image.toDrive({
  image: ref, description: 'tanganyika_reference_layers',
  folder: FOLDER, region: AOI, scale: SCALE, crs: 'EPSG:32735', maxPixels: 1e10, fileFormat: 'GeoTIFF'
});

Map.centerObject(AOI, 10);
Map.addLayer(yearly.select('y2020').selfMask(), {min: 0, max: 100, palette: ['white', 'blue']}, 'S1 water frequency 2020');
print('Now open the Tasks tab (top right) and press Run on the 3 export tasks.');
