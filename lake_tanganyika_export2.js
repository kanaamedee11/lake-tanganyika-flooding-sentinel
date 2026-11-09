// Lake Tanganyika study, second export: independent checks for the paper's main finding.
// Paste into https://code.earthengine.google.com , press Run, then press "Run" on the 4 tasks in the Tasks tab.
// No API keys are needed. Output goes to the same Drive folder as before.

var AOI = ee.Geometry.Rectangle([29.20, -3.55, 29.55, -3.15]);        // same box as the first export
var AOI_W = ee.Geometry.Rectangle([29.00, -3.45, 29.30, -3.25]);      // extension to the west (Rusizi delta beyond the old edge)
var START_YEAR = 2017, END_YEAR = 2025, SCALE = 30, FOLDER = 'TanganyikaStudy';
var MIN_CLEAR_S2 = 3;     // a pixel needs at least this many clear Sentinel-2 observations in a year

var slope = ee.Terrain.slope(ee.Image('USGS/SRTMGL1_003'));
var flat = slope.lt(5);
var s1all = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .select('VV');
var s2all = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');

function yearRange(y) { y = ee.Number(y); return [ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y, 12, 31).advance(1, 'day')]; }
function stackBands(list) {
  return ee.Image(ee.List(list).iterate(function (cur, acc) { return ee.Image(acc).addBands(ee.Image(cur)); }, ee.Image([])));
}
var YEARS = ee.List.sequence(START_YEAR, END_YEAR);

// Sentinel-1: yearly percent of passes classed as water below a threshold (dB), optionally one orbit direction only
function s1Freq(geom, thr, pass, prefix) {
  return stackBands(YEARS.map(function (y) {
    var r = yearRange(y);
    var col = s1all.filterBounds(geom).filterDate(r[0], r[1]);
    if (pass) col = col.filter(ee.Filter.eq('orbitProperties_pass', pass));
    var w = col.map(function (im) { return im.lt(thr).rename('w'); });
    return w.mean().multiply(100).updateMask(flat).unmask(255).toUint8()
      .rename(ee.String(prefix).cat('_y').cat(ee.Number(y).format('%d')));
  })).clip(geom);
}
// Sentinel-2: yearly percent of CLEAR observations classed as water (MNDWI > 0); 255 where fewer than MIN_CLEAR_S2 clear observations
function s2Freq(geom, prefix) {
  return stackBands(YEARS.map(function (y) {
    var r = yearRange(y);
    var col = s2all.filterBounds(geom).filterDate(r[0], r[1]).map(function (im) {
      var scl = im.select('SCL');
      var clear = scl.eq(4).or(scl.eq(5)).or(scl.eq(6));
      return im.normalizedDifference(['B3', 'B11']).gt(0).rename('w').updateMask(clear);
    });
    var n = col.count();
    return col.mean().multiply(100).updateMask(n.gte(MIN_CLEAR_S2)).unmask(255).toUint8()
      .rename(ee.String(prefix).cat('_y').cat(ee.Number(y).format('%d')));
  })).clip(geom);
}

// 1. Sentinel-2 yearly water frequency, same box as before (independent confirmation of the local change)
Export.image.toDrive({image: s2Freq(AOI, 's2'), description: 'tanganyika_s2_frequency_by_year',
  folder: FOLDER, region: AOI, scale: SCALE, crs: 'EPSG:32735', maxPixels: 1e10, fileFormat: 'GeoTIFF'});

// 2. Sentinel-1 threshold sensitivity: -14 dB and -18 dB (the first export used -16 dB)
Export.image.toDrive({image: s1Freq(AOI, -14, null, 't14').addBands(s1Freq(AOI, -18, null, 't18')),
  description: 'tanganyika_s1_threshold_sensitivity',
  folder: FOLDER, region: AOI, scale: SCALE, crs: 'EPSG:32735', maxPixels: 1e10, fileFormat: 'GeoTIFF'});

// 3. Sentinel-1 split by orbit direction at -16 dB (tests the changing acquisition mix)
Export.image.toDrive({image: s1Freq(AOI, -16, 'ASCENDING', 'asc').addBands(s1Freq(AOI, -16, 'DESCENDING', 'desc')),
  description: 'tanganyika_s1_by_orbit_direction',
  folder: FOLDER, region: AOI, scale: SCALE, crs: 'EPSG:32735', maxPixels: 1e10, fileFormat: 'GeoTIFF'});

// 4. West extension of the Rusizi delta: both sensors over lon 29.00 to 29.30
Export.image.toDrive({image: s1Freq(AOI_W, -16, null, 's1').addBands(s2Freq(AOI_W, 's2')),
  description: 'tanganyika_delta_west_extension',
  folder: FOLDER, region: AOI_W, scale: SCALE, crs: 'EPSG:32735', maxPixels: 1e10, fileFormat: 'GeoTIFF'});

Map.centerObject(AOI_W, 10);
print('Open the Tasks tab (top right) and press Run on the 4 export tasks. Value 255 = no data.');
