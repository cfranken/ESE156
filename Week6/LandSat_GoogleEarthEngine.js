// Plot Landsat 8 band value 
var imageCollection = ee.ImageCollection("LANDSAT/LC08/C01/T1_RT_TOA")
var lat = 39.9596//° N, 121.6219
var lon = -121.4219 // Dec 9 2017
// var lat = 18.333 // August 25 vs Sept. 10
// var lon = -64.9

var dLat = 0.055
var dLon = 0.055

var re =
    ee.Geometry.Rectangle(lon-dLon, lat-dLat, lon+dLon, lat+dLat);

var bands = ['B2','B4','B3','B5','B6','B7','B9','B10'];//'B10', 'B11'
var landsat8Toa = imageCollection
    .filterDate('2018-10-01', '2018-11-20')
    .select(bands);

// Create an image time series chart.
var chart = ui.Chart.image.series({
  imageCollection: landsat8Toa,
  region: re,
  reducer: ee.Reducer.mean(),
  scale: 100
});

// Add the chart to the map.
chart.style().set({
  position: 'bottom-right',
  width: '500px',
  height: '300px'
});
Map.add(chart);

// Outline and center San Francisco on the map.
var sfLayer = ui.Map.Layer(re, {color: 'FF0000'}, 'SF');
Map.layers().add(sfLayer);
Map.setCenter(lon, lat, 10);

// Create a label on the map.
var label = ui.Label('Click a point on the chart to show the image for that date.');
Map.add(label);

// When the chart is clicked, update the map and label.
chart.onClick(function(xValue, yValue, seriesName) {
  if (!xValue) return;  // Selection was cleared.

  // Show the image for the clicked date.
  var equalDate = ee.Filter.equals('system:time_start', xValue);
  var image = ee.Image(landsat8Toa.filter(equalDate).first());
  // Compute the EVI using an expression.
  var Te = image.expression(
    '1*NIR', {
      'NIR': image.select('B10')
    });
  var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('B5'),
      'RED': image.select('B4'),
      'BLUE': image.select('B2')
    });
  var visParams = {bands: ['B4', 'B3', 'B2'], max: 0.3};

  var NIRv = image.expression(
    '((NIR - RED) / (NIR + RED)) * NIR', {
      'NIR': image.select('B5'),
      'RED': image.select('B4'),
      'BLUE': image.select('B2')
});
  var ndvi = image.normalizedDifference(['B5', 'B4']);
  var l8Layer = ui.Map.Layer(NIRv, {
  //  gamma: 1.3,
    min: 0.0,
    max: 0.6,palette: ['white', 'green']
//    bands: ['B5', 'B5', 'B5']
  }, 'NIRv');
  var l8Layer2 = ui.Map.Layer(evi, {
  //  gamma: 1.3,
    min: 0.0,
    max: 0.8,palette: ['white', 'green']
//    bands: ['B5', 'B5', 'B5']
  }, 'EVI');
  
var TMasked = Te.updateMask(Te.gte(320));
//var TParams = {bands: ['B10'], max: 420, min:320};

  var l8LayerT = ui.Map.Layer(TMasked, {
   // gamma: 1.3,
    min: 300,
    max: 400,
    palette: ['green', 'red']
  }, 'Thermal');
  
  
  
  //Map.layers().reset([l8Layer2, l8Layer,visParams,l8LayerT]);
  Map.layers().add(l8Layer2);
  //ap.layers().add(l8LayerT);
  Map.addLayer(image, visParams, 'true-color composite');
  Map.layers().add(l8LayerT);
  
  //Map.addLayer(image, TParams, 'Thermal');
 // Map.addLayer(image, l8LayerT );
  // Show a label with the date on the map.
  label.setValue((new Date(xValue)).toUTCString());
});
