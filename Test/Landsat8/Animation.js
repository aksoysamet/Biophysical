/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"),
    geometry = /* color: #d63000 */ee.Geometry.Point([39.879248523300284, 36.80380062897793]);
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var Biophysical = require("users/aksoysamett/Biophysical:Biophysical");

l8 = l8.filterBounds(geometry).filterDate("2021-3-1","2021-6-1").map(applyScaleFactors)
l8 =  Biophysical.addBiophysical(l8, "LAI")
print(l8)
var animation = require('users/gena/packages:animation')
animation.animate(l8, {maxFrames: l8.size(), timeStep: 500, vis:{min:0, max:4, bands:["LAI"], palette:["blue","cyan","yellow","green"]}});
Map.centerObject(geometry, 12)

function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}