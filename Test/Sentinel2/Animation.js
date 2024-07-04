/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED"),
    geometry = /* color: #d63000 */ee.Geometry.Point([39.879248523300284, 36.80380062897793]);
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var Biophysical = require("users/aksoysamett/Biophysical:Biophysical");

s2 = s2.filterBounds(geometry).filterDate("2021-3-1","2021-6-1").map(applyScaleFactors)
s2 =  Biophysical.addBiophysical(s2, "LAI")
print(s2)
var animation = require('users/gena/packages:animation')
animation.animate(s2, {maxFrames: s2.size(), timeStep: 500, vis:{min:0, max:4, bands:["LAI"], palette:["blue","cyan","yellow","green"]}});
Map.centerObject(geometry, 12)

function applyScaleFactors(image) {
  return image.addBands(image.select('B.*').multiply(0.0001), null, true)
}