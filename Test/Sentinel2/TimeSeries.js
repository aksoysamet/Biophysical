/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED"),
    geometry = /* color: #d63000 */ee.Geometry.Point([39.879248523300284, 36.80380062897793]);
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var Biophysical = require("users/aksoysamett/Biophysical:Biophysical");

var index = "LAI" //FAPAR, FVC, LAI, CAB, CWC
var use10mBandsOnly = false
var useXMLAngles = false
var gridCheck = false

s2 = s2.filterBounds(geometry).filterDate("2021-4-1","2021-6-1").map(applyScaleFactors)
s2 = Biophysical.addBiophysical(s2, index, use10mBandsOnly, useXMLAngles, gridCheck)
if(gridCheck)
  s2 = s2.map(function(i){return i.updateMask(i.select("Flag").eq(0))})
print(s2)
print(ui.Chart.image.series(s2.select(index), geometry, ee.Reducer.first(), 10))

function applyScaleFactors(image) {
  return image.addBands(image.select('B.*').multiply(0.0001), null, true)
}