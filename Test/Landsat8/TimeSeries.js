/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var geometry = /* color: #d63000 */ee.Geometry.Point([39.879248523300284, 36.80380062897793]),
    l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var Biophysical = require("users/aksoysamett/Biophysical:Biophysical");

var index = "LAI" //FAPAR, FVC, LAI, CAB, CWC
var use10mBandsOnly = false
var useXMLAngles = false
var gridCheck = false

l8 = l8.filterBounds(geometry).filterDate("2021-4-1","2021-6-1").map(applyScaleFactors)
l8 = Biophysical.addBiophysical(l8, index, use10mBandsOnly, useXMLAngles, gridCheck)
if(gridCheck)
  l8 = l8.map(function(i){return i.updateMask(i.select("Flag").eq(0))})
print(l8)
print(ui.Chart.image.series(l8.select(index), geometry, ee.Reducer.first(), 10))


function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}