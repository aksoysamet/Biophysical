/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED"),
    cilingoz_forest = 
    /* color: #0b4a8b */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[28.13198186114467, 41.49778174748733],
          [28.13198186114467, 41.47127480490167],
          [28.15219457453176, 41.47127480490167],
          [28.15219457453176, 41.49778174748733]]]),
    iznik_olive_grove = 
    /* color: #ffc82d */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[29.39577140340616, 40.41522902860147],
          [29.39577140340616, 40.39960649212002],
          [29.42844906525518, 40.39960649212002],
          [29.42844906525518, 40.41522902860147]]]),
    ceylanpinar_agricultural = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[39.55938319427734, 36.76146884083449],
          [39.55938319427734, 36.73570129785626],
          [39.587454173876, 36.73570129785626],
          [39.587454173876, 36.76146884083449]]]),
    snapCollection = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/TestCases/SNAP");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var sAOI = "Çilingoz" //Çilingoz, İznik, Ceylanpinar
var index = "LAI"; // FAPAR, FVC, LAI, Cab, CWC
var is10m = false, //When true, uses only 10 m bands
  fromXML = true, //When true, sun and view angles are extracted from GCP XML
  gridCheck = false; //When true, input flags contains grid check


var Biophysical = require("users/aksoysamett/Biophysical:Biophysical");
var aois = {"Çilingoz":cilingoz_forest,"İznik":iznik_olive_grove, "Ceylanpinar":ceylanpinar_agricultural}
var imgs = {"Çilingoz":"L2A_T35TNF_A031464_20210701T090332", "İznik":"L2A_T35TPE_A022727_20210713T085013", "Ceylanpinar":"L2A_T37SEA_A021554_20210422T081230"}

var aoi = aois[sAOI];
var img = s2.filter(ee.Filter.eq("GRANULE_ID", imgs[sAOI])).first();
var snap = snapCollection.filter(ee.Filter.eq("GRANULE_ID", imgs[sAOI] + (is10m ? "_10m": ""))).first()
//print(img)
print(snap)
Map.centerObject(aoi, 14)
Map.addLayer(img, {min:0, max:3000, bands:["B4","B3","B2"]}, "RGB", false)
img = applyScaleFactors(img)
var bio =  Biophysical.addBiophysical(ee.ImageCollection([img]), index, is10m, fromXML, gridCheck)
print(bio)
bio = bio.first()
snap = snap.clip(aoi)
bio = bio.clip(aoi)
bio.select(index).reduceRegion(ee.Reducer.minMax(), aoi, 20).evaluate(function(minMax){
  Map.addLayer(snap.select([index, index+"_flags"]), {min:minMax[index + "_min"], max:minMax[index + "_max"], bands:[index]}, "Snap")
  Map.addLayer(bio.select([index, "Flag"]), {bands:[index], min:minMax[index + "_min"], max:minMax[index + "_max"]}, "GEE")
})
var error = bio.select(index).subtract(snap.select(index)).rename("Error")
var hist = ee.Dictionary(error.reduceRegion(ee.Reducer.histogram(50).unweighted(), aoi, 10).get("Error"))
var histFC = ee.FeatureCollection(ee.List(hist.get("bucketMeans")).zip(hist.get("histogram")).map(function(i){
  i = ee.List(i)
  return ee.Feature(null, {"bucketMeans":i.get(0), "histogram":i.get(1)})
}))
print(histFC)
Export.table.toDrive({collection:histFC, description:sAOI.replace("Ç", "C").replace("İ", "I")+"_"+index+ (is10m ? "_10m" : "") + (fromXML ? "_XML" : ""),selectors:["bucketMeans","histogram"]})
print("Histogram", ui.Chart.image.histogram(error, aoi, 10, 50))

error.reduceRegion(ee.Reducer.frequencyHistogram(), aoi, 10).evaluate(function(freq){
  if(freq["Error"].length <= 5000)
    print("Frequency of errors", freq["Error"]);
  var arr = Object.keys(freq["Error"]);
  var min = Math.min.apply(null, arr),
  max = Math.max.apply(null, arr);
  error = error.mask(error.lt(0)).visualize({min:min, max:0, palette:["Blue","Black"]})
    .blend(error.mask(error.eq(0)).visualize({palette:["Black"]}))
    .blend(error.mask(error.gt(0)).visualize({min:0, max:max, palette:["Black","Red"]}))
  Map.addLayer(error.reproject("EPSG:32635", null, 10), {}, "Errors")
  //Legend %50 0 - min, %50 0 - max
  Map.add(createLegend(min, max))
})

function applyScaleFactors(image) {
  return image.addBands(image.select('B.*').multiply(0.0001), null, true)
}

function createLegend(min, max) {
    var legend = ui.Panel({
    style: {
      position: 'bottom-left',
      padding: '8px 15px'
    }
  })

  // Create legend title
  var legendTitle = ui.Label({
    value: index+' Errors',
    style: {
      fontWeight: 'bold',
      fontSize: '18px',
      margin: '0 0 4px 0',
      padding: '0'
      }
  });
  
   // Add the title to the panel
  legend.add(legendTitle); 
  
  var gtPanel = ui.Panel({style: {width: '120px', height:'200px'}});
  gtPanel.setLayout("absolute");
  var lon = ee.Image.pixelLonLat().select('latitude');
  var gradient, legendImage;
  if (min === 0 || max === 0){
    gradient = lon.multiply((max-min)/100.0).add(min);
    legendImage = gradient.visualize({min:min, max:max, palette:(min ? ["Blue","Black"] : ["Black","Red"])})
  } else {
    gradient = lon.mask(lon.lt(50)).multiply((0-min)/50.0).add(min);
    var gradient2 = lon.mask(lon.gte(50)).subtract(50).multiply((max-0)/50.0).add(0);
    legendImage = gradient.visualize({min:min, max:0, palette:["Blue","Black"]})
      .blend(gradient2.visualize({min:0, max:max, palette:["Black","Red"]}))
  }
  var thumbnail = ui.Thumbnail({
    image: legendImage, 
    params: {bbox:'0,0,10,100', dimensions:'10x200'},  
    style: {padding: '1px', position: 'bottom-left'}
  });
  gtPanel.add(thumbnail);
  gtPanel.add(ui.Label(getNumberString(max),{position: 'top-right', width: '60px', margin: '-20px 20px 0 0'}))
  if (!(min === 0 || max === 0))gtPanel.add(ui.Label("0.0000", {position: 'top-right', margin: '50% 20px 0 0', width: '60px'}))
  gtPanel.add(ui.Label(getNumberString(min),{position: 'bottom-right', width: '60px', margin: '0 20px -3px 0'}))

  legend.add(gtPanel);
  return legend
}

function getNumberString(number){
  if(Math.abs(number) > 0.0009)
    return number.toFixed(4)
  else if(number === 0)
   return "0.0000"
  return number.toExponential(1)
}