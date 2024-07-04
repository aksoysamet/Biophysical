/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var DefinitionDomainMinMax = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/DefinitionDomainMinMax"),
    Denormalisation = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/Denormalisation"),
    ExtremeCases = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/ExtremeCases"),
    Normalisation = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/Normalisation"),
    Weights_Layer1_Neurons = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/Weights_Layer1_Neurons"),
    Weights_Layer1_Bias = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/Weights_Layer1_Bias"),
    Weights_Layer2_Neurons = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/Weights_Layer2_Neurons"),
    Weights_Layer2_Bias = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/Weights_Layer2_Bias");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var satellite = "S2B" //S2A, S2B, S2A_10m, S2B_10m
var index = "Cab" //FAPAR, FVC, LAI, Cab

var testCases = ee.FeatureCollection("users/aksoysamett/BiophysicalAuxData/TestCases/"+satellite+"_"+index)

var RO = (Math.PI/180.0);
var BandNames = ["B3","B4","B5","B6","B7","B8A", "B11","B12"];
var BandNames10m = ["B3","B4","B8"];
var definitionDomainMinMax = DefinitionDomainMinMax.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
var denormalisation = Denormalisation.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
var extremeCases = ExtremeCases.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
var normalisation = Normalisation.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
var wL1 = Weights_Layer1_Neurons.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
var bL1 = Weights_Layer1_Bias.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
var wL2 = Weights_Layer2_Neurons.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
var bL2 = Weights_Layer2_Bias.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
var bandNames = satellite.split("_").pop() === "10m" ? BandNames10m : BandNames

/*****************Test cases to image collection**********************/
var collection = ee.ImageCollection(testCases.map(function(f){
  return ee.Image.constant(f.toArray(bandNames.concat(["VZ","SZ","RA"])).toList()).set(index, f.get(index))
    .rename(bandNames.concat(["VZ","SZ","RA"]))
}))
/*********************************************************************/

var indexCol = collection.map(function(image){
  //var prepared = image.multiply(0.0001) // No need for test cases
  var prepared = image.select(bandNames) // Selects appropriate bands to check flag
  var arrI = prepared.toArray()
  /*prepared = prepared
    .addBands(ee.Image(image.getNumber("MEAN_INCIDENCE_ZENITH_ANGLE_B4").multiply(RO).cos()).rename("VZ"))
    .addBands(ee.Image(image.getNumber("MEAN_SOLAR_ZENITH_ANGLE").multiply(RO).cos()).rename("SZ"))
    .addBands(ee.Image(image.getNumber("MEAN_SOLAR_AZIMUTH_ANGLE").subtract(image.getNumber("MEAN_INCIDENCE_AZIMUTH_ANGLE_B4")).multiply(RO).cos()).rename("RA"))
    No need for adding angle bands, it is already in test image
  */
  prepared = image //Select all bands to calculate index
  var arrayImage = prepared.toArray()
  var normalizedImage = ee.Image(2).multiply(arrayImage.subtract(normalisation.select("min"))).divide(normalisation.select("max").subtract(normalisation.select("min"))).subtract(1).toArray(1);
  var layer1 = wL1.matrixMultiply(normalizedImage).add(bL1);
  var tansig = ee.Image(2).divide(ee.Image(1).add(ee.Image(-2).multiply(layer1).exp())).subtract(1);
  var output = wL2.matrixMultiply(tansig).arrayProject([0]).add(bL2);
  var nnout  = ee.Image(0.5).multiply(output.add(1)).multiply(denormalisation.select("max").subtract(denormalisation.select("min"))).add(denormalisation.select("min")).arrayFlatten([[index]]);
  
  var emin = extremeCases.select("min")
  var emax = extremeCases.select("max")
  var tolerance = extremeCases.select("tolerance")
  var outMin = nnout.gt(emin.subtract(tolerance));
  var outMax = nnout.lt(emax.add(tolerance));
  var outMinThr = outMin.and(nnout.lt(emin));
  var outMaxThr = outMax.and(nnout.gt(emax));
  var outLow = outMin.not();
  var outHig = outMax.not();
  nnout = nnout.where(outMinThr,emin)
    .where(outMaxThr,emax)
  return image.addBands(nnout)
})

var compare = indexCol.select(index).map(function(i){
  var pre = i.sample(ee.Geometry.Point([0,0]), 111000).first().toDictionary().getNumber(index);
  var tru = i.getNumber(index);
  var err = pre.subtract(tru)
  return ee.Feature(null,{Predicted:pre, GroundTruth:tru, err:err});
})
//print(compare.sort("dif", false).limit(10))
Export.table.toDrive(compare.sort("err", false), satellite+"_"+index, "TestBiophysical")
