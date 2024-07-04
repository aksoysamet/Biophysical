/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var DefinitionDomainGrid = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/DefinitionDomainGrid"),
    DefinitionDomainMinMax = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/DefinitionDomainMinMax"),
    Denormalisation = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/Denormalisation"),
    ExtremeCases = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/ExtremeCases"),
    Normalisation = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/Normalisation"),
    Weights_Layer1_Neurons = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/Weights_Layer1_Neurons"),
    Weights_Layer1_Bias = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/Weights_Layer1_Bias"),
    Weights_Layer2_Neurons = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/Weights_Layer2_Neurons"),
    Weights_Layer2_Bias = ee.ImageCollection("users/aksoysamett/BiophysicalAuxData/Weights_Layer2_Bias");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var XMLAngles = require("users/aksoysamett/Biophysical:XMLAngles");
var RO = (Math.PI/180.0);
var S2BandNames = ["B3","B4","B5","B6","B7","B8A","B11","B12"];
var S2BandNames10m = ["B3","B4","B8"];
var L8BandNames = [".*B3",".*B4",".*B5",".*B6",".*B7"];

function addBiophysicalInternal(collection, index, satellite, angleFunction, gridCheck){
  var definitionDomainGrid = DefinitionDomainGrid.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
  var definitionDomainMinMax = DefinitionDomainMinMax.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
  var denormalisation = Denormalisation.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
  var extremeCases = ExtremeCases.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
  var normalisation = Normalisation.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
  var wL1 = Weights_Layer1_Neurons.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
  var bL1 = Weights_Layer1_Bias.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
  var wL2 = Weights_Layer2_Neurons.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
  var bL2 = Weights_Layer2_Bias.filter(ee.Filter.eq("satellite",  satellite)).filter(ee.Filter.eq("index", index)).first()
  var bandNames = satellite === "L8" ? L8BandNames : (satellite.split("_").pop() === "10m" ? S2BandNames10m : S2BandNames)
  var initialFlagFunction = gridCheck ? minMaxGridCheck.bind(null, bandNames, definitionDomainMinMax, definitionDomainGrid) : minMaxCheck.bind(null, definitionDomainMinMax)

  return collection.map(function(image){
    var prepared = image.select(bandNames)
    var arrI = prepared.toArray();
    var flag = initialFlagFunction(arrI);
    prepared = prepared.addBands(angleFunction(image));
    var arrayImage = prepared.toArray();
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
    flag = flag.bitwiseOr(outMinThr.multiply(2)).bitwiseOr(outMaxThr.multiply(4)).bitwiseOr(outLow.multiply(8)).bitwiseOr(outHig.multiply(16)).rename("Flag");
    nnout = nnout.where(outMinThr,emin)
      .where(outMaxThr,emax)
    return image.addBands(nnout).addBands(flag)
  })
}

function minMaxCheck(definitionDomainMinMax, arrI){
  return arrI.lt(definitionDomainMinMax.select("min")).arrayReduce(ee.Reducer.min(), [0]).and(arrI.gt(definitionDomainMinMax.select("max")).arrayReduce(ee.Reducer.min(), [0])).arrayFlatten([["flag"]]);
}

function minMaxGridCheck(bandNames, definitionDomainMinMax, definitionDomainGrid, arrI){
  var arrImage = ee.Image(10).multiply(arrI.subtract(definitionDomainMinMax.select("min")).divide(definitionDomainMinMax.select("max").subtract(definitionDomainMinMax.select("min")))).add(1).floor().uint8();
  return definitionDomainGrid.eq(arrImage.arrayReshape(ee.Array([1,bandNames.length]),2).arrayRepeat(0, definitionDomainGrid.arrayLength(0))).arrayReduce(ee.Reducer.min(), [1]).arrayReduce(ee.Reducer.max(), [0]).arrayProject([0]).arrayFlatten([["flag"]]).not()
  .or(minMaxCheck(definitionDomainMinMax, arrI))
}

function addAnglesFromMetada(image){
  return ee.Image(image.getNumber("MEAN_INCIDENCE_ZENITH_ANGLE_B4").multiply(RO).cos()).rename("VZ")
    .addBands(ee.Image(image.getNumber("MEAN_SOLAR_ZENITH_ANGLE").multiply(RO).cos()).rename("SZ"))
    .addBands(ee.Image(image.getNumber("MEAN_SOLAR_AZIMUTH_ANGLE").subtract(image.getNumber("MEAN_INCIDENCE_AZIMUTH_ANGLE_B4")).multiply(RO).cos()).rename("RA"))
}
function addAnglesFromXML(image){
  image = XMLAngles.addImageAngles(image)
  return convertAngles(image)
}

function convertAngles(image){
  return image.select("viewZenithMean").multiply(RO).cos().rename("VZ")
    .addBands(image.select("sunZenith").multiply(RO).cos().rename("SZ"))
    .addBands(image.select("sunAzimuth").subtract(image.select("viewAzimuthMean")).multiply(RO).cos().rename("RA"))
}

function addAnglesCollection(collection){
  return ee.ImageCollection(ee.Join.saveFirst("matched").apply(collection, ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA").select(["SAA","SZA","VAA","VZA"],["sunAzimuth","sunZenith","viewAzimuthMean","viewZenithMean"]), ee.Filter.equals({leftField: 'system:time_start', rightField: 'system:time_start'}))
    .map(function(i){return ee.Image(i).addBands(ee.Image(i.get('matched')).divide(100), null, true)}))
}

function addBiophysical(collection, index, is10m, fromXML, gridCheck){
  is10m = is10m || false;
  fromXML = fromXML || false;
  gridCheck = gridCheck || false;
  var angleFunction = fromXML ? addAnglesFromXML : addAnglesFromMetada
  var collectionS2A = collection.filter(ee.Filter.eq("SPACECRAFT_NAME","Sentinel-2A"))
  var collectionS2B = collection.filter(ee.Filter.eq("SPACECRAFT_NAME","Sentinel-2B"))
  var collectionL8 = collection.filter(ee.Filter.eq("SPACECRAFT_ID","LANDSAT_8"))
  return addBiophysicalInternal(collectionS2A, index, is10m ? "S2A_10m" : "S2A", angleFunction, gridCheck).merge(
    addBiophysicalInternal(collectionS2B, index, is10m ? "S2B_10m" : "S2B", angleFunction, gridCheck)).merge(
    addBiophysicalInternal(addAnglesCollection(collectionL8), index, "L8", convertAngles, gridCheck))
}
exports.addBiophysical = addBiophysical