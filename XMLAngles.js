var names = ee.List.sequence(0,21).map(function(i){return ee.String(ee.Number(i).int())});
function getImageXMLURL(img){
  return ee.String("gs://gcp-public-data-sentinel-2/")
  .cat(ee.Algorithms.If(img.getString("GRANULE_ID").index("L2").neq(-1), "L2/", ""))
  .cat("tiles/")
  .cat(img.getString("MGRS_TILE").replace("(\\d+)(\\w)(\\w+)", "$1/$2/$3/"))
  .cat(img.getString("PRODUCT_ID"))
  .cat(".SAFE/GRANULE/")
  .cat(img.getString("GRANULE_ID"))
  .cat("/MTD_TL.xml");
}
exports.getImageXMLURL = getImageXMLURL;

function getImageXMLContent(img){
  return ee.Blob(getImageXMLURL(img)).string();
}
exports.getImageXMLContent = getImageXMLContent;

function addImageSunAngles(img){
  var xml = getImageXMLContent(img);
  var proj = img.select(0).projection().atScale(5000);
  var pixelCoords = ee.Image.pixelCoordinates(proj);
  var mask = pixelCoords.gte(0).and(pixelCoords.lte(22))
  mask = mask.select("x").and(mask.select("y"))
  pixelCoords = pixelCoords.updateMask(mask).int().select(["y","x"])
  return addImageSunAnglesInternal(img, xml, pixelCoords);
}
exports.addImageSunAngles = addImageSunAngles;

function addImageViewAngles(img){
  var xml = getImageXMLContent(img);
  var proj = img.select(0).projection().atScale(5000);
  var pixelCoords = ee.Image.pixelCoordinates(proj);
  var mask = pixelCoords.gte(0).and(pixelCoords.lte(22))
  mask = mask.select("x").and(mask.select("y"))
  pixelCoords = pixelCoords.updateMask(mask).int().select(["y","x"])
  return addImageViewAnglesInternal(img, xml, pixelCoords);
}
exports.addImageViewAngles = addImageViewAngles;

function addImageAngles(img){
  var xml = getImageXMLContent(img);
  var proj = img.select(0).projection().atScale(5000);
  var pixelCoords = ee.Image.pixelCoordinates(proj);
  var mask = pixelCoords.gte(0).and(pixelCoords.lte(22))
  mask = mask.select("x").and(mask.select("y"))
  pixelCoords = pixelCoords.updateMask(mask).int().select(["y","x"])
  var sungAngles = addImageSunAnglesInternal(img, xml, pixelCoords);
  return addImageViewAnglesInternal(sungAngles, xml, pixelCoords);
}
exports.addImageAngles = addImageAngles;

function addImageSunAnglesInternal(img, xml, pixelCoords){
  var sunAnglesGrid = xml.slice(xml.index("<Sun_Angles_Grid>"),xml.index("</Sun_Angles_Grid>"));
  var sunZenith = sunAnglesGrid.slice(sunAnglesGrid.index("<Zenith>"),sunAnglesGrid.index("</Zenith>")).match("<VALUES>.*</VALUES>","g").slice(0,22).map(angleListfromString);
  var sunAzimuth = sunAnglesGrid.slice(sunAnglesGrid.index("<Azimuth>"),sunAnglesGrid.index("</Azimuth>")).match("<VALUES>.*</VALUES>","g").slice(0,22).map(angleListfromString);
  var sunZAnglesImage = arrayToImage(sunZenith, pixelCoords).rename("sunZenith");
  var sunAAnglesImage = arrayToImage(sunAzimuth, pixelCoords).rename("sunAzimuth");
  return img.addBands(sunZAnglesImage).addBands(sunAAnglesImage);
}

function addImageViewAnglesInternal(img, xml, pixelCoords){
  var viewingAnglesGrid = xml.slice(xml.index("<Viewing_Incidence_Angles_Grids"),xml.rindex("</Viewing_Incidence_Angles_Grids>").add(34)).split("</Viewing_Incidence_Angles_Grids>", "g").slice(0,-1);
  var viewingAngles = ee.FeatureCollection(viewingAnglesGrid.map(function(str){
    str = ee.String(str);
    var bandDetector = str.match("(?:<Viewing_Incidence_Angles_Grids bandId=\")(\\d*)(?:\" detectorId=\")(\\d*)(?:\">)");
    var bandId = ee.Number.parse(bandDetector.get(1));
    var detectorId = ee.Number.parse(bandDetector.get(2));
    var zenithAngles = angleFeaturesFromlist(str.slice(str.index("<Zenith>"),str.index("</Zenith>")).match("<VALUES>.*</VALUES>","g"), bandId, detectorId, "Zenith");
    var azimuthAngles = angleFeaturesFromlist(str.slice(str.index("<Azimuth>"),str.index("</Azimuth>")).match("<VALUES>.*</VALUES>","g"), bandId, detectorId, "Azimuth");
    return [zenithAngles,azimuthAngles];
  }).flatten());
  var azimuthAngles = mosaicDetector(viewingAngles.filter(ee.Filter.eq("Type","Azimuth")));
  var zenithAngles = mosaicDetector(viewingAngles.filter(ee.Filter.eq("Type","Zenith")));

  var viewAnglesImages = ee.ImageCollection(ee.List.sequence(0, 12).map(function(band){
    var va = arrayToImage(changeNaNFromListofList(azimuthAngles.get(band), 500), pixelCoords).rename("va")
    var vz = arrayToImage(changeNaNFromListofList(zenithAngles.get(band), 500), pixelCoords).rename("vz")
    return va.updateMask(va.neq(500)).addBands(vz.updateMask(va.neq(500)))
  }))
  var reducer = ee.Reducer.mean().combine({
    reducer2:ee.Reducer.count(),
    sharedInputs:true
  });
  var viewZenithMean = viewAnglesImages.select("vz").reduce(reducer);
  viewZenithMean = viewZenithMean.select([0],["viewZenithMean"]).mask(viewZenithMean.select(1).eq(13));
  var viewAzimuthMean = viewAnglesImages.select("va").reduce(reducer);
  viewAzimuthMean = viewAzimuthMean.select([0],["viewAzimuthMean"]).mask(viewAzimuthMean.select(1).eq(13));
  return img.addBands(viewZenithMean).addBands(viewAzimuthMean);
}

function arrayToImage(array, pixelCoords){
  var arrayImage = ee.Image(ee.Array(array));
  return arrayImage.arrayGet(pixelCoords);
}
exports.arrayToImage = arrayToImage;

function changeNaNFromListofList(listOfList, newValue){
  return ee.List(listOfList).map(function(list){
    return ee.List(list).map(function(number){
      return ee.Algorithms.If(number, number, newValue)
    })
  })
}

function angleFeaturesFromlist(angleList, bandId, detectorId, type){
  return ee.List.sequence(0, 21).map(function(i){
    var angles = angleListfromString(angleList.get(i));
    return ee.Feature(null, ee.Dictionary.fromLists(names, angles))
    .set("Type",type)
    .set("bandId",bandId)
    .set("idx",i)
    .set("detectorId",detectorId);
  });
}

function mosaicDetector(angles)
{
  return ee.List.sequence(0, 12).map(function(band){
    var bandAngles = angles.filter(ee.Filter.eq("bandId",band));
    return ee.List.sequence(0, 21).map(function(i){
      return bandAngles.filter(ee.Filter.eq("idx",i))
      .sort("detectorId")
      .reduceColumns(ee.Reducer.lastNonNull().repeat(22), names).get("last");
    });
  });
}

function angleListfromString(angleString) {
  return ee.String(angleString)
  .slice(8,-9)
  .split("\\s")
  .slice(0, -1)
  .map(function(k){
    return ee.Algorithms.If(ee.String(k).compareTo("NaN"), ee.Number.parse(k), null);
  });
}