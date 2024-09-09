# The Biophysical Engine: A Global Biophysical Processor in the Google Earth Engine & Assessment of Leaf Area Index
This paper introduces Sentinel-2 Biophysical Engine implemented on the Google Earth Engine (GEE) platform, addressing the critical need for high-resolution, near-real-time monitoring of Earth's surfaces in response to population growth and climate change. Leveraging Copernicus Ground-Based Observations for Validation (GBOV) data, the tool utilizes neural network models to derive essential biophysical variables such as Leaf Area Index (LAI). Integrating reflectance values, acquisition geometry parameters, and auxiliary data, the GEE implementation demonstrates accuracy through in-situ validation. Results show promising Root Mean Square Error (RMSE), Mean Absolute Error (MAE), correlation, and bias values. With applications in land surface monitoring, agriculture, and hydrological modeling, this tool has strong potential to contribute sustainable resource management and climate change mitigation. Future developments involve expanding the tool to incorporate Landsat-8 data, further enhancing its global applicability.
## How to use
```javascript
//Dataset definition
var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED"),
    geometry = ee.Geometry.Point([39.88, 36.80]);

//Import library
var Biophysical = require("users/aksoysamett/Biophysical:Biophysical");

//Preparing and filtering satellite images
s2 = s2.filterBounds(geometry).filterDate("2021-4-1","2021-6-1").map(applyScaleFactors)

//Biophysical parameter estimation
var index = "LAI" //FAPAR, FVC, LAI, CAB, CWC
s2 = Biophysical.addBiophysical(s2, index)

//Visualization
Map.addLayer(s2, {min:0, max:4, bands:["LAI"], palette:["blue","cyan","yellow","green"]})
Map.centerObject(geometry, 12)


function applyScaleFactors(image) {
  return image.addBands(image.select('B.*').multiply(0.0001), null, true)
}
```
## How to cite
S. Aksoy, H. Akcay and E. Sertel, "The Biophysical Engine: A Global Biophysical Processor in the Google Earth Engine & Assessment of Leaf Area Index," IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium, Athens, Greece, 2024, pp. 5200-5203, doi: 10.1109/IGARSS53475.2024.10641262. keywords: {Earth;XML;Land surface;Metadata;Internet;Indexes;Engines;Leaf Area Index (LAI);Google Earth Engine;Sentinel 2;In-situ validation;Biophysical variables},
