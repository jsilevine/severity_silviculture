////////////////////////////    WARNING: Use at your own risk     //////////////////////////////////
//
// Script to develop and export maps of Composite Burn Index (CBI) predictions, and spectral imagery that
// supports the CBI predictions, for fires in North America.
// For updates or questions on this code, please contact:
//        Sean Parks,     sean.parks@usda.gov


//////////////////////     IMPORTANT - PLEASE READ       ////////////////////////////////
//  1) Create a folder called 'cbi' in your Google Drive account before running script.
//  2) Click the ‘Run’ button (above); then under "Tasks" tab (left), click ‘Run’ buttons to export each image
//  3) Output files are named by concatenating the Fire_ID (see description below) and the relevant Earth Engine output.
//     For example: MT4751811328120110719_CBI or MT4751811328120110719_CBI_bc
//  4) If there is no imagery available for pre or post-fire period (e.g. early 1980s may lack data)
//     the resulting CBI, CBI_bc and other imagery will have NO DATA.
//  5) If you use this script or any of the products produced with the script, please cite the following paper:
//            Parks SA, Holsinger LM, Koontz MJ, Collins, Whitman E, et al. 2019. Giving ecological meaning to satellite-derived
//            fire severity metrics across North American forests within Google Earth Engine. Remote Sensing 2019, 11.


////////////////            UPDATE 4/8/2021             ///////////////////////////////////
// The original Random Forest classifier used in this study and associated code was deprecated on 3/1/2021.
// We have updated the code to use the new Random Forest classifier (smileRandomForest).
// We re-ran the 5-fold cross-validation using the new classifier and the resulting statistics originally
// reported are identical (cross-validated R2 = 0.72). Consequently, we have high confidence that the gridded
// fire severity maps produced with this updated code are not affected.
// Thanks to Zack Steel for alerting us about this issue and suggesting an alternative.
//
// Also, we found an error in the RdNBR equation.
// Instead of using pre-fire NBR^0.5 in the denominator, we used pre-fire NBR^0.25. This is now fixed.
// RdNBR was a bonus feature of this code and it was not used to produce CBI or even described in the associated paper,
// Consequently, this error does not affect the findings of the paper.
// Be aware of this change if you used the previous version to produce RdNBR.
// Thanks to Alison Paulson for alerting us about this issue.
///////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
///////////////            UPDATE 12/16/2021             /////////////////////////////
//  We updated the Landsat imagery to import Collection 2, Level 2 (SR), Tier 1
//  and to use new rescaling and mask parameters

///////////////            UPDATE 10/13/2022             /////////////////////////////
//  We updated the Landsat imagery to import Landsat 9 from Collection 2, Level 2 (SR), Tier 1

/////////////////////////////////////////////////////////////////////////////////////
//                        INPUTS                                                   //
/////////////////////////////////////////////////////////////////////////////////////
//--------------------       FIRE PERIMETER INPUT       ---------------------------//
// Import shapefile with fire polygons. The shapefile must have the following attributes with these exact names:
//       NAME           DESCRIPTION
//       Fire_ID        unique identifier for each fire
//       Fire_Year      year of fire
//       Start_Day      start day of fire season in julian days, e.g. June 15 = 166
//       End_Day        end day of fire season in julian days
var fires = ee.FeatureCollection("users/jl104/gee_perimeters");

//----------        IMAGERY SELECTION        ---------------------------//
// Select the imagery to export from the suite of available indices below.
// Add the VARIABLE NAMES (as desired) to brackets below, using quotes, e.g.  ['CBI', 'CBI_bc', 'dnbr', 'rbr']

//    VARIABLE NAME     DESCRIPTION
//    CBI               Composite Burn Index
//    CBI_bc            Bias-corrected Composite Burn Index
//    dnbr              delta normalized burn ratio
//    rbr               relativized burn ratio
//    rdnbr             relativized delta normalized burn ratio
///   dndvi             delta normalized differenced vegetation index
//    devi              delta enhanced vegetation index
//    dndmi             delta normalized difference moisture index
//    dmirbi            delta mid-infrared bi-spectral index
//    post_nbr          post-fire normalized burn ratio
//    post_mirbi        post-fire mid-infrared bi-spectral index
var bandsToExport      = ['CBI_bc']


/////////////////////////////////////////////////////////////////////////////////////
//                        END  OF INPUTS                                           //
/////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
//----------------  RANDOM FOREST SPECIFICATIONS---------------------//
// Bands for the random forest classification
var rf_bands = ['def', 'lat',  'rbr', 'dmirbi', 'dndvi', 'post_mirbi'];

//Load training data for Random Forest classification
var cbi = ee.FeatureCollection("users/grumpyunclesean/CBI_predictions/data_for_ee_model");

//Load climatic water deficit variable (def) for random forest classification
var def = ee.Image("users/grumpyunclesean/CBI_predictions/def").rename('def').toInt();

//Create latitude image for random forest classification
var lat = ee.Image.pixelLonLat().select('latitude').rename('lat').round().toInt();

// Parameters for random forest classification
var nrow_training_fold = cbi.size();    // number of training observations
var minLeafPopulation  = nrow_training_fold.divide(75).divide(6).round();

// Random forest classifier
var fsev_classifier = ee.Classifier.smileRandomForest(
    {numberOfTrees: 500,
    minLeafPopulation: minLeafPopulation,
    seed: 123
    })
  .train(cbi, 'CBI', rf_bands)
  .setOutputMode('REGRESSION');

////////////////////////////////////////////////////////////////////////////////////
//--------------------     PROCESSING     ----------------------------//
//-------- Initialize variables for fire perimeters  -----------------//
// create list with fire IDs
var fireID    = ee.List(fires.aggregate_array('Fire_ID')).getInfo();
var nFires = fireID.length;

//  Suite of spectral indices available for export.
var bandList      = ['dnbr', 'rbr', 'rdnbr', 'dndvi', 'devi', 'dndmi', 'dmirbi', 'post_nbr', 'post_mirbi', 'CBI', 'CBI_bc'];

//////////////////////////////////
//  GET LANDSAT COLLECTIONS
//////////////////////////////////
// Landsat 5, 7, 8 and 9 Surface Reflectance (Level 2) Tier 1 Collection 2
var ls9SR = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'),
    ls8SR = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'),
    ls7SR = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'),
    ls5SR = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2'),
    ls4SR = ee.ImageCollection('LANDSAT/LT04/C02/T1_L2');

/////////////////////////////////////////
// FUNCTIONS TO CREATE SPECTRAL INDICES
/////////////////////////////////////////
// Apply scaling factors.
var applyScaleFactors=function(lsImage) {
  var opticalBands=lsImage.select(["SR_B."]).multiply(0.0000275).add(-0.2);
  return lsImage.addBands(opticalBands, null, true)
}

// Returns vegetation indices for LS8 and LS9
var ls8_9_Indices = function(lsImage){
  var nbr  = lsImage.normalizedDifference(['SR_B5', 'SR_B7']).toFloat();
  var ndvi = lsImage.normalizedDifference(['SR_B5', 'SR_B4']).toFloat();
  var ndmi = lsImage.normalizedDifference(['SR_B5', 'SR_B6']).toFloat();
  var evi  = lsImage.expression(
              '2.5 * ((SR_B5 - SR_B4) / (SR_B5 + 6 * SR_B4 - 7.5 * SR_B2 + 1))',
              {'SR_B5': lsImage.select('SR_B5'),
              'SR_B4': lsImage.select('SR_B4'),
              'SR_B2': lsImage.select('SR_B2')}).toFloat();
  var mirbi = lsImage.expression(
              '((10 * SR_B6) - (9.8 * SR_B7) + 2)',
              {'SR_B6': lsImage.select('SR_B6'),
               'SR_B7': lsImage.select('SR_B7'),
              }).toFloat();
  var qa = lsImage.select(['QA_PIXEL']);
  return nbr.addBands([ndvi,ndmi,evi,mirbi,qa])
          .select([0,1,2,3,4,5], ['nbr','ndvi','ndmi','evi','mirbi','QA_PIXEL'])
          .copyProperties(lsImage, ['system:time_start']);

  };

// Returns indices for LS4, LS5 and LS7
var ls4_7_Indices = function(lsImage){
  var nbr  = lsImage.normalizedDifference(['SR_B4', 'SR_B7']).toFloat();
  var ndvi = lsImage.normalizedDifference(['SR_B4', 'SR_B3']).toFloat();
  var ndmi = lsImage.normalizedDifference(['SR_B4', 'SR_B5']).toFloat();
  var evi = lsImage.expression(
              '2.5 * ((SR_B4 - SR_B3) / (SR_B4 + 6 * SR_B3 - 7.5 * SR_B1 + 1))',
              {'SR_B4': lsImage.select('SR_B4'),
               'SR_B3': lsImage.select('SR_B3'),
               'SR_B1': lsImage.select('SR_B1')}).toFloat();
  var mirbi = lsImage.expression(
              '((10 * SR_B5) - (9.8 * SR_B7) + 2)',
              {'SR_B5': lsImage.select('SR_B5'),
               'SR_B7': lsImage.select('SR_B7'),
              }).toFloat();
  var qa = lsImage.select(['QA_PIXEL']);
  return nbr.addBands([ndvi,ndmi,evi,mirbi,qa])
          .select([0,1,2,3,4,5], ['nbr','ndvi','ndmi','evi','mirbi','QA_PIXEL'])
          .copyProperties(lsImage, ['system:time_start']);
  };

/////////////////////////////////////////////////
// FUNCTION TO MASK CLOUD, WATER, SNOW, ETC.
////////////////////////////////////////////////
var lsCfmask = function(lsImage) {
	// Bits 3,4,5,7: cloud,cloud-shadow,snow,water respectively.
	var cloudsBitMask      = (1 << 3);
	var cloudShadowBitMask = (1 << 4);
	var snowBitMask        = (1 << 5);
	var waterBitMask       = (1 << 7);
	// Get the pixel QA band.
	var qa = lsImage.select('QA_PIXEL');
	// Flags should be set to zero, indicating clear conditions.
  var clear =  qa.bitwiseAnd(cloudsBitMask).eq(0)
              .and(qa.bitwiseAnd(cloudShadowBitMask).eq(0))
              .and(qa.bitwiseAnd(snowBitMask).eq(0))
              .and(qa.bitwiseAnd(waterBitMask).eq(0));
	return lsImage.updateMask(clear).select([0,1,2,3,4])
              .copyProperties(lsImage,["system:time_start"]);
}

// Create water mask from Hansen's Global Forest Change to use in processing function
var waterMask = ee.Image('UMD/hansen/global_forest_change_2015').select(['datamask']).eq(1);

/////////////////////////////////////////////////
// RUN FUNCTIONS ON LANDSAT COLLECTION
////////////////////////////////////////////////
var ls9 = ls9SR.map(applyScaleFactors)
                .map(ls8_9_Indices)
                .map(lsCfmask);
var ls8 = ls8SR.map(applyScaleFactors)
               .map(ls8_9_Indices)
               .map(lsCfmask);
var ls7 = ls7SR.map(applyScaleFactors)
               .map(ls4_7_Indices)
               .map(lsCfmask);
var ls5 = ls5SR.map(applyScaleFactors)
               .map(ls4_7_Indices)
               .map(lsCfmask);
var ls4 = ls4SR.map(applyScaleFactors)
               .map(ls4_7_Indices)
               .map(lsCfmask);

// Merge Landsat Collections
var lsCol = ee.ImageCollection(ls9.merge(ls8).merge(ls7).merge(ls5).merge(ls4));

///////////////////////////////////////////////////////////////////////////////////////////////
// ------------------ Create Spectral Imagery for each fire -----------------//
var indices = ee.ImageCollection(fires.map(function(ft){
  // use 'Fire_ID' as unique identifier
  var fName    = ft.get("Fire_ID");

  // select fire
  var fire = ft;
  var fireBounds = ft.geometry().bounds();

  // create pre- and post-fire imagery
  var fireYear = ee.Date.parse('YYYY', fire.get('Fire_Year'));
  var startday = ee.Number.parse(fire.get('Start_Day'));
  var endday   = ee.Number.parse(fire.get('End_Day'));

  // Pre-Imagery
  var preFireYear = fireYear.advance(-1, 'year');
  //    check if imagery is available for time window; if so, get means across pixels; otherwise, output imagery will be masked
  var preFireIndices  = ee.Algorithms.If(lsCol.filterBounds(fireBounds)
                          .filterDate(preFireYear, fireYear)
                          .filter(ee.Filter.dayOfYear(startday, endday))
                          .size(),
                          lsCol.filterBounds(fireBounds)
                            .filterDate(preFireYear, fireYear)
                            .filter(ee.Filter.dayOfYear(startday, endday))
                            .mean()
                            .select([0,1,2,3,4], ['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi']),
                          ee.Image.cat(ee.Image(),ee.Image(),ee.Image(),ee.Image(),ee.Image())
                            .rename(['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi'])
                          )

  //    if any pixels within fire have only one 'scene' or less, add additional year backward to fill in behind
  //    as above, check if imagery is available for time window; if so, get means across pixels; otherwise, output imagery will be masked
  var preFireYear2 = fireYear.advance(-2, 'year');
  var preFireIndices2  = ee.Algorithms.If(lsCol.filterBounds(fireBounds)
                          .filterDate(preFireYear2, fireYear)
                          .filter(ee.Filter.dayOfYear(startday, endday))
                          .size(),
                          lsCol.filterBounds(fireBounds)
                            .filterDate(preFireYear2, fireYear)
                            .filter(ee.Filter.dayOfYear(startday, endday))
                            .mean()
                            .select([0,1,2,3,4], ['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi']),
                          ee.Image.cat(ee.Image(),ee.Image(),ee.Image(),ee.Image(),ee.Image())
                            .rename(['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi'])
                          )

  var pre_filled=ee.Image(preFireIndices).unmask(preFireIndices2);

  // Post-Imagery
  var postFireYear = fireYear.advance(1, 'year');

  //    check if imagery is available for time window; if so, get means across pixels; otherwise, output imagery will be masked
  var postFireIndices  = ee.Algorithms.If(lsCol.filterBounds(fireBounds)
                           .filterDate(postFireYear, fireYear.advance(2, 'year'))
                          .filter(ee.Filter.dayOfYear(startday, endday))
                          .size(),
                          lsCol.filterBounds(fireBounds)
                            .filterDate(postFireYear, fireYear.advance(2, 'year'))
                            .filter(ee.Filter.dayOfYear(startday, endday))
                            .mean()
                            .select([0,1,2,3,4], ['post_nbr','post_ndvi', 'post_ndmi', 'post_evi','post_mirbi']),
                          ee.Image.cat(ee.Image(),ee.Image(),ee.Image(),ee.Image(),ee.Image())
                            .rename(['post_nbr','post_ndvi', 'post_ndmi', 'post_evi','post_mirbi'])
                          )
  //    if any pixels within fire have only one 'scene' or less, add additional year forward to fill in behind
  //    as above, check if imagery is available for time window; if so, get means across pixels; otherwise, output imagery will be masked
  var postFireIndices2  = ee.Algorithms.If(lsCol.filterBounds(fireBounds)
                          .filterDate(postFireYear, fireYear.advance(3, 'year'))
                          .filter(ee.Filter.dayOfYear(startday, endday))
                          .size(),
                          lsCol.filterBounds(fireBounds)
                            .filterDate(postFireYear, fireYear.advance(3, 'year'))
                            .filter(ee.Filter.dayOfYear(startday, endday))
                            .mean()
                            .select([0,1,2,3,4], ['post_nbr','post_ndvi', 'post_ndmi', 'post_evi','post_mirbi']),
                          ee.Image.cat(ee.Image(),ee.Image(),ee.Image(),ee.Image(),ee.Image())
                            .rename(['post_nbr','post_ndvi', 'post_ndmi', 'post_evi','post_mirbi'])
                          )
  var post_filled = ee.Image(postFireIndices).unmask(postFireIndices2)
  var fireIndices = pre_filled.addBands(post_filled);

  // calculate dNBR
  var burnIndices = fireIndices.expression(
              "(b('pre_nbr') - b('post_nbr')) * 1000")
              .rename('dnbr').toInt().addBands(fireIndices);

  // calculate RBR
  var burnIndices2 = burnIndices.expression(
            "b('dnbr') / (b('pre_nbr') + 1.001)")
            .rename('rbr').toInt().addBands(burnIndices);

  // calculate RdNBR
   var burnIndices3 = burnIndices2.expression(
            "abs(b('pre_nbr')) < 0.001 ? 0.001" +
            ": b('pre_nbr')")
            .abs().sqrt().rename('pre_nbr2').toFloat().addBands(burnIndices2);

  var burnIndices4 = burnIndices3.expression(
            "b('dnbr') / b('pre_nbr2')")
            .rename('rdnbr').toInt().addBands(burnIndices3);

  // calculate dNDVI
  var burnIndices5 = burnIndices4.expression(
              "(b('pre_ndvi') - b('post_ndvi')) * 1000")
              .rename('dndvi').toInt().addBands(burnIndices4);

  // calculate dEVI
  var burnIndices6 = burnIndices5.expression(
              "(b('pre_evi') - b('post_evi')) * 1000")
              .rename('devi').toInt().addBands(burnIndices5);

   // calculate dNDMI
  var burnIndices7 = burnIndices6.expression(
              "(b('pre_ndmi') - b('post_ndmi')) * 1000")
              .rename('dndmi').toInt().addBands(burnIndices6);

   // calculate dMIRBI
  var burnIndices8 = burnIndices7.expression(
              "(b('pre_mirbi') - b('post_mirbi')) * 1000")
              .rename('dmirbi').toInt().addBands(burnIndices7);

  // Multiply post_mirbi band by 1000 to put it on the same scale as CBI plot extractions
  var post_mirbi_1000 = burnIndices8.select("post_mirbi").multiply(1000).toInt();
  burnIndices8 = burnIndices8.addBands(post_mirbi_1000, null, true); // null to copy all bands; true to overwrite original post_mirbi

 //  add in climatic water deficit variable, i.e. def
 var burnIndices9 = burnIndices8.addBands(def);

  //  add in latitude
 var burnIndices10 = burnIndices9.addBands(lat);

 // Classify the image with the same bands used to train the Random Forest classifier.
  var cbi_rf = burnIndices10.select(rf_bands).classify(fsev_classifier).
            rename('CBI').toFloat()
            .multiply(Math.pow(10,2)).floor().divide(Math.pow(10,2));  // set precision to two decimal places

  var burnIndices11 = cbi_rf.addBands(burnIndices10)

   // Create bias corrected CBI
   var bias_correct = function(bandName) {
      var cbi_lo = bandName.expression("((b('CBI') - 1.5) * 1.3)  + 1.5");
      var cbi_hi = bandName.expression("((b('CBI') - 1.5) * 1.175) + 1.5");
      var cbi_mg = bandName.where(bandName.lte(1.5),cbi_lo).where(bandName.gt(1.5),cbi_hi);
      return cbi_mg.where(cbi_mg.lt(0),0).where(cbi_mg.gt(3),3)
                  .multiply(Math.pow(10,2)).floor().divide(Math.pow(10,2)) // set precision to two decimal places
                  .rename('CBI_bc');
  };

  var validMask   = pre_filled.select('pre_nbr').add(10).add(post_filled.select('post_nbr')).add(10) // Adding 10 ensures resulting image doesn't have 0 values that would become masked in final output bands
  var mask = waterMask.updateMask(validMask)
  var burnIndices12 = bias_correct(burnIndices11.select('CBI')).addBands(burnIndices11)
  burnIndices12 = burnIndices12.select(bandList,bandList)
              .updateMask(mask)
              .clip(fire);

  return burnIndices12.set({
                        'fireID' : ft.get('Fire_ID'),
                        'fireYear' : ft.get('Fire_Year')
  });
}));

///////////////////////////////////////////////////////////////////////////////////////////
// ----------------------          Visualize Results               ----------------------//
Map.centerObject(fires,3);
var paletteclass = ["#0000CD","#6B8E23", "#FFFF00","#FFA500","#FF0000" ];
Map.addLayer(fires, {color: 'Red'}, "Fire perimeters");
Map.addLayer(indices.select('CBI_bc'), {min:0, max:3, palette:paletteclass}, 'CBI bias corrected');
Map.addLayer(indices.select('CBI'),    {min:0, max:3, palette:paletteclass}, 'CBI');

///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
// ----------------------   Export CBI and Other Spectral Imagery      ----------------------//
var nBands = bandsToExport.length

for (var j = 0; j < nFires; j++){
  var id   = fireID[j];
  var Name = id;
  var fireExport = ee.Image(indices.filterMetadata('fireID', 'equals', id).first());
  var fireBounds = ee.Feature(fires.filterMetadata('Fire_ID', 'equals', id).first()).geometry().bounds();

  for (var i = 0; i < nBands; i++) {
    var bandExport = bandsToExport[i];
    var exportImg = fireExport.select(bandExport);
    Export.image.toDrive({
      image: exportImg.toFloat(),  //casting to float maintains NA values in masked pixels
      description: Name + '_' + bandExport,
      fileNamePrefix: Name + '_' + bandExport,
      maxPixels: 1e13,
      scale: 30,
      crs: "EPSG:4326",
      folder: 'cbi',
      region: fireBounds
  });
}}

//////////
// DONE //
//////////
