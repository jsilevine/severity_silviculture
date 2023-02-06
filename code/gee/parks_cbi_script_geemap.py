import ee 
import math
from ee_plugin import Map

##############    WARNING: Use at your own risk     #################
#
# Script to develop and export maps of Composite Burn Index (CBI) predictions, and spectral imagery that
# supports the CBI predictions, for fires in North America.
# For updates or questions on this code, please can contact:
#        Lisa Holsinger, lisa.holsinger@usda.gov
#        Sean Parks,     sean.parks@usda.gov


###########     IMPORTANT - PLEASE READ       ################
#  1) Create a folder called 'cbi' in your Google Drive account before running script.
#  2) Click the ‘Run’ button (above); then under "Tasks" tab (left), click ‘Run’ buttons to export each image
#  3) Output files are named by concatenating the Fire_ID (see description below) and the relevant Earth Engine output.
#     For example: MT4751811328120110719_CBI or MT4751811328120110719_CBI_bc
#  4) If there is no imagery available for pre or post-fire period (e.g. early 1980s may lack data)
#     the resulting CBI, CBI_bc and other imagery will have NO DATA.
#  5) If you use this script or any of the products produced with the script, please cite the following paper:
#            Parks SA, Holsinger LM, Koontz MJ, Collins, Whitman E, et al. 2019. Giving ecological meaning to satellite-derived
#            fire severity metrics across North American forests within Google Earth Engine. Remote Sensing 2019, 11.


########            UPDATE 4/8/2021             #################/
# The original Random Forest classifier used in this study and associated code was deprecated on 3/1/2021.
# We have updated the code to use the new Random Forest classifier (smileRandomForest).
# We re-ran the 5-fold cross-validation using the new classifier and the resulting statistics originally
# reported are identical (cross-validated R2 = 0.72). Consequently, we have high confidence that the gridded
# fire severity maps produced with this updated code are not affected.
# Thanks to Zack Steel for alerting us about this issue and suggesting an alternative.
#
# Also, we found an error in the RdNBR equation.
# Instead of using pre-fire NBR^0.5 in the denominator, we used pre-fire NBR^0.25. This is now fixed.
# RdNBR was a bonus feature of this code and it was not used to produce CBI or even described in the associated paper,
# Consequently, this error does not affect the findings of the paper.
# Be aware of this change if you used the previous version to produce RdNBR.
# Thanks to Alison Paulson for alerting us about this issue.
###################################################/


##########################################/
#                        INPUTS                                                   #
##########################################/
#--------------------       FIRE PERIMETER INPUT       ---------------------------#
# Import shapefile with fire polygons. The shapefile must have the following attributes with these exact names:
#       NAME           DESCRIPTION
#       Fire_ID        unique identifier for each fire
#       Fire_Year      year of fire
#       Start_Day      start day of fire season in julian days, e.g. June 15 = 166
#       End_Day        end day of fire season in julian days
fires = ee.FeatureCollection("~/")

#----------        IMAGERY SELECTION        ---------------------------#
# Select the imagery to export from the suite of available indices below.
# Add the VARIABLE NAMES (as desired) to brackets below, using quotes, e.g.  ['CBI', 'CBI_bc','dnbr', 'rbr']

#    VARIABLE NAME     DESCRIPTION
#    CBI               Composite Burn Index
#    CBI_bc            Bias-corrected Composite Burn Index
#    dnbr              delta normalized burn ratio
#    rbr               relativized burn ratio
#    rdnbr             relativized delta normalized burn ratio
#/   dndvi             delta normalized differenced vegetation index
#    devi              delta enhanced vegetation index
#    dndmi             delta normalized difference moisture index
#    dmirbi            delta mid-infrared bi-spectral index
#    post_nbr          post-fire normalized burn ratio
#    post_mirbi        post-fire mid-infrared bi-spectral index
bandsToExport      = ['CBI', 'CBI_bc', 'rbr']

##########################################/
#                        END  OF INPUTS                                           #
##########################################/

##########################################
#----------------  RANDOM FOREST SPECIFICATIONS---------------------#
# Bands for the random forest classification
rf_bands = ['def', 'lat',  'rbr', 'dmirbi', 'dndvi', 'post_mirbi']

#Load training data for Random Forest classification
cbi = ee.FeatureCollection("users/grumpyunclesean/CBI_predictions/data_for_ee_model")

#Load climatic water deficit variable (def) for random forest classification
def = ee.Image("users/grumpyunclesean/CBI_predictions/def").rename('def').toInt()

#Create latitude image for random forest classification
lat = ee.Image.pixelLonLat().select('latitude').rename('lat').round().toInt()

# Parameters for random forest classification
nrow_training_fold = cbi.size();    # number of training observations
minLeafPopulation  = nrow_training_fold.divide(75).divide(6).round()

# Random forest classifier
fsev_classifier = ee.Classifier.smileRandomForest(
    {'numberOfTrees': 500,
    'minLeafPopulation': minLeafPopulation,
    'seed': 123
    }) \
  .train(cbi, 'CBI', rf_bands) \
  .setOutputMode('REGRESSION')

##########################################
#--------------------     PROCESSING     ----------------------------#
#-------- Initialize variables for fire perimeters  -----------------#
# create list with fire IDs
fireID    = ee.List(fires.aggregate_array('Fire_ID')).getInfo()
nFires = fireID.length

#  Suite of spectral indices available for export.
bandList      = ['dnbr', 'rbr', 'rdnbr', 'dndvi', 'devi', 'dndmi', 'dmirbi', 'post_nbr', 'post_mirbi', 'CBI', 'CBI_bc']

#################
#  GET LANDSAT COLLECTIONS
#################
# Landsat 4, 5, 7, and 8 Surface Reflectance Tier 1 collections
ls8SR = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR'),
    ls7SR = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR'),
    ls5SR = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR'),
    ls4SR = ee.ImageCollection('LANDSAT/LT04/C01/T1_SR')

####################/
# FUNCTIONS TO CREATE SPECTRAL INDICES
####################/
# Returns indices for LS8
def ls8_Indices(lsImage):
  nbr  = lsImage.normalizedDifference(['B5', 'B7']).toFloat()
  ndvi = lsImage.normalizedDifference(['B5', 'B4']).toFloat()
  ndmi = lsImage.normalizedDifference(['B5', 'B6']).toFloat()
  evi  = lsImage.expression(
              '(2.5 * ((B5 - B4) / (B5 + (6 * B4) - (7.5 * B2) + 1)))',
              {'B5': lsImage.select('B5').multiply(0.0001),
              'B4': lsImage.select('B4').multiply(0.0001),
              'B2': lsImage.select('B2').multiply(0.0001),
              }).toFloat()
  mirbi = lsImage.expression(
              '((10 * B6) - (9.8 * B7) + 2)',
              {'B6': lsImage.select('B6').multiply(0.0001),
              'B7': lsImage.select('B7').multiply(0.0001),
              }).toFloat()
  qa = lsImage.select(['pixel_qa'])
  return nbr.addBands([ndvi,ndmi,evi,mirbi,qa]) \
          .select([0,1,2,3,4,5], ['nbr','ndvi','ndmi','evi','mirbi','pixel_qa']) \
          .copyProperties(lsImage, ['system:time_start'])

  

# Returns indices for LS4, LS5 and LS7
def ls4_7_Indices(lsImage):
  nbr = lsImage.normalizedDifference(['B4', 'B7']).toFloat()
  ndvi = lsImage.normalizedDifference(['B4', 'B3']).toFloat()
  ndmi = lsImage.normalizedDifference(['B4', 'B5']).toFloat()
  evi = lsImage.expression(
              '(2.5 * ((B4 - B3) / (B4 + (6 * B3) - (7.5 * B1) + 1)))',
              {'B4': lsImage.select('B4').multiply(0.0001),
              'B3': lsImage.select('B3').multiply(0.0001),
              'B1': lsImage.select('B1').multiply(0.0001),
              }).toFloat()
  mirbi = lsImage.expression(
              '((10 * B5) - (9.8 * B7) + 2)',
              {'B5': lsImage.select('B5').multiply(0.0001),
              'B7': lsImage.select('B7').multiply(0.0001),
              }).toFloat()
  qa = lsImage.select(['pixel_qa'])
  return nbr.addBands([ndvi,ndmi,evi,mirbi,qa]) \
          .select([0,1,2,3,4,5], ['nbr','ndvi','ndmi','evi','mirbi','pixel_qa']) \
          .copyProperties(lsImage, ['system:time_start'])
  

########################/
# FUNCTION TO MASK CLOUD, WATER, SNOW, ETC.
########################
def lsCfmask(lsImg):
  quality =lsImg.select(['pixel_qa'])
  clear = quality.bitwiseAnd(8).eq(0)             # cloud shadow \
                .And(quality.bitwiseAnd(32).eq(0) \
                .And(quality.bitwiseAnd(4).eq(0) \
                .And(quality.bitwiseAnd(16).eq(0)))); 
  return lsImg.updateMask(clear).select([0,1,2,3,4]) \
            .copyProperties(lsImg, ['system:time_start'])


# Create water mask from Hansen's Global Forest Change to use in processing function
waterMask = ee.Image('UMD/hansen/global_forest_change_2015').select(['datamask']).eq(1)

########################/
# RUN FUNCTIONS ON LANDSAT COLLECTION
########################
ls8 = ls8SR.map(ls8_Indices) \
                .map(lsCfmask)
ls7 = ls7SR.map(ls4_7_Indices) \
                .map(lsCfmask)
ls5 = ls5SR.map(ls4_7_Indices) \
                .map(lsCfmask)
ls4 = ls4SR.map(ls4_7_Indices) \
                .map(lsCfmask)

# Merge Landsat Collections
lsCol = ee.ImageCollection(ls8.merge(ls7).merge(ls5).merge(ls4))

###############################################/
# ------------------ Create Spectral Imagery for each fire -----------------#

def func_fxf(ft):
  # use 'Fire_ID' as unique identifier
  fName    = ft.get("Fire_ID")

  # select fire
  fire = ft
  fireBounds = ft.geometry().bounds()

  # create pre- and post-fire imagery
  fireYear = ee.Date.parse('YYYY', fire.get('Fire_Year'))
  startday = ee.Number.parse(fire.get('Start_Day'))
  endday   = ee.Number.parse(fire.get('End_Day'))

  # Pre-Imagery
  preFireYear = fireYear.advance(-1, 'year')
  #    check if imagery is available for time window; if so, get means across pixels; otherwise, output imagery will be masked
  preFireIndices  = ee.Algorithms.If(lsCol.filterBounds(fireBounds) \
                          .filterDate(preFireYear, fireYear) \
                          .filter(ee.Filter.dayOfYear(startday, endday)) \
                          .size(),
                          lsCol.filterBounds(fireBounds) \
                            .filterDate(preFireYear, fireYear) \
                            .filter(ee.Filter.dayOfYear(startday, endday)) \
                            .mean() \
                            .select([0,1,2,3,4], ['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi']),
                          ee.Image.cat(ee.Image(),ee.Image(),ee.Image(),ee.Image(),ee.Image()) \
                            .rename(['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi'])
                          )

  #    if any pixels within fire have only one 'scene' or less, add additional year backward to fill in behind
  #    as above, check if imagery is available for time window; if so, get means across pixels; otherwise, output imagery will be masked
  preFireYear2 = fireYear.advance(-2, 'year')
  preFireIndices2  = ee.Algorithms.If(lsCol.filterBounds(fireBounds) \
                          .filterDate(preFireYear2, fireYear) \
                          .filter(ee.Filter.dayOfYear(startday, endday)) \
                          .size(),
                          lsCol.filterBounds(fireBounds) \
                            .filterDate(preFireYear2, fireYear) \
                            .filter(ee.Filter.dayOfYear(startday, endday)) \
                            .mean() \
                            .select([0,1,2,3,4], ['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi']),
                          ee.Image.cat(ee.Image(),ee.Image(),ee.Image(),ee.Image(),ee.Image()) \
                            .rename(['pre_nbr','pre_ndvi', 'pre_ndmi', 'pre_evi','pre_mirbi'])
                          )

  pre_filled=ee.Image(preFireIndices).unmask(preFireIndices2)

  # Post-Imagery
  postFireYear = fireYear.advance(1, 'year')

  #    check if imagery is available for time window; if so, get means across pixels; otherwise, output imagery will be masked
  postFireIndices  = ee.Algorithms.If(lsCol.filterBounds(fireBounds) \
                           .filterDate(postFireYear, fireYear.advance(2, 'year')) \
                          .filter(ee.Filter.dayOfYear(startday, endday)) \
                          .size(),
                          lsCol.filterBounds(fireBounds) \
                            .filterDate(postFireYear, fireYear.advance(2, 'year')) \
                            .filter(ee.Filter.dayOfYear(startday, endday)) \
                            .mean() \
                            .select([0,1,2,3,4], ['post_nbr','post_ndvi', 'post_ndmi', 'post_evi','post_mirbi']),
                          ee.Image.cat(ee.Image(),ee.Image(),ee.Image(),ee.Image(),ee.Image()) \
                            .rename(['post_nbr','post_ndvi', 'post_ndmi', 'post_evi','post_mirbi'])
                          )
  #    if any pixels within fire have only one 'scene' or less, add additional year forward to fill in behind
  #    as above, check if imagery is available for time window; if so, get means across pixels; otherwise, output imagery will be masked
  postFireIndices2  = ee.Algorithms.If(lsCol.filterBounds(fireBounds) \
                          .filterDate(postFireYear, fireYear.advance(3, 'year')) \
                          .filter(ee.Filter.dayOfYear(startday, endday)) \
                          .size(),
                          lsCol.filterBounds(fireBounds) \
                            .filterDate(postFireYear, fireYear.advance(3, 'year')) \
                            .filter(ee.Filter.dayOfYear(startday, endday)) \
                            .mean() \
                            .select([0,1,2,3,4], ['post_nbr','post_ndvi', 'post_ndmi', 'post_evi','post_mirbi']),
                          ee.Image.cat(ee.Image(),ee.Image(),ee.Image(),ee.Image(),ee.Image()) \
                            .rename(['post_nbr','post_ndvi', 'post_ndmi', 'post_evi','post_mirbi'])
                          )
  post_filled = ee.Image(postFireIndices).unmask(postFireIndices2)
  fireIndices = pre_filled.addBands(post_filled)

  # calculate dNBR
  burnIndices = fireIndices.expression(
              "(b('pre_nbr') - b('post_nbr')) * 1000") \
              .rename('dnbr').toInt().addBands(fireIndices)

  # calculate RBR
  burnIndices2 = burnIndices.expression(
            "b('dnbr') / (b('pre_nbr') + 1.001)") \
            .rename('rbr').toInt().addBands(burnIndices)

  # calculate RdNBR
   burnIndices3 = burnIndices2.expression(
            "abs(b('pre_nbr')) < 0.001 ? 0.001" + \
            ": b('pre_nbr')") \
            .abs().sqrt().rename('pre_nbr2').toFloat().addBands(burnIndices2)

  burnIndices4 = burnIndices3.expression(
            "b('dnbr') / b('pre_nbr2')") \
            .rename('rdnbr').toInt().addBands(burnIndices3)

  # calculate dNDVI
  burnIndices5 = burnIndices4.expression(
              "(b('pre_ndvi') - b('post_ndvi')) * 1000") \
              .rename('dndvi').toInt().addBands(burnIndices4)

  # calculate dEVI
  burnIndices6 = burnIndices5.expression(
              "(b('pre_evi') - b('post_evi')) * 1000") \
              .rename('devi').toInt().addBands(burnIndices5)

   # calculate dNDMI
  burnIndices7 = burnIndices6.expression(
              "(b('pre_ndmi') - b('post_ndmi')) * 1000") \
              .rename('dndmi').toInt().addBands(burnIndices6)

   # calculate dMIRBI
  burnIndices8 = burnIndices7.expression(
              "(b('pre_mirbi') - b('post_mirbi')) * 1000") \
              .rename('dmirbi').toInt().addBands(burnIndices7)

  # Multiply post_mirbi band by 1000 to put it on the same scale as CBI plot extractions
  post_mirbi_1000 = burnIndices8.select("post_mirbi").multiply(1000).toInt()
  burnIndices8 = burnIndices8.addBands(post_mirbi_1000, None, True); # None to copy all bands; True to overwrite original post_mirbi

 #  add in climatic water deficit variable, i.e. def
 burnIndices9 = burnIndices8.addBands(def)

  #  add in latitude
 burnIndices10 = burnIndices9.addBands(lat)

 # Classify the image with the same bands used to train the Random Forest classifier.
  cbi_rf = burnIndices10.select(rf_bands).classify(fsev_classifier).
            rename('CBI').toFloat() \
            .multiply(math.pow(10,2)).floor().divide(math.pow(10,2));  

  burnIndices11 = cbi_rf.addBands(burnIndices10)

   # Create bias corrected CBI
   def bias_correct(bandName):
      cbi_lo = bandName.expression("((b('CBI') - 1.5) * 1.3)  + 1.5")
      cbi_hi = bandName.expression("((b('CBI') - 1.5) * 1.175) + 1.5")
      cbi_mg = bandName.where(bandName.lte(1.5),cbi_lo).where(bandName.gt(1.5),cbi_hi)
      return cbi_mg.where(cbi_mg.lt(0),0).where(cbi_mg.gt(3),3) \
                  .multiply(math.pow(10,2)).floor().divide(math.pow(10,2)) \
                  .rename('CBI_bc')
  

  validMask   = pre_filled.select('pre_nbr').add(10).add(post_filled.select('post_nbr')).add(10) # Adding 10 ensures resulting image doesn't have 0 values that would become masked in final output bands
  mask = waterMask.updateMask(validMask)
  burnIndices12 = bias_correct(burnIndices11.select('CBI')).addBands(burnIndices11)
  burnIndices12 = burnIndices12.select(bandList,bandList) \
              .updateMask(mask) \
              .clip(fire)

  return burnIndices12.set({
                        'fireID' : ft.get('Fire_ID'),
                        'fireYear' : ft.get('Fire_Year')
  })

indices = ee.ImageCollection(fires.map(func_fxf
))





























































































































































))

#############################################/
# ----------------------          Visualize Results               ----------------------#
Map.centerObject(fires,3)
paletteclass = ["#0000CD","#6B8E23", "#FFFF00","#FFA500","#FF0000" ]
Map.addLayer(fires, {'color': 'Red'}, "Fire perimeters")
Map.addLayer(indices.select('CBI_bc'), {'min':0, 'max':3, 'palette':paletteclass}, 'CBI bias corrected')
Map.addLayer(indices.select('CBI'),    {'min':0, 'max':3, 'palette':paletteclass}, 'CBI')

###############################################/

###############################################/
# ----------------------   Export CBI and Other Spectral Imagery      ----------------------#
nBands = bandsToExport.length

for j in range(0, nFires, 1):
  id   = fireID[j]
  Name = id
  fireExport = ee.Image(indices.filterMetadata('fireID', 'equals', id).first())
  fireBounds = ee.Feature(fires.filterMetadata('Fire_ID', 'equals', id).first()).geometry().bounds()

  for i in range(0, nBands, 1):
    bandExport = bandsToExport[i]
    exportImg = fireExport.select(bandExport)
    Export.image.toDrive({
      'image': exportImg.toFloat(),  #casting to float maintains NA values in masked pixels
      'description': Name + '_' + bandExport,
      'fileNamePrefix': Name + '_' + bandExport,
      'maxPixels': 1e13,
      'scale': 30,
      'crs': "EPSG:4326",
      'folder': 'cbi',
      'region': fireBounds
  })


###################################################
# DONE #
#####
Map