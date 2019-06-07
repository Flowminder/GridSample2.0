# -*- coding: utf-8 -*-
"""

#!/usr/bin/env python3

@author: icw
"""


import os
import osr

from osgeo import gdal, ogr, gdalconst
import psycopg2

import numpy as np

import sample_builder_utils as sbu

import subprocess

import csv


import matplotlib
# set matplotlib to not use Xwindows backend
matplotlib.use('Agg')
import matplotlib.pyplot as plt


import matplotlib.cm as cm
enable_debug = False

########################################################

# The functions here are principally for manipulating rasters and shapefiles
# They are a mixture of manipulations in "code" and calls to command line gdal functions

########################################################

def createColouredImageWithScale(inputTif,outputImage):

	if enable_debug:
		print '## Entering createColouredTifWithScale'
		print inputTif
		print outputImage
		
	#import matplotlib.pyplot as plt
	#import matplotlib.cm as cm
	
	im_r = gdal.Open(inputTif)
	im_r_band = im_r.GetRasterBand(1)
	img = np.array(im_r_band.ReadAsArray())
	
	my_cmap = cm.coolwarm
	my_cmap.set_under('w', alpha=0)
	
	figz  = plt.figure()

	fig = plt.imshow(img,cmap = my_cmap,clim=[0, None])
	
	cb =figz.colorbar(fig)
	
	plt.axis('off')
	
	figz.savefig(outputImage, bbox_inches='tight')

def convertRasterToPNGWithTransparentNDV(input_raster,output_png):

	if enable_debug:
		print '## Entering convertRasterToPNGWithTransparentNDV '
		print input_raster
		print output_png
		
	# return 0 for sucess, other for error
	
	if not os.path.exists(input_raster):
		return 1
	
	if os.path.exists(output_png):
		os.remove(output_png)
	
	# get original raster array and save array
	driver = gdal.GetDriverByName('GTiff')
	
	# get no data value
	ir = gdal.Open(input_raster,gdal.GA_Update)
	ndv= ir.GetRasterBand(1).GetNoDataValue()
	
	ir = None
	intermediate_raster_name = input_raster.replace(".tif","_intermediate.tif")
	
	# first add alpha channel for transparency
	gdal_cmd_as_string = 'gdalwarp ' + input_raster + " -dstalpha -srcnodata '" + str(ndv) + "' -co compress=LZW -overwrite " + intermediate_raster_name

	#gdal_cmd_as_string = 'gdalwarp -to SRC_METHOD=NO_GEOTRANSFORM ' + input_raster + " -dstalpha -srcnodata '" + str(ndv) + "' -co compress=LZW -overwrite " + intermediate_raster_name

	
	if enable_debug:
		print gdal_cmd_as_string
	
	ret = os.system(gdal_cmd_as_string)

	if ret == 0:
		# now convert to png
		gdal_cmd_as_string = 'gdal_translate -of PNG -scale ' + intermediate_raster_name + ' ' + output_png
		
		if enable_debug:
			print gdal_cmd_as_string
			
		ret = os.system(gdal_cmd_as_string)

	# if it exists remove the intermediate raster
	if os.path.exists(intermediate_raster_name):
		os.remove(intermediate_raster_name)
	
	return ret


def make_raster_copy(src_rast,new_filename,initalise_empty = False,initial_value= None,ndv=None,integer_pixels = False, inital_array= None):
	# makes a copy of a raster
	# it can be initialised empty
	# the no data value can be changed
	# an initial value to set any cells that have a value can be defined	

	if enable_debug:
		print '## Entering make_raster_copy ' 
		print 'src_rast' + src_rast
		print 'new_filename' + new_filename
		


	driver = gdal.GetDriverByName('GTiff')
	
	template_raster = gdal.Open(src_rast)
	
	cols = template_raster.RasterXSize
	rows = template_raster.RasterYSize
	projection = template_raster.GetProjection()
	geotransform = template_raster.GetGeoTransform()
	bands = template_raster.RasterCount

	pixel_type = gdal.GDT_Float32
	if integer_pixels == True:
		pixel_type = gdal.GDT_Int32
		
	new_raster = driver.Create(str(new_filename), cols, rows, bands, pixel_type, options = [ 'COMPRESS=LZW' ])
	new_raster.SetProjection(projection)
	new_raster.SetGeoTransform(geotransform)

	for i in range(bands):
		band = new_raster.GetRasterBand(i + 1)
		
		template_band=template_raster.GetRasterBand(i + 1)
		template_array = np.array(template_band.ReadAsArray())
		
		# either leave the no data value as it is or set it as required
		if ndv is None:
			ndv = template_band.GetNoDataValue()
		else:
			template_array[(template_array == template_band.GetNoDataValue() )] = ndv
		band.SetNoDataValue(ndv)
		
		if not initial_value is None:
			template_array[(template_array != ndv )] = initial_value
		
		if initalise_empty:
			band.Fill(ndv)
			
		if not inital_array is None:
			band.WriteArray(inital_array)
		else:
			band.WriteArray(template_array)
	
	new_raster.FlushCache()
	
	# close 
	template_raster = None
	new_raster = None

def makeShapefileFromDatabaseLayers(layer_ids,dest_dir,output_file_name='coverage_file'):	
	
	if enable_debug:
		print '## Entering makeShapefileFromDatabaseLayers ' 
		
	# using pgsql2shp to create a shapefile
	try:
		shapefile_name=output_file_name
		
		if enable_debug:
			print 'layer_ids:' + str(layer_ids)
		
		# ids are json encoded so remove  any square brackets
		layer_ids=str(layer_ids).replace('[','')
		layer_ids=layer_ids.replace(']','')
		
		pgsql2shp ='pgsql2shp -f %s -h db.worldpop.org.uk -p 5432 -u wpadmin -P wp-southville wpinputs "select * from grid_sample.gadm_all_phase2 where id in (%s)"' % ( shapefile_name , layer_ids )
		
		#print pgsql2shp
		
		os.chdir(dest_dir)
		os.system(pgsql2shp)
	except Exception as error:
		print str(error) 
		shapefile_name=None
		
	if enable_debug:
		print '## Leaving makeShapefileFromDatabaseLayers ' 
		
	return shapefile_name
			
	
def clipRasterWithRasterMask(mask_file,file_to_cut_from,output_file):

	if enable_debug:
		print '## Entering clipRasterWithRasterMask '
		
	# clip one raster with another - this assumes one lies within another and they have both the same spatial ref and have the pixels aligned

	# open mask file and get it's boundaries
	data = gdal.Open(mask_file)
	geoTransform = data.GetGeoTransform()
	minx = geoTransform[0]
	maxy = geoTransform[3]
	maxx = minx + geoTransform[1] * data.RasterXSize
	miny = maxy + geoTransform[5] * data.RasterYSize
	
	outsize = ' -outsize ' + str(data.RasterXSize) + ' ' + str(data.RasterYSize)
	
	# create a gdal_translate string to cut to a window from the source
	gdal_cmd_as_string = 'gdal_translate -projwin ' + ' '.join([str(x) for x in [minx, maxy, maxx, miny]]) + outsize + ' -co compress=LZW -of GTiff ' + file_to_cut_from + ' ' + output_file
	
	
	if enable_debug:
		print gdal_cmd_as_string
	
	return os.system(gdal_cmd_as_string)
	
def clipRasterWithRasterMask2(mask_file,file_to_cut_from,output_file):

	# this is just a bit of a hack to get aligned rasters
	if enable_debug:
		print '## Entering clipRasterWithRasterMask2 '
		
	# clip raster normally
	intermediate_output_file = output_file.replace(".tif","_intermediate.tif")
	ret = clipRasterWithRasterMask(mask_file,file_to_cut_from,intermediate_output_file)

	if ret == 0:
		print '##########################################'
		of = gdal.Open(intermediate_output_file,gdal.GA_Update)
		of_band = of.GetRasterBand(1)
		of_ndv = of_band.GetNoDataValue()
		of_arr = np.array(of_band.ReadAsArray())	
		
		# create new raster using
		make_raster_copy(mask_file,output_file,integer_pixels = True, inital_array= of_arr)
		
		#if os.path.exists(intermediate_output_file):
		#		os.remove(intermediate_output_file)
		
	
	return ret
	#return os.system(gdal_cmd_as_string)
	
def cropRasterWithShapefileMask(mask_file,file_to_cut_from,output_file,pixel_size = 0.000833333,other='',tap=True):

	if enable_debug:
		print '## Entering cropRasterWithShapefileMask '
		
	# crop one raster with a shapefile - this assumes one lies within another and they have both the same spatial ref and have the pixels aligned

	# create a gdal_translate string to cut to a window form the source
	#gdal_cmd_as_string = 'gdalwarp -tr ' + str(pixel_size) + ' -' + str(pixel_size) 
	
	
	#if tap:
	#	gdal_cmd_as_string +=  ' -tap  ' 
		
	#gdal_cmd_as_string  +=  ' -overwrite ' + other + ' -co compress=LZW -wo CUTLINE_ALL_TOUCHED=TRUE -crop_to_cutline -cutline ' + mask_file + ' ' + file_to_cut_from + ' ' + output_file
	gdal_cmd_as_string = 'gdalwarp '
	gdal_cmd_as_string  +=  ' -overwrite ' + other + ' -co compress=LZW -wo CUTLINE_ALL_TOUCHED=TRUE  -cutline ' + mask_file + ' ' + file_to_cut_from + ' ' + output_file
	
	if enable_debug:
		print gdal_cmd_as_string
	
	return os.system(gdal_cmd_as_string)

def reprojectToWGS84DegreesTargetAlignedPixel(file_to_reproject,reprojected_file_name,pixel_size = 0.000833333,tap=True,no_data_value = None):

	if enable_debug:
		print '## Entering reprojectToWGS84DegreesTargetAlignedPixel'
		
	# reproject raster to wgs 84, in degrees, targetting aligned pixels
	gdal_cmd_as_string = 'gdalwarp -overwrite -tr ' + str(pixel_size) + ' -' + str(pixel_size) 
	
	if tap:
		gdal_cmd_as_string += ' -tap '
	
	if not no_data_value is None:
		gdal_cmd_as_string +=  ' -dstnodata ' + str(no_data_value) + ' '
		
	gdal_cmd_as_string += ' -co compress=LZW -t_srs epsg:4326 ' + file_to_reproject + ' ' +  reprojected_file_name

	if enable_debug:
		print gdal_cmd_as_string
	
	return os.system(gdal_cmd_as_string)

def mergeShapefiles(first_file,second_file):
	
	if enable_debug:
		print '## Entering mergeShapefiles '
		
	# merges 2 shapefiles by adding the second s to the first - so the first file is changed 
	gdal_cmd_as_string = 'ogr2ogr -f "ESRI Shapefile" -update -append ' + first_file + ' ' +  second_file
	#gdal_cmd_as_string = 'ogr2ogr -f "ESRI Shapefile" -nlt MULTIPOLYGON -update -append ' + first_file + ' ' +  second_file
	
	if enable_debug:
		print gdal_cmd_as_string
	
	return os.system(gdal_cmd_as_string)
	

	
def reprojectShapefileToWGS84(file_to_reproject,reprojected_file_name):

	if enable_debug:
		print ("Entering reprojectShapefileToWGS84")	
	
	# reproject 
	gdal_cmd_as_string = 'ogr2ogr -f "ESRI Shapefile" -t_srs EPSG:4326 ' + reprojected_file_name + ' ' +  file_to_reproject
	#gdal_cmd_as_string = 'ogr2ogr -f "ESRI Shapefile"  ' + reprojected_file_name + ' ' +  file_to_reproject

	if enable_debug:
		print gdal_cmd_as_string
	
	return os.system(gdal_cmd_as_string)
	
def clipShapefilewithShapefile(file_to_clip,output_file_name,clip_src_file):

	# this seems to fail sometimes when the polys being interesected are at the edge of the shape 
	# clipShapefilewithShapefileUsingIntersection seems to have resolved this 
	if enable_debug:
		print ("Entering clipShapefilewithShapefile")	
	
	# clip 
	gdal_cmd_as_string = 'ogr2ogr -f "ESRI Shapefile" -clipsrc ' + clip_src_file + ' '  + output_file_name + ' '  + file_to_clip + ' -nlt MULTIPOLYGON -skipfailures'


	if enable_debug:
		print gdal_cmd_as_string
	
	# this generates error messages which we ignore as it still creates the shapefile so ignore the status - but what if it fails completely!
	os.system(gdal_cmd_as_string)
	
	#return os.system(gdal_cmd_as_string)
	return 0
	
def clipShapefilewithShapefileUsingIntersection(file_to_clip,output_file_name,clip_src_file,id_field_name,name_field_name):

	if enable_debug:
		print ("Entering clipShapefilewithShapefileUsingIntersection")	
		
	# get filenames without .shp
	clip_src_file_name = os.path.basename(clip_src_file)
	clip_src_file_name = clip_src_file_name.replace(".shp","")
	file_to_clip_name = os.path.basename(file_to_clip)
	file_to_clip_name = file_to_clip_name.replace(".shp","")
	
	# clip 
	gdal_cmd_as_string = 'ogr2ogr -f "ESRI Shapefile"  ' + output_file_name  + ' '  + file_to_clip + ' -dialect sqlite -sql "SELECT p.' + id_field_name + ',p.' + name_field_name + ',ST_Intersection(p.Geometry, t.Geometry) FROM '  + file_to_clip_name + ' p, \'' + clip_src_file + '\'.' + clip_src_file_name + ' t WHERE ST_Intersects(p.Geometry, t.Geometry) "  -nlt MULTIPOLYGON -skipfailures'

	# ogr2ogr -f 'ESRI Shapefile' /data/generatedFiles/4/output_clip_a25.shp /data/generatedFiles/4/reprojected_strata_poly.shp -dialect sqlite -sql "SELECT p.ID,p.NAME_2,ST_Intersection(p.Geometry, t.Geometry) FROM reprojected_strata_poly p, '/data/generatedFiles/4/coverage_poly_reprojected.shp'.coverage_poly_reprojected t WHERE ST_Intersects(p.Geometry, t.Geometry) " -nlt MULTIPOLYGON -skipfailures


	if enable_debug:
		print gdal_cmd_as_string
	
	# this generates error messages which we ignore as it still creates the shapefile so ignore the status - but what if it fails completely!
	os.system(gdal_cmd_as_string)
	
	#return os.system(gdal_cmd_as_string)
	return 0
	
def clipShapefilewithShapefilebyAttribute(file_to_clip,output_file_name,clip_src_file,where):

	if enable_debug:
		print ("Entering clipShapefilewithShapefilebyAttribute")	
	
	# clip 
	gdal_cmd_as_string = 'ogr2ogr -f "ESRI Shapefile" -t_srs EPSG:4326 -clipsrc ' + clip_src_file + ' -clipsrcwhere "' + where + '" '  + output_file_name + ' '  + file_to_clip + ' -nlt MULTIPOLYGON -skipfailures'

	if enable_debug:
		print gdal_cmd_as_string
	
	# this generates error messages which we ignore as it still creates the shapefile so ignore the status - but what if it fails completely!
	os.system(gdal_cmd_as_string)
	
	#return os.system(gdal_cmd_as_string)
	return 0
	
def resampleWorldPopRasterByScaleFactorTargetAlignedPixel_PROBABLY_NOT_IN_USE(file_to_reproject,reprojected_file_name):

	if enable_debug:
		print '## Entering resampleWorldPopRasterByScaleFactorTargetAlignedPixel_PROBABLY_NOT_IN_USE'
		
	# reproject raster to wgs 84, in degrees, targetting aligned pixles
	gdal_cmd_as_string = 'gdalwarp -overwrite -tr 0.000833333 -0.000833333 -tap -co compress=LZW -t_srs epsg:4326 ' + file_to_reproject + ' ' +  reprojected_file_name

	if enable_debug:
		print gdal_cmd_as_string
	
	return os.system(gdal_cmd_as_string)
	
def runVirusCheckOnFile(file_to_check):
	# this connects to a remote clamAV container
	
	# 1 contains virus or can't scan
	ret = 1
		
	cmd ='clamdscan -c /data/clamd.remote.conf --fdpass --stream ' + file_to_check
	
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	(output, err) = p.communicate()
 
	## Wait for check to terminate. Get return returncode  and any message
	p_status = p.wait()
	
	#print "Command output : ", output
	#print "Command exit status/return code : ", p_status
	
	# check it all worked Ok by looking for  a string made up of the file name and OK in hte output
	sucess_string = file_to_check + ": OK"
	if not p_status and sucess_string in output:
		# no virus
		ret = 0
	
	#print "ret" + str(ret)
	
	return ret

# this create a boundary polyline then creates a buffer polygon
def get_shapefile_info(source_file):

	if enable_debug:
		print "Entering get_shapefile_info"
		
	# get the input layer and get its spatial reference
	driver = ogr.GetDriverByName("ESRI Shapefile")

	# 0 is read only
	sf = driver.Open(source_file,0)

	bnd = None
	
	results = []
	if sf:
		layer = sf.GetLayer()
		srs = layer.GetSpatialRef()
		
		featureCount = layer.GetFeatureCount()
		print "Number of features in shapefile %d" % (featureCount)
	  
		layerDefinition = layer.GetLayerDefn()

		# get attribute names
		for i in range(layerDefinition.GetFieldCount()):
			fieldName = layerDefinition.GetFieldDefn(i).GetName()
			print layerDefinition.GetFieldDefn(i).GetName()
			results.append()

		feat = layer.GetNextFeature()

		while feat is not None:
			print "\nID:" + str(feat.GetField("ID"))
			print "NAME_1:" + feat.GetField("NAME_1")
			print "SHAPE_AREA:" + str(feat.GetField("SHAPE_AREA"))
				
			feat = layer.GetNextFeature()
			print "feat:" + str(feat)
	
	return results

# get the names of hte shapefile attributes
def get_shapefile_attributes(source_file):

	if enable_debug:
		print "Entering get_shapefile_attributes"
		
	# get the input layer and get its spatial reference
	driver = ogr.GetDriverByName("ESRI Shapefile")

	# 0 is read only
	sf = driver.Open(source_file,0)

	bnd = None
	
	results = []
	if sf:
		layer = sf.GetLayer()
		srs = layer.GetSpatialRef()
		
		featureCount = layer.GetFeatureCount()
		#print "Number of features in shapefile %d" % (featureCount)
	  
		layerDefinition = layer.GetLayerDefn()

		# get attribute names
		for i in range(layerDefinition.GetFieldCount()):
			fieldName = layerDefinition.GetFieldDefn(i).GetName()
			results.append(fieldName)	
	
	return results
	
def get_shapefile_values(source_file,arr_field_names= None, id_field_name = None):
	
	if enable_debug:
		print '## Entering get_shapefile_values '
		
		print 'source_file:' + source_file
		
		if id_field_name:
			print 'id_field_name:##' + id_field_name + '##'
		
	# return a list of the values found for each featureCount
	# if no field names given then return all
	# else return just hte required field values
	
	# get the input layer and get its spatial reference
	driver = ogr.GetDriverByName("ESRI Shapefile")

	# 0 is read only
	sf = driver.Open(source_file,0)
 		
	bnd = None
	
	# if an id_field_name is given then we return a dictionary whose key is the id field value
	# otherwise a list of dictionaries - effectively an associative array
	results = []
	if id_field_name:
		results = {}
		
	if sf:
		layer = sf.GetLayer()
	
		layerDefinition = layer.GetLayerDefn()
	
		feat = layer.GetNextFeature()
		
		# for each feature
		while feat is not None:
			# for each attribute
			row = {}
			for i in range(layerDefinition.GetFieldCount()):
				attribute_name =layerDefinition.GetFieldDefn(i).GetName()
				
				#print "\n" + attribute_name + "\n" + str(feat.GetField(attribute_name))
				
				# if there is a list of required field names and attibute is in this list then add to the results
				# if there is no list then include all
				if arr_field_names and attribute_name in arr_field_names:
					row[attribute_name] = feat.GetField(attribute_name)
				elif arr_field_names is None:
					row[attribute_name] = feat.GetField(attribute_name)
			
			if id_field_name:
				results[row[id_field_name]] = row
				#results[feat.GetField(id_field_name)] = row
			else:
				results.append(row)
				
			feat = layer.GetNextFeature()
			#print "feat:" + str(feat)
			
	return results
	
def get_shapefile_values_as_list(source_file,field_name ):
	
	if enable_debug:
		print '##***** Entering get_shapefile_values_as_list '
		
		print 'source_file:' + source_file
		print 'field_name:' + field_name + "@@"
		

	# get the input layer and get its spatial reference
	driver = ogr.GetDriverByName("ESRI Shapefile")

	# 0 is read only
	sf = driver.Open(source_file,0)
 		
	bnd = None
	
	# if an id_field_name is given then we return a dictionary whose key is the id field value
	# otherwise a list of dictionaries - effectively an associative array
	results = []

	if sf:
		layer = sf.GetLayer()
	
		layerDefinition = layer.GetLayerDefn()
	
		feat = layer.GetNextFeature()
		
		# for each feature
		while feat is not None:
			# for each attribute
			#print '#####' + field_name + '##'
			results.append(feat.GetField(str(field_name)))
					
			feat = layer.GetNextFeature()
			#print "feat:" + str(feat)
			
	return results
	
	
	
def rasterize_shapefile_with_gdal(shapefile_to_rasterise,destination_raster_name,attribute_name,ndv=-99,keep_orig_raster_no_data_cells = False,pixel_size = 0.0008333,ot="Int32"):
		
	# this method  guarantees the same number of pixels in the output raster as in the input one - but it doesn't guarantee a pixel size!!
	
	if enable_debug:
		print '## Entering rasterize_shapefile with gdal '
		print 'shapefile_to_rasterise:' + shapefile_to_rasterise 
		print 'destination_raster_name:' + destination_raster_name 
		print 'attribute_name:' + attribute_name 
		print 'pixel_size:' + str(pixel_size) 
	

	# get original raster array and save array
	driver = gdal.GetDriverByName('GTiff')
		
	template_raster_orig = gdal.Open(destination_raster_name)
	template_band_orig = template_raster_orig.GetRasterBand(1)
	template_array_orig = np.array(template_band_orig.ReadAsArray())
	template_ndv_orig= template_band_orig.GetNoDataValue()
	
	
	geoTransform = template_raster_orig.GetGeoTransform()
	minx = geoTransform[0]
	maxy = geoTransform[3]
	maxx = minx + geoTransform[1] * template_raster_orig.RasterXSize
	miny = maxy + geoTransform[5] * template_raster_orig.RasterYSize
	
	extents = ' -te ' + ' '.join([str(x) for x in [minx, miny, maxx, maxy]]) + ' -tr 0.0008333 -0.0008333 '
	
	print '###EXTENTS:' + extents + ' -tr 0.0008333 -0.0008333 '
	
		
		
	template_raster_orig = None
	
	width = template_array_orig.shape[1]
	height = template_array_orig.shape[0]
	

	
	# create a gdal_translate string to cut to a window from the source
	#gdal_cmd_as_string = 'gdal_rasterize -a "'  + attribute_name + '" ' + ' -ot ' + ot  + ' -ts ' + str(width) + '  ' + str(height) + ' -a_nodata '  + str(ndv)+ ' -co compress=LZW "' + shapefile_to_rasterise + '" "'+ destination_raster_name +'"'
	
	gdal_cmd_as_string = 'gdal_rasterize -a "'  + attribute_name + '" ' + ' -ot ' + ot  + extents + ' -a_nodata '  + str(ndv)+ ' -co compress=LZW "' + shapefile_to_rasterise + '" "'+ destination_raster_name +'"'
	
	if enable_debug:
		print gdal_cmd_as_string
	
	ret = os.system(gdal_cmd_as_string)

	if keep_orig_raster_no_data_cells:		
		
		# open new raster made form shapefile and set any ndv in orig to ndv in new
		template_raster_new = gdal.Open(destination_raster_name,gdal.GA_Update)
		template_band_new = template_raster_new.GetRasterBand(1)
		template_array_new = np.array(template_band_new.ReadAsArray())
		template_ndv_new= template_band_new.GetNoDataValue()

		# set any ndvs in orginal array to ndv in new
		template_array_new[template_array_orig == template_ndv_orig] = template_ndv_new
	
		template_band_new.WriteArray(template_array_new)

		template_raster_new = None

	return ret
	
def rasterize_shapefile_with_gdal_20190321(shapefile_to_rasterise,destination_raster_name,attribute_name,ndv=-99,keep_orig_raster_no_data_cells = False,pixel_size = 0.0008333,ot="Int32"):
		
	# this method  guarantees the same number of pixels in the output raster as in the input one - but it doesn't guarantee a pixel size!!
	
	if enable_debug:
		print '## Entering rasterize_shapefile with gdal '
		print 'shapefile_to_rasterise:' + shapefile_to_rasterise 
		print 'destination_raster_name:' + destination_raster_name 
		print 'attribute_name:' + attribute_name 
		print 'pixel_size:' + str(pixel_size) 
	

	# get original raster array and save array
	driver = gdal.GetDriverByName('GTiff')
		
	template_raster_orig = gdal.Open(destination_raster_name)
	template_band_orig = template_raster_orig.GetRasterBand(1)
	template_array_orig = np.array(template_band_orig.ReadAsArray())
	template_ndv_orig= template_band_orig.GetNoDataValue()
		
	template_raster_orig = None
	
	width = template_array_orig.shape[1]
	height = template_array_orig.shape[0]
	
	
	# create a gdal_translate string to cut to a window from the source
	gdal_cmd_as_string = 'gdal_rasterize -a "'  + attribute_name + '" ' + ' -ot ' + ot  + ' -ts ' + str(width) + '  ' + str(height) + ' -a_nodata '  + str(ndv)+ ' -co compress=LZW "' + shapefile_to_rasterise + '" "'+ destination_raster_name +'"'
	
	if enable_debug:
		print gdal_cmd_as_string
	
	ret = os.system(gdal_cmd_as_string)

	if keep_orig_raster_no_data_cells:		
		
		# open new raster made form shapefile and set any ndv in orig to ndv in new
		template_raster_new = gdal.Open(destination_raster_name,gdal.GA_Update)
		template_band_new = template_raster_new.GetRasterBand(1)
		template_array_new = np.array(template_band_new.ReadAsArray())
		template_ndv_new= template_band_new.GetNoDataValue()

		# set any ndvs in orginal array to ndv in new
		template_array_new[template_array_orig == template_ndv_orig] = template_ndv_new
	
		template_band_new.WriteArray(template_array_new)

		template_raster_new = None

	return ret

	
def rasterize_shapefile_with_gdal_orig(shapefile_to_rasterise,destination_raster_name,attribute_name,ndv=-99,keep_orig_raster_no_data_cells = False,pixel_size = 0.000833333,ot="Int32"):
	
	# this method doesn't guarantee the same number of pixels in the output raster as in the input one - it guarantees a pixel size!!
	
	if enable_debug:
		print '## Entering rasterize_shapefile with gdal '

	if keep_orig_raster_no_data_cells:
		
		# get original raster array and save array
		driver = gdal.GetDriverByName('GTiff')
		
		template_raster_orig = gdal.Open(destination_raster_name)
		template_band_orig = template_raster_orig.GetRasterBand(1)
		template_array_orig = np.array(template_band_orig.ReadAsArray())
		template_ndv_orig= template_band_orig.GetNoDataValue()
		
		template_raster_orig = None
	
	# create a gdal_translate string to cut to a window from the source
	gdal_cmd_as_string = 'gdal_rasterize -a "'  + attribute_name + '" ' + ' -ot ' + ot  + ' -tr ' + str(pixel_size) + ' -' + str(pixel_size) + ' -tap -a_nodata '  + str(ndv)+ ' -co compress=LZW  "' + shapefile_to_rasterise + '" "'+ destination_raster_name +'"'
	
	if enable_debug:
		print gdal_cmd_as_string
	
	ret = os.system(gdal_cmd_as_string)

	if keep_orig_raster_no_data_cells:		
		
		# open new raster made form shapefile and set any ndv in orig to ndv in new
		template_raster_new = gdal.Open(destination_raster_name,gdal.GA_Update)
		template_band_new = template_raster_new.GetRasterBand(1)
		template_array_new = np.array(template_band_new.ReadAsArray())
		template_ndv_new= template_band_new.GetNoDataValue()

		# set any ndvs in orginal array to ndv in new
		template_array_new[template_array_orig == template_ndv_orig] = template_ndv_new

		template_band_new.WriteArray(template_array_new)

		template_raster_new = None

	return ret
def rescale_raster_mode(src_raster_name,destination_raster_name,scale,orig_pixel_size = 0.000833333):
	
	new_pixel_size = orig_pixel_size*scale
	
	# resample the orignal raster to the required scale - this gives the mode of the no data cells
	gdal_cmd_as_string = 'gdalwarp -overwrite -tr ' + str(new_pixel_size) + '  ' + str(new_pixel_size) 
	gdal_cmd_as_string += ' -tap -co compress=LZW -t_srs epsg:4326 -r mode ' + src_raster_name + ' ' + destination_raster_name

	if enable_debug:
		print gdal_cmd_as_string
	
	return os.system(gdal_cmd_as_string)
	
def rescale_raster_mode_by_width_and_height(src_raster_name,destination_raster_name,width,height):
		
	# resample the orignal raster to the required scale - this gives the mode of the no data cells
	gdal_cmd_as_string = 'gdalwarp -overwrite -ts ' + str(width) + '  ' + str(height) 
	gdal_cmd_as_string += ' -co compress=LZW -t_srs epsg:4326 -r mode ' + src_raster_name + ' ' + destination_raster_name

	if enable_debug:
		print gdal_cmd_as_string
	
	return os.system(gdal_cmd_as_string)
	
def rescale_raster_sum(src_raster_name,destination_raster_name,scale,orig_pixel_size = 0.000833333):
	# this ought to be easy but there is no "sum" option for gdal warp so we have to messa round with "average"
	# assume orig is 0.000833333 pixel size

	new_pixel_size = orig_pixel_size*scale
	
	dest_dir = os.path.dirname(destination_raster_name)
	pop_raster_ones = dest_dir +  "/pop_raster_ones.tif"
	pop_raster_ones_resampled = dest_dir + "/pop_raster_ones_resampled.tif"
	
	# resample the orignal raster to the required scale - this gives the average of the no data cells
	gdal_cmd_as_string = 'gdalwarp -overwrite -tr ' + str(new_pixel_size) + '  ' + str(new_pixel_size) 
	gdal_cmd_as_string += ' -tap -co compress=LZW -t_srs epsg:4326 -r average ' + src_raster_name + ' ' + destination_raster_name

	if enable_debug:
		print gdal_cmd_as_string
		
	ret = os.system(gdal_cmd_as_string)
	
	if not ret:
		# we don't want the average we want the sum, so we need to calculate the number of no cells that went into the average
		# to then ge the sum by multiplying the average by that number
		# set cell to 1 where it has a value: ndv  to 0
		gdal_cmd_as_string = 'gdal_calc.py --overwrite -A ' + src_raster_name + ' --outfile=' + pop_raster_ones + ' --calc="1*(A>=0)" --NoDataValue=0 --co="COMPRESS=LZW"'
		
		if enable_debug:
			print gdal_cmd_as_string
		
		ret = os.system(gdal_cmd_as_string)
	
		if not ret:
			# reset the no-data value so that for further processing we have a raster of 1s and 0 with none with ndv
			gdal_cmd_as_string = 'gdal_edit.py ' + pop_raster_ones + ' -a_nodata -99'
			
			if enable_debug:
				print gdal_cmd_as_string
		
			ret = os.system(gdal_cmd_as_string)
		
		if not ret:
			# now we resample the raster of ones and zeros using average to get e.g 1,0.25,0.75
			# the no data value will be 0 as this indicates where none had a value so we aren't interested in them
			gdal_cmd_as_string = 'gdalwarp -overwrite -tr ' + str(new_pixel_size) + '  ' + str(new_pixel_size) 
			gdal_cmd_as_string += ' -tap -co compress=LZW -t_srs epsg:4326 -r average -dstnodata 0 ' + pop_raster_ones + ' ' + pop_raster_ones_resampled 
			
			if enable_debug:
				print gdal_cmd_as_string
				
			ret = os.system(gdal_cmd_as_string)
		
		if not ret:
			gdal_cmd_as_string = 'gdal_calc.py --overwrite -A '  + destination_raster_name + ' -B ' + ' ' + pop_raster_ones_resampled 
			gdal_cmd_as_string += ' --outfile=' + destination_raster_name + ' --calc="A*B*' + str(scale*scale) + '" --co="COMPRESS=LZW"'
			
			if enable_debug:
				print gdal_cmd_as_string
				
			ret = os.system(gdal_cmd_as_string)
	
	# finally if they exsit remove the intermediate files
	if os.path.exists(pop_raster_ones):
		os.remove(pop_raster_ones)
		
	if os.path.exists(pop_raster_ones_resampled):
		os.remove(pop_raster_ones_resampled)	
		
	return ret
	
############################################################################
########################### functions for calcualte border pixels


# this create a boundary polyline then creates a buffer polygon
def make_boundary_shapefile(source_file,new_shapefile_name,buffer_distance = 0.000001,feature_name=None,feature_value=None):

	if enable_debug:
		print "entering make_boundary_shapefile"
		
	# get the input layer and get its spatial reference
	driver = ogr.GetDriverByName("ESRI Shapefile")

	# 0 is read only
	sf = driver.Open(source_file,0)

	bnd = None
	
	if sf:
		layer = sf.GetLayer()
		srs = layer.GetSpatialRef()
		
		featureCount = layer.GetFeatureCount()
		#print "Number of features in shapefile %d" % (featureCount)
	  
		layerDefinition = layer.GetLayerDefn()

		# get attribute names
		#for i in range(layerDefinition.GetFieldCount()):
		#	print layerDefinition.GetFieldDefn(i).GetName()
		

		feat = layer.GetNextFeature()
		unionz = ogr.Geometry(ogr.wkbPolygon)
		while feat is not None:
			#print "\nID:" + str(feat.GetField("ID"))
			#print "NAME_1:" + feat.GetField("NAME_1")
			#print "SHAPE_AREA:" + str(feat.GetField("SHAPE_AREA"))
			
			if feature_name and feature_value:
				#print "### " + feature_name + " " + str(feature_value) + " " + str(feat.GetField(feature_name))
				if feat.GetField(str(feature_name)) == feature_value:
					geom = feat.GetGeometryRef()
					unionz = unionz.Union(geom)
					
			else:
				geom = feat.GetGeometryRef()
				unionz = unionz.Union(geom)

			
			#unionz = unionz.Union(geom)
			
	
			feat = layer.GetNextFeature()
			
				
				
		bnd = unionz.Boundary()
		
		#print "geom name" + str(bnd.GetGeometryName())
		#print  "geom type" +  str(bnd.GetGeometryType())
				
		# this should be the length of the diagonal across a pixel
		#buffer_distance = 0.001178464
		# ICW: depends what we are doing with it
		# currently just want a very thin band as are going to do a rasterize poly using hte touched by option as a means 
		# of getting a "crosses" set of pixels
		buffer = None
		if buffer_distance is None:
			# convert line to poly
			# convex hull works
			#buffer = bnd.ConvexHull()
			buffer = bnd.Polygonize()

		else:

			buffer= bnd.Buffer(buffer_distance)
		 		
	

	# create the output layer to fill with reprojected polygons
	
	outputShapefile = new_shapefile_name
	if os.path.exists(outputShapefile):
		driver.DeleteDataSource(outputShapefile)
	outDataSet = driver.CreateDataSource(outputShapefile)
	
	#outLayer = outDataSet.CreateLayer("boundary", srs,geom_type=ogr.wkbMultiLineString)
	

	outLayer = outDataSet.CreateLayer("buffer", srs,geom_type=ogr.wkbMultiPolygon)
	outLayerDefn = outLayer.GetLayerDefn()
	

	# if we have a boundary - add buffer to shapefile
	if bnd is not None:
		outFeature = ogr.Feature(outLayerDefn)
		# set the geometry and attribute
		#outFeature.SetGeometry(bnd)

		outFeature.SetGeometry(buffer)
	
		# ICW : is this range correct?
		# add the feature to the shapefile
		outLayer.CreateFeature(outFeature)
	
	outDataSet = None
	
#create a raster based on a template raster. This raster will be empty and will need to be filled later using a different function
def new_raster_from_base(template_raster, outputURI, format, nodata, datatype):

	if enable_debug:
		print '## Entering new_raster_from_base'
		
	#Creates a new raster with a specific format and nodata value
	#Equal in resolution & projection & extent to template_raster
	cols = template_raster.RasterXSize
	rows = template_raster.RasterYSize
	projection = template_raster.GetProjection()
	geotransform = template_raster.GetGeoTransform()
	bands = template_raster.RasterCount

	driver = gdal.GetDriverByName(format)

	new_raster = driver.Create(str(outputURI), cols, rows, bands, datatype)
	new_raster.SetProjection(projection)
	new_raster.SetGeoTransform(geotransform)

	for i in range(bands):
		new_raster.GetRasterBand(i + 1).SetNoDataValue(nodata)
		new_raster.GetRasterBand(i + 1).Fill(nodata)
			
			
	return new_raster
	

	
def indexrast_create(src_rast,filename,maintain_src_ndv = True):

	if enable_debug:
		print '## Entering indexrast_create '

	sr = gdal.Open(src_rast)
	sr_band = sr.GetRasterBand(1)
	sr_ndv = sr_band.GetNoDataValue()
	sr_array = np.array(sr_band.ReadAsArray())
	
	#This function creates a raster labeled 1:n, where n is the number of cells in the raster.
	raster_out = filename
	index_rast = new_raster_from_base(sr, raster_out, 'GTiff',-99, gdal.GDT_Int32)
	index_rast_array = np.array(index_rast.GetRasterBand(1).ReadAsArray())
	
	index = 1
	# iterate over the raster, giving each cell a unique value starting from 1
	for i in range(index_rast_array.shape[0]):
		for j in range(index_rast_array.shape[1]):
			index_rast_array[i][j]=index
			index += 1
			
	if maintain_src_ndv:
		index_rast_array[sr_array == sr_ndv] = -99
		
	outBand = index_rast.GetRasterBand(1)
	#save the raster and reload
	outBand.WriteArray(index_rast_array, 0, 0)
	outBand.FlushCache()
	outBand.SetNoDataValue(-99)

	index_rast = None
	index_rast_reload = gdal.Open(filename)
	 
 
	 
	return index_rast_reload

	
def rasterize_shapefile_ORIG(shapefile_to_rasterise,template_raster_name,new_filename,ndv = -99):
	
	#open raster
	template_raster = gdal.Open(template_raster_name)
	
	# open shapefile , 0 is read only
	driver = ogr.GetDriverByName("ESRI Shapefile")
	poly = driver.Open(shapefile_to_rasterise,0) 
	poly_layer = poly.GetLayer()

    # create a new raster
	poly_raster = new_raster_from_base(template_raster, new_filename, 'GTiff',
								ndv, gdal.GDT_Int32)
	band = poly_raster.GetRasterBand(1)
	nodata = band.GetNoDataValue()
	band.Fill(nodata)
	
	# now rasterise the layer
	#ALL_TOUCHED = TRUE specifies that the rasterization function will be inclusive, including all cells that touch a polygon, regardless of whether it covers the cell completely
	gdal.RasterizeLayer(poly_raster, [1], poly_layer, options = ["ALL_TOUCHED=TRUE"])

	band.FlushCache()
	band.SetNoDataValue(ndv)

	poly_raster = None

def rasterize_shapefile(shapefile_to_rasterise,template_raster_name,new_filename,ndv = -99):

	if enable_debug:
		print '## Entering rasterize_shapefile' 

	# this returns a raster with either value 255 in a pixl or 0 for no data value
	#open raster
	template_raster = gdal.Open(template_raster_name)
	
	
	# open shapefile , 0 is read only
	driver = ogr.GetDriverByName("ESRI Shapefile")
	poly = driver.Open(shapefile_to_rasterise,0) 
	poly_layer = poly.GetLayer()

    # create a new raster
	poly_raster = new_raster_from_base(template_raster, new_filename, 'GTiff',
								ndv, gdal.GDT_Int32)
	band = poly_raster.GetRasterBand(1)
	nodata = band.GetNoDataValue()
	band.Fill(nodata)
	
	# now rasterise the layer
	#ALL_TOUCHED = TRUE specifies that the rasterization function will be inclusive, including all cells that touch a polygon, regardless of whether it covers the cell completely
	gdal.RasterizeLayer(poly_raster, [1], poly_layer, options = ["ALL_TOUCHED=TRUE"])
	
	# now try changing the 255 values in the new raster to the values from the template
	template_rast_array = np.array(template_raster.GetRasterBand(1).ReadAsArray())    
	new_rast_array = np.array(band.ReadAsArray())
	
	bordervals = template_rast_array[(new_rast_array == 255)]
	new_rast_array[(new_rast_array == 255)] = bordervals
	band.WriteArray(new_rast_array)

	band.FlushCache()
	band.SetNoDataValue(ndv)

	poly_raster = None	
	
	
	
def polygonise_raster (raster_to_polygonise,new_shapefile_name):	
	
	if enable_debug:
		print '## Entering polygonise_raster ' 
		
	## now try polygonizing it as we want the original pixels rasterized !!

	rtp_ds = gdal.Open( raster_to_polygonise )
	rtp_band = rtp_ds.GetRasterBand(1)
	
	srs = osr.SpatialReference()
	srs.ImportFromWkt(rtp_ds.GetProjectionRef())

	layername = "border"
	drv = ogr.GetDriverByName("ESRI Shapefile")
	ds = drv.CreateDataSource( new_shapefile_name )
	layer = ds.CreateLayer(layername, srs = srs )
	
	newField = ogr.FieldDefn('INDEX_ID', ogr.OFTInteger)
	layer.CreateField(newField)
	
	# second parameter is an optional mask band. All pixels in the mask band with a value other than zero will be considered suitable for collection as polygon
	# as we have set the no data value to 0 this works fine
	# 4th parameter is 	the attribute field index indicating the feature attribute into which the pixel value of the polygon should be written.
	gdal.Polygonize( rtp_band,  rtp_band, layer, 0, [], callback=None )
	
	if enable_debug:
		print '## Leaving polygonise_raster ' 
		print
	
def polygonise_open_raster (opened_raster,new_shapefile_name,layername = "my_layer",field_name='INDEX_ID'):	
	
	if enable_debug:
		print '## Entering polygonise_open_raster' 
	
		
	raster_band = opened_raster.GetRasterBand(1)

	srs = osr.SpatialReference()
	srs.ImportFromWkt(opened_raster.GetProjectionRef())
	
	#print "proj ref" + str(opened_raster.GetProjectionRef())
	
	drv = ogr.GetDriverByName("ESRI Shapefile")
	ds = drv.CreateDataSource( new_shapefile_name )
	layer = ds.CreateLayer(layername, srs = srs )
	
	newField = ogr.FieldDefn(field_name, ogr.OFTInteger)
	layer.CreateField(newField)
	
	# second parameter is an optional mask band. All pixels in the mask band with a value other than zero will be considered suitable for collection as polygon
	# as we have set the no data value to 0 this works fine
	# 4th parameter is 	the attribute field index indicating the feature attribute into which the pixel value of the polygon should be written.
	gdal.Polygonize( raster_band,  raster_band, layer, 0, [], callback=None )

	
	
def polygonise_raster_and_add_attributes (index_raster,new_shapefile_name,arr_attribute_rasters = None,arr_attribute_names = None,arr_attribute_types = None):	

	if enable_debug:
		print '## Entering polygonise_raster_and_add_attributes' 
		
	## now try polygonizing it as we want the original pixels rasterized !!

	ir_ds = gdal.Open( index_raster )
	ir_band = ir_ds.GetRasterBand(1)
	
	srs = osr.SpatialReference()
	srs.ImportFromWkt(ir_ds.GetProjectionRef())

	layername = "border"
	drv = ogr.GetDriverByName("ESRI Shapefile")
	ds = drv.CreateDataSource( new_shapefile_name )
	layer = ds.CreateLayer(layername, srs = srs )
	
	newField = ogr.FieldDefn('INDEX_ID', ogr.OFTInteger)
	layer.CreateField(newField)
	
	# second parameter is an optional mask band. All pixels in the mask band with a value other than zero will be considered suitable for collection as polygon
	# as we have set the no data value to 0 this works fine
	# 4th parameter is 	the attribute field index indicating the feature attribute into which the pixel value of the polygon should be written.
	gdal.Polygonize( ir_band,  ir_band, layer, 0, [], callback=None )
	
	ds = None
	

	if arr_attribute_rasters  and arr_attribute_names and arr_attribute_types and len(arr_attribute_rasters) == len(arr_attribute_names):
		
		# pen shapefile
		driver = ogr.GetDriverByName("ESRI Shapefile")
		sf = driver.Open(new_shapefile_name,1)
		
		if sf:
			# doing each raste in turn takes up less memory
			for i in range(len(arr_attribute_rasters)):
				
				attribute_name = arr_attribute_names[i]
					
				ds = gdal.Open(arr_attribute_rasters[i])
				raster_values = np.array(ds.GetRasterBand(1).ReadAsArray())
				
				# get dimensions of raster (they should all have the same dimensions)
				array_rows = raster_values.shape[0]
				array_columns = raster_values.shape[1]		

				# get layer and add a new field (attribute)
				layer = sf.GetLayer()
					
				newField = ogr.FieldDefn(attribute_name, arr_attribute_types[i])
				layer.CreateField(newField)
					
				# for each feature set the new attribute value
				feat = layer.GetNextFeature()

				while feat is not None:
					
					# the index tells us hte array position
					array_index = feat.GetField("INDEX_ID")
				
					# index starts from 1 not 0
					array_col = (array_index % array_columns) -1
					array_row = (array_index // array_columns) 
					
					val = float(raster_values[array_row][array_col])
				
					feat.SetField(attribute_name, val)
					
					layer.SetFeature(feat)
						
					feat = layer.GetNextFeature()
					
				# to start reading features for the beginning again	
				layer.ResetReading()
			
			# close data source
			sf.Destroy()

def add_attributes_to_shapefile (shapefile_name,arr_attribute_values = None,arr_attribute_names = None,arr_attribute_types = None):	
	if enable_debug:
		print ("Entering add_attributes_to_shapefile")
		
	if arr_attribute_values  and arr_attribute_names and arr_attribute_types and len(arr_attribute_values) == len(arr_attribute_names):
		
		# pen shapefile
		driver = ogr.GetDriverByName("ESRI Shapefile")
		sf = driver.Open(shapefile_name,1)
		
		if sf:
			# doing each raste in turn takes up less memory
			for i in range(len(arr_attribute_values)):
				
				attribute_name = arr_attribute_names[i]
					
				attribute_value= arr_attribute_values[i]	

				# get layer and add a new field (attribute)
				layer = sf.GetLayer()
					
				newField = ogr.FieldDefn(attribute_name, arr_attribute_types[i])
				layer.CreateField(newField)
					
				# for each feature set the new attribute value
				feat = layer.GetNextFeature()

				while feat is not None:				
				
					feat.SetField(attribute_name, attribute_value)
					
					layer.SetFeature(feat)
						
					feat = layer.GetNextFeature()
					
				# to start reading features from the beginning again	
				layer.ResetReading()
			
			# close data source
			sf.Destroy()

			
def	add_strata_id_to_results_shapefile(shapefile_name,PSU_id_strata_id_info):

	if enable_debug:
		print ("Entering add_strata_id_to_results_shapefile")
		print PSU_id_strata_id_info
	
	# the keys of PSU_id_strata_id_info are the cluster ids


	# open shapefile
	driver = ogr.GetDriverByName("ESRI Shapefile")
	sf = driver.Open(shapefile_name,1)
		
	if sf:
	
		layer = sf.GetLayer()
		
		# add the new attributes - just the strata id at the moment
		
		newField = ogr.FieldDefn("str_id", ogr.OFTString)
		layer.CreateField(newField)
		
		newField = ogr.FieldDefn("cl_type", ogr.OFTString)
		layer.CreateField(newField)
		
		newField = ogr.FieldDefn("s_cl_st", ogr.OFTString)
		layer.CreateField(newField)
		
		newField = ogr.FieldDefn("s_cl_tot", ogr.OFTString)
		layer.CreateField(newField)
		
		newField = ogr.FieldDefn("s_cl_pop", ogr.OFTString)
		layer.CreateField(newField)
		
		newField = ogr.FieldDefn("s_cl_hh", ogr.OFTString)
		layer.CreateField(newField)
							
		feat = layer.GetNextFeature()
		
		# for each feature
		while feat is not None:
			cl_id = feat.GetField('cl_id')
				
			# if there is info for this key then add to 
			if( cl_id in PSU_id_strata_id_info.keys() ):
				
				feat.SetField("str_id", str(PSU_id_strata_id_info[cl_id]['strata_id']) )
				feat.SetField("cl_type", str(PSU_id_strata_id_info[cl_id]['cl_type']) )
				feat.SetField("s_cl_st", str(PSU_id_strata_id_info[cl_id]['s_cl_st']) )
				feat.SetField("s_cl_tot", str(PSU_id_strata_id_info[cl_id]['s_cl_tot']) )
				feat.SetField("s_cl_pop", str(PSU_id_strata_id_info[cl_id]['s_cl_pop']) )
				feat.SetField("s_cl_hh", str(PSU_id_strata_id_info[cl_id]['s_cl_hh']) )
									
				layer.SetFeature(feat)
				
			feat = layer.GetNextFeature()
		
		# close data source
		sf.Destroy()		
		
def add_strata_attributes_to_results_shapefile (shapefile_name,strata_info_csv):	
	
	if enable_debug:
		print ("Entering add_strata_attributes_to_shapefile")
	
	# create a dictionary whose keys are strata ids
	strata_info={}
	with open(strata_info_csv) as csvfile:
		reader = csv.DictReader(csvfile)
		for row in reader:
			strata_info[row['id']] = row

			
	#print strata_info
	#print strata_info.keys()
	
	# open shapefile
	driver = ogr.GetDriverByName("ESRI Shapefile")
	sf = driver.Open(shapefile_name,1)
	
	id = 1
		
	if sf:
	
		layer = sf.GetLayer()
		
		# add the new attributes
		newField = ogr.FieldDefn("u_n", ogr.OFTString)
		layer.CreateField(newField)
		
		newField = ogr.FieldDefn("u_medpop", ogr.OFTString)
		layer.CreateField(newField)
		
		newField = ogr.FieldDefn("st_hhsiz", ogr.OFTString)
		layer.CreateField(newField)

		newField = ogr.FieldDefn("st_pop_n", ogr.OFTString)
		layer.CreateField(newField)
		
		newField = ogr.FieldDefn("st_pop_p", ogr.OFTString)
		layer.CreateField(newField)
		
		newField = ogr.FieldDefn("st_hh_n", ogr.OFTString)
		layer.CreateField(newField)
		
		newField = ogr.FieldDefn("st_hh_p", ogr.OFTString)
		layer.CreateField(newField)
		
		newField = ogr.FieldDefn("st_name", ogr.OFTString)
		layer.CreateField(newField)
		
		# cl_id needs to be incremental 1:n
		# we will keep the orig id as 
		newField = ogr.FieldDefn("ori_id", ogr.OFTInteger)
		layer.CreateField(newField)
		
		#newField = ogr.FieldDefn("tot_psu", ogr.OFTString)
		#layer.CreateField(newField)
		
		feat = layer.GetNextFeature()
	
		# for each feature
		while feat is not None:
			
			st_id = feat.GetField('str_id')
			
			# move cluster id to original id and set cl_id to be incremental integer starting form 1
			feat.SetField("ori_id", int(float(feat.GetField('cl_id'))) )
			feat.SetField("cl_id", id)
			id+=1
			#print st_id
			
			# if there is info for this key then add to 
			if(str(st_id) in strata_info.keys()):
			
				#print strata_info[str(st_id)]
				
				feat.SetField("u_n", strata_info[str(st_id)]['total_PSUs_in_coverage_area'])
				feat.SetField("u_medpop", strata_info[str(st_id)]['median_pop_per_psu'])
				feat.SetField("st_hhsiz", strata_info[str(st_id)]['str_hh_size'])
				feat.SetField("st_pop_n", strata_info[str(st_id)]['population'])
				feat.SetField("st_pop_p", strata_info[str(st_id)]['population_percent'])
				feat.SetField("st_hh_n", strata_info[str(st_id)]['str_households'])
				feat.SetField("st_hh_p", strata_info[str(st_id)]['households_percent'])

				
				feat.SetField("st_name", strata_info[str(st_id)]['strata_name'])
				#feat.SetField("tot_psu", strata_info[str(st_id)]['total_PSUs'])
				
				## we do a bit of a hack here : because we were selecting PSUs in terms of households not populations
				## the values passed through here are in terms of HH not pop so need to re-calculate
				feat.SetField("s_cl_pop",round( float(feat.GetField('s_cl_pop'))*float(strata_info[str(st_id)]['str_hh_size']),5))
				feat.SetField("s_cl_hh",round(float(feat.GetField('s_cl_hh'))*float(strata_info[str(st_id)]['str_hh_size']),5))
					
				layer.SetFeature(feat)
				
			feat = layer.GetNextFeature()
		
		# close data source
		sf.Destroy()
	
	if enable_debug:
		print ("Leaving add_strata_attributes_to_shapefile")		
	
def calculate_border_pixels_area_inside_strata(border_pixels_shapefile,strata_file,pop_raster,index_raster,apply_scaling=False):	
	# for each border pixel polygon get the interesction with the strata polygon
	# from the intersection area we know what % of the pixel is in the strata whihc is a scaling factor
	# for the pixel id we know which raster pixel it corresponds to and can adjust the population accordingly

	# the final output is the array pop_file_scaling which can be multiplied by the pop_array to get the scaled pop array
	# which allows for the pixels only partly covered by the shape of interest 
	
	# NOTE: currently not doing anything with the results but the logic all works
	# eventually with either return a scaling raster or a copy of a scaled popn raster
	
	# get an array for the population file values
	# make a similar sized array with all elements set to 1 that will be used store the scaling values the pop file
	
	if enable_debug:
		print ("Entering calculate_border_pixels_area_inside_strata")	
	
	pop_file_raster = gdal.Open(pop_raster,gdal.GA_Update)
	pop_file_raster_band = pop_file_raster.GetRasterBand(1)
	pop_file_raster_ndv = pop_file_raster_band.GetNoDataValue()
	
	pop_file_raster_array = np.array(pop_file_raster_band.ReadAsArray())
	

	pop_file_scaling = np.array(pop_file_raster_band.ReadAsArray())
	pop_file_scaling.fill(1)


	# index_file
	index_file_raster = gdal.Open(index_raster)
	index_file_raster_band = index_file_raster.GetRasterBand(1)
	index_file_raster_array = np.array(index_file_raster_band.ReadAsArray())

	array_rows = pop_file_scaling.shape[0]
	array_columns = pop_file_scaling.shape[1]

	# open the shape files
	driver = ogr.GetDriverByName("ESRI Shapefile")
	
	pr = driver.Open(border_pixels_shapefile,0)
	sf = driver.Open(strata_file,0)

	srs=None
	
	number_of_overlapping_cells_found = 0;
	
	if pr and sf:
		pr_layer = pr.GetLayer()	
		pr_featureCount = pr_layer.GetFeatureCount()
		#print "Number of features in polygonised raster %d" % (pr_featureCount)
		
		sf_layer = sf.GetLayer()	
		sf_featureCount = sf_layer.GetFeatureCount()
		#print "Number of features in strata  %d" % (sf_featureCount)
					
		## these need to have the same projection!!	
		pr_layer_s_ref = pr_layer.GetSpatialRef()
		sf_layer_s_ref = sf_layer.GetSpatialRef()

		coordTrans = osr.CoordinateTransformation(sf_layer_s_ref,pr_layer_s_ref)
		
		for pr_feat in pr_layer:
			
			pr_geom = pr_feat.GetGeometryRef()
			pixel_area = pr_geom.GetArea()
			
					
			for sf_feat in sf_layer:
				# we need to handle what actually happens when there is more than one intersection
				sf_geom = sf_feat.GetGeometryRef()
				sf_geom.Transform(coordTrans)
					
				if  pr_geom.Overlaps(sf_geom):
				
					number_of_overlapping_cells_found += 1
				
					# at this point  with the area and the id we know which pixel we need to adjust the population count for
					intersection =  pr_geom.Intersection(sf_geom)
					#print "\nDN:" + str(pr_feat.GetField("INDEX_ID"))	
					array_index = pr_feat.GetField("INDEX_ID");
					
					# we can calcualte the array position from the index 
					# index starts from 1 not 0
					array_col = (array_index % array_columns) -1
					array_row = (array_index // array_columns) 
					
					# assuming 100m*100m cells - we can't always assume htis as user may rescale but OK for now
					cell_percentage_covered = round(intersection.GetArea()/pixel_area,3)
					
					# if 1 then this is the first section of intersection for this pixel
					# otherwise it is additional coverage so just add them
					if pop_file_scaling[array_row][array_col]  == 1:
						pop_file_scaling[array_row][array_col] = cell_percentage_covered
					else :
						pop_file_scaling[array_row][array_col] = pop_file_scaling[array_row][array_col] +  cell_percentage_covered
						
					#print "index raster value for this cell:" + str(index_file_raster_array[array_row][array_col])
					#print "intersection found: areaz =" + str(cell_percentage_covered)

			# reset for next loop
			sf_layer.ResetReading() 
			
		if apply_scaling:
			#print "applying scaling"
			#print np.sum(pop_file_scaling[pop_file_scaling != 1])
			#print np.sum(pop_file_raster_array[pop_file_raster_array != pop_file_raster_ndv])
			#print np.sum(pop_file_scaling[pop_file_raster_array != pop_file_raster_ndv])
			# do the scaling
			#pop_file_raster_array = pop_file_raster_array*pop_file_scaling
			# don't scale ndv
			pop_file_raster_array[pop_file_raster_array != pop_file_raster_ndv] = pop_file_raster_array[pop_file_raster_array != pop_file_raster_ndv]*pop_file_scaling[pop_file_raster_array != pop_file_raster_ndv]
			pop_file_raster_band.WriteArray(pop_file_raster_array)
			pop_file_raster_band.FlushCache()
			pop_file_raster = None
	
			
		#print "number_of_overlapping_cells_found:" + str(number_of_overlapping_cells_found)
				
	else:
		print "Couldn't open shapefile"
	
	# apply scaling to the population raster array and return it	
	pop_file_raster_array[pop_file_raster_array != pop_file_raster_ndv] = pop_file_raster_array[pop_file_raster_array != pop_file_raster_ndv]*pop_file_scaling[pop_file_raster_array != pop_file_raster_ndv]
	
	return pop_file_raster_array

def	create_histogram_from_values(vals,image_name,title = None,info = None):
	# TO DO - check the location exists
	if enable_debug:
		print ("Entering create_histogram_from_values")

	
	plt.hist(vals,bins=100)
	
	plt.xlabel('Population')
	
	if title:
		plt.title(title)
	
	if info:	
		left, right = plt.xlim()
		bottom,top = plt.ylim()
		width = right - left
		xpos= 0.6*right
		ypos = 0.8*top
		
		#bbox=dict(facecolor='red',width=100,height=75, alpha=0.2)
		plt.text( xpos,ypos, info)

	plt.savefig(image_name)
	
	plt.close()


def extractPopulationValuesForStrata(pop_raster_name,strata_raster_name=None,multicell_raster_name = None):
	# if no strata raster is given then treat the situation as if there is a single strata with id 1
	# if multicell then each PSU is a number of pixels rather than a single pixel
	if enable_debug:
		print '## Entering extractPopulationValuesForStrata '
		print pop_raster_name
		print strata_raster_name
		print multicell_raster_name

	# this returns a dictionary with keys the strata ids and values the populations
	pop_raster = gdal.Open(pop_raster_name)
	pop_raster_band = pop_raster.GetRasterBand(1)
	pop_raster_array = np.array(pop_raster_band.ReadAsArray())
	pop_ndv= pop_raster_band.GetNoDataValue()
	
	multicell_raster_array = None
	
	#  slow - might want to disable for testing
	if  multicell_raster_name:
		multicell_raster = gdal.Open(multicell_raster_name)
		multicell_raster_band = multicell_raster.GetRasterBand(1)
		multicell_raster_array = np.array(multicell_raster_band.ReadAsArray())
		multicell_ndv= multicell_raster_band.GetNoDataValue()
	
	strata_populations = {}
	
	if strata_raster_name:
		strata_raster = gdal.Open(strata_raster_name)
		strata_raster_band = strata_raster.GetRasterBand(1)
		strata_raster_array = np.array(strata_raster_band.ReadAsArray())
		strata_ndv= strata_raster_band.GetNoDataValue()
		
		strata_ids = np.unique(strata_raster_array[strata_raster_array !=strata_ndv])
		
		#print strata_ids

		
		for s_id in strata_ids:
			
			vals = pop_raster_array[ (strata_raster_array == s_id) & (pop_raster_array != pop_ndv)]
			population = int(round(sum(vals)))
			median_pop = np.median(vals)
			total_psus = len(vals)
			
			# if multicell we need to do a bit of work to get the values for each multicell PSU
			# unfortunately this is slow as it iterates each multicell id
			if not multicell_raster_array is None:
				vals = []
				mcs = np.unique(multicell_raster_array[ (strata_raster_array == s_id) & (pop_raster_array != pop_ndv)])
				total_psus = len(mcs)
				for mc in mcs:
					mc_pop = np.sum(pop_raster_array[ (strata_raster_array == s_id) & (pop_raster_array != pop_ndv) & (multicell_raster_array == mc)])
					vals.append(mc_pop)

				median_pop = np.median(vals)
			
			
			# make the dictionary value a tuple so we can add more info
			strata_populations[int(s_id)] = ( population,median_pop,total_psus )
			#print "total for strata " + str(int(s_id)) + ":" + str(int(round(sum(vals)))) 
			
		strata_raster = None
	else:
		vals = pop_raster_array[ pop_raster_array != pop_ndv]
		population = int(round(sum(vals)))
		median_pop = np.median(vals)
		total_psus = len(vals)
		
		# if multicell we need to do a bit of work to get the values for each multicell PSU
		# unfortunately this is slow as it iterates each multicell id
		if not multicell_raster_array is None:
			vals = []
			mcs = np.unique(multicell_raster_array[ (strata_raster_array == s_id) & (pop_raster_array != pop_ndv)])
			total_psus = len(mcs)
			for mc in mcs:
				mc_pop = np.sum(pop_raster_array[ (strata_raster_array == s_id) & (pop_raster_array != pop_ndv) & (multicell_raster_array == mc)])
				vals.append(mc_pop)

			median_pop = np.median(vals)

			
		# make the dictionary value a tuple so we can add more info
		strata_populations[1] = ( population,median_pop,total_psus )
	
	pop_raster = None
	


	return strata_populations
	
def dictionary_to_csv(csv_name,dict,required_fields= None):

	if enable_debug:
		print '## Entering dictionary_to_csv '
		
		print csv_name
		print dict
		
	# this expects a dictionary of dictionaries
	# extract values from dictionary and write to CSV
	# if field names given just those fields else all fields
	
	#print "entering dictionary_to_csv"
	
	if  os.path.isdir(os.path.dirname(csv_name) ):
		csv_file = open(csv_name,'w+')

		# create header and write header line
		# then add all fields
		
		keys_level_1 = dict.keys()
		keys_level_2 = dict[keys_level_1[0]].keys()
		
		header = ''
		header_fields = []
		# if no list specified then use them all
		if required_fields is None:
			header = 'id,'
			header += ','.join(keys_level_2)
			header_fields = keys_level_2
			#print header
		else:
			header = 'id'
			for rf in required_fields:
				if rf in keys_level_2:

					header += ',' + rf 
					header_fields.append(rf)			
		
		#print header
		#print header_fields
				
		csv_file.write(header + "\n")
		
		# now go through each row in dictionary and add to csv
		for id in keys_level_1:
			row = str(id) 
			
			for hf in header_fields:
				row += ',' + str(dict[id][hf])
		
			csv_file.write(row + "\n")
		
		csv_file.close()

		csv_file = None	
		
def getMulticellPopValues(pop_raster,multicell_id_raster):
	# for each id in the multicell_id raster (whihc corresponds to an area)
	# calcuate the population in that area and add to array
	
	if enable_debug:
		print "entering getMulticellPopValues"
	
	mcr = gdal.Open(multicell_id_raster)
	mcr_band = mcr.GetRasterBand(1)
	mcr_arr_vals = np.array(mcr_band.ReadAsArray())
	mcr_ndv = mcr_band.GetNoDataValue()
	
	pr = gdal.Open(pop_raster)
	pr_band = pr.GetRasterBand(1)
	pr_arr_vals = np.array(pr_band.ReadAsArray())
	pr_ndv = pr_band.GetNoDataValue()
	
	valuez = []

	#print "mcr_ndv " + str(mcr_ndv) 
	
#
	# this is slow: is there a faster way?
	for multicell_id in np.unique(mcr_arr_vals[mcr_arr_vals != mcr_ndv]):
		# we are accidentally including the ndv becasue of rounding so try and exclude
		if multicell_id < 0:
			continue
			
		v = np.sum( pr_arr_vals[(mcr_arr_vals==multicell_id) & (pr_arr_vals !=pr_ndv)] )
					
		#print str(multicell_id) + '--' + str(v)

		valuez.append(  round(v,1)   )
	
		# for testing
		#if multicell_id > 1000:
		#	break
		
	return valuez
		