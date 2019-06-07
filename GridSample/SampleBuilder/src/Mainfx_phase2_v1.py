# -*- coding: utf-8 -*-
"""
Created on 2018/2019
#!/usr/bin/env python3

@author: Ian Waldock
"""


from osgeo import gdal, ogr, gdalconst

import math
import numpy as np
import random

import string

import ogr2ogr

import os
import osr

import glob

import imghdr
import sample_builder_utils as sbu

import pdfkit

import csv

import datetime

import pickle

# the gridEA code is written in R which can be run from python using this module
import rpy2	
import rpy2.robjects as robjects

r=robjects.r
# we need to set the path to the gridEA file
r.source("/src/gridEA_code.R")
icw_gridEA =r['icw_gridEA']

# print out a bunch of debug
enable_debug = False

# for debug purposes we can choose not to delete the intermediate output files
delete_intermediate_outputs = True

# this function is intended for use as part of the GridSample web app
# It creates a summary pdf (via an intermediate html page) based on hte outputs of the gridsample function
def create_html_and_pdf_summary_page( 
				output_dir , # base directory for output
				uniqueID, # Unique user/project ID
				email , 
				proj_name,
				country,
				subnational_survey_options ,		
				pop_file_name,
				frame_type, 
				frame_type_param1, # if frame type:single -Single-cell: Input cell size (m X m),multi - gridEA: Cluster size,own-Own: Shapefile name 
				frame_type_param2,  # if frame type:single -Single-cell: Exclusion,multi - gridEA: Exclusion-Own: Shapefile name 
				stages,
				stratification_yn,
				oversample_yn,
				random_seed,
				oversample_size ,
				target_population_label,
				av_targ_pop_membs_per_hh,
				hh_size,
				coverage_area_names,
				
				selectedSampledHHsPerRuralPSU,
				selectedSampledHHsPerUrbanPSU,
				selectedAllocationPerStrata,
			
				coverage_shapefile_name = 'NA',
				strata_description = '',
				strata_description_long = ''
				
):
	if enable_debug:
		print "## entering function:create_html_and_pdf_summary_page"
		print output_dir
		print uniqueID
			
	status = 1	
	
	# Create a path for outputs to go in
	working_dir = output_dir + str(uniqueID)
	
	html_summary_file_name = working_dir + '/gridsample_summary.html'
	pdf_summary_file_name = working_dir + '/gridsample_summary.pdf'
	

	# top of page and css
	top = '<!DOCTYPE html><html>\
			<head>\
				<title>Page Title</title>\
				<style>\
					h1 {color: maroon;margin-left: 40px;} \
					table {width:100%;} \
					table, th, td {border-collapse: collapse;border: 1px solid black;} \
					td {padding-left:20px;} \
					img {width:50%;} \
				</style>\
			</head>\
		<body>'
		
	bottom = '</body></html>'
	
	try:
		# create the html summary file
		sf = open(html_summary_file_name,'w') 
 
		sf.write(top)
		
		# first page
		job_details = '<b>GridSample Job ID:' + str(uniqueID) + '</b><br/>'
		d = datetime.datetime.now()

		job_details += '<b>Date:' + d.strftime("%c") + '</b><br/>'
		job_details += '<b>Sample Name:' + proj_name + '</b><br/>'
		job_details += '<b>User:' + email + '</b><br/>'
		
		job_details = '<p>' + job_details + '</p>'
		
		job_details += '<p>' + 'This report summarizes the parameters and input datasets used to select clusters for the '
		
		job_details += '"' + proj_name + '" gridded population survey on ' + d.strftime("%c") + ' by ' + email + '.'  +'</p>'

		job_details += '<p>' + 'If your sample was not successful or did not produce usable cluster boundaries, you can update your survey parameters and/or input datasets and resubmit.' + '</p>'

		job_details += '<p>' + 'Contact the GridSample team with any questions:gridsample@flowminder.org' + '</p>'
	
		job_details += '<p style="page-break-after: always;">&nbsp;</p>'
		
		sf.write(job_details)
		
		
		
		# second page = table of parameters chosen
		
		params_chosen = '<h2>Parameters and input datasets<h2>'
		
		params_chosen += '<table>'
		
		params_chosen += '<tr><th>Tab</th><th>Field</th><th>Value</th></tr>'
		
		# coverage 
		params_chosen += '<tr><td>Coverage</td><td>Sample name</td><td>' + proj_name + '</td></tr>'
		
		params_chosen += '<tr><td>Coverage</td><td>Country</td><td>' + country + '</td></tr>'
		
		params_chosen += '<tr><td>Coverage</td><td>Subnational survey options</td><td>' + subnational_survey_options + '</td></tr>'
		
		params_chosen += '<tr><td>Coverage</td><td>Shapefile name</td><td>' + coverage_shapefile_name + '</td></tr>'
		
		# frame
		params_chosen += '<tr><td>Frame</td><td>WorldPop dataset</td><td>' + pop_file_name + '</td></tr>'
		
		params_chosen += '<tr><td>Frame</td><td>Frame type</td><td>' + frame_type + '</td></tr>'

		
		if frame_type == 'single':

			params_chosen += '<tr><td>Frame</td><td>Single-cell: Input cell size (m X m)</td><td>' + frame_type_param1 + '</td></tr>'
			params_chosen += '<tr><td>Frame</td><td>Single-cell: Exclusion</td><td>' + frame_type_param2 + '</td></tr>'
		elif frame_type == 'multi':

			params_chosen += '<tr><td>Frame</td><td>gridEA: Cluster size</td><td>' + frame_type_param1 + '</td></tr>'
			params_chosen += '<tr><td>Frame</td><td>gridEA: Exclusion</td><td>' + frame_type_param2 + '</td></tr>'
		else:   # frame_type == 'own':

			params_chosen += '<tr><td>Frame</td><td>Own: Shapefile name</td><td>' + frame_type_param1 + '</td></tr>'

			params_chosen += '<tr><td>Frame</td><td>Own: Unit ID</td><td>' + frame_type_param2 + '</td></tr>'

		# design
		params_chosen += '<tr><td>Design</td><td>Stages</td><td>' + stages + '</td></tr>'
		
		params_chosen += '<tr><td>Design</td><td>Stratification</td><td>' + stratification_yn + '</td></tr>'
		
		params_chosen += '<tr><td>Design</td><td>Spatial oversample</td><td>' + oversample_yn + '</td></tr>'
		
		params_chosen += '<tr><td>Design</td><td>Random Number</td><td>' + random_seed + '</td></tr>'
		
		
		
		#strata
		if stratification_yn == "Yes":
			params_chosen += '<tr><td>Strata</td><td>Admin area or shapefile name</td><td>' + strata_description + '</td></tr>'
		
		# spatial
		if oversample_yn == "Yes":
			params_chosen += '<tr><td>Spatial</td><td>Area (km X km)</td><td>' + str(oversample_size) + '</td></tr>'
		
		# target
		params_chosen += '<tr><td>Target</td><td>Target population name</td><td>' + target_population_label + '</td></tr>'
		
		params_chosen += '<tr><td>Target</td><td>Target population per household</td><td>' + str(av_targ_pop_membs_per_hh) + '</td></tr>'
		
		params_chosen += '<tr><td>Target</td><td>Average household size</td><td>' + str(hh_size) + '</td></tr>'
				
		params_chosen += '</table>'
		
		params_chosen += '<p style="page-break-after: always;">&nbsp;</p>'
		
		sf.write(params_chosen)

		
		# third page page = table of sample design parameters parameters chosen
		
		sample_design = '<h2>Sample design<h2>'
		
		strata_info_csv = working_dir + '/strata_areas_population_info.csv'
		
		if not os.path.exists(strata_info_csv):
			sample_design += '<p>Can\'t display sample design information</p>'
		else:
			sample_design += '<table>'
			
			sample_design += '<tr><th>Strata ID</th><th>Strata Name</th><th>Number of clusters</th><th>HHs sampled per urban cluster</th><th>HHs sampled per rural cluster</th><th>Total HH sample size</th></tr>'
				
			grand_total_clusters = 0
			grand_total_hhs = 0
			
			# we need to read strata ids and names from CSV and match strata ids with info from front-end to make the table
			with open(strata_info_csv) as csvfile:
				reader = csv.DictReader(csvfile)
				for row in reader:
					
					number_of_clusters = 0
					if isinstance(selectedAllocationPerStrata,dict):
						
						if isinstance(selectedAllocationPerStrata.keys()[0], int):
							
							number_of_clusters = int(selectedAllocationPerStrata[int(row['id'])])
						else:
							
							number_of_clusters = int(selectedAllocationPerStrata[str(row['id'])])
					else: 
						
						number_of_clusters = int(selectedAllocationPerStrata)					
					
					sample_design += '<tr><td>' + str(row['id']) + '</td>'
					sample_design += '<td>' + str(row['strata_name']) + '</td>'
				
					sample_design += '<td>' + str(number_of_clusters) + '</td>'
				
					sample_design += '<td>' + str(selectedSampledHHsPerUrbanPSU) + '</td>'
					sample_design += '<td>' + str(selectedSampledHHsPerRuralPSU) + '</td>'
	
					total_hh_sample_size = round(number_of_clusters*( selectedSampledHHsPerRuralPSU + selectedSampledHHsPerUrbanPSU )/2)
					

					sample_design += '<td>' + str(total_hh_sample_size) + '</td></tr>'
					
					grand_total_clusters += number_of_clusters
					grand_total_hhs += total_hh_sample_size

			sample_design += '<tr><td>Total</td><td></td><td>'+ str(grand_total_clusters) + '</td><td></td><td></td><td>' + str(grand_total_hhs) + '</td></tr>'
					
			
			sample_design += '</table>'
		
		
		
		sample_design += '<p style="page-break-after: always;">&nbsp;</p>'
	
		sf.write(sample_design)
		
		# strata page
				
		strata_page = '<h2>Sample coverage and strata boundaries</h2>'
		
		strata_page += '<p>The map below shows the strata units restricted to the coverage area boundary. The strata dataset is visualized as a raster with cells that reflect the size of the underlying gridded population dataset used to generate the sample frame.</p>'
		
		strata_page += '<h3>Coverage</h3>'
		
		if subnational_survey_options == 'custom_upload':
			text = 'A custom survey coverage area was defined in ' + country + ' using a shapefile called ' + coverage_shapefile_name
			strata_page += '<p>' + text + '</p>'
		elif subnational_survey_options == "urbanOnly":
			text = 'The survey is restricted to an area defined by urban only boundaries (value 3) in  ' + country + ' as defined by GHS-SMOD ' 
			text += 'and available at <a href="https://ghsl.jrc.ec.europa.eu/ghs_smod.php">https://ghsl.jrc.ec.europa.eu/ghs_smod.php</a>'
			strata_page += '<p>' + text + '</p>'
		elif  subnational_survey_options == "ruralOnly":
			text = 'The survey is restricted to an area defined by rural only boundaries (values 0,1,2) in  ' + country + ' as defined by GHS-SMOD ' 
			text += 'and available at <a href="https://ghsl.jrc.ec.europa.eu/ghs_smod.php">GHS-SMOD</a>'
			strata_page += '<p>' + text + '</p>'
		elif subnational_survey_options == "National":
			text = 'The survey is restricted to an area defined by national boundaries  in  ' + country + ' as defined by GADM ' 
			text += 'and available at <a href="https://gadm.org/download_country_v3.html">https://gadm.org/download_country_v3.html</a>'
			strata_page += '<p>' + text + '</p>'
		else:
			text = 'The survey is restricted to an area defined by ' + subnational_survey_options + ' boundaries  (' + coverage_area_names + ') in  ' + country + ' as defined by GADM ' 
			text += 'and available at <a href="https://gadm.org/download_country_v3.html">https://gadm.org/download_country_v3.html</a>'
			strata_page += '<p>' + text + '</p>'
		
		strata_page += '<h3>Strata</h3>'		
		
		if stratification_yn == "No":
			text = 'No strata were chosen for this sample so it was treated as having one strata over the coverage area'
			strata_page += '<p>' + text + '</p>'
		else:
			strata_page += '<p>' + strata_description_long + '</p>'
			
		# there should be a strata_raster.tif file in the working dir
		# if there is then convert to a png and add to the document
		strata_raster = working_dir + "/strata_raster.tif"
		strata_raster_png = working_dir + "/strata_raster.png"
		
		if not os.path.exists(strata_raster):
			strata_page += '<p>Can\'t display strata raster</p>'
		else:
			if not sbu.convertRasterToPNGWithTransparentNDV(strata_raster,strata_raster_png):
				strata_page += '<img src="' + strata_raster_png + '" alt="strata raster" >'
				
		strata_page += '<p style="page-break-after: always;">&nbsp;</p>'
				
		sf.write(strata_page)
		
		## sample frame page
		sample_frame_page = '<h2>Sample frame & exclusion</h2>'
		
		sample_frame_page += '<p>The map below shows the final sample frame units used to select the sample of PSUs with excluded areas. The histogram below shows the distribution of sample frame unit population values and, if applicable, the population values in excluded cells (note that cells might be of a different size than the rest of the sample frame if you have used the gridEA algorithm or uploaded your own shapefile to define sample frame units)</p>'
		
		sample_frame_page += '<h3>Frame</h3>'
		
		if frame_type == 'single':
			text = 'The final sample frame is from WorldPop and called pop_raster_cut.tif.' 	
			text += ' Cells were aggregated from approximately 100mX100m to ' + str(frame_type_param1) + 'm X ' + str(frame_type_param1) + 'm'
			sample_frame_page += '<p>' + text + '</p>'		
		elif frame_type == 'multi':
			text = 'The sample frame was generated with the gridEA algorithm which groups gridded population datasets into small areas with a target population and maximum area. The following datasets were passed to the gridEA algorithm:' 	
			text += '<ul>'
			text += '<li>WorldPop,'  + pop_file_name +' available at <a href="https://www.worldpop.org/project/categories?id=3">WorldPop Population Data</a></li>'
			text += '<li>GADM: The lowest-level administrative units available for ' + country + ', available at: <a href="https://gadm.org/download_country_v3.html">GADM</a></li>'
			text += '<li>GHS-SMOD: High density urban (value 3), low density urban (value 2), rural (value 1), and unsettled (value 0), available at:<a href="https://ghsl.jrc.ec.europa.eu/ghs_smod.php">GHS-SMOD</li>'

			text += '</ul>'
			sample_frame_page += '<p>' + text + '</p>'
		else:   # frame_type == 'own':
			text = 'The sample frame reflects the boundaries uploaded in ' + coverage_shapefile_name 	
			text += 'and population totals for each sample frame unit aggregated from the ' + pop_file_name + ' WorldPop dataset.'
			text += 'WorldPop data are available at <a href="https://www.worldpop.org/project/categories?id=3">WorldPop Population Data</a>'
			sample_frame_page += '<p>' + text + '</p>'

		sample_frame_page += '<h3>Exclusions</h3>'
		
		if frame_type == 'single':	
			if frame_type_param2 == 'None':
				text = 'No population was excluded from the sample frame.'
			elif frame_type_param2 == 'Exclude GHS-SMOD unsettled areas':
				text = 'Areas defined by GHS-SMOD to be unsettled (value 0) were excluded from the sample frame.' 
			else:
				text = 'Cells ' + str(frame_type_param1) + 'm X ' + str(frame_type_param1) + 'm with fewer than ' + str(frame_type_param2) + ' people are excluded from the sample frame'		
		elif frame_type == 'multi':
			if frame_type_param2 == 'Exclude GHS-SMOD unsettled areas':
				text = 'Areas defined by GHS-SMOD to be unsettled (value 0) were excluded from the sample frame.' 
			else:
				text = 'No population was excluded from the sample frame.'
		else:
			text = 'No population was excluded from the sample frame.'
		
		sample_frame_page += '<p>' + text + '</p>'
		
		# add images
		sample_frame_raster = working_dir
		sample_frame_raster_colorized = working_dir +"/sample_frame_raster_colourized.tif"	
		sample_frame_png = working_dir + "/sample_frame_raster.png"
		
		if frame_type == 'single':
			#sample_frame_raster +=  "/index_raster.tif"
			sample_frame_raster +=  "/pop_raster_cut.tif"
		elif frame_type == 'multi':	
			sample_frame_raster +=  "/gridEA/EA_raster_master1.tif"
		else:
			sample_frame_raster += '/multicell_psu_ids.tif'

		if not os.path.exists(sample_frame_raster):
			sample_frame_page += '<p>Can\'t display sample frame raster</p>'
		else:
			sbu.createColouredImageWithScale(sample_frame_raster,sample_frame_png)	
			if os.path.exists(sample_frame_png):
				sample_frame_page += '<p><img src="' + sample_frame_png + '" alt="sample frame raster" ></p>'		
			#sbu.createColouredTifWithScale(sample_frame_raster,sample_frame_raster_colorized)	
			#if not sbu.convertRasterToPNGWithTransparentNDV(sample_frame_raster_colorized,sample_frame_png):
			#	sample_frame_page += '<p><img src="' + sample_frame_png + '" alt="sample frame raster" ></p>'
				
		# add histogram images if htey exist
		histogram_image_above_cut_off =  working_dir +  '/pop_histogram_cut_off_above.png'
		histogram_image_below_cut_off =  working_dir +  '/pop_histogram_cut_off_below.png'
		
		if os.path.exists(histogram_image_above_cut_off):
			sample_frame_page += '<img src="' + histogram_image_above_cut_off + '" alt="histogram_image_above_cut_off" >'

		if os.path.exists(histogram_image_below_cut_off):
			sample_frame_page += '<img src="' + histogram_image_below_cut_off + '" alt="histogram_image_above_cut_off" >'
			
		sample_frame_page += '<p style="page-break-after: always;">&nbsp;</p>'
		
		sf.write(sample_frame_page)
		
		# spatial oversampling
		if oversample_yn == "Yes":
			spatial_oversampling_page = '<h2>Spatial oversampling</h2>'
		
			text = 'The map below shows the dataset generated and used for spatial oversampling, clipped to the coverage area boundary. The super-grid cells created for spatial oversampling were approximately'
			text += str(oversample_size) + 'km X ' + str(oversample_size) + 'km.<br/>'
			text += 'The GridSample algorithm ensured that at least one sampling unit (cluster) was selected per super-grid cell. This means that if a cluster was not selected inside a given super-grid cell during the main (PPS) sampling procedure, a cluster within the super-grid cell was selected at random during the spatial oversampling step and added to the sample.'
		
			spatial_oversampling_page += '<p>' + text + '</p>'

			# add image
			spatial_oversampling_raster = working_dir + '/spatial_oversample_raster.tif'
			spatial_oversampling_png = working_dir + "/spatial_oversample_raster.png"
			
			if not os.path.exists(spatial_oversampling_raster):
				spatial_oversampling_page += '<p>Can\'t display spatial oversampling raster</p>'
			else:
				if not sbu.convertRasterToPNGWithTransparentNDV(spatial_oversampling_raster,spatial_oversampling_png):
					spatial_oversampling_page += '<img src="' + spatial_oversampling_png + '" alt="sample frame raster" >'
						
			spatial_oversampling_page += '<p style="page-break-after: always;">&nbsp;</p>'
		
			sf.write(spatial_oversampling_page)
		
		# implicit stratification
		implicit_stratification_page = '<h2>Implicit Stratification</h2>'
		
		implicit_stratification_page += '<p>All samples are implicitly stratified using the GHS-SMOD dataset. This means that sampling units were ordered first by strata, second by GHS-SMOD class (high-density urban, low-density urban, rural, and unsettled) and third by geography from north-to-south, west-to-east before the PPS sampling procedure was performed.</p>'
		
		implicit_stratification_page += '<p>GHS-SMOD is available at <a href="https://ghsl.jrc.ec.europa.eu/ghs_smod.php">https://ghsl.jrc.ec.europa.eu/ghs_smod.php</a></p>'
		
		# add image
		ghs_mod_raster = working_dir + '/ghs_mod_cut.tif'
		ghs_mod_raster_colourized = working_dir + "/ghs_mod_cut_colourized.tif"
		ghs_mod_raster_png = working_dir + "/ghs_mod_cut.png"
			
		if not os.path.exists(ghs_mod_raster):
			spatial_oversampling_page += '<p>Can\'t display GHS-SMOD raster</p>'
		else:
			sbu.createColouredImageWithScale(ghs_mod_raster,ghs_mod_raster_png)		
			if os.path.exists(ghs_mod_raster_png):
				implicit_stratification_page += '<img src="' + ghs_mod_raster_png + '" alt="sample frame raster" >'	
			#sbu.createColouredTifWithScale(ghs_mod_raster,ghs_mod_raster_colourized)		
			#if not sbu.convertRasterToPNGWithTransparentNDV(ghs_mod_raster_colourized,ghs_mod_raster_png):
			#	implicit_stratification_page += '<img src="' + ghs_mod_raster_png + '" alt="sample frame raster" >'		
		
		implicit_stratification_page += '<p style="page-break-after: always;">&nbsp;</p>'
		
		sf.write(implicit_stratification_page)	
		
		sf.write(bottom)
 
		sf.close()
		
		# convert the summary file to pdf
		pdfkit.from_file(html_summary_file_name, pdf_summary_file_name)		
		
	
	except Exception as error:
		print "Exception in create_html_and_pdf_summary_page code:" + str(error) 
		status=0
		message = str(error)

			
	return status
	
# because the gridEA algorithm takes a long we don't want to run it unless the parameters have changed
# this function compares the parameters to the previous ones that have been stored in a file (upating the file if they have changed 	
def gridEADoesNotNeedRemaking(working_dir,pop_raster,grid_ea_coverage_info,cfg_exclude_ghssmod0_bool,cfg_multi_cell_cluster_size):

	if enable_debug:
		print "##entering function:gridEADoesNotNeedRemaking"
	

	doesNotNeedReMaking = True
	
	multicell_psu_raster = working_dir + '/gridEA/EA_raster_master1.tif'
	gridEA_last_run_parameters_file = working_dir + '/gridEA_last_run_parameters.txt'

	# string of parameters that we will compare with what's stored form hte last run
	gridEA_parameters_string = pop_raster + str(grid_ea_coverage_info) + str(cfg_exclude_ghssmod0_bool) + str(cfg_multi_cell_cluster_size)
	
	if not os.path.exists(multicell_psu_raster):
		# not multicell has ben made yet
		doesNotNeedReMaking = False		
	elif not os.path.exists(gridEA_last_run_parameters_file):
		# the parameters from the last run haven't been stored
		doesNotNeedReMaking = False
	else:
		f= open(gridEA_last_run_parameters_file,'r') 
		old_parameters = f.read()
		f.close()
		
		if old_parameters  !=  gridEA_parameters_string:
			doesNotNeedReMaking = False
	
	# reset the old parameters file	
	if doesNotNeedReMaking is False:

		f=open(gridEA_last_run_parameters_file,'w') 
		f.write(gridEA_parameters_string)
		f.close()
				
					
					
	return doesNotNeedReMaking

## this function originally was just to create the grid sample
## because there is now additional functionality that requires the pre-preparation of the data then an action option has been added
## options are : "create_sample" (the default) creates the grid sample
## "create_histogram": makes two images of histograms showing the distribution of the PSU cells populations
## "get_strata_pop_values": calculates the population in each strata and saves as a CSV
def gridsample( 
				output_dir, # Directory for output
				uniqueID, # Unique ID for this build/project -output will be put in a subdir of results_dir with this id
		
				ghs_mod_raster , # the location of the GHS raster - this is mandatory as it is used to do an u/r stratification behind the scenes
				
				pop_raster, # location of the population raster 

				# Either an integer which is the  psu required per strata and will be used as the psu required for each strata
				# or a dictionary with keys the strata id and value the required number of psu for that strata
				cfg_psu_per_strata  ,
				
				PSU_filename = "gridsample_output.shp", # Desired filename of output shapefile
										
				coverage_poly_filename = None, # filename of coverage polygon
				coverage_ur_option = None, # urbanOnly, ruralOnly : uses the national boundary but only includes u/r according to GHS-MOD
				
				strata_poly_filename = None, # strata shapefile object
				strata_ID_field = None, # string: name of variable that defines unique strata (ex. "ID")  (in strata shape file)
				strata_name_field = None, # string: name of variable that defines strata name (in strata shape file)
			
				  
				cfg_stratification_method = None, # None for no stratification , else one of urbanRural,adminArea,custom_upload
				
				# Size of each household. This value will be used to divide the population raster, so that all calculations are being done in terms of households
				# either an integer if all strata have the same value or a dictionary of strata id and values
				# if not set we just work in terms of population
				# mandatory when making a grid sample : not required for histogram or strata population calls
				cfg_hh_size = None, 
				
				
				cfg_exclude_pop_per_cell = None, # Minimum population in a cell that will be used to exclude that cell entirely
				cfg_exclude_ghssmod0_bool = False, # whether to exclude cells that GHS_MOD say are unpopulated
				
				cfg_random_number = None, # Random number seed for code
				
					
				## whether using a standard population raster with each cell a PSU, a shapefile to indicate PSU, or a population raster and a multicell raster to indicate the PSUs 
				# one of single,own, or multi
				cfg_frame_type = 'single',

				cfg_resample_size = None ,#number of cells to aggregate for resampling all rasters (None means don't resample) -- based on 100m x 100m, so 3 = aggregate to 300m x 300m, 100 = aggregate to 10km x 10km 
				
				action = "create_sample" , # options are create_sample, create_histogram, get_strata_pop_values
				
				cfg_oversample_grid_spatial_scale = None, # Size of cell to be used to make an oversized grid that we will ensure at least one cell of each has been selected from
				
				# if cfg_frame_type is "own" then the location of the shapefle containing the PSUs.
				# if cfg_frame_type is "multi" then either the location of the raster containing the cell ids or nothing in which case cfg_multi_cell_cluster_size must be set and the file will be geenrated by running the gridEA algorithm
				cfg_own_frame_file = None,
				
				# if cfg_frame_type is "own" then the name of the shapefile id columns for the shapefile containing the PSUs.
				cfg_own_frame_id = None,
				
				# if cfg_frame_type is "multi" then one of small,medium large
				cfg_multi_cell_cluster_size = None,
				
				# if user chooses multi-cell (GridEA) then they can optionally give a shapefile indicating the strata to use
				# if none given the entire country is treated as one strata (which might causes the algorithm problems)
				grid_ea_strata_file = None	,
				# if a grid_ea_strata_file the the name of the id field
				grid_ea_strata_file_id_field = None,
				# used to help decide if parameters have changed to decide whether it is necessary to rebuild gridEA as it's slow
				grid_ea_coverage_info = None
				
				):
				
				
		
	if enable_debug:
		print "##entering function:gridsample"
			
	status = 1	 
	message = ''
	
	# we will use a file to indicate whether histograms are currently being made
	# we set it when they are under construction and remove it when created
	working_dir = output_dir + str(uniqueID)
	pop_histogram_currently_being_made_file = working_dir +  '/pop_histogram_currently_being_made.txt'
	
	try:
	
		#################################################################
		############### initial set-up
							 
		default_pixel_size = 0.000833333
		
		# for testing to avoid re-creating all the intermediate rasters/files to make things quicker
		recreate_existing_raster = True
		
		# Create a path for outputs to go in
		if not os.path.exists(working_dir):
			os.makedirs(working_dir)
						
			
		###############################################################################
		########################## check parameters make sense
		
		## Different sets of parameters are checked  depending on the chosen action
		
		if strata_ID_field:
			strata_ID_field = str(strata_ID_field)
			
		if strata_name_field:
			strata_name_field = str(strata_name_field)
		
	
		if action == "create_sample" or action == "create_histogram" or action == "get_strata_pop_values":
			if not PSU_filename:
				status = 0
				message += "PSU_filename not given\n"
			elif PSU_filename[-4:] != ".shp":
				status = 0
				message += "PSU_filename must end in .shp\n"
			
			# frame type must be single :e.g using a pop raster to be able to excldue cells based on population or to resample
			#if (cfg_exclude_pop_per_cell or cfg_exclude_ghssmod0_bool or cfg_resample_size) and not cfg_frame_type == "single":
			#	status = 0
			#	message += "cfg_frame_type must be 'single if any of cfg_exclude_pop_per_cell or cfg_exclude_ghssmod0_bool or cfg_resample_size are chosen' \n"
			
			if	cfg_frame_type == "own" and ( not cfg_own_frame_id or not cfg_own_frame_file):		
				status = 0
				message += "cfg_own_frame_id and cfg_own_frame_file must both be set when cfg_frame_type is own' \n"
			
			if	cfg_frame_type == "multi" and ( not cfg_own_frame_file and not cfg_multi_cell_cluster_size):		
				status = 0
				message += "either cfg_own_frame_file or cfg_multi_cell_cluster_size must  be set when cfg_frame_type is multi' \n"
				
			if cfg_resample_size:
				# this must be a postive integer
				try:
					cfg_resample_size = int(cfg_resample_size)
					
					if cfg_resample_size < 1:
						status = 0
						message += "cfg_resample_size must be a positive integer 1"
				except:
					status = 0
					message += "cfg_resample_size must be a positive integer 2"
					

			
			# only one of these two options is allowed
			if cfg_exclude_pop_per_cell and cfg_exclude_ghssmod0_bool:
				status = 0
				message += "Values for both arr_cfg_exclude_pop_per_cell and arr_cfg_exclude_ghssmod0_bool have been given. Only one should be chosen\n"			
		
			# if ghs mod no pop cells are to be excluded then a ghs mod raster must be provided
			# a superfluous check as currently ghs mod is mandatory as it is used to stratify by urban rural		
			if cfg_exclude_ghssmod0_bool and not ghs_mod_raster:
				status = 0
				message += "cfg_exclude_ghssmod0_bool has been chosen but there is no value for ghs_mod_raster\n"			
				
			## this is the optional minimum population per cell
			if cfg_exclude_pop_per_cell:
				try:
					cfg_exclude_pop_per_cell = float(cfg_exclude_pop_per_cell)
				except:
					status = 0
					message += "cfg_exclude_pop_per_cell is not a number\n"
					
				if  status:			
					if cfg_exclude_pop_per_cell <= 0:
						status = 0
						message += "cfg_exclude_pop_per_cell must be greater than 0\n"

		if action == "create_sample" or action == "get_strata_pop_values":					
			if cfg_oversample_grid_spatial_scale:
				# this must be a postive integer
				try:
					cfg_oversample_grid_spatial_scale = int(cfg_oversample_grid_spatial_scale)
					
					if cfg_oversample_grid_spatial_scale < 1:
						status = 0
						message += "cfg_oversample_grid_spatial_scale must be a positive integer"
				except:
					status = 0
					message += "cfg_oversample_grid_spatial_scale must be a positive integer"
					
			if cfg_resample_size and cfg_oversample_grid_spatial_scale and status:
				#  cfg_oversample_grid_spatial_scale should be a multiple of  cfg_resample_size
				cfg_oversample_grid_spatial_scale = int(cfg_oversample_grid_spatial_scale/cfg_resample_size)
				if cfg_oversample_grid_spatial_scale < 2:
					status = 0
					message += "cfg_oversample_grid_spatial_scale should be a multiple of  cfg_resample_size"		
			
			# if stratification then one of 3 stratification options must be given 
			# a shapefile with stratification method  plus id and name field are required
			if  cfg_stratification_method and cfg_stratification_method !="urbanRural" and (not strata_poly_filename or not strata_ID_field or not strata_name_field):
				status = 0
				message += "stratification has been chosen but at least one of the 3 required paramters is missing (strata_poly_filename or strata_ID_field or strata_name_field) \n"			
		
		if action == "create_sample":
			if  not cfg_psu_per_strata:
				status = 0
				message += "No parameter cfg_psu_per_strata has been given\n"
						
			if not cfg_hh_size:  
				status = 0
				message += "cfg_hh_size is mandatory when creating a sample\n"	

			if cfg_random_number:  
				# this is not mandatory but if it exists then it should be an integer
				if not isinstance(cfg_random_number, int) :
					try:
						cfg_random_number = int(cfg_random_number)
					except:
						status = 0
						message += "cfg_random_number must be an integer\n"	

		
		#### if there was an error with the parameters go no further	
		if not status:
			if os.path.exists(pop_histogram_currently_being_made_file):
				os.remove(pop_histogram_currently_being_made_file)
			return [status,message]

		################################### all the parameters are OK	
		################### prepare files for analysis
		### as there are 3 possible different outputs which don't require all steps we return at different places
	
		##  check pop file exists, and is a tif, and if so reproject it targeting aligned pixels
		## we should probably use the file name of the raster itself so we don't need to recreate it each time
		## NOTE: the code to target aligned pixels has been commented out due to a requirement not to change the source population raster
		## It can be re-enabled if required
		pop_raster_tapped = working_dir + "/pop_raster_tapped.tif"
		pop_raster_final = pop_raster_tapped
		if not os.path.exists(pop_raster):
			status = 0
			message = "Can't find pop_raster in location given\n" + pop_raster
		elif ( imghdr.what(pop_raster) != 'tiff'):
			status = 0
			message += "pop_raster in not a tiff file\n"
		elif recreate_existing_raster or not os.path.exists(pop_raster_tapped):
			# NOW DON'T CHANGE ORIGIN OF wp RASTER
			#ret = sbu.reprojectToWGS84DegreesTargetAlignedPixel(pop_raster,pop_raster_tapped,tap=True)
			#if ret:
			#	status = 0
			#	message += " Could not reproject raster\n"
			# so we just set the final raster to point to the original one and don't manipulate it
			pop_raster_final = pop_raster		
		if not status:
			if os.path.exists(pop_histogram_currently_being_made_file):
				os.remove(pop_histogram_currently_being_made_file)
			return [status,message]
			
			
		## check if there is a coverage file, and if so that it is a shapefile, then reproject as required
		## if there is no coverage file then it is national coverage (in which case maybe we should pass in the national boundary ?)	
		coverage_poly_reprojected = working_dir + "/coverage_poly_reprojected.shp"
		coverage_poly_boundary = working_dir + "/coverage_poly_boundary.shp"
		coverage_poly_boundary_test = working_dir + "/coverage_poly_boundary_test.shp"
		driver = ogr.GetDriverByName("ESRI Shapefile")
		if coverage_poly_filename:
			if not os.path.exists(coverage_poly_filename):
				status = 0
				message = "Can't find coverage poly in location given\n" + str(coverage_poly_filename)
			elif not driver.Open(coverage_poly_filename):
				status = 0
				message = "Couldn't open coverage shapefile\n"
			elif recreate_existing_raster or not os.path.exists(coverage_poly_reprojected):
				ret = sbu.reprojectShapefileToWGS84(coverage_poly_filename,coverage_poly_reprojected)
				if ret:
					status = 0
					message += " Could not reproject coverage shapefile\n"
				else:
					# assume this always works - having problems with making a boundary shapefile from coverage shapefile
					#sbu.make_boundary_shapefile(coverage_poly_reprojected,coverage_poly_boundary,None)
					#sbu.make_boundary_shapefile(coverage_poly_reprojected,coverage_poly_boundary_test,0.01)
					#coverage_poly_reprojected = coverage_poly_boundary
					pass
					
		if not status:
			if os.path.exists(pop_histogram_currently_being_made_file):
				os.remove(pop_histogram_currently_being_made_file)
			return [status,message]
			
		## if there is a coverage shapefile then cut the population raster to the coverage shapefile
		if coverage_poly_filename:
			pop_raster_cut = working_dir + "/pop_raster_cut.tif"
			# commented out as rather than use target aligned pixels option we align to the population raster
			#ret = sbu.cropRasterWithShapefileMask(coverage_poly_reprojected,pop_raster_tapped,pop_raster_cut,tap=True)
			ret = sbu.cropRasterWithShapefileMask(coverage_poly_reprojected,pop_raster,pop_raster_cut,other=' -dstnodata -99 ')
			if recreate_existing_raster or not os.path.exists(pop_raster_cut):
				if ret:
					status = 0
					message += " Could not cut raster\n"
				else:
					pop_raster_final = pop_raster_cut
			else:
				pop_raster_final = pop_raster_cut

		if not status:
			if os.path.exists(pop_histogram_currently_being_made_file):
				os.remove(pop_histogram_currently_being_made_file)
			return [status,message]
			
		## if there is a coverage shapefile then cut GHS-MOD raster to the coverage shapefile
		## otherwise to the population raster
		ghs_mod_raster_cut = working_dir + "/ghs_mod_cut.tif"
		if not ghs_mod_raster:
			# no ghs mod raster
			status = 0
			message += " No location for the GHS MOD raster was passed in\n"
		elif not os.path.exists(ghs_mod_raster):
			# no ghs mod raster
			status = 0
			message += " Could not find GHS MOD raster\n"
		elif coverage_poly_filename :
			# a coverage poly exists
			if recreate_existing_raster or not os.path.exists(ghs_mod_raster_cut):
				ghs_mod_raster_cut_intermediate = working_dir + "/ghs_mod_cut_intermediate.tif"
				if sbu.clipRasterWithRasterMask2(pop_raster_final,ghs_mod_raster,ghs_mod_raster_cut_intermediate):
					status = 0
					message += " Could not cut raster\n"
			
				# set ndv to -99 and change pixel type to int16
				#ret = sbu.cropRasterWithShapefileMask(coverage_poly_reprojected,ghs_mod_raster,ghs_mod_raster_cut,pixel_size = 0.0008333,other=' -ot Int16 -dstnodata -99')
				ret = sbu.cropRasterWithShapefileMask(coverage_poly_reprojected,ghs_mod_raster_cut_intermediate,ghs_mod_raster_cut,pixel_size = 0.0008333,other=' -ot Int16 -dstnodata -99')		
				
				if ret:
					status = 0
					message += " Could not cut raster\n"
		else:
			# cut to pop raster
			if recreate_existing_raster or not os.path.exists(ghs_mod_raster_cut):
				# note this is cutting to the bounds of the mask it can therefore include ghs-mod pixels from another country
				# just be careful not to include them in later calcs - or set to nodata 
				if sbu.clipRasterWithRasterMask2(pop_raster_final,ghs_mod_raster,ghs_mod_raster_cut):
					status = 0
					message += " Could not cut raster\n"			
		
		if not status:
			if os.path.exists(pop_histogram_currently_being_made_file):
				os.remove(pop_histogram_currently_being_made_file)
			return [status,message]
		

		if coverage_ur_option:
			# we need to restrict to either just urban or just rural areas		
			driver = gdal.GetDriverByName('GTiff')
	
			p_r = gdal.Open(pop_raster_final,gdal.GA_Update)
			p_r_band = p_r.GetRasterBand(1)
			p_r_ndv = p_r_band.GetNoDataValue()
			p_r_arr = np.array(p_r_band.ReadAsArray())
			
			ghs_r = gdal.Open(ghs_mod_raster_cut)
			ghs_r_band = ghs_r.GetRasterBand(1)
			ghs_r_ndv = ghs_r_band.GetNoDataValue()
			ghs_r_arr = np.array(ghs_r_band.ReadAsArray())
			
			# set the rural or urban values to ndv depending on option chosen
			if coverage_ur_option == "urbanOnly":
				# set rural to ndv
				p_r_arr[ (ghs_r_arr == 0)  | (ghs_r_arr == 1) | (ghs_r_arr == 2) ] = p_r_ndv
			else :
				# set urban to ndv
				p_r_arr[ (ghs_r_arr == 3) ] = p_r_ndv
				
			p_r_band.WriteArray(p_r_arr)
			p_r_band.FlushCache()
			p_r = None
		
			
		if cfg_resample_size and cfg_resample_size > 1:
			# resample the pop and ghs rasters: all others are made based on these so no more resampling should be necessary!	
			default_pixel_size = cfg_resample_size*default_pixel_size
			
			# pop raster
			pop_raster_final_resampled = working_dir + "/pop_raster_cut_resampled.tif"
			
			if recreate_existing_raster or not os.path.exists(pop_raster_final_resampled):
				ret = sbu.rescale_raster_sum(pop_raster_final,pop_raster_final_resampled,cfg_resample_size)
				if ret:
					status = 0
					message += " Could not scale raster\n"
				else:
					pop_raster_final = pop_raster_final_resampled
			else:
				pop_raster_final = pop_raster_final_resampled
			
			if not status:
				if os.path.exists(pop_histogram_currently_being_made_file):
					os.remove(pop_histogram_currently_being_made_file)
				return [status,message]
			
			# ghs raster
			ghs_mod_raster_cut_resampled = working_dir + "/ghs_mod_cut_resampled.tif"	
			
			if recreate_existing_raster or not os.path.exists(ghs_mod_raster_cut_resampled):
				# specifying pixel size doesn't seem to guarantee the same output raster size
				#ret = sbu.rescale_raster_mode(ghs_mod_raster_cut,ghs_mod_raster_cut_resampled,cfg_resample_size)
				# so do this alternate way to guarantee raster size
				pop_r = gdal.Open(pop_raster_final)
				pop_r_band = pop_r.GetRasterBand(1)
				pop_r_arr = np.array(pop_r_band.ReadAsArray())
				h = pop_r_arr.shape[0]
				w = pop_r_arr.shape[1]
				
				ret = sbu.rescale_raster_mode_by_width_and_height(ghs_mod_raster_cut,ghs_mod_raster_cut_resampled,w,h)
				if ret:
					status = 0
					message += " Could not scale raster\n"
				else:
					ghs_mod_raster_cut = ghs_mod_raster_cut_resampled
			else:
				ghs_mod_raster_cut = ghs_mod_raster_cut_resampled
			
			if not status:
				if os.path.exists(pop_histogram_currently_being_made_file):
					os.remove(pop_histogram_currently_being_made_file)
				return [status,message]				
		
		# If we have spatial oversampling then it should be on the entire coverage area without excluding cells which happens next
		pop_raster_without_exclusions = pop_raster_final
		cfg_oversample_grid_spatial_scale
		
		# we may exclude cells below some threshold from the population raster 
		# if we do we keep a raster of the excluded cells for use in creating a histogram later
		pop_raster_excluded_cells = working_dir + "/pop_raster_excluded.tif"	
		if cfg_exclude_ghssmod0_bool:
			# we need to set the pop_raster cells that ghs mod says are empty to ndv - or 0 
			# we are also setting the ghs mod empty cells to ndv though it probably isn't necessary as it should be handled automatically later
			# it is handled later so we don't change GHS_mod
			if cfg_oversample_grid_spatial_scale:
				# make a copy of the raster without exclusion for use in spatial oversampling work
				pop_raster_without_exclusions = working_dir + "/pop_raster_no_exclusions.tif"
				sbu.make_raster_copy(pop_raster_final,pop_raster_without_exclusions)
			
			# make a copy of the final raster to keep excluded cells values
			if recreate_existing_raster or not os.path.exists(pop_raster_excluded_cells):
				sbu.make_raster_copy(pop_raster_final,pop_raster_excluded_cells)
			
				# we know the files exist - possibly they are already open from previous stages
				ghs_r = gdal.Open(ghs_mod_raster_cut,gdal.GA_Update)
				pop_r = gdal.Open(pop_raster_final,gdal.GA_Update)
				pop_r_excl = gdal.Open(pop_raster_excluded_cells,gdal.GA_Update)
				
				ghs_r_band = ghs_r.GetRasterBand(1)
				ghs_r_arr = np.array(ghs_r_band.ReadAsArray())
				ghs_r_ndv = ghs_r_band.GetNoDataValue()
				
				pop_r_band = pop_r.GetRasterBand(1)
				pop_r_arr = np.array(pop_r_band.ReadAsArray())
				pop_r_ndv = pop_r_band.GetNoDataValue()
				
				pop_r_excl_band = pop_r_excl.GetRasterBand(1)
				pop_r_excl_arr = np.array(pop_r_excl_band.ReadAsArray())
				pop_r_excl_ndv = pop_r_excl_band.GetNoDataValue()
				
				# set pop raster values in un-pop areas to ndv
				pop_r_arr[ghs_r_arr == 0] = pop_r_ndv
				pop_r_band.WriteArray(pop_r_arr)	
				pop_r_band.FlushCache()
				
				# check if all have been set to ndv seems unlikely but could be if strange areas chosen
				if not len (pop_r_arr[pop_r_arr != pop_r_ndv]):
					status = 0
					message += " The entire population is excluded when excluding GHS-SMOD unsettled areas"

				# set excluded raster values in GHS populated areas to ndv
				pop_r_excl_arr[ghs_r_arr > 0] = pop_r_excl_ndv
				pop_r_excl_band.WriteArray(pop_r_excl_arr)	
				pop_r_excl_band.FlushCache()

				## close everything to save resources and make sure actually written to disc before next stage
				pop_r_arr = None
				ghs_r_arr = None
				pop_r_excl_arr = None
				ghs_r = None
				pop_r = None
				pop_r_excl = None
			
		if cfg_exclude_pop_per_cell:
			# set the pop raster cells that are less than this value to ndv			
			if cfg_oversample_grid_spatial_scale:
				# make a copy of the raster without exclusion for use in spatial oversampling work
				pop_raster_without_exclusions = working_dir + "/pop_raster_no_exclusions.tif"
				sbu.make_raster_copy(pop_raster_final,pop_raster_without_exclusions)
			
			# make a copy of the final raster to keep excluded cells values
			if recreate_existing_raster or not os.path.exists(pop_raster_excluded_cells):
				sbu.make_raster_copy(pop_raster_final,pop_raster_excluded_cells)
				
				# we know the files exist - possibly there are already open from previous stages
				pop_r = gdal.Open(pop_raster_final,gdal.GA_Update)
				pop_r_excl = gdal.Open(pop_raster_excluded_cells,gdal.GA_Update)
				
				pop_r_band = pop_r.GetRasterBand(1)
				pop_r_arr = np.array(pop_r_band.ReadAsArray())
				pop_r_ndv = pop_r_band.GetNoDataValue()
				
				pop_r_excl_band = pop_r_excl.GetRasterBand(1)
				pop_r_excl_arr = np.array(pop_r_excl_band.ReadAsArray())
				pop_r_excl_ndv = pop_r_excl_band.GetNoDataValue()
				
				# set pop raster values in cells less than the cut-off to ndv
				pop_r_arr[pop_r_arr < cfg_exclude_pop_per_cell] = pop_r_ndv	
				pop_r_band.WriteArray(pop_r_arr)
				pop_r_band.FlushCache()
				
				# check if all have been set to ndv - which could happen if user has chosen a very high value!!!
				if not len (pop_r_arr[pop_r_arr != pop_r_ndv]):
					status = 0
					message += " The entire population is excluded when the value of cfg_exclude_pop_per_cell is set to " + str(cfg_exclude_pop_per_cell)
			
				# set pop raster values  in cells greater than or equal the cutoff to ndv
				pop_r_excl_arr[pop_r_arr >= cfg_exclude_pop_per_cell] = pop_r_excl_ndv
				pop_r_excl_band.WriteArray(pop_r_excl_arr)	
				pop_r_excl_band.FlushCache()

				## close everything to save resources and make sure actually written to disc before next stage
				pop_r_arr = None
				pop_r_excl_arr = None
				pop_r = None
				pop_r_excl = None

		if not status:
			if os.path.exists(pop_histogram_currently_being_made_file):
				os.remove(pop_histogram_currently_being_made_file)
			return [status,message]
			
		# if user has defined their own PSU via a shapefile we need to make a multicell type raster indicating which cells are which PSU
		# so both multicell and own shapefile reuslt in 2 input rasters the pop raster and an id raster indicating the PSU groups
		multicell_psu_raster = None
		if cfg_frame_type == 'own':
			if enable_debug:
				print '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ own make multicell'
			# user has uploaded their own shapefile
			reprojected_cfg_own_frame_file = working_dir + "/own_PSUs_reprojected.shp"
			cut_cfg_own_frame_file = working_dir + "/own_PSUs_cut.shp"
			# this is deliberately the same as the above as a default 
			own_PSU_shapefile = working_dir + "/own_PSUs_cut.shp"
			multicell_psu_raster = working_dir + "/multicell_psu_ids.tif"
			
			# check shapefile exists
			if not os.path.exists(cfg_own_frame_file):
				status = 0
				message += " Could not find 'own' PSU shapefile at location given"
			else:
				if enable_debug:
					print "recreate_existing_raster" + str(recreate_existing_raster)
					print "multicell_psu_raster" + str(multicell_psu_raster)
				if recreate_existing_raster or not os.path.exists(multicell_psu_raster): 
					# check the id attribute exists
					shapefile_attributes = sbu.get_shapefile_attributes(cfg_own_frame_file)

					if not cfg_own_frame_id in shapefile_attributes:
						status = 0
						message += " cfg_own_frame_id does not exist as an attribute in the PSU shapefile"				
					else:
						# after reprojecting rasterize shapefile based on id field - to the size of the final cut pop file
						# if we do this prior to resampling then we should include all areas by using 100m*100m scale - unless there are some very small shapes 
						if not sbu.reprojectShapefileToWGS84(cfg_own_frame_file,reprojected_cfg_own_frame_file):
							# we expect there to be a coverage poly but it isn't mandatory
							if coverage_poly_filename:
								if not sbu.clipShapefilewithShapefile(reprojected_cfg_own_frame_file,cut_cfg_own_frame_file,coverage_poly_reprojected):
									own_PSU_shapefile = cut_cfg_own_frame_file
								else:
									status = 0
									message += " Problem clipping own PSU poly "

							else:
								own_PSU_shapefile= reprojected_cfg_own_frame_file
								
							if own_PSU_shapefile:
								sbu.make_raster_copy(pop_raster_final,multicell_psu_raster,False,1,-99,True)
								sbu.rasterize_shapefile_with_gdal(own_PSU_shapefile,multicell_psu_raster,cfg_own_frame_id,-99,True,default_pixel_size)
						else:
							status = 0
							message += " Problem reprojecting own PSU poly"		
		elif cfg_frame_type == 'multi':
			if enable_debug:
				print '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ doing multicell stuff'
			
			# user either gives their own file as input or we create one
			if cfg_own_frame_file:
				multicell_psu_raster = working_dir + "/multicell_psu_ids.tif"
				multicell_psu_raster_tapped = working_dir + "/multicell_psu_ids_tapped.tif"
				
				# as we aren't sure this will be aligned to coverage area etc, tap and cut to pop raster size
				if recreate_existing_raster or not os.path.exists(multicell_psu_raster_tapped):
					ret = sbu.reprojectToWGS84DegreesTargetAlignedPixel(cfg_own_frame_file,multicell_psu_raster_tapped,no_data_value = -99)
					if ret:
						status = 0
						message += " Could not reproject multicell raster\n"
				
				if status and (recreate_existing_raster or not os.path.exists(multicell_psu_raster)):
					if sbu.clipRasterWithRasterMask2(pop_raster_final,multicell_psu_raster_tapped,multicell_psu_raster):
						status = 0
						message += " Could not cut multicell raster\n"	
			elif cfg_multi_cell_cluster_size:
				if enable_debug:
					print '@@@@@@ About to run grid EA algorithm'
				
				# although the algorithm accepts different strata and runs faster we won't be using strata so need to create a strata file of 1s
				strata_raster_for_multicell = working_dir +  '/strata_raster_for_multicell.tif'
						
				reprojected_grid_ea_strata_file = working_dir + "/reprojected_grid_ea_strata_file.shp"
				cut_grid_ea_strata_file = working_dir + "/cut_grid_ea_strata_file.shp"
				
				# a shapefile of admin areas to create a strata raster from might have been given 
				# otherwise we create a single strata raster of 1's for the entire country
				if gridEADoesNotNeedRemaking(working_dir,pop_raster,grid_ea_coverage_info,cfg_exclude_ghssmod0_bool,cfg_multi_cell_cluster_size):
					# it's a slow process so don't remake if unnecessary
					multicell_psu_raster = working_dir + '/gridEA/EA_raster_master1.tif'
				elif grid_ea_strata_file and grid_ea_strata_file_id_field:
					if os.path.exists(grid_ea_strata_file):
						if not sbu.reprojectShapefileToWGS84(grid_ea_strata_file,reprojected_grid_ea_strata_file):	
							# we expect there to be a coverage poly but it isn't mandatory
							if not sbu.clipShapefilewithShapefile(reprojected_grid_ea_strata_file,cut_grid_ea_strata_file,coverage_poly_reprojected):
								sbu.make_raster_copy(pop_raster_final,strata_raster_for_multicell,False,1,-99,True)
								sbu.rasterize_shapefile_with_gdal(cut_grid_ea_strata_file,strata_raster_for_multicell,'ID',-99,True,default_pixel_size)
													
								if enable_debug:
									print "################## GridEA parameters: ##########################"
									print "pop_raster_final:" + pop_raster_final
									print "ghs_mod_raster_cut:" + ghs_mod_raster_cut
									print "strata_raster_for_multicell:" + strata_raster_for_multicell
									print "cfg_multi_cell_cluster_size:" + str(cfg_multi_cell_cluster_size)
									print "cfg_exclude_ghssmod0_bool:" + str(cfg_exclude_ghssmod0_bool)
								
								res = icw_gridEA(pop_raster_final, ghs_mod_raster_cut, strata_raster_for_multicell, 
									working_dir + '/gridEA',cfg_multi_cell_cluster_size,cfg_exclude_ghssmod0_bool)
									
								res =  res.r_repr()
				
								if res == '0':
									status = 0
									message += " problem running Grid EA algorithm\n"	
								else:
									multicell_psu_raster = working_dir + '/gridEA/EA_raster_master1.tif'
							else:
								status = 0
								message += " Problem clipping grid_ea_strata_file to coverage area "	
						else:
							status = 0
							message += " Problem reprojecting grid_ea_strata_file "						
					else:
						status = 0
						message += " multicell admin area shapefile not found in location given\n"
					
				elif recreate_existing_raster or not os.path.exists(strata_raster_for_multicell):
					
					sbu.make_raster_copy(pop_raster_final,strata_raster_for_multicell,False,1,-99,True)
								
					res = icw_gridEA(pop_raster_final, ghs_mod_raster_cut, strata_raster_for_multicell, 
						working_dir + '/gridEA',cfg_multi_cell_cluster_size,cfg_exclude_ghssmod0_bool)
						
					res =  res.r_repr()
	
					if res == '0':
						status = 0
						message += " problem running Grid EA algorithm\n"	
					else:
						multicell_psu_raster = working_dir + '/gridEA/EA_raster_master1.tif'
				else:
					multicell_psu_raster = working_dir + '/gridEA/EA_raster_master1.tif'
			else:
				status = 0
				message += " multicell chosen but neither of a frame file or multicell cluster size were given\n"	
			
		if not status:
			if os.path.exists(pop_histogram_currently_being_made_file):
				os.remove(pop_histogram_currently_being_made_file)
			return [status,message]
				
		############ if creating a histogram we don't need to do anymore preparation so just get on and do it
		if action  == "create_histogram":
			# we don't do much error checking here so put in a try except section
			if enable_debug:
				print '################# creating histogram'
			try:
				histogram_image_above_cut_off =  working_dir +  '/pop_histogram_cut_off_above.png'
				histogram_image_below_cut_off =  working_dir +  '/pop_histogram_cut_off_below.png'
				
				# if they exist then delete them first
				if os.path.exists(histogram_image_above_cut_off):
					os.remove(histogram_image_above_cut_off)
					
				if os.path.exists(histogram_image_below_cut_off):
					os.remove(histogram_image_below_cut_off)
			
				vals = None
				
				if cfg_frame_type == 'single':
					pop_raster = gdal.Open(pop_raster_final)
					pop_raster_band = pop_raster.GetRasterBand(1)
					pop_raster_array = np.array(pop_raster_band.ReadAsArray())
					pop_ndv= pop_raster_band.GetNoDataValue()
			
					# above cut-off
					vals = pop_raster_array[ pop_raster_array != pop_ndv]
				elif cfg_frame_type == 'multi' or cfg_frame_type == 'own':
					vals = sbu.getMulticellPopValues(pop_raster_final,multicell_psu_raster)
				
				
				median_value = 0
				# there may not be any included cells
				if len(vals):	
					if enable_debug:
						print "Making Histogram:included cells"
					m = np.median(vals)
					median_value = m
					info = "Median: " + str(m) + "\n"
					info += "Min: " + str(np.min(vals)) + "\n"
					info += "Max: " + str(np.max(vals)) + "\n"
					info += "Total Pop : " + str(np.sum(vals)) + "\n"
					info += "Total Units: " + str(len(vals)) 
					
					if cfg_exclude_ghssmod0_bool:
						title = "Histogram of sample frame unit population in GHS-SMOD settled areas"
					elif cfg_exclude_pop_per_cell:
						title = "Histogram of sample frame unit population above cut-off of " + str(cfg_exclude_pop_per_cell)
					else:
						title = "Histogram of sample frame unit population" 
					 
					sbu.create_histogram_from_values(vals,histogram_image_above_cut_off,title,info)

				
				# create excluded cells histogram 	
				if cfg_exclude_ghssmod0_bool or cfg_exclude_pop_per_cell:
							
					# open the excluded cells raster
					pop_raster = gdal.Open(pop_raster_excluded_cells)
					pop_raster_band = pop_raster.GetRasterBand(1)
					pop_raster_array = np.array(pop_raster_band.ReadAsArray())
					pop_ndv= pop_raster_band.GetNoDataValue()
			
					# below cut-off
					vals = pop_raster_array[ pop_raster_array != pop_ndv]

					# there may not be any excluded cells
					if len(vals):
						if enable_debug:
							print "Making Histogram:excluded cells"
						m = np.median(vals)
						if not median_value:
							median_value = m
						info = "Median: " + str(m) + "\n"
						info += "Min: " + str(np.min(vals)) + "\n"
						info += "Max: " + str(np.max(vals)) + "\n"
						info += "Total Pop: " + str(np.sum(vals)) + "\n"
						info += "Total Unitss: " + str(len(vals)) 
 
						if cfg_exclude_ghssmod0_bool:
							title = "Histogram of sample frame unit population in GHS-SMOD unsettled areas"
						elif cfg_exclude_pop_per_cell:
							title = "Histogram of sample frame unit population below cut-off of " + str(cfg_exclude_pop_per_cell)
							
						sbu.create_histogram_from_values(vals,histogram_image_below_cut_off,title,info)

				# we need to save the median value in a file to retrive when making the histograms urls
				# if there is an excluded cells value we use that
							
				median_dict = { 1: {'median_pop':median_value} }
				median_info_csv = working_dir + "/median_population.csv"
					
				if os.path.exists(median_info_csv):
					os.remove(median_info_csv)	
				sbu.dictionary_to_csv(median_info_csv,median_dict)

			except Exception as error:
				print "Exception in gridsample code - creating histogram:" + str(error) 
				status=0
				message = str(error)
				
			# we will use a file to indicate whether histograms are currently being made
			# we set it when they are under construction and remove it when created			
			if os.path.exists(pop_histogram_currently_being_made_file):
				os.remove(pop_histogram_currently_being_made_file)
				
			# we always return here whether successful or not as we are only interested in making the histograms
			return [status,message]	
			
		if not status:
			return [status,message]
			
		## create a raster for all the strata 
		## if there is no stratification then create a raster with one strata id for the entire coverage area as there is effectively one strata
		strata_raster = working_dir + "/strata_raster.tif"
		# if there is a shapefile passed in we have to reproject it as an intermediate step
		reprojected_strata_poly_file_name = working_dir + "/reprojected_strata_poly.shp"
		# and  cut it to the coverage area
		cut_strata_poly_file_name = working_dir + "/cut_strata_poly.shp"
		
		strata_poly_file_name = None

	
		if not cfg_stratification_method:	
			if recreate_existing_raster or not os.path.exists(strata_raster):
				# integer pixels
				sbu.make_raster_copy(pop_raster_final,strata_raster,False,1,-99,True)
		elif cfg_stratification_method == "urbanRural":	
			# use ghs  raster
			sbu.make_raster_copy(ghs_mod_raster_cut,strata_raster,False,None,-99,True)
			
			# we need to convert this so that 3 is urban,0,1,2 rural
			# 1 = rural, 3=urban (leave as it is)
			driver = gdal.GetDriverByName('GTiff')
	
			ur_str_r = gdal.Open(strata_raster,gdal.GA_Update)
			ur_str_r_band = ur_str_r.GetRasterBand(1)
			ur_str_r_arr = np.array(ur_str_r_band.ReadAsArray())
				
			ur_str_r_arr[ (ur_str_r_arr == 0)  | (ur_str_r_arr == 1) | (ur_str_r_arr == 2) ] = 1
				
			ur_str_r_band.WriteArray(ur_str_r_arr)
			ur_str_r_band.FlushCache()
			ur_str_r = None
			
		else:
			# either using shapefile made from ghs-admin aras or user has uploaded their own
			# we have checked above that parameters exist
			# check shapefile exists
			if not os.path.exists(strata_poly_filename):
				status = 0
				message += " Could not find strata shapefile at location given"
			else:
				# get list of ids and names
				strata_attributes = sbu.get_shapefile_values(strata_poly_filename,[strata_ID_field,strata_name_field], strata_ID_field)
				if enable_debug:
					print strata_attributes
				
				# check id field and name field exist as attributes
				# this is a dictionary of dictionaries - get a top level key, then get the keys for its sub-dictionary
				keys = strata_attributes[strata_attributes.keys()[0]].keys()

				if not strata_ID_field in keys or not strata_name_field in keys:
					status = 0
					message += " One or both of strata_ID_field and strata_name_field do not exists as attributes in the strata shapefile"				
				else:

					# check that the id field is actually unique by getting all the values
					strata_field_vals = sbu.get_shapefile_values_as_list(strata_poly_filename,strata_ID_field)
									
					if len(strata_field_vals) > len(set(strata_field_vals)):
						status = 0
						message += " The strata ID field chosen does not identify units uniquely"
					else:
						# after reprojecting rasterize shapefile based on id field - to the size of the final cut  pop file
						# if we do this prior to resampling then we should include all areas by using 100m*100m scale - unless there are some very small shapes 
						if not sbu.reprojectShapefileToWGS84(strata_poly_filename,reprojected_strata_poly_file_name):
							
							# we expect there to be a coverage poly but it isn't mandatory
							# if we define strata then this defines hte new boundaries
							if coverage_poly_filename:
								if not sbu.clipShapefilewithShapefileUsingIntersection(reprojected_strata_poly_file_name,cut_strata_poly_file_name,coverage_poly_reprojected,strata_ID_field,strata_name_field):
									strata_poly_file_name = cut_strata_poly_file_name
								else:
									status = 0
									message += " Problem clipping strata poly"	
							else:
								strata_poly_file_name = reprojected_strata_poly_file_name
								
							if strata_poly_file_name:	
								sbu.make_raster_copy(pop_raster_final,strata_raster,False,1,-99,True)
								sbu.rasterize_shapefile_with_gdal(strata_poly_file_name,strata_raster,strata_ID_field,-99,True,default_pixel_size)
				
						else:
							status = 0
							message += " Problem reprojecting strata poly"						
					
		if not status:
			return [status,message]
			
		#strata_info_csv = None
		strata_info_csv = working_dir + "/strata_areas_population_info.csv"
		if  ( recreate_existing_raster or not os.path.exists(strata_info_csv) ) and ( action == "get_strata_pop_values" or action == "create_sample") :
			try:
				
				strata_info_csv = working_dir + "/strata_areas_population_info.csv"
				strata_pop_info_pkl = working_dir + "/strata_areas_population_info.pkl"
				
				if enable_debug:
					print str(strata_poly_filename) + str(strata_ID_field)  + str(strata_name_field)
				
				if os.path.exists(strata_info_csv):
					os.remove(strata_info_csv)
				# first get id,name info from shapefile - this is a list of ids and field names
				shapefile_info = None
				if strata_poly_filename:
					# either ghs admin areas or user defined areas
					shapefile_info = sbu.get_shapefile_values(strata_poly_filename,[strata_ID_field,strata_name_field], strata_ID_field)
				elif cfg_stratification_method == "urbanRural":
					# there are 2 possible strata
					shapefile_info = {}
					shapefile_info[0] = {}
					shapefile_info[1] = {}
				else:
					# there are no strata at all
					shapefile_info = {}
					shapefile_info[1] = {}
				# this returns the strata id and population in that strata
				strata_population_info = None
				
				if recreate_existing_raster or not os.path.exists(strata_pop_info_pkl):
					strata_population_info = sbu.extractPopulationValuesForStrata(pop_raster_final,strata_raster,multicell_psu_raster)
					
					# save the strata pop info - we can possibly then avoid recreating it as it takes ages
					afile = open(strata_pop_info_pkl, 'wb')
					pickle.dump(strata_population_info, afile)
					afile.close()
				else:
					# reload from file
					file2 = open(strata_pop_info_pkl, 'rb')
					strata_population_info = pickle.load(file2)
					file2.close()			
				


				total_population = 0
				total_PSUs = 0
				for key in strata_population_info.keys()  :
					total_population +=  strata_population_info[key][0]
					total_PSUs +=  strata_population_info[key][2]
												
				# now we can add the population info to the shapefile info previously obtained
				# the key is the strata id
				if enable_debug:
					print 'shapefile_info'
					print shapefile_info
				
				for key in shapefile_info.keys()  :
					# if we aren't interested in this area then remove it from our dict
					if not key in strata_population_info.keys():
						shapefile_info.pop(key)
						continue
					
					household_size_weighting = 1
					
					strata_id = key
					if cfg_hh_size and isinstance(cfg_hh_size,dict):
						if cfg_hh_size and str(strata_id) in cfg_hh_size.keys():
							household_size_weighting =	float(cfg_hh_size[str(strata_id)])
						elif cfg_hh_size and int(strata_id) in cfg_hh_size.keys():
							household_size_weighting =	float(cfg_hh_size[int(strata_id)])
					elif cfg_hh_size and isinstance(cfg_hh_size,int):
						household_size_weighting = cfg_hh_size
					elif cfg_hh_size and  isinstance(cfg_hh_size,str):
						household_size_weighting = float(cfg_hh_size)
					elif cfg_hh_size :
						household_size_weighting = float(cfg_hh_size)
		
	
					if key in strata_population_info.keys():
						shapefile_info[key]['population'] = strata_population_info[key][0]
						shapefile_info[key]['population_percent'] = round(float(strata_population_info[key][0])*100/total_population,3)
						shapefile_info[key]['median_pop_per_psu'] = strata_population_info[key][1]
						shapefile_info[key]['str_hh_size'] = household_size_weighting
						shapefile_info[key]['str_households'] = round(float(strata_population_info[key][0])/household_size_weighting,1)
						shapefile_info[key]['households_percent'] = shapefile_info[key]['population_percent']
						shapefile_info[key]['total_PSUs'] = strata_population_info[key][2]
						shapefile_info[key]['total_PSUs_in_coverage_area'] = total_PSUs
						
						# if there is a strata name field
						if strata_name_field in shapefile_info[key].keys():
							shapefile_info[key]['strata_name'] = shapefile_info[key][strata_name_field]
						else:
							shapefile_info[key]['strata_name'] =  ''
							

					else:
						shapefile_info[key]['population'] = '0'
						shapefile_info[key]['population_percent'] = '0'
						shapefile_info[key]['median_pop_per_psu'] = '0'
						shapefile_info[key]['str_hh_size'] = '0'
						shapefile_info[key]['str_households'] = '0'
						shapefile_info[key]['households_percent'] = '0'
						shapefile_info[key]['total_PSUs'] = '0'
						shapefile_info[key]['total_PSUs_in_coverage_area'] = 0
						
						shapefile_info[key]['strata_name'] =  ''
							
				# finally save as CSV
				try:
					sbu.dictionary_to_csv(strata_info_csv,shapefile_info)
				except Exception as error:	
					status=0
					message = "Could not extract strata population info. If using your own coverage file check the correct id and name field selected"
				
				if not os.path.exists(strata_info_csv):
					status=0
					message = "Problem creating strata info file"
				
			except Exception as error:
				print "Exception in gridsample code - get_strata_pop_values:" + str(error) 
				status=0
				message = str(error)

			if  action == "get_strata_pop_values":
				return [status,message]
		elif action == "get_strata_pop_values" or action == "create_sample":
			strata_info_csv = working_dir + "/strata_areas_population_info.csv"
		else:
			strata_info_csv = None
			
		if not status:
			return [status,message]	
	

		## now we create an index raster the same size as the final pop raster numbered 1,2,3 ... across rows
		## so each  pixel has a different value. 
		## e.g.
		## 1   2  3  4
		## 5   6 ND ND
		## 7   8  9 ND
		## 10 11 12 13
		index_raster = working_dir + "/index_raster.tif"
	
		if recreate_existing_raster or not os.path.exists(index_raster ):
			sbu.indexrast_create(pop_raster_final,index_raster )		
						
		## create a raster to contain the results - this will be the chosen pixels
		## it should be the same size as the pop_raster_final  - the pop raster cut to the coverage area
		## originally no data values will be -99 and all other cells set to 0
		results_raster = working_dir + "/results_raster.tif"
		sbu.make_raster_copy(pop_raster_final,results_raster,False,0,-99)
		
		# check if we need to make an oversized coverage grid
		spatial_raster = None
		
		# create spatial raster
		# we need to make an oversize grid raster to ensure that at least one of these grid members is filled in
		# it is based on scaling up the pop _raster
		# we will create a copy of pop raster but will fill as follows assuming scaling factor of 2
		# 1 1 2 2 3 3 4
		# 1 1 2 2 3 3 4
		# 5 5 6 6 7 7 8
		# 5 5 6 6 7 7 8		
		if not cfg_oversample_grid_spatial_scale:
			# user doesn't want spatial overgrid
			pass
		elif cfg_oversample_grid_spatial_scale < 0:
				status = 0
				message += " cfg_oversample_grid_spatial_scale either not given or is negative\n"	
		else: 
			# create the spatial raster 
			spatial_raster = working_dir + "/spatial_oversample_raster.tif"
				
			if recreate_existing_raster or not os.path.exists(spatial_raster ):
				sbu.make_raster_copy(pop_raster_final,spatial_raster,False,1,-99,True)
				# this was here to include all areas without exclusions - but that isn't what is required
				#sbu.make_raster_copy(pop_raster_without_exclusions,spatial_raster,False,1,-99,True)
								
				spatial_raster_scale = int(cfg_oversample_grid_spatial_scale)
				
				sr = gdal.Open(spatial_raster,gdal.GA_Update)
				sr_band = sr.GetRasterBand(1)
				spatial_raster_arr = np.array(sr_band.ReadAsArray())
				spatial_raster_arr_orig = np.array(sr_band.ReadAsArray())
			
				spatial_raster_cols = sr.RasterXSize
				spatial_raster_rows = sr.RasterYSize
					
				spatial_index = 1
				for row in range(0,spatial_raster_rows,spatial_raster_scale):
					for col in range(0,spatial_raster_cols,spatial_raster_scale):
						spatial_raster_arr[row:row + spatial_raster_scale,col:col + spatial_raster_scale] = spatial_index
						spatial_index += 1
					
				# reset the no-data	cells		
				spatial_raster_arr[spatial_raster_arr_orig == -99] = -99		
				sr_band.WriteArray(spatial_raster_arr)
				sr_band.FlushCache()

				## close everything to save resources and make sure actually written to disc before next stage
				spatial_raster_arr = None
				spatial_raster_arr_orig = None
				sr = None					

		if not status:
			return [status,message]			

		
		## work out PSU_allocation to strata
		if action == "create_sample":
			# Now create raster of selected PSUs 
			results_shapefile = working_dir + "/" + PSU_filename
			
			PSU_id_strata_id_info = None
			
			if multicell_psu_raster:
				strata_file_1 =  None
				strata_id_field_1  = None
				if strata_poly_file_name:
					strata_file_1 = strata_poly_file_name  
					strata_id_field_1  = strata_ID_field
		
				PSU_id_strata_id_info = selectPSUs_MultiCellPSUs(results_raster,pop_raster_final,ghs_mod_raster_cut,strata_raster,index_raster,cfg_psu_per_strata,multicell_psu_raster,results_shapefile,strata_file_1,strata_id_field_1,cfg_hh_size,spatial_raster)
			else:
				strata_file_1 =  None
				strata_id_field_1  = None
				if strata_poly_file_name:
					strata_file_1 = strata_poly_file_name  
					strata_id_field_1  = strata_ID_field
				PSU_id_strata_id_info = selectPSUs(results_raster,pop_raster_final,ghs_mod_raster_cut,strata_raster,index_raster,cfg_psu_per_strata,results_shapefile,strata_file_1,strata_id_field_1,cfg_hh_size,spatial_raster)
						
			# now add extra attributes to the the results shapefile 
			if multicell_psu_raster:
				# The results shapefile is created within selectPSUs_MultiCellPSUs
				sbu.add_strata_id_to_results_shapefile (results_shapefile,PSU_id_strata_id_info)
				sbu.add_strata_attributes_to_results_shapefile (results_shapefile,strata_info_csv)

			else:
				# The results shapefile is  created within selectPSUs
				sbu.add_strata_id_to_results_shapefile (results_shapefile,PSU_id_strata_id_info)
				sbu.add_strata_attributes_to_results_shapefile (results_shapefile,strata_info_csv)
				
		message = " gridsample completed"
	
		# finally delete the intermediate rasters
		# we may need some of these for the html summary page
		if delete_intermediate_outputs:
			delete_intermediate_files(working_dir)
			
	except Exception as error:
		print "Exception in gridsample code:" + str(error) 
		status=0
		message = str(error)

	
	if os.path.exists(pop_histogram_currently_being_made_file):
		os.remove(pop_histogram_currently_being_made_file)
				
	if enable_debug:
		print "## Leaving function:gridsample"
			
	return [status,message]

		
def getHouseholdSizeWeighting(cfg_hh_size,strata_id):
	if enable_debug:
		print "## Entering function:getHouseholdSizeWeighting"
	
	household_size_weighting = 1
	if cfg_hh_size and isinstance(cfg_hh_size,dict):
		if cfg_hh_size and str(strata_id) in cfg_hh_size.keys():
			household_size_weighting =	int(cfg_hh_size[str(strata_id)])
	elif cfg_hh_size and isinstance(cfg_hh_size,int):
		household_size_weighting = cfg_hh_size
	elif cfg_hh_size and  isinstance(cfg_hh_size,str):
		household_size_weighting = int(cfg_hh_size)
		
	return household_size_weighting
	
# This function will create a new raster and a shapefile indicating selected PSUs 
def selectPSUs(results_raster,pop_raster_final,ghs_mod_raster,strata_raster,index_raster,cfg_psu_per_strata,results_shapefile,strata_file = None,strata_id_field = None,cfg_hh_size = None,spatial_raster = None,cfg_random_number = None  ):
	# this gets a specific number of PSU per strata
	if enable_debug:
		print "## Entering function:selectPSUs"
			
		print pop_raster_final
		print ghs_mod_raster
		print strata_raster	
	
		print "results_shapefile" + results_shapefile

	### we may prefer to pass through references to already open raster but for the moment open them here
	## and get values as arrays
	driver = gdal.GetDriverByName('GTiff')
	
	rr = gdal.Open(results_raster,gdal.GA_Update)
	results_raster_band = rr.GetRasterBand(1)
	results_raster_ndv = results_raster_band.GetNoDataValue()
	results_raster_arr = np.array(results_raster_band.ReadAsArray())
	
	pr = gdal.Open(pop_raster_final)
	pop_raster_final_arr = np.array(pr.GetRasterBand(1).ReadAsArray())
	pr_ndv = pr.GetRasterBand(1).GetNoDataValue()


	ghs = gdal.Open(ghs_mod_raster)
	ghs_raster_arr = np.array(ghs.GetRasterBand(1).ReadAsArray())
	ghs_ndv = ghs.GetRasterBand(1).GetNoDataValue()

	sr = gdal.Open(strata_raster)
	s_band = sr.GetRasterBand(1)
	s_ndv = s_band.GetNoDataValue()
	strata_raster_arr = np.array(s_band.ReadAsArray())
	
	ir = gdal.Open(index_raster,gdal.GA_Update)
	ir_band = ir.GetRasterBand(1)
	index_raster_arr = np.array(ir_band.ReadAsArray())	
	index_raster_arr_orig = np.array(ir_band.ReadAsArray())
	ir_ndv = ir_band.GetNoDataValue()
	
	## get values we are going to loop through - exclude no data value	
	strata_ids = np.unique(strata_raster_arr[strata_raster_arr != s_ndv])
	
	if enable_debug:
		print "strata_ids"
	
	# we are going to try creating the output shapefiles here rather than externally as it may be easier
	# the idea is to create one per strata then merge with the main results shapefile
	working_dir = os.path.dirname(results_shapefile)
	if os.path.exists(results_shapefile):
		os.remove(results_shapefile)
	
	## there should be 0,1,2,3 - but it's possible there aren't any of certain values within specific areas
	## we may have manipulated the ghs and pop rasters before this stage to exclude 0 values by setting them to ndv
	ghs_values = np.unique(ghs_raster_arr[ghs_raster_arr!=ghs_ndv])
		
	# this will keep track of how much population has been allocated
	total_pop_allocated = 0
	# the number of cells/PSUs allocated
	cells_allocated = 0
	
	# we need to save quite a lot of info about number of cells per stratum,population per cluster etc
	# vocab has changed somewhat during development but PSU,cell and cluster all mean the same thing
	# a dictionary of dictionaries holding the id of each PSU selected and other info such as its strata_id
	selected_PSU_info = {}
	selected_strata_info = {}
	
	if cfg_random_number:
		random.seed(cfg_random_number)
	
	## we will go through each strata then each GHS value
	for strata_id in strata_ids:
		if enable_debug:
			print("strata_id: " + str(strata_id))
		
		# make a new clean copy for this new strata
		strata_raster_arr_copy = np.array(s_band.ReadAsArray())
		index_raster_arr_copy = np.array(ir_band.ReadAsArray())
				
		
		###################### pixel edge recalculations 
		## as boundary lines cross pixels the population value of each bounadry pixel is re-calculated
		# allowing for what percentage falls inside the boundary
		# if we are doing the edge pixel recalculations then make a new clean copy for this new strata 
		if 1 and strata_file and strata_id_field:
			if enable_debug:
				print "### About to start pixel edge recalculations"
				
			# firstly create a small buffer polygon around the border of the strata of interest
			buffered_boundary_shapefile = working_dir + "/buffered_boundary_shapefile_" + str(strata_id)  + ".shp"
			
			sbu.make_boundary_shapefile(strata_file,buffered_boundary_shapefile,buffer_distance = 0.000001,feature_name=strata_id_field,feature_value=strata_id)
			
			# now create a raster for the shapefile - this is essentially ths intersection of the shapefile with the index raster
			# so we have the indexes of the boundary cells
			boundary_raster =  working_dir + "/buffered_boundary_raster_" + str(strata_id)  + ".tif"
			sbu.rasterize_shapefile(buffered_boundary_shapefile,index_raster,boundary_raster,0)
			
			# now convert pixels to polygons so we can calculate intersections between rasters
			# each polygonised pixel has the value of the index raster cell
			polygonized_boundary_pixels = working_dir + "/boundary_pixels_" + str(strata_id)  + ".shp"
			sbu.polygonise_raster (boundary_raster,polygonized_boundary_pixels)
			
			pop_raster_final_arr = sbu.calculate_border_pixels_area_inside_strata(polygonized_boundary_pixels,strata_file,pop_raster_final,index_raster)
			
			# remember to delete the intermediate files
			if delete_intermediate_outputs:
				for f in glob.glob(working_dir + "/buffered_boundary_shapefile_*"):
					os.remove(f)			
					
				for f in glob.glob(working_dir + "/boundary_pixels_*"):
					os.remove(f)	

			# always remove this otherwise we make loads and run out of space
			if os.path.exists(boundary_raster):
				os.remove(boundary_raster)
					
			if enable_debug:			
				print "pixel edge calculations done"
		#####################################################################
	
		## Note for future work maybe scale everything to make integers - it might speed things up
		# allowing for ndv in pop raster - e.g when stratifying by ur
		pop_raster_for_this_strata = pop_raster_final_arr[ (strata_raster_arr == strata_id) & (pop_raster_final_arr != pr_ndv)  ]
		
		strata_pop = sum(pop_raster_for_this_strata)
		
		if enable_debug:
			print("*** strata_pop: " + str(strata_pop))
			print "cfg_psu_per_strata" + str(cfg_psu_per_strata)
		
		## calculate the number of cells (units we can select from) 
		
		# this could happen when user excludes by say GHS mod and a strata is full of NDV
		# if no strata then cfg_psu_per_strata can just be an int
		# key could be a string or an int so allow for either
		if isinstance(cfg_psu_per_strata,dict) and not str(strata_id) in cfg_psu_per_strata.keys() and not strata_id in cfg_psu_per_strata.keys():
			if enable_debug:
				print '#####' + str(strata_id)
				print str(cfg_psu_per_strata.keys())
			continue
			
		# we calculate a weighting to convert population to households
		# if not given or no values exists for this id set weighting to 1 and we are effectively working on population basis
		#household_size_weighting = 1
		household_size_weighting = getHouseholdSizeWeighting(cfg_hh_size,strata_id)
		
		psu_per_strata = None
		if isinstance(cfg_psu_per_strata,dict):		
			if isinstance(cfg_psu_per_strata.keys()[0], int):
				psu_per_strata = int(cfg_psu_per_strata[int(strata_id)])
			else:
				psu_per_strata = int(cfg_psu_per_strata[str(strata_id)])
		else: 
			psu_per_strata = cfg_psu_per_strata
		
		if enable_debug:		
			print '###############' + str(strata_id) + '---' + str(psu_per_strata)
		
		## sampling point works out the cumulative poulation increment to get the required number of psus in this strata
		sampling_interval = strata_pop/(psu_per_strata*household_size_weighting)	
		
		if enable_debug:
			print("sampling_interval: " + str(sampling_interval))
			
		starting_point = random.random()*sampling_interval	
		
		if enable_debug:
			print("starting_point: " + str(starting_point))
		
		# because the next section works out the cumulative total population/households for a separate subset we need to add the previous total to it  
		overall_cum_total = 0
		# similarly we need to update the interval start point
		number_of_sampling_intervals_start = 0 
		
		
		for ghs_value in ghs_values:
			if enable_debug:
				print("ghs_value" + str(ghs_value))

			# get the population values for this strata and ghs value 
			# now allowing for ndv in pop raster - e.g when stratifying by ur
			arr_locations_required = (strata_raster_arr == strata_id)  & (ghs_raster_arr == ghs_value) & (pop_raster_final_arr != pr_ndv)

			pop_raster_for_this_strata_and_ghs_value = pop_raster_final_arr[ arr_locations_required ]
			
			# change the pop raster to allow for the hh size weighting 
			pop_raster_for_this_strata_and_ghs_value = pop_raster_for_this_strata_and_ghs_value/household_size_weighting

			
			if not len(pop_raster_for_this_strata_and_ghs_value):
				if enable_debug:
					print "No populated cells were found for this strata and ghs value " + str(strata_id) + " " + str(ghs_value)
				continue
			
			# we need to get the proportional number of psus to select from this ghs value
			# use the cumulative values to select 
			cum_values_arr = np.cumsum(pop_raster_for_this_strata_and_ghs_value) + overall_cum_total
			
			if enable_debug:
				print(cum_values_arr)
			
			# update the overall cumulative total - the last value in the cum values arr is the total
			overall_cum_total = cum_values_arr[cum_values_arr.size-1]
			
			# create an array the same size as the selected pop array and intialize to zero
			# as a location is chosen it will be set to 1
			selected_locations_for_this_strata_and_ghs_value = np.zeros(cum_values_arr.size)
					
			for j in range(number_of_sampling_intervals_start,psu_per_strata):
				
				# this value is the cumulative population value that we want to find in array of cumulative population
				value_to_find_interval_for = starting_point +  j*sampling_interval
				
				if enable_debug:
					print("j:" + str(j) + "\n")
					print("number_of_sampling_intervals_start:" + str(number_of_sampling_intervals_start) + "\n")
					print("psu_per_strata:" + str(psu_per_strata) + "\n")
				
				indx = np.searchsorted(cum_values_arr, value_to_find_interval_for)
				
				## don't go over the end of the array
				if(indx < cum_values_arr.size):
					cells_allocated += 1
					total_pop_allocated += pop_raster_for_this_strata_and_ghs_value[indx]
					selected_locations_for_this_strata_and_ghs_value[indx] = 1
					
					# save the number of sampled clusters in this stratum
					if strata_id in selected_strata_info.keys():
						selected_strata_info[strata_id]['s_cl_st'] = selected_strata_info[strata_id]['s_cl_st'] + 1
					else:
						# this is the first PSU selected 
						row = {'s_cl_st':1}
						selected_strata_info[strata_id] = row 
					
				else:
					#  update the start point for the next time through 
					number_of_sampling_intervals_start = j

					break
			
			if enable_debug:			
				print("# cells_allocated:" + str(cells_allocated))	
				print("# total_pop_allocated:" + str(total_pop_allocated))
			
			# update the results raster with locations selected for this strata and ghs level
			results_raster_arr[ arr_locations_required ] = selected_locations_for_this_strata_and_ghs_value
			
			if strata_id in selected_strata_info.keys() and selected_strata_info[strata_id]['s_cl_st'] >= psu_per_strata:
				break
	
		if enable_debug:
			print "About to make shapefile"
		
		# make the results shapefile strata by strata
		# and combine them into one final shapefile of selected PSUs
		if not os.path.exists(results_shapefile):
			
			# this is the first chosen PSU

			index_raster_arr_copy[strata_raster_arr != strata_id] = ir_ndv
			index_raster_arr_copy[ results_raster_arr != 1] = ir_ndv	

			# only make a shapefile if there is something to put in it
			if len(index_raster_arr_copy[ index_raster_arr_copy != ir_ndv]):
				# create a shapefile from the results
				ir_band.WriteArray(index_raster_arr_copy)		
				sbu.polygonise_open_raster(ir,results_shapefile,field_name="cl_id")
				
				# reset raster band
				ir_band.WriteArray(index_raster_arr)				
				
				# because we have selected raster cells some go over the strata boundaries so cut to the strata boundary 
				if strata_file and strata_id_field:
					where = str(strata_id_field) + '=' + str(strata_id)
					sbu.clipShapefilewithShapefilebyAttribute(results_shapefile,results_shapefile,strata_file,where)			
		else:
			# make a shapefile for this strata and merge with main results shapefile
			temp_shapefile_for_this_strata = working_dir + "/" + "temp_results_for_strata_" + str(strata_id) +".shp"
					
			index_raster_arr_copy[strata_raster_arr != strata_id] = ir_ndv
			index_raster_arr_copy[ results_raster_arr != 1] = ir_ndv	

			# only make a shapefile if there is something to put in it
			if len(index_raster_arr_copy[ index_raster_arr_copy != ir_ndv]):
				# create a shapefile from the results
				ir_band.WriteArray(index_raster_arr_copy)		
				sbu.polygonise_open_raster(ir,temp_shapefile_for_this_strata,field_name="cl_id")
				
				# reset raster band
				ir_band.WriteArray(index_raster_arr)				
				
				# because we have selected raster cells some go over the strata boundaries so cut to the strata boundary 
				if strata_file and strata_id_field:
					where = str(strata_id_field) + '=' + str(strata_id)
					sbu.clipShapefilewithShapefilebyAttribute(temp_shapefile_for_this_strata,temp_shapefile_for_this_strata,strata_file,where)			
				
				# now merge with results shapefile - we are building up the results shapefile bit by bit
				sbu.mergeShapefiles(results_shapefile,temp_shapefile_for_this_strata)									

	if enable_debug:		
		print("@ cells_allocated:" + str(cells_allocated))	
		print("@ total_pop_allocated:" + str(total_pop_allocated))
	
	# we need to record the strata id,and selection section for each id
	# we add to a dictionary which will be used afterwards to add the strata id to the shapefile to each cell
	selected_cell_ids = index_raster_arr[results_raster_arr == 1]
	selected_strata_ids = strata_raster_arr[results_raster_arr == 1]
	selected_pop = pop_raster_final_arr[results_raster_arr == 1]
	
	for i in range(len(selected_cell_ids)):

		strata_id = selected_strata_ids[i]
		household_size_weighting = getHouseholdSizeWeighting(cfg_hh_size,strata_id)
		pop_in_cluster= selected_pop[i]
		households_in_cluster = round(pop_in_cluster/household_size_weighting,2)
		
		# to indicate that this cell was selected via the main algorithm
		row = {'strata_id':strata_id,'cl_type':'main','s_cl_pop':pop_in_cluster,'s_cl_hh':households_in_cluster}
		selected_PSU_info[selected_cell_ids[i]] = row 

	
	## the spatial raster is a larger scale grid we want to make sure that at least one cell in each larger grid is selected
	# we pick a random empty cell for each one that isn't selected
	# this can take a long time if the spatial scale is small
	# we don't worry about the total allocated
	if spatial_raster:	
		if enable_debug:
			print "doing spatial oversample - this can be time consuming"
		
		sp_r = gdal.Open(spatial_raster)
		sp_r_band = sp_r.GetRasterBand(1)
		sp_r_ndv = sp_r_band.GetNoDataValue()
		sp_r_arr = np.array(sp_r_band.ReadAsArray())
		
		spatial_ids =  np.unique(sp_r_arr[sp_r_arr!= -99])
		
		# make new clean copies for the spatial raster
		strata_raster_arr_copy = np.array(s_band.ReadAsArray())
		results_raster_arr_copy = np.array(results_raster_band.ReadAsArray())
		
		index_raster_arr_copy = np.array(ir_band.ReadAsArray())
		
		results_raster_locations_with_values = (results_raster_arr !=results_raster_ndv)
		
		for sp_id in spatial_ids:
		
			if enable_debug:
				print "spatial_ids:" + str(sp_id)
			
			# the locations that have this spatial id that aren't ndv
			# do just once here to speed things up
			locations_of_interest = (sp_r_arr == sp_id) & results_raster_locations_with_values
			
			
			# get results raster values for this spatial id
			rr_vals = results_raster_arr[ locations_of_interest ]	
						
			# check if any have been chosen - at least one cell will have the value 1 - otherwise all the values will be 0	
			if len(rr_vals) and np.sum(rr_vals) == 0:
				if enable_debug:
					print ' @@@@@@@@@@@@ Choosing a spatial oversample cell randomly'
				
				# get the oversample area cells and the results raster cells for this oversample id
				oversample_area_cells_of_interest = ( sp_r_arr == sp_id ) & (strata_raster_arr != s_ndv) & (pop_raster_final_arr != pr_ndv)
				
				this_oversample_area = sp_r_arr[oversample_area_cells_of_interest]
				rr_vals = results_raster_arr[ oversample_area_cells_of_interest]	
				strata_id_vals = strata_raster_arr	[ oversample_area_cells_of_interest ]	

				if len(this_oversample_area):				
		
					index_chosen = 0					
					if len(this_oversample_area) > 1:					
						# pick a random empty cell
						index_chosen = random.randint(0,len(this_oversample_area)-1)
				
					# set it to chosen
					rr_vals[index_chosen] = 1
					
					# update the results raster - this also includes the locations previously chosen
					results_raster_arr[ oversample_area_cells_of_interest ] = rr_vals
					# this is just the locations chosen via this section of code
					
					results_raster_arr_copy[ oversample_area_cells_of_interest ] = rr_vals	

					cells_allocated += 1	
					
					strata_id = strata_id_vals[index_chosen]
					# save the number of sampled clusters in this stratum
					if strata_id in selected_strata_info.keys():
						selected_strata_info[strata_id]['s_cl_st'] = selected_strata_info[strata_id]['s_cl_st'] + 1
					else:
						# this is the first PSU selected 
						row = {'s_cl_st':1}
						selected_strata_info[strata_id] = row 				

		# if any new locations have been selected update the results shapefile
		if np.sum(results_raster_arr_copy[results_raster_arr_copy!=results_raster_ndv]) > 0:
			
			# make a shapefile for this strata and merge with main results shapefile
			temp_shapefile_for_this_strata = working_dir + "/" + "temp_results_for_strata_spatial_oversample.shp"
					
			# set values that haven't been selected to ndv
			index_raster_arr_copy[ results_raster_arr_copy != 1] = ir_ndv
			
			ir_band.WriteArray(index_raster_arr_copy)		
			sbu.polygonise_open_raster(ir,temp_shapefile_for_this_strata,field_name="cl_id")
				
			# reset raster band
			ir_band.WriteArray(index_raster_arr)	
						
			# now merge with results shapefile - we are building up the results shapefile bit by bit
			sbu.mergeShapefiles(results_shapefile,temp_shapefile_for_this_strata)
			
	
	# finally overwrite the original results raster
	results_raster_band.WriteArray(results_raster_arr)
	results_raster_band.FlushCache()
	
	# remove unwanted files
	if delete_intermediate_outputs:
		for f in glob.glob(working_dir + "/" + "temp_results_for_strata_*"):
			os.remove(f)	
		
	# reset - it seems to be gettng unset somewhere
	ir_band.WriteArray(index_raster_arr_orig)
	ir_band.FlushCache()
	 	
	
	rr = None
	ir = None
	
	# add oversample info to dictionary
	# we need to get the strata id matching each index
	# we create a dictionary which will be used afterwards to add the strata id to the shapefile to each cell			
	selected_cell_ids = index_raster_arr[results_raster_arr == 1]  
	selected_strata_ids = strata_raster_arr[results_raster_arr == 1]
	selected_pop = pop_raster_final_arr[results_raster_arr == 1]
	
	for i in range(len(selected_cell_ids)):
		# if id not in list then it was selected as part of the oversample algorithm
		if not selected_cell_ids[i] in selected_PSU_info.keys():
			# to indicate that this cell was selected via the oversample algorithm
			strata_id = selected_strata_ids[i]
			household_size_weighting = getHouseholdSizeWeighting(cfg_hh_size,strata_id)
			pop_in_cluster= selected_pop[i]
			households_in_cluster = round(pop_in_cluster/household_size_weighting,2)
			
			row = {'strata_id':strata_id,'cl_type':'oversample','s_cl_pop':pop_in_cluster,'s_cl_hh':households_in_cluster}

			selected_PSU_info[selected_cell_ids[i]] = row 
		
		# add number of sampled clusters in stratum info 
		# and total number of selected clusters	
		strata_id = selected_PSU_info[selected_cell_ids[i]]['strata_id']
		# this should always be the case but check anyway
		if strata_id in selected_strata_info.keys():
			selected_PSU_info[selected_cell_ids[i]]['s_cl_st'] = selected_strata_info[strata_id]['s_cl_st']
			selected_PSU_info[selected_cell_ids[i]]['s_cl_tot'] = cells_allocated
		else:
			selected_PSU_info[selected_cell_ids[i]]['s_cl_st'] = 0
			selected_PSU_info[selected_cell_ids[i]]['s_cl_tot'] = cells_allocated
		
	return selected_PSU_info


	
# This function will create a new raster and shapefile made of selected PSU seeds 
# The complication here is that a PSU is not a raster cell, rather it is a selection of cells as given by an input raster multicell_psu_raster
def selectPSUs_MultiCellPSUs(results_raster,pop_raster_final,ghs_mod_raster,strata_raster,index_raster,cfg_psu_per_strata,multicell_psu_raster,results_shapefile,strata_file = None,strata_id_field = None,cfg_hh_size = None,spatial_raster = None,cfg_random_number = None  ):
	# this gets a specific number of PSU per strata
	if enable_debug:
		print "#### Entering function:selectPSUs_MultiCellPSUs"
				

	### we may prefer to pass through references to already open raster but for the moment open them here
	## and get values as arrays
	driver = gdal.GetDriverByName('GTiff')
	
	rr = gdal.Open(results_raster,gdal.GA_Update)
	results_raster_band = rr.GetRasterBand(1)
	results_raster_ndv = results_raster_band.GetNoDataValue()
	results_raster_arr = np.array(results_raster_band.ReadAsArray())
	
	pr = gdal.Open(pop_raster_final)
	pop_raster_final_arr = np.array(pr.GetRasterBand(1).ReadAsArray())
	pr_ndv = pr.GetRasterBand(1).GetNoDataValue()
	
	ghs = gdal.Open(ghs_mod_raster)
	ghs_raster_arr = np.array(ghs.GetRasterBand(1).ReadAsArray())
	ghs_ndv = ghs.GetRasterBand(1).GetNoDataValue()
	
	sr = gdal.Open(strata_raster)
	s_band = sr.GetRasterBand(1)
	strata_raster_arr = np.array(s_band.ReadAsArray())

	ir = gdal.Open(index_raster,gdal.GA_Update)
	ir_band = ir.GetRasterBand(1)
	index_raster_arr = np.array(ir_band.ReadAsArray())	
	
	ir_ndv = ir_band.GetNoDataValue()
	
	mcr = gdal.Open(multicell_psu_raster)
	mcr_band = mcr.GetRasterBand(1)
	# keep an original copy
	multicell_psu_raster_arr_orig = np.array(mcr_band.ReadAsArray())
	mcr_ndv = mcr_band.GetNoDataValue()
	
	## get values we are going to loop through - exclude no data value
	s_ndv = s_band.GetNoDataValue()
	strata_ids = np.unique(strata_raster_arr[strata_raster_arr != s_ndv])
	
	# we are going to try creating the output shapefiles here rather than externally as it may be easier
	# the idea is to create one per strata then merge with the main results shapefile
	working_dir = os.path.dirname(results_shapefile)
	if os.path.exists(results_shapefile):
	#	os.remove(results_shapefile)
		for f in glob.glob( results_shapefile.replace('.shp','*')  ):
			os.remove(f)
	
	## there should be 0,1,2,3 - but it's possible there aren't any of certain values within specific areas
	## we may have manipulated the ghs and pop rasters before this stage to exclude 0 values by setting them to ndv
	ghs_values = np.unique(ghs_raster_arr[ghs_raster_arr!=ghs_ndv])
	
	# this is an array holding the total value of each multi-cell in just one cell per multi-cell id
	pop_raster_grouped_values = np.array(pr.GetRasterBand(1).ReadAsArray())
	pop_raster_grouped_values.fill(pr_ndv)
	

	
	# we need to save quite a lot of info about number of cells per stratum,population per cluster etc
	# vocab has changed somewhat during development but PSU,cell and cluster all mean the same thing
	# a dictionary of dictionies holding hte id of each PSU selected and other info such as its strata_id
	selected_PSU_info = {}
	selected_strata_info = {}
	
	# this will keep track of how much population has been allocated
	total_pop_allocated = 0
	
	cells_allocated = 0
	
	if cfg_random_number:
		random.seed(cfg_random_number)
	
	## we will go through each strata then each GHS value
	for strata_id in strata_ids:
		if enable_debug:
			print("strata_id: " + str(strata_id))
			
		# we need to be careful because var_2_arr = var_1_arr is just a pointer to car_1_arr not a copy
		# make a new clean copy for this new strata
		multicell_psu_raster_arr = np.array(mcr_band.ReadAsArray())
		
		# make a new clean copy for this new strata
		strata_raster_arr_copy = np.array(s_band.ReadAsArray())
	
		pop_raster_for_this_strata = pop_raster_final_arr[ (strata_raster_arr == strata_id) & (pop_raster_final_arr != pr_ndv) ]
		
		strata_pop = sum(pop_raster_for_this_strata)
		
		if enable_debug:
			print("strata_pop: " + str(strata_pop))
		
		
		## calculate the number of cells (units we can select from) 

		# this could happen when user excludes by say GHS mod and a strata is full of NDV
		# if no strata then cfg_psu_per_strata can just be an int
		if isinstance(cfg_psu_per_strata,dict) and not str(strata_id) in cfg_psu_per_strata.keys():
			continue
			
		# we calculate a weighting to convert population to households
		# if not given or no values exists for this id set weighting to 1 and we are effectively working on population basis
		household_size_weighting = 1
		if cfg_hh_size and isinstance(cfg_hh_size,dict):
			if cfg_hh_size and str(strata_id) in cfg_hh_size.keys():
				household_size_weighting =	int(cfg_hh_size[str(strata_id)])
		elif cfg_hh_size and isinstance(cfg_hh_size,int):
			household_size_weighting = cfg_hh_size
		elif cfg_hh_size and  isinstance(cfg_hh_size,str):
			household_size_weighting = int(cfg_hh_size)
		
		if enable_debug:
			print "household_size_weighting:" + str(household_size_weighting) 
				
		psu_per_strata = None
		if isinstance(cfg_psu_per_strata,dict):		
			psu_per_strata = int(cfg_psu_per_strata[str(strata_id)])
		else: 
			psu_per_strata = cfg_psu_per_strata
			
		if enable_debug:
			print '###############' + str(strata_id) + '---' + str(psu_per_strata)
		
		## sampling point works out the cumulative poplation increment to get the required number of psus in this strata
		sampling_interval = strata_pop/(psu_per_strata*household_size_weighting)

		if enable_debug:		
			print("sampling_interval: " + str(sampling_interval))
			
		starting_point = random.random()*sampling_interval	
		
		if enable_debug:
			print("starting_point: " + str(starting_point))
		
		# because the next section works out the cumulative total population/households for a separate subset we need to add the previous total to it  
		overall_cum_total = 0
		# similarly we need to update the interval start point
		number_of_sampling_intervals_start = 0 
		
		# Have removed the GHS bit as it means we are selecting parts of PSUs and probably the whole PSU is required
		# for ghs_value in ghs_values:
		if 1:
			# now allowing for ndv in pop raster - e.g when stratifying by ur
			arr_locations_required = (strata_raster_arr == strata_id)   & (pop_raster_final_arr != pr_ndv)			
			
			pop_raster_for_this_strata_and_ghs_value = pop_raster_final_arr[arr_locations_required]
			
			multicell_raster_for_this_strata_and_ghs_value = multicell_psu_raster_arr[arr_locations_required]
					
			## convert the pop raster so that the first cell of a multi-cell PSU has the total value in it and the remaining cells for that PSU have the value 0
			# we are sort of converting the multicell raster to a single cell raster
			jj=0
			for multicell_id in np.unique(multicell_raster_for_this_strata_and_ghs_value):
				jj+=1
				multicell_locations_of_interest = (multicell_raster_for_this_strata_and_ghs_value==multicell_id)
				
				temp_arr= pop_raster_for_this_strata_and_ghs_value[multicell_locations_of_interest]
				temp_arr[0]=np.sum(temp_arr)
				
				temp_arr[1:]=0

				pop_raster_for_this_strata_and_ghs_value[multicell_locations_of_interest] = temp_arr
				
			# this should hold them all
			pop_raster_grouped_values[arr_locations_required] = pop_raster_for_this_strata_and_ghs_value
			
			# change the pop raster to allow for the hh size weighting 
			pop_raster_for_this_strata_and_ghs_value = pop_raster_for_this_strata_and_ghs_value/household_size_weighting
			
			if not len(pop_raster_for_this_strata_and_ghs_value):
				if enable_debug:
					print "No populated cells were found for this strata and ghs value " + str(strata_id) 
				continue

			
			# we need to get the proportional number of psus to select from this ghs value
			# use the cumulative values to select 
			cum_values_arr = np.cumsum(pop_raster_for_this_strata_and_ghs_value) + overall_cum_total
			
			# update the overall cumulative total - the last value in the cum values arr is the total
			overall_cum_total = cum_values_arr[cum_values_arr.size-1]
			
			# create an array the same size as the selected pop array and intialize to zero
			# as a location is chosen it will be set to 1
			selected_locations_for_this_strata_and_ghs_value = np.zeros(cum_values_arr.size)
			
			for j in range(number_of_sampling_intervals_start,psu_per_strata):
				# this value is the cumulative population value that we want to find in array of cumulative population
				value_to_find_interval_for = starting_point +  j*sampling_interval

				indx =  np.searchsorted(cum_values_arr, value_to_find_interval_for) 
					
				## don't go over the end of the array
				if( indx < cum_values_arr.size):
					
					cells_allocated += 1
					total_pop_allocated += pop_raster_for_this_strata_and_ghs_value[indx]

					# set all the locations with the index value for this PSU to 1
					selected_locations_for_this_strata_and_ghs_value[multicell_raster_for_this_strata_and_ghs_value == multicell_raster_for_this_strata_and_ghs_value[indx]] = 1
									
					# save the number of sampled clusters in this stratum
					if strata_id in selected_strata_info.keys():
						selected_strata_info[strata_id]['s_cl_st'] = selected_strata_info[strata_id]['s_cl_st'] + 1
					else:
						# this is the first PSU selected 
						row = {'s_cl_st':1}
						selected_strata_info[strata_id] = row 
						
				else:
					#  update the start point for the next time through 
					number_of_sampling_intervals_start = j

					break
			
			if enable_debug:			
				print("# cells_allocated:" + str(cells_allocated))	
				print("# total_pop_allocated:" + str(total_pop_allocated))
			
			# update the results raster with locations selected for this strata and ghs level
			results_raster_arr[ arr_locations_required ] = selected_locations_for_this_strata_and_ghs_value
			
			# this is only needed if the ghs loop is in play
			#if strata_id in selected_strata_info.keys() and selected_strata_info[strata_id]['s_cl_st'] >= psu_per_strata:
			#	break
		# make the results shapefile strata by strata
		# and combine them into one final shapefile of selected PSUs
		
		if 1:
			# this is giving the shapes the multicell id not the strata id
			# keep for the time being while we make sure what we really want - it appears we do want it after all!
			if not os.path.exists(results_shapefile):
				# this is the first chosen PSU
				
				# set values that haven't been selected to ndv
				multicell_psu_raster_arr[strata_raster_arr != strata_id] = mcr_ndv
				
				multicell_psu_raster_arr[ results_raster_arr != 1] = mcr_ndv
				
				
				# if strata id is zero it causes a problem with the polygonise raster so we have to change it to another value
				if strata_id == 0:
					multicell_psu_raster_arr[multicell_psu_raster_arr == 0] = 1
				
				# only make a shapefile if there is something to put in it
				if len(multicell_psu_raster_arr[ multicell_psu_raster_arr != mcr_ndv]):	
					
					# create a shapefile from the results
					mcr_band.WriteArray(multicell_psu_raster_arr)		
					sbu.polygonise_open_raster(mcr,results_shapefile,field_name='cl_id')			
					# reset raster band
					mcr_band.WriteArray(multicell_psu_raster_arr_orig)
					
					# because we have selected raster cells some go over the strata boundaries so cut to the strata boundary 
					if strata_file and strata_id_field:
						
						where = str(strata_id_field) + '=' + str(strata_id)
						sbu.clipShapefilewithShapefilebyAttribute(results_shapefile,results_shapefile,strata_file,where)
					
			else:
				# make a shapefile for this strata and merge with main results shapefile
				temp_shapefile_for_this_strata = working_dir + "/" + "temp_results_for_strata_" + str(strata_id) +".shp"
								
				# set values that haven't been selected to ndv
				multicell_psu_raster_arr[strata_raster_arr != strata_id] = mcr_ndv
				multicell_psu_raster_arr[ results_raster_arr != 1] = mcr_ndv
				
				# if strata id is zero it causes a problem with the polygonise raster so we have to change it to another value
				if strata_id == 0:
					multicell_psu_raster_arr[multicell_psu_raster_arr == 0] = 1
					
				# only make a shapefile if there is something to put in it
				if len(multicell_psu_raster_arr[ multicell_psu_raster_arr != mcr_ndv]):	
					# create a shapefile from the results
					mcr_band.WriteArray(multicell_psu_raster_arr)		
					sbu.polygonise_open_raster(mcr,temp_shapefile_for_this_strata,field_name='cl_id')				
					# reset raster band
					mcr_band.WriteArray(multicell_psu_raster_arr_orig)
					
					# because we have selected raster cells some go over the strata boundaries so cut to the strata boundary
					if strata_file and strata_id_field:
						where = str(strata_id_field) + '=' + str(strata_id)
						sbu.clipShapefilewithShapefilebyAttribute(temp_shapefile_for_this_strata,temp_shapefile_for_this_strata,strata_file,where)
					
					# now merge with results shapefile - we are building up the results shapefile bit by bit
					sbu.mergeShapefiles(results_shapefile,temp_shapefile_for_this_strata)

	# we need to get the strata id matching each index
	# which we save in a dictionary
	selected_cell_ids = multicell_psu_raster_arr_orig[(pop_raster_grouped_values > 0) & (multicell_psu_raster_arr_orig !=mcr_ndv)]  
	selected_strata_ids = strata_raster_arr[(pop_raster_grouped_values > 0) & (multicell_psu_raster_arr_orig !=mcr_ndv)]
	selected_pop = pop_raster_grouped_values[(pop_raster_grouped_values > 0) & (multicell_psu_raster_arr_orig !=mcr_ndv)]
	previous_cell_id = -100
	for i in range(len(selected_cell_ids)):
				
		if selected_cell_ids[i] == previous_cell_id:
			continue
		
		previous_cell_id = selected_cell_ids[i]
		
		strata_id = selected_strata_ids[i]
		household_size_weighting = getHouseholdSizeWeighting(cfg_hh_size,strata_id)
		pop_in_cluster= selected_pop[i]
			
		households_in_cluster = round(pop_in_cluster/household_size_weighting,2)
		
		# to indicate that this cell was selected via the main algorithm
		row = {'strata_id':strata_id,'cl_type':'main','s_cl_pop':pop_in_cluster,'s_cl_hh':households_in_cluster}
		selected_PSU_info[selected_cell_ids[i]] = row 
		
	## the spatial raster is a larger scale grid we want to make sure that at least one cell in each larger grid is selected
	# we pick a random empty cell for each one that isn't selected
	# this can take a long time if the spatial scale is small				
	if spatial_raster:	
		if enable_debug:
			print "doing spatial oversample - this can be time consuming"

		# spatial rasters is a bunch of repeated ids forming square blocks depending on the resolution chosen for the oversample
		sp_r = gdal.Open(spatial_raster)
		sp_r_band = sp_r.GetRasterBand(1)
		sp_r_ndv = sp_r_band.GetNoDataValue()
		sp_r_arr = np.array(sp_r_band.ReadAsArray())
		
		# get a list of each oversample id: will will want to check there is a selection within each id (group) 
		spatial_ids =  np.unique(sp_r_arr[sp_r_arr!= -99])
		
		# make new clean copies for the spatial raster
		strata_raster_arr_copy = np.array(s_band.ReadAsArray())
		results_raster_arr_copy = np.array(results_raster_band.ReadAsArray())
			
		results_raster_locations_with_values = (results_raster_arr !=results_raster_ndv)
		
		for sp_id in spatial_ids:
			
			if enable_debug:
				print "sp_id" + str(sp_id)
			# the locations that have this spatial id that aren't ndv
			# do just once here to speed things up
			locations_of_interest = (sp_r_arr == sp_id) & results_raster_locations_with_values
						
			# get results raster values for this spatial id
			rr_vals = results_raster_arr[ locations_of_interest ]	
			
			# check if any have been chosen - at least one cell will have the value 1 - otherwise all the values will be 0
			# I think this works OK for spatial oversample as we don't have to worry about exclusions
			if len(rr_vals) and np.sum(rr_vals) == 0:
			
				# pick a random empty PSU
				multicell_ids_to_choose_from = np.unique(multicell_psu_raster_arr[locations_of_interest])
				
				if len(multicell_ids_to_choose_from):
					
					index_chosen = 0
					if len(multicell_ids_to_choose_from) > 1:
						index_chosen = random.randint(0,len(multicell_ids_to_choose_from)-1)
					
					multicell_id_chosen = multicell_ids_to_choose_from[index_chosen]
					
					# get all of the locations of interest with that multicell id
					locations_of_interest = locations_of_interest  & (multicell_psu_raster_arr == multicell_id_chosen)
					
					# update the results raster - by setting them to 1 where selected
					results_raster_arr[ locations_of_interest ] = 1
					# this is just the locations chosen via this section of code
					results_raster_arr_copy[ locations_of_interest ] = 1
					
					cells_allocated += 1	
			
					strata_id_vals = strata_raster_arr[ locations_of_interest ]
					
					# a multicell can cross multiple strata
					# for each strata crossed increment the clusters in strata - we could try incrementing just the one most covered
					#strata_id = strata_id_vals[index_chosen]
					for strata_id in np.unique(strata_id_vals):						
						# save the number of sampled clusters in this stratum
						if strata_id in selected_strata_info.keys():						
							selected_strata_info[strata_id]['s_cl_st'] = selected_strata_info[strata_id]['s_cl_st'] + 1
						else:					
							# this is the first PSU selected 
							row = {'s_cl_st':1}
							selected_strata_info[strata_id] = row 
				
				
		# if any new locations have been selected update the results shapefile
		if np.sum(results_raster_arr_copy[results_raster_arr_copy!=results_raster_ndv]) > 0:
			
			# make a shapefile for this strata and merge with main results shapefile
			temp_shapefile_for_this_strata = working_dir + "/" + "temp_results_for_strata_spatial_oversample.shp"
						
			multicell_psu_raster_arr_copy = np.array(mcr_band.ReadAsArray())
			
			# set values that haven't been selected to ndv
			multicell_psu_raster_arr_copy[ results_raster_arr_copy != 1] = mcr_ndv
			
			# create a shapefile from the results
			mcr_band.WriteArray(multicell_psu_raster_arr_copy)		
			sbu.polygonise_open_raster(mcr,temp_shapefile_for_this_strata,field_name="cl_id")				
			# reset raster band
			mcr_band.WriteArray(multicell_psu_raster_arr_orig)
			
			
			# now merge with results shapefile - we are building up the results shapefile bit by bit
			sbu.mergeShapefiles(results_shapefile,temp_shapefile_for_this_strata)
					


	# finally overwrite the original results raster
	results_raster_band.WriteArray(results_raster_arr)
	results_raster_band.FlushCache()

	# remove unwanted files
	if delete_intermediate_outputs:
		for f in glob.glob(working_dir + "/" + "temp_results_for_strata_*"):
			os.remove(f)
		

	rr = None
	ir = None
	
	
	# we need to get the strata id matching each index
	selected_cell_ids = multicell_psu_raster_arr_orig[(pop_raster_grouped_values !=pr_ndv) & (multicell_psu_raster_arr_orig !=mcr_ndv)]  
	selected_strata_ids = strata_raster_arr[(pop_raster_grouped_values !=pr_ndv) & (multicell_psu_raster_arr_orig !=mcr_ndv)]
	selected_pop = pop_raster_grouped_values[(pop_raster_grouped_values !=pr_ndv) & (multicell_psu_raster_arr_orig !=mcr_ndv)]
			
	####################################
	
	previous_cell_id = -100
	for i in range(len(selected_cell_ids)):
		
		if selected_cell_ids[i] == previous_cell_id:
			continue
			
		previous_cell_id = selected_cell_ids[i]
		
		# if id not in list then it was selected as part of the oversample algorithm
		if not selected_cell_ids[i] in selected_PSU_info.keys():
			# to indicate that this cell was selected via the oversample algorithm
			strata_id = selected_strata_ids[i]
			household_size_weighting = getHouseholdSizeWeighting(cfg_hh_size,strata_id)
			pop_in_cluster= selected_pop[i]
			
			households_in_cluster = round(pop_in_cluster/household_size_weighting,2)
			
			row = {'strata_id':strata_id,'cl_type':'oversample','s_cl_pop':pop_in_cluster,'s_cl_hh':households_in_cluster}

			selected_PSU_info[selected_cell_ids[i]] = row 
		
		# add number of sampled clusters in stratum info 
		# and total number of selected clusters	
		strata_id = selected_PSU_info[selected_cell_ids[i]]['strata_id']
		# this should always be the case but check anyway
		if strata_id in selected_strata_info.keys():
			selected_PSU_info[selected_cell_ids[i]]['s_cl_st'] = selected_strata_info[strata_id]['s_cl_st']
			selected_PSU_info[selected_cell_ids[i]]['s_cl_tot'] = cells_allocated
		else:
			selected_PSU_info[selected_cell_ids[i]]['s_cl_st'] = 0
			selected_PSU_info[selected_cell_ids[i]]['s_cl_tot'] = cells_allocated
	
		
	return selected_PSU_info
	

def delete_intermediate_files(working_dir):

	# the user might choose not to delete all intermediate outputs for debug, or further use such as creating the summary page
	# this function should delete them all
	
	# It doesn't remove gridEA files as they take ages to create and we might want to keep them for further runs
	
	for f in glob.glob(working_dir + "/ghs_mod_cut*"):
		os.remove(f)
				
	for f in glob.glob(working_dir + "/index_raster*"):
		os.remove(f)

	for f in glob.glob(working_dir + "/pop_raster*"):
		os.remove(f)

	for f in glob.glob(working_dir + "/results_raster*"):
		os.remove(f)
				
	for f in glob.glob(working_dir + "/spatial_oversample*"):
		os.remove(f)
				
	for f in glob.glob(working_dir + "/strata_raster*"):
		os.remove(f)			

	for f in glob.glob(working_dir + "/buffered_boundary_shapefile_*"):
		os.remove(f)			
					
	for f in glob.glob(working_dir + "/boundary_pixels_*"):
		os.remove(f)	
					
	for f in glob.glob(working_dir + "/" + "temp_results_for_strata_*"):
		os.remove(f)	

	for f in glob.glob(working_dir + "/" + "coverage_poly_*"):
		os.remove(f)

	for f in glob.glob(working_dir + "/" + "*strata_poly*"):
		os.remove(f)		

