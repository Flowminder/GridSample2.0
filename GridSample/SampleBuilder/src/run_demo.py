

import Mainfx_phase2_v1 as mfx

###############################################
# This is an example of creating a grid sample for a stratified, sub-national area, based on a WorldPop population raster
# excluding cells with a population below 0.5 and with an oversample grid of 10km*10km

## You will need to set the two paths below if running outside of a docker environment

# directory where the sample rasters and shapefiles are found
src_files_dir = '/test_files/'

# Directory for output
results_dir = '/results/'

# Unique ID for this build/project - output will be put in a subdir of results_dir with this id
userID = 1

# the location of the GHS raster - this is mandatory as it is used to do an u/r stratification behind the scenes
ghs_mod_raster = src_files_dir + "GHS_SMOD_tap_WGS84.tif"

# the population raster location
pop_file = src_files_dir + 'PHL_ppp_v2b_2015.tif'

# filename of coverage shapefile
coverage_poly_filename = src_files_dir + 'phl_8_areas.shp'

# Desired filename for output shapefile
output_file_name = 'gridsample_output.shp'

# In this case a dictionary with keys the strata id and value the required number of psu for that strata
cfg_psu_per_strata = {8:100,9:100,19:100,20:100,21:100,22:100,23:100,24:100}

# If no coverage polygon is given the coverage area can be set to urban or rural
# for this test run it is set to None as there is a coverage shapefile
#  options: urbanOnly, ruralOnly : uses the national boundary but only includes u/r according to GHS-MOD
coverage_ur_option = None

# If stratification is required then the method to be used
# None for no stratification , else one of urbanRural,adminArea,custom_upload
# In this case a shapefile is being used
cfg_stratification_method = "custom_upload"

## if stratifying using a shapefile the next 3 parameters must be set
# strata shapefile name
strata_poly_filename = src_files_dir + "phl_8_areas.shp"
# name of variable (shapefile attribute) that defines unique strata (ex. "ID")
strata_ID_field = "ID_1" 
#  name of variable (shapefile attribute)  that defines unique strata name
strata_name_field = "NAME_1" 

# Size of each household. This value will be used to divide the population raster, so that all calculations are being done in terms of households
# either an integer if all strata have the same value or a dictionary of strata id and values
# if not set we just work in terms of population
# mandatory when making a grid sample : not required for histogram or strata population calls
cfg_hh_size = 4.5 

# one and only one of the following options may be chosen to exclude cells in the population raster with low values

# Minimum population in a cell that will be used to exclude that cell entirely
cfg_exclude_pop_per_cell = 0.5
# If true excldue cells that GHS-MOD says are un-populated
cfg_exclude_ghssmod0_bool = False
 
# an optional Random number seed for code 
cfg_random_number = 1234

# indicate whether using a standard population raster with each cell a PSU
# a shapefile to indicate PSU
#or a population raster and a multicell raster to indicate the PSUs 
# one of single,own, or multi
cfg_frame_type = 'single'


# optional size in sq km for resampling all rasters (None means don't resample)
# assumes population raster is comprised of 100m*100m pixels
cfg_resample_size = None 
	
# The gridsample function has 3 possible outputs and options which are create_sample, create_histogram, get_strata_pop_values
# create_sample: make a grid sample
# create_histogram: create up to two histograms shwing the PSU population distributions
# get_strata_pop_values: creates a CSV of information about each strata containing information such as the population in the strata
action = 'create_sample'

# Size of cell to be used to make an oversized grid that we will ensure at least one cell of each has been selected from
# An integer multiple of the population raster pixel size
cfg_oversample_grid_spatial_scale = 100

# The following 3 options are used if the PSUs are obtained from a shapefile, or via the gridEA algorithm

# if cfg_frame_type is "own" then the location of the shapefle containing the PSUs.
# if cfg_frame_type is "multi" then either the location of the raster containing the cell ids or nothing in which case cfg_multi_cell_cluster_size must be set and the file will be generated by running the gridEA algorithm
cfg_own_frame_file = None
				
# if cfg_frame_type is "own" then the name of the shapefile id columns for the shapefile containing the PSUs.
cfg_own_frame_id = None
				
# if cfg_frame_type is "multi" then one of small,medium large - the gridEA algorithm option
cfg_multi_cell_cluster_size = None

# The following 3 options are optional when gridEA is chosen

# if user chooses multi-cell (GridEA) then they can optionally give a shapefile indicating the strata to use
# if None given then the entire country is treated as one strata (which might causes the algorithm problems/ make it very slow)
grid_ea_strata_file = None
# if a grid_ea_strata_file is given then the the name of the id field
grid_ea_strata_file_id_field = None

# Used to help decide if parameters values  have changed to decide whether it is necessary to rebuild gridEA as it's slow
# The grid EA parameters are optional stored in a file. If the parameters passed in have changed from those used in the last run then
# the algorithm is not re-run
# gridEA_parameters_string = pop_raster + str(grid_ea_coverage_info) + str(cfg_exclude_ghssmod0_bool) + str(cfg_multi_cell_cluster_size)
grid_ea_coverage_info = None

print "Grid Sample Demo"

# set to 1 to enable this section of code, or to 0 to disable it
if 1:
	print "About to run test"

	testout = mfx.gridsample(
						
		output_dir = results_dir,
		uniqueID = userID,
							
		ghs_mod_raster = ghs_mod_raster,
						
		pop_raster = pop_file, 
						
		cfg_psu_per_strata = cfg_psu_per_strata,
						
		PSU_filename = output_file_name,
														
		coverage_poly_filename = coverage_poly_filename,
						
		coverage_ur_option = coverage_ur_option ,
						
		strata_poly_filename = strata_poly_filename, 
		strata_ID_field = strata_ID_field, 
		strata_name_field = strata_name_field,
							
		cfg_stratification_method = cfg_stratification_method,
					  
		cfg_hh_size = cfg_hh_size,
						
		
		cfg_exclude_pop_per_cell = cfg_exclude_pop_per_cell,
		cfg_exclude_ghssmod0_bool = cfg_exclude_ghssmod0_bool,
				  
		cfg_random_number = cfg_random_number,
									
		cfg_frame_type = cfg_frame_type,
		
		cfg_resample_size = cfg_resample_size, 
		
		action = action,			
						
		cfg_oversample_grid_spatial_scale = cfg_oversample_grid_spatial_scale,
						
		cfg_own_frame_file = cfg_own_frame_file,
		cfg_own_frame_id = cfg_own_frame_id,
		cfg_multi_cell_cluster_size = cfg_multi_cell_cluster_size,

		grid_ea_strata_file =grid_ea_strata_file	,

		grid_ea_strata_file_id_field = grid_ea_strata_file_id_field,

		grid_ea_coverage_info = grid_ea_coverage_info
		
						)
						
	# an array [status,message]
	# status =1 for success, 0 for failure
	# if a file then an error message is the second array parameter	
	print "Test completed"	
	print testout

# If the delete_intermediate_outputs flag is set to false to check for debug, or to make a summary pdf
# they can be deleted with this function
mfx.delete_intermediate_files( results_dir+ str(userID) )