rm(list=ls())

#session info:
#R version 3.4.2 (2017-09-28)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows >= 8 x64 (build 9200)

library(sp)           #version: sp_1.2-7
library(raster)       #version: raster_2.6-7
library(rgdal)        #version: rgdal_1.2-18
library(parallel)     #base package
#memory.limit(5000000)

ncores <- detectCores() - 1

### gridEA produces a raster of EAs 

#population_raster: Raster* layer
#settlement_raster: Raster* layer. Values for different settlement types
#strata_raster: Raster* layer. Values for different strata types
#exclude_unsettled. logical. Default is FALSE.          ...Alternatively, unsettled cells could already be NA in the settlement raster, in which case exclude_unsettled should be set to FALSE 
#unsettled_ID. numeric. Needs to be defined if exclude_unsettled = TRUE. A single value corresponding to the value for unsettled cells in the settlement_raster. Default is NA. Not that all cells in settlement_raster will automatically be excluded from the region for which EAs are generated.
#using_ghs_smod_pop2015. logical. Default is FALSE
#predefined_EA_size: logical. Default is TRUE
#EA_target_size: If predefined_EA_size is TRUE. Categorical. "small", "medium", "large". Default is "medium"
#EA_by_hh: logical. Default is FALSE. 
#target_hh_per_EA: If EA_by_hh is TRUE. Element or data frame. Element if a single household size per EA is to be used across all strata and settlement types. Data frame if different houshold sizes per EA across strata and settlement types are needed. Data frame should contain 4 columns: "strata_ID_number", "settlement_type_ID_number", "target_number_hh_per_EA", "pop_number_per_hh". The "strata_ID_number" and "settlement_type_ID_number" must correspond to the sett/strata values. Default = 100
#pop_per_hh: If EA_by_hh is TRUE and target_hh_per_EA is a single number. Element. A single population size per hh is to be used. Default = 5
#EA_by_pop: logical. Default is FALSE. 
#target_pop_per_EA: If EA_by_pop is TRUE. Element or data frame. Element if a single population size per EA is to be used across all strata and settlement types. Data frame if different population per EA across strata and settlement types are needed. Data frame should contain 3 columns: "strata_ID_number", "settlement_type_ID_number", "target_pop_number_per_EA". The "strata_ID_number" and "settlement_type_ID_number" must correspond to the sett/strata values. Default = 500
#max_cells_per_EA: If predefined_EA_size = FALSE. Numeric. Default = 900 which is approx. 3km by 3 km when using worldpop 100m x 100m population rasters, depending on location   
#run_ID: character. This will appear at the end of output file. Allows easy distinction between resultant sampling frames from different runs of the code, e.g. for different countries or for different initial parameters. Default = "_run1"


gridEA <- function(population_raster, settlement_raster, strata_raster, 
                   exclude_unsettled = FALSE, unsettled_ID = NA, 
                   using_ghs_smod_pop2015 = FALSE,
                   predefined_EA_size = TRUE, EA_target_size = "medium", 
                   EA_by_hh = FALSE, target_hh_per_EA = 100, pop_per_hh = 5,
                   EA_by_pop = FALSE, target_pop_per_EA = 500, max_cells_per_EA = 900,
                   output_path, run_ID = "_run1"){
  if(file.exists(paste(output_path,"/EA_raster_master", run_ID, ".tif",sep=""))){
    stop(paste("Desired output filename already exists; delete file or change run_ID so that output file is given a different name. File location and name is:", output_path,"/EA_raster_master", run_ID, ".tif",sep=""))
  }
  if(file.exists(paste(output_path,"/EA_pop_raster_master", run_ID, ".tif",sep=""))){
    stop(paste("Desired output filename already exists; delete file or change run_ID so that output file is given a different name. File location and name is:", output_path,"/EA_pop_raster_master", run_ID, ".tif",sep=""))
  }
  if (sp::proj4string(population_raster) != "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") {
    population_raster <- projectRaster(population_raster, crs = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  }
  if (sp::proj4string(settlement_raster) != "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") {
    settlement_raster <- projectRaster(settlement_raster, crs = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  }
  if (sp::proj4string(strata_raster) != "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") {
    strata_raster <- projectRaster(strata_raster, crs = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  }
  # need extra check here to make sure extents match 
  if(exclude_unsettled == TRUE){
    if(is.numeric(unsettled_ID) == FALSE | !length(unsettled_ID) == 1){
      stop("unsettled_ID needs to be a single number")
    }
    settlement_raster[settlement_raster[] == unsettled_ID] <- NA
  }  
  if(using_ghs_smod_pop2015 == TRUE){
    if(exclude_unsettled == FALSE){
      settlement_raster[settlement_raster[] == 0] <- 1
    }
    settlement_raster[settlement_raster[] == 2] <- 1
  }
  settlement_raster[which(is.na(population_raster[]))] <- NA
  strata_raster[which(is.na(population_raster[]))] <- NA
  settlement_raster[which(is.na(strata_raster[]))] <- NA
  strata_raster[which(is.na(settlement_raster[]))] <- NA
  sett_vals <- sort(unique(getValues(settlement_raster))) # including sort excludes NA from the vector
  strat_vals <- sort(unique(getValues(strata_raster)))
  if(predefined_EA_size == TRUE){
    if(EA_by_hh == TRUE | EA_by_pop == TRUE){
      stop("Need to specify predefined_EA_size = FALSE if you want to specify a target number of households or people per Enumeration Zone")
    }
    if(EA_target_size == "small"){
      target_pop_per_EA <- 75
      max_cells_per_EA <- 100
    }
    if(EA_target_size == "medium"){
      target_pop_per_EA <- 500
      max_cells_per_EA <- 900
    }
    if(EA_target_size == "large"){
      target_pop_per_EA <- 1200
      max_cells_per_EA <- 2500
    }
    print(paste("Creating EZs based on predefined EZ size criteria: target population per EZ = ", target_pop_per_EA, " & max number of cells per EA = ", max_cells_per_EA))
  }else{
    if(EA_by_hh == TRUE & EA_by_pop == TRUE){
      stop("Both EA_by_hh & EA_by_pop can't be TRUE. Set just one to be TRUE")
    }
    if(EA_by_hh == TRUE){ 
      if(class(target_hh_per_EA) == "numeric"){
        ### add check in here for pop_per_hh being numeric too
        target_pop_per_EA <- target_hh_per_EA * pop_per_hh
      }else{
        if(!class(target_hh_per_EA) == "data.frame"){
          stop(paste("target_hh_per_EA has class", class(target_hh_per_EA), "; target_hh_per_EA needs to be a single number or a data.frame object, see ?class"))
        }
        if(!ncol(target_hh_per_EA) == 4){
          stop(paste("target_hh_per_EA has", ncol(target_hh_per_EA), " columns; target_hh_per_EA should have 4 columns - 'strata_ID_number', 'settlement_type_ID_number', 'target_number_hh_per_EA', 'pop_number_per_hh'"))
        }
        if(FALSE %in% (c("strata_ID_number", "settlement_type_ID_number", "target_number_hh_per_EA", "pop_number_per_hh") %in% names(target_hh_per_EA))){
          stop("target_hh_per_EA column names need to be 'strata_ID_number', 'settlement_type_ID_number', 'target_number_hh_per_EA', 'pop_number_per_hh', see names(target_hh_per_EA) to check")
        }
        if(exclude_unsettled == TRUE){     # leave here rather than further down ifs in case user puts NAs for the unsettled entries
          target_hh_per_EA <- target_hh_per_EA[-(target_hh_per_EA$settlement_type_ID_number == unsettled_ID),]
        } 
        if(TRUE %in% is.na(unlist(target_hh_per_EA))){
          stop("target_hh_per_EA contains NAs, all entries should be numeric with 'strata_ID_number' and 'settlement_type_ID_number' entries corresponding to the values in strata_raster and settlement_raster, respectively")
        }   
        if(!is.numeric(unlist(target_hh_per_EA))){
          stop("target_hh_per_EA contains non-numerics, see is.numeric(unlist(target_hh_per_EA)). All entries should be numeric with 'strata_ID_number' and 'settlement_type_ID_number' entries corresponding to the values in strata_raster and settlement_raster, respectively")
        }   
        if(FALSE %in% (sort(unique(getValues(strata_raster))) %in% target_hh_per_EA$strata_ID_number)){
          stop("incorrect strata_ID_number values in target_hh_per_EA; all values in strata_raster need to be included in target_hh_per_EA")
        }
        if(FALSE %in% (sort(unique(getValues(settlement_raster))) %in% target_hh_per_EA$settlement_type_ID_number)){
          stop("incorrect settlement_type_ID_number values in target_hh_per_EA; all values settlement_raster need to be included in target_hh_per_EA")
        }
        # because not all sett types may be in a stratum, user may not specify (number of strata x number of settlement) rows, if there are missing combos this is picked up later in code
        target_pop_number_per_EA <- target_hh_per_EA$target_number_hh_per_EA * target_hh_per_EA$pop_number_per_hh
        target_pop_per_EA <- as.data.frame(cbind(target_hh_per_EA$strata_ID_number,target_hh_per_EA$settlement_type_ID_number, target_pop_number_per_EA))
      }
      print("EAs will be created based on user's specified number of households per EA")
    }
    if(EA_by_pop == TRUE){
      if(!class(target_pop_per_EA) == "numeric"){
        if(!class(target_pop_per_EA) == "data.frame"){
          stop(paste("target_pop_per_EA has class", class(target_pop_per_EA), "; target_pop_per_EA needs to be a single number or a data.frame object, see ?class"))
        }
        if(!ncol(target_hh_per_EA) == 3){
          stop(paste("target_pop_per_EA has", ncol(target_pop_per_EA), " columns; target_pop_per_EA should have 3 columns - 'strata_ID_number', 'settlement_type_ID_number', 'target_pop_number_per_EA'"))
        }
        if(FALSE %in% (c("strata_ID_number", "settlement_type_ID_number", "target_pop_number_per_EA") %in% names(target_pop_per_EA))){
          stop("target_pop_per_EA column names need to be 'strata_ID_number', 'settlement_type_ID_number', 'target_pop_number_per_EA', see names(target_pop_per_EA) to check")
        }
        if(exclude_unsettled == TRUE){     # leave here rather than further down ifs in case user puts NAs for the unsettled entries
          target_pop_per_EA <- target_pop_per_EA[-(target_pop_per_EA$settlement_type_ID_number == unsettled_ID),]
        } 
        if(TRUE %in% is.na(unlist(target_pop_per_EA))){
          stop("target_pop_per_EA contains NAs, all entries should be numeric with 'strata_ID_number' and 'settlement_type_ID_number' entries corresponding to the values in strata_raster and settlement_raster, respectively")
        }   
        if(!is.numeric(unlist(target_pop_per_EA))){
          stop("target_pop_per_EA contains non-numerics, see is.numeric(unlist(target_pop_per_EA)). All entries should be numeric with 'strata_ID_number' and 'settlement_type_ID_number' entries corresponding to the values in strata_raster and settlement_raster, respectively")
        }   
        if(FALSE %in% (sort(unique(getValues(strata_raster))) %in% target_pop_per_EA$strata_ID_number)){
          stop("incorrect strata_ID_number values in target_pop_per_EA; all values in strata_raster need to be included in target_pop_per_EA")
        }
        if(FALSE %in% (sort(unique(getValues(settlement_raster))) %in% target_pop_per_EA$settlement_type_ID_number)){
          stop("incorrect settlement_type_ID_number values in target_pop_per_EA; all values settlement_raster need to be included in target_pop_per_EA")
        }
        # because not all sett types may be in a stratum, user may not specify (number of strata x number of settlement) rows, if there are missing combos this is picked up later in code
      }
      print("EAs will be created based on user's specified population per EA")
    }    
    if(is.numeric(max_cells_per_EA) == FALSE | !length(max_cells_per_EA) == 1){
      stop("max_cells_per_EA needs to be a single number")
    }
    # might want to add a minimum value for max_cells_per_EA here
  }
  population_mat <- as.matrix(population_raster)
  n_rows <- dim(population_raster)[1]
  n_cols <- dim(population_raster)[2]
  strat_sett_df <- cbind(rep(strat_vals, each = length(sett_vals)), rep(sett_vals, length(strat_vals))) 
  print(paste("creating", " /temp_folder", run_ID, sep=""))       #keep this up here, well before trying to save to the folder
  dir.create(paste(output_path,"/temp_folder", run_ID, sep=""))   
  clump_IDs_mat <- matrix(NA, nrow = n_rows, ncol = n_cols)
  clumps_df <- NULL
  for(ssID in 1:nrow(strat_sett_df)){     
    i <- strat_sett_df[ssID, 1] #stratum ID
    j <- strat_sett_df[ssID, 2] #sett_type ID
    if(is.data.frame(target_pop_per_EA) == TRUE){
      if(!length(target_pop_per_EA$target_pop_number_per_EA[target_pop_per_EA$strata_ID_number == i & target_pop_per_EA$settlement_type_ID_number == j] == 1)){
        stop(paste("no row in user specifications for stratum", i, "and settlement type", j, "in either target_hh_per_EA or target_pop_per_EA data frame"))
      }
      ss_target_pop_per_EA <- target_pop_per_EA$target_pop_number_per_EA[target_pop_per_EA$strata_ID_number == i & target_pop_per_EA$settlement_type_ID_number == j] 
    }else{
      ss_target_pop_per_EA <- target_pop_per_EA
    }
    # Create raster of contiguous regions by settlement type so that an EA consists only of one settlement type
    sett_clump_raster <- setValues(settlement_raster, NA)
    sett_clump_raster[strata_raster[] == i & settlement_raster[] == j] <- 1
    if(length(sort(unique(getValues(sett_clump_raster)))) == 0) next   # if no cells of j sett type are in i state
    clump_IDs_raster <- clump(sett_clump_raster, directions = 4, gaps = FALSE)
    clump_IDs <- sort(unique(getValues(clump_IDs_raster)))
    if(class(clumps_df)=="matrix"){     # true after first loop (NULL in first loop)
      clump_IDs_raster <- clump_IDs_raster + max(clumps_df[,1]) 
    }
    clump_IDs <- sort(unique(getValues(clump_IDs_raster)))
    clumps_df <- rbind(clumps_df, cbind(clump_IDs, rep(ss_target_pop_per_EA, length(clump_IDs))))  
    temp_clump_mat <- as.matrix(clump_IDs_raster)
    clump_IDs_mat[!is.na(temp_clump_mat)] <- temp_clump_mat[!is.na(temp_clump_mat)]
  }
  print(Sys.time())
  print("done with clumps designation")
  clump_info <- list()
  for(clump_row in 1:nrow(clumps_df)){
    clump_ID <- clumps_df[clump_row, 1]
    cell_IDs <- which(clump_IDs_mat[] == clump_ID, arr.ind = TRUE)
    min_clump_row <- min(cell_IDs[,1])
    max_clump_row <- max(cell_IDs[,1])
    min_clump_col <- min(cell_IDs[,2])
    max_clump_col <- max(cell_IDs[,2])
    subclump_IDs_mat <- clump_IDs_mat[min_clump_row:max_clump_row, min_clump_col:max_clump_col]
    subclump_IDs_mat[which(!subclump_IDs_mat[] == clump_ID)] <- NA
    subpop_mat <- population_mat[min_clump_row:max_clump_row, min_clump_col:max_clump_col]
    subpop_mat[which(is.na(subclump_IDs_mat))] <- NA
    ss_target_pop_per_EA <- clumps_df[clump_row, 2]
    clump_info[[clump_row]] <- list(subpop_mat, c(min_clump_row, max_clump_row, min_clump_col,max_clump_col),ss_target_pop_per_EA)
  }
  r_xmn <- extent(population_raster)[1]
  r_xmx <- extent(population_raster)[2]
  r_ymn <- extent(population_raster)[3]
  r_ymx <- extent(population_raster)[4]
  r_crs <- crs(population_raster)
  # does deleting rasters stop them being duplicated across nodes, if they're not used does the node 'take' them???
  rm(population_raster); rm(settlement_raster); rm(strata_raster); rm(clump_IDs_raster); rm(sett_clump_raster)
  cl <- makeCluster(ncores)  
  print("starting cluster")
  clusterExport(cl, split(ls("package:raster"),1:length(ls("package:raster")))) #, envir = .GlobalEnv
  EA_gen_fn <- function(clump_x, max_cells_per_EA, output_path, run_ID){
    EA_cell_IDs <- list()
    ss_target_pop_per_EA <- clump_x[[3]]
    min_pop_per_EA <- ss_target_pop_per_EA * 0.75  # 25%
    max_pop_per_EA <- ss_target_pop_per_EA * 1.25  # 25%
    min_clump_row <- clump_x[[2]][1]
    max_clump_row <- clump_x[[2]][2]
    min_clump_col <- clump_x[[2]][3]
    max_clump_col <- clump_x[[2]][4]
    subpop_mat <- clump_x[[1]]
    clump_x[[1]] <- 0 # to reduce memory use
    if(class(subpop_mat) == "matrix"){                       #i.e. more than 1 cell (almost always)
      cell_IDs <- which(!is.na(subpop_mat), arr.ind = TRUE)
    }else{                                                      #i.e. just 1 cell (occasionally is other cells of smod block have NA pop or strata)
      cell_IDs <- matrix(c(1,1), nrow=1)
    }
    clump_pop_size <- sum(subpop_mat[cell_IDs])
    #if total pop of clump is smaller than min pop per EA AND total number of cells in clump (given by nrow(cell_IDs)) is smaller than max cells per EA, then don't split into multiple EAs
    if((clump_pop_size <= max_pop_per_EA) & (nrow(cell_IDs) <= max_cells_per_EA * 1.5)){
      cell_IDs[,1] <- cell_IDs[,1] + min_clump_row - 1
      cell_IDs[,2] <- cell_IDs[,2] + min_clump_col - 1
      EA_cell_IDs <- list(cell_IDs)                     # needs to be a list because other files will have cell_IDs for many EZs and so reading in code is based on list format
    }else{
      max_cell_pop_size <- max(subpop_mat[!(is.na(subpop_mat))]) 
      med_cell_pop_size <- median(subpop_mat[!(is.na(subpop_mat))]) 
      min_sett_row <- min_sett_col <- 1
      max_sett_row <- max_clump_row - (min_clump_row - 1)   
      max_sett_col <- max_clump_col - (min_clump_col - 1)   
      #if cell pops are to small across whole of clump, create EAs that are the max geographic size, else calc appropriate initial block size
      if((max_cell_pop_size * max_cells_per_EA) < min_pop_per_EA){
        block_horiz <- block_vert <- round(sqrt(max_cells_per_EA))
      }else{
        val1 <- round(sqrt(ss_target_pop_per_EA/med_cell_pop_size))
        if(val1 > round(sqrt(max_cells_per_EA))){
          block_horiz <- block_vert <- round(sqrt(max_cells_per_EA))
        }else{
          if(val1 %in% 1:2){ 
            block_vert <- block_horiz <- val1
          }else{
            if(val1 < 5){val2 <- val1 - 1}
            if(val1 >= 5){val2 <- val1 - 2}
            row_range <- max_sett_row - min_sett_row 
            col_range <- max_sett_col - min_sett_col 
            if(row_range >= col_range){
              block_vert <- val1
              block_horiz <- val2
            }else{
              block_vert <- val2
              block_horiz <- val1
            }
          }
        }
      }
      # ID location of max pop cell so that we can put in centre of a block and build block frame from there
      ID_max_pop_cell <- (which(subpop_mat == max_cell_pop_size, arr.ind = TRUE))[1,] 
      y_start <- min((ID_max_pop_cell[1] + floor(block_vert/2)), max_sett_row) #/2 so that max cell is middle of block
      x_start <- min((ID_max_pop_cell[2] + floor(block_horiz/2)), max_sett_col)
      seq_vert <- y_start; seq_horiz <- x_start   # so that these seqs include at least the initial val
      # build grid for blocks. Lots of ifs so that use alternate values if close to the edge of the raster section
      if((y_start - block_vert) > min_sett_row){seq_vert <- rev(seq(y_start, min_sett_row, -block_vert))}
      if((y_start + block_vert) < max_sett_row){seq_vert <- c(seq_vert, seq((y_start + block_vert), max_sett_row, block_vert))}
      if(!seq_vert[1] == min_sett_row){seq_vert <- c(min_sett_row, seq_vert)}
      if(!seq_vert[length(seq_vert)] == max_sett_row){seq_vert <- c(seq_vert, max_sett_row)}
      if((x_start - block_horiz) > min_sett_col){seq_horiz <- rev(seq(x_start, min_sett_col, -block_horiz))}
      if((x_start + block_horiz) < max_sett_col){seq_horiz <- c(seq_horiz, seq((x_start + block_horiz), max_sett_col, block_horiz))}
      if(!seq_horiz[1] == min_sett_col){seq_horiz <- c(min_sett_col, seq_horiz)}
      if(!seq_horiz[length(seq_horiz)] == max_sett_col){seq_horiz <- c(seq_horiz, max_sett_col)}
      block_ref_vert <- cbind(seq_vert[1:(length(seq_vert)-1)], c((seq_vert[2:(length(seq_vert)-1)]-1), max_sett_row))
      block_ref_horiz <- cbind(seq_horiz[1:(length(seq_horiz)-1)], c((seq_horiz[2:(length(seq_horiz)-1)]-1), max_sett_col))
      block_ref_df <- cbind(block_ref_vert[rep(seq_len(nrow(block_ref_vert)), nrow(block_ref_horiz)),], block_ref_horiz[rep(seq_len(nrow(block_ref_horiz)), each = nrow(block_ref_vert)),])
      # assign temporary EA IDs to all blocks within extent of clump. Extent defined by min/max_sett_row/col
      EA_mat <- matrix(NA, nrow = max_sett_row, ncol = max_sett_col)  
      for(block_ID in 1:nrow(block_ref_df)){    # loop instead of function or apply because want vals in one matrix to be updated
        EA_mat[unique(block_ref_df[block_ID,1]:block_ref_df[block_ID,2]), unique(block_ref_df[block_ID,3]:block_ref_df[block_ID,4])] <- block_ID
      }
      # mask locations in EA matrix where cells are not in clump
      EA_mat[which(is.na(subpop_mat))] <- NA    
      EA_ID_vec <- unique(EA_mat[which(!is.na(EA_mat))])  
      EA_pop <- sapply(EA_ID_vec, function(EA_index){sum(subpop_mat[which(EA_mat == EA_index)])})
      # identify EAs with pop below min pop/EA
      small_pop_EAs <- EA_ID_vec[which(EA_pop < min_pop_per_EA)]
      # of these small EAs, exclude from small EA list those that meet geog size 
      cells_count_test <- table(EA_mat[EA_mat %in% small_pop_EAs]) > (max_cells_per_EA * 0.8) #have  80% of max limit here, if an EA is almost at the lim don't want to add to another block 
      if("TRUE" %in% cells_count_test){
        max_geog_in_small_pop_EAs <- as.numeric(names(cells_count_test[cells_count_test == "TRUE"]))
        small_pop_EAs <- small_pop_EAs[!small_pop_EAs %in% max_geog_in_small_pop_EAs]
      }
      # re-order small_pop_EAs by geog size 
      small_pop_EAs <- as.numeric(names(sort(table(EA_mat[EA_mat %in% small_pop_EAs])))) 
      if(length(small_pop_EAs)>0){
        print("sorting small EAs")
        #plot(raster(EA_mat), col=col1)
        repeat{
          EA_ID <- small_pop_EAs[1]
          # neighbour EAs:
          EA_cell_IDs <- which(EA_mat == EA_ID, arr.ind = TRUE)
          nb_EAs <- cell_nbs <- NULL
          for(df_row in 1:nrow(EA_cell_IDs)){
            loc_vert <- EA_cell_IDs[df_row,1]
            loc_horiz <- EA_cell_IDs[df_row,2]
            adj_cells <- cbind(c((loc_vert + 1), (loc_vert - 1), loc_vert, loc_vert), c(loc_horiz, loc_horiz, (loc_horiz + 1), (loc_horiz - 1)))
            if(TRUE %in% (!adj_cells[,1] %in% 1:max_sett_row)){adj_cells <- adj_cells[-which((!adj_cells[,1] %in% 1:max_sett_row)),]}
            if(TRUE %in% (!adj_cells[,2] %in% 1:max_sett_col)){adj_cells <- adj_cells[-which((!adj_cells[,2] %in% 1:max_sett_col)),]}
            cell_nbs <- EA_mat[adj_cells]
            nb_EAs <- c(nb_EAs, cell_nbs)
          }
          nb_EAs <- unique(nb_EAs)
          nb_EAs <- nb_EAs[!is.na(nb_EAs) & !nb_EAs == EA_ID]
          if(EA_ID %in% nb_EAs){nb_EAs <- nb_EAs[-which(nb_EAs == EA_ID)]}  ### WHY IS THIS LINE NEEDED???? DOES THE LINE ABOVE NOT DO THIS???
          # measure compactness for each of the possible combinations of blocks/EAs
          compact_vals <- NULL
          for(new_EA in nb_EAs){
            test_EA_locs <- which(EA_mat == EA_ID | EA_mat == new_EA ,arr.ind = TRUE)
            edge_peri <- 0
            for(df_row in 1:nrow(test_EA_locs)){
              loc_vert <- test_EA_locs[df_row,1]
              loc_horiz <- test_EA_locs[df_row,2]
              adj_cells <- cbind(c((loc_vert + 1), (loc_vert - 1), loc_vert, loc_vert), c(loc_horiz, loc_horiz, (loc_horiz + 1), (loc_horiz - 1)))
              d_test <- duplicated(rbind(test_EA_locs,adj_cells))
              d_test <- d_test[(nrow(test_EA_locs)+1):length(d_test)] #FALSES for EA_cell_IDs that are not duplicated in cell_IDs for one side of split
              edge_count <- length(which(d_test == FALSE))
              edge_peri <- sum(edge_peri, edge_count)
            }
            compact_vals <- c(compact_vals, (edge_peri^2)/nrow(test_EA_locs)) 
          }
          EA_ID_to_add <- (nb_EAs[compact_vals == min(compact_vals)])[1] 
          EA_mat[EA_mat == EA_ID_to_add] <- EA_ID   #ie give new EA original small EA's ID
          #plot(raster(EA_mat), col=col1)
          if(EA_ID_to_add %in% small_pop_EAs){small_pop_EAs <- small_pop_EAs[-which(small_pop_EAs == EA_ID_to_add)]} #if block being added to the small pop EA (EA_ID) is on small_pop_EA list, it need to be taken off
          if(sum(subpop_mat[which(EA_mat == EA_ID)]) > min_pop_per_EA){ #if newly constructed EA has bigger pop than min/EA take off small_EA list
            small_pop_EAs <- small_pop_EAs[-which(small_pop_EAs == EA_ID)]
          }else{  #if newly constructed EA does not have bigger pop than min/EA, check if it has bigger geog than geog max/EA, if yes then take off small_EA list
            if(length(which(EA_mat[] == EA_ID)) > (max_cells_per_EA * 0.8)){
              small_pop_EAs <- small_pop_EAs[-which(small_pop_EAs == EA_ID)]
            }
          } #if newly constructed EA doesn't meet min pop or max geog per EA then leave it on small EA list so will repeat adding block to it 
          if(!length(small_pop_EAs)>0) break
          small_pop_EAs <- as.numeric(names(sort(table(EA_mat[EA_mat %in% small_pop_EAs])))) 
        }
        EA_ID_vec <- unique(EA_mat[which(!is.na(EA_mat))]) 
        EA_pop <- sapply(EA_ID_vec, function(EA_index){sum(subpop_mat[which(EA_mat == EA_index)])})
      }
      vbig_pop_EAs <- EA_ID_vec[which(EA_pop >= (2 * ss_target_pop_per_EA))] #identify very big pop EAs,if over 2*target pop (very big) then can split in two, along axis with most cells
      if(length(vbig_pop_EAs)>0){
        print("sorting vbig EAs") 
        EA_IDs_new <- (max(EA_ID_vec) + 1):(max(EA_ID_vec) + 2)
        repeat{
          EA_ID <- vbig_pop_EAs[1]
          EA_cell_IDs <- which(EA_mat == EA_ID, arr.ind = TRUE)
          tot_pop_size <- sum(subpop_mat[EA_cell_IDs])
          #create list, an element for each potential split, element contains matrix of cell row & col ids one EA of the potential split, do this because can't assume EAs are only four edged shapes
          row_IDs <- sort(unique(EA_cell_IDs[,1]))
          col_IDs <- sort(unique(EA_cell_IDs[,2]))
          if(length(row_IDs) >= length(col_IDs)){
            side_IDs <- row_IDs
            min_side_IDs <- min(row_IDs)
            row_or_col <- 1
          }else{
            side_IDs <- col_IDs
            min_side_IDs <- min(col_IDs)
            row_or_col <- 2
          }
          cell_IDs_ls <- lapply(side_IDs[1:(length(side_IDs) - 1)], function(x){matrix(EA_cell_IDs[EA_cell_IDs[,row_or_col] %in% (min_side_IDs:x),], ncol=2)})    
          pop_sizes <- sapply(cell_IDs_ls, function(x){sum(subpop_mat[x])})
    #diff_in_pops <- abs(tot_pop_size - (2*pop_sizes))  # i.e. diff in the two potential EAs' pop sizes. Want |C - B| and C = A - B, so |C - B| = |A - B - B| 
    #test_splits <- order(diff_in_pops)
          pop_sizes_second <- tot_pop_size - pop_sizes
          keepers1 <- which(pop_sizes >= min_pop_per_EA)
          keepers2 <- which(pop_sizes_second >= min_pop_per_EA)
          test_splits <- keepers1[which(keepers1 %in% keepers2)]
          if(!length(test_splits) > 0){
            vbig_pop_EAs <- vbig_pop_EAs[-(which(vbig_pop_EAs == EA_ID))]
            if(!length(vbig_pop_EAs)>0){ 
              break
            }else{
              next    
            }
          }
          compact_vals <- NULL
          for(test_EA_split in test_splits){
            for(x in 1:2){
              if(x==1){test_EA_locs <- cell_IDs_ls[[test_EA_split]]}
              if(x==2){
                d1 <- duplicated(rbind(cell_IDs_ls[[test_EA_split]],EA_cell_IDs))
                d1 <- d1[(nrow(cell_IDs_ls[[test_EA_split]])+1):length(d1)] #FALSES for EA_cell_IDs that are not duplicated in cell_IDs for one side of split
                test_EA_locs <- EA_cell_IDs[which(d1 == FALSE),]
              }
              if(!class(test_EA_locs) == "matrix"){test_EA_locs <- matrix(test_EA_locs, nrow=1)}
              edge_peri <- 0
              for(df_row in 1:nrow(test_EA_locs)){
                loc_vert <- test_EA_locs[df_row,1]
                loc_horiz <- test_EA_locs[df_row,2]
                adj_cells <- cbind(c((loc_vert + 1), (loc_vert - 1), loc_vert, loc_vert), c(loc_horiz, loc_horiz, (loc_horiz + 1), (loc_horiz - 1)))
                d_test <- duplicated(rbind(test_EA_locs,adj_cells))
                d_test <- d_test[(nrow(test_EA_locs)+1):length(d_test)] #FALSES for EA_cell_IDs that are not duplicated in cell_IDs for one side of split
                edge_count <- length(which(d_test == FALSE))
                edge_peri <- sum(edge_peri, edge_count)
              }
              compact_vals <- c(compact_vals, (edge_peri^2)/nrow(test_EA_locs)) 
            }
          }
          compact_vals <- matrix(compact_vals, ncol=2, byrow=TRUE)
          compact_sum <- compact_vals[,1] + compact_vals[,2]
          compact_threshold <- min(compact_sum) * 1.25
          selected_split_ID <- test_splits[which(compact_sum <= compact_threshold)[1]]
          EA_mat[cell_IDs_ls[[selected_split_ID]]] <- EA_IDs_new[1]  # use a brand new EA_ID, taken from beginning of vec of new IDs
          EA_ID_pop <- sum(subpop_mat[which(EA_mat == EA_ID)])
          if(EA_ID_pop < (2 * ss_target_pop_per_EA)){vbig_pop_EAs <- vbig_pop_EAs[-(which(vbig_pop_EAs == EA_ID))]}
          EA_ID_new_pop <- sum(subpop_mat[which(EA_mat == EA_IDs_new[1])])
          if(EA_ID_new_pop >= (2 * ss_target_pop_per_EA)){vbig_pop_EAs <- c(vbig_pop_EAs, EA_IDs_new[1])}
          EA_IDs_new <- EA_IDs_new[-1]; EA_IDs_new <- c(EA_IDs_new, (max(EA_IDs_new) + 1))
          if(!length(vbig_pop_EAs)>0) break
        }
        EA_ID_vec <- unique(EA_mat[which(!is.na(EA_mat))]) 
        EA_pop <- sapply(EA_ID_vec, function(EA_index){sum(subpop_mat[which(EA_mat == EA_index)])})
      }
      big_pop_EAs <- EA_ID_vec[which(EA_pop > max_pop_per_EA)] #remaining big pop EAs (exceeding max pop limit), add to near block, then divide by 3
      cells_count_test <- table(EA_mat[EA_mat %in% big_pop_EAs]) >= 2  #at least 2 cells needed for splitting
      if("TRUE" %in% cells_count_test){
        big_pop_EAs <- as.numeric(names(cells_count_test[cells_count_test == "TRUE"]))
      }
      if(length(big_pop_EAs)>0){
        EA_IDs_new <- (max(EA_ID_vec) + 1):(max(EA_ID_vec) + 3)   
        repeat{
          EA_ID <- big_pop_EAs[1]
          EA_cell_IDs <- which(EA_mat == EA_ID, arr.ind = TRUE)
          #neighbour EAs:
          nb_EAs <- cell_nbs <- NULL
          for(df_row in 1:nrow(EA_cell_IDs)){
            loc_vert <- EA_cell_IDs[df_row,1]
            loc_horiz <- EA_cell_IDs[df_row,2]
            adj_cells <- cbind(c((loc_vert + 1), (loc_vert - 1), loc_vert, loc_vert), c(loc_horiz, loc_horiz, (loc_horiz + 1), (loc_horiz - 1)))
            if(TRUE %in% (!adj_cells[,1] %in% 1:max_sett_row)){adj_cells <- adj_cells[-which((!adj_cells[,1] %in% 1:max_sett_row)),]}
            if(TRUE %in% (!adj_cells[,2] %in% 1:max_sett_col)){adj_cells <- adj_cells[-which((!adj_cells[,2] %in% 1:max_sett_col)),]}
            cell_nbs <- EA_mat[adj_cells]
            nb_EAs <- c(nb_EAs, cell_nbs)
          }
          nb_EAs <- unique(nb_EAs)
          nb_EAs <- nb_EAs[!is.na(nb_EAs) & !nb_EAs == EA_ID]
          if(length(nb_EAs) > 0){         #i.e. add to neighbour only if neighbour exists
            EA_pop_size <- sum(subpop_mat[EA_cell_IDs])
            nb_EA_pop_sizes <- NULL
            for(nb_EA in nb_EAs){
              nb_EA_cell_IDs <- which(EA_mat == nb_EA, arr.ind = TRUE)
              nb_EA_pop_size <- sum(subpop_mat[nb_EA_cell_IDs])
              nb_EA_pop_sizes <- c(nb_EA_pop_sizes, nb_EA_pop_size)
            }
            new_pop_sizes <- (nb_EA_pop_sizes + EA_pop_size)/2
            if(max(new_pop_sizes) < min_pop_per_EA){
                big_pop_EAs <- big_pop_EAs[-which(big_pop_EAs == EA_ID)]
                if(!length(big_pop_EAs)>0){ 
                  break
                }else{
                  next   #don't do anything if splitting will result in new EZs with too small pop 
                }
            }
            #combine with appropriate neighbouring EZ 
            nb_EA_ID <- which(new_pop_sizes >= min_pop_per_EA)
            if(length(nb_EA_ID) > 1){
              compact_vals <- NULL
              nb_EAs <- nb_EAs[nb_EA_ID]
              for(new_EA in nb_EAs){
                test_EA_locs <- which(EA_mat == EA_ID | EA_mat == new_EA ,arr.ind = TRUE)
                edge_peri <- 0
                for(df_row in 1:nrow(test_EA_locs)){
                  loc_vert <- test_EA_locs[df_row,1]
                  loc_horiz <- test_EA_locs[df_row,2]
                  adj_cells <- cbind(c((loc_vert + 1), (loc_vert - 1), loc_vert, loc_vert), c(loc_horiz, loc_horiz, (loc_horiz + 1), (loc_horiz - 1)))
                  d_test <- duplicated(rbind(test_EA_locs,adj_cells))
                  d_test <- d_test[(nrow(test_EA_locs)+1):length(d_test)] #FALSES for EA_cell_IDs that are not duplicated in cell_IDs for one side of split
                  edge_count <- length(which(d_test == FALSE))
                  edge_peri <- sum(edge_peri, edge_count)
                }
                compact_vals <- c(compact_vals, (edge_peri^2)/nrow(test_EA_locs)) 
              }
              nb_EA <- (nb_EAs[compact_vals == min(compact_vals)])[1] 
            }else{
              nb_EA <- nb_EAs[nb_EA_ID]
            }
            #create single big EA to then split into three or two
            EA_mat[EA_mat == nb_EA] <- EA_ID
            if(nb_EA %in% big_pop_EAs){big_pop_EAs <- big_pop_EAs[-which(big_pop_EAs == nb_EA)]}
          }
          EA_cell_IDs <- which(EA_mat == EA_ID, arr.ind = TRUE)
          tot_pop_size <- sum(subpop_mat[EA_cell_IDs])
          if(((tot_pop_size/3) > min_pop_per_EA) & (nrow(EA_cell_IDs) >= 3)){ 
            split_types <- c(1,2)     # want to split twice. When split # is 1 going to split according to pop of 1 EZ, when it's 2 split according to pop diff between the two EZs
          }else{
            split_types <- 2          # want to split once
          }
          #split shape into most compact shapes
          for(split_type in split_types){
            if(split_type == 2 & length(split_types) == 2){     #need to update if doing 2nd split
              EA_cell_IDs <- which(EA_mat == EA_ID, arr.ind = TRUE)
              tot_pop_size <- sum(subpop_mat[EA_cell_IDs])
              if(nrow(EA_cell_IDs) == 1) break
            }
            #create list, an element for each potential split, element contains matrix of cell row & col ids one EA of the potential split, do this because can't assume EAs are only four edged shapes
            row_IDs <- sort(unique(EA_cell_IDs[,1]))
            col_IDs <- sort(unique(EA_cell_IDs[,2]))
            if(length(row_IDs) >= length(col_IDs)){
              side_IDs <- row_IDs
              min_side_IDs <- min(row_IDs)
              row_or_col <- 1
            }else{
              side_IDs <- col_IDs
              min_side_IDs <- min(col_IDs)
              row_or_col <- 2
            }
            cell_IDs_ls <- lapply(side_IDs[1:(length(side_IDs) - 1)], function(x){matrix(EA_cell_IDs[EA_cell_IDs[,row_or_col] %in% (min_side_IDs:x),], ncol=2)})    
            pop_sizes <- sapply(cell_IDs_ls, function(x){sum(subpop_mat[x])})
            if(split_type == 1){
              diff_to_ideal_pop <- abs((tot_pop_size/3) - pop_sizes)
            }
            if(split_type == 2){
              diff_to_ideal_pop <- abs(tot_pop_size - (2*pop_sizes))  # A = B + C, this gives is |C - B| i.e. diff in the two potential EAs' pop sizes
            }
            test_splits <- order(diff_to_ideal_pop)
            if(length(test_splits) >= 10){test_splits <- test_splits[1:10]}
            compact_vals <- NULL
            for(test_EA_split in test_splits){
              for(x in 1:2){
                if(x==1){test_EA_locs <- cell_IDs_ls[[test_EA_split]]}
                if(x==2){
                  d1 <- duplicated(rbind(cell_IDs_ls[[test_EA_split]],EA_cell_IDs))
                  d1 <- d1[(nrow(cell_IDs_ls[[test_EA_split]])+1):length(d1)] #FALSES for EA_cell_IDs that are not duplicated in cell_IDs for one side of split
                  test_EA_locs <- EA_cell_IDs[which(d1 == FALSE),]
                }
                if(!class(test_EA_locs) == "matrix"){test_EA_locs <- matrix(test_EA_locs, nrow=1)}
                edge_peri <- 0
                for(df_row in 1:nrow(test_EA_locs)){
                  loc_vert <- test_EA_locs[df_row,1]
                  loc_horiz <- test_EA_locs[df_row,2]
                  adj_cells <- cbind(c((loc_vert + 1), (loc_vert - 1), loc_vert, loc_vert), c(loc_horiz, loc_horiz, (loc_horiz + 1), (loc_horiz - 1)))
                  d_test <- duplicated(rbind(test_EA_locs,adj_cells))
                  d_test <- d_test[(nrow(test_EA_locs)+1):length(d_test)] #FALSES for EA_cell_IDs that are not duplicated in cell_IDs for one side of split
                  edge_count <- length(which(d_test == FALSE))
                  edge_peri <- sum(edge_peri, edge_count)
                }
                compact_vals <- c(compact_vals, (edge_peri^2)/nrow(test_EA_locs)) 
              }
            }
            compact_vals <- matrix(compact_vals, ncol=2, byrow=TRUE)
            compact_sum <- compact_vals[,1] + compact_vals[,2]
            compact_threshold <- min(compact_sum) * 1.25
            selected_split_ID <- test_splits[which(compact_sum <= compact_threshold)[1]]
            EA_mat[cell_IDs_ls[[selected_split_ID]]] <- EA_IDs_new[1]  # use a brand new EA_ID, taken from beginning of vec of new IDs
            EA_IDs_new <- EA_IDs_new[-1]; EA_IDs_new <- c(EA_IDs_new, (max(EA_IDs_new) + 1))
          }
          big_pop_EAs <- big_pop_EAs[-which(big_pop_EAs == EA_ID)]
          if(!length(big_pop_EAs)>0) break
        } 
        EA_ID_vec <- unique(EA_mat[which(!is.na(EA_mat))]) 
      }
      EA_geog <- sapply(EA_ID_vec, function(EA_index){length(which(EA_mat == EA_index))})
      big_geog_EAs <- EA_ID_vec[which(EA_geog > (1.5 * max_cells_per_EA))]
      if(length(big_geog_EAs)>0){
        EA_IDs_new <- (max(EA_ID_vec) + 1):(max(EA_ID_vec) + length(big_geog_EAs))
        for(EA_ID in big_geog_EAs){  #all EAs with big pops have already been split, so these are low pop EAs with big geog, therefore split by geog rather than pop so that within target range of geog size
          EA_cell_IDs <- which(EA_mat == EA_ID, arr.ind = TRUE)
          row_IDs <- sort(unique(EA_cell_IDs[,1]))
          col_IDs <- sort(unique(EA_cell_IDs[,2]))
          if(length(row_IDs) >= length(col_IDs)){
            side_IDs <- row_IDs
            min_side_IDs <- min(row_IDs)
            row_or_col <- 1
          }else{
            side_IDs <- col_IDs
            min_side_IDs <- min(col_IDs)
            row_or_col <- 2
          }
          cell_IDs_ls <- lapply(side_IDs[1:(length(side_IDs) - 1)], function(x){matrix(EA_cell_IDs[EA_cell_IDs[,row_or_col] %in% (min_side_IDs:x),], ncol=2)})    
          cell_counts <- sapply(cell_IDs_ls, function(x){nrow(x)})
          target_geog <- EA_geog[which(EA_ID_vec == EA_ID)]/2
          selected_split_ID <- which(abs(target_geog - cell_counts) == min(abs(target_geog - cell_counts)))[1]
          EA_mat[cell_IDs_ls[[selected_split_ID]]] <- EA_IDs_new[1]
          EA_IDs_new <- EA_IDs_new[-1]
        }
        EA_ID_vec <- unique(EA_mat[which(!is.na(EA_mat))]) 
      }
      EA_ID_max <- max(EA_ID_vec)
      for(EA_ID in EA_ID_vec){     
        cell_IDs <- which(EA_mat[] == EA_ID, arr.ind = TRUE)
        min_EA_row <- min(cell_IDs[,1])
        max_EA_row <- max(cell_IDs[,1])
        min_EA_col <- min(cell_IDs[,2])
        max_EA_col <- max(cell_IDs[,2])
        subEA_IDs_mat <- EA_mat[min_EA_row:max_EA_row, min_EA_col:max_EA_col]
        subEA_IDs_mat[which(!subEA_IDs_mat[] == EA_ID)] <- NA
        if(class(subEA_IDs_mat) == "matrix"){
          EA_clump_raster <- raster::clump(raster::raster(subEA_IDs_mat), directions = 4, gaps = FALSE)
          clump_ls <- sort(unique(getValues(EA_clump_raster)))
          if(length(clump_ls) > 1){
            new_EA_IDs_vec <- c(EA_ID, seq((EA_ID_max + 1), (EA_ID_max + length(clump_ls) - 1), length.out = (length(clump_ls) - 1)))
            EA_clump_mat <- as.matrix(EA_clump_raster)
            rm(EA_clump_raster)
            for(new_clump in clump_ls[-1]){   #the first clump can keep the original EA_ID so start with updating the second clump
              new_clump_sub_IDs <- which(EA_clump_mat[] == new_clump, arr.ind = TRUE)
              new_clump_IDs <- cbind(new_clump_sub_IDs[,1] + min_EA_row - 1, new_clump_sub_IDs[,2] + min_EA_col - 1)
              EA_mat[new_clump_IDs] <- new_EA_IDs_vec[which(clump_ls == new_clump)]
            }
            EA_ID_max <- max(new_EA_IDs_vec) + 1
            new_EA_pop <- sapply(new_EA_IDs_vec, function(EA_index){sum(subpop_mat[which(EA_mat == EA_index)])})
            new_EA_size <- sapply(new_EA_IDs_vec, function(EA_index){length(which(EA_mat == EA_index))})
            small_new_EAs <- new_EA_IDs_vec[which(new_EA_size < (0.75 * max_cells_per_EA) & (new_EA_pop < min_pop_per_EA))]
            for(small_EA in small_new_EAs){
              # neighbour EAs:
              EA_cell_IDs <- which(EA_mat == small_EA, arr.ind = TRUE)
              nb_EAs <- cell_nbs <- NULL
              for(df_row in 1:nrow(EA_cell_IDs)){
                loc_vert <- EA_cell_IDs[df_row,1]
                loc_horiz <- EA_cell_IDs[df_row,2]
                adj_cells <- cbind(c((loc_vert + 1), (loc_vert - 1), loc_vert, loc_vert), c(loc_horiz, loc_horiz, (loc_horiz + 1), (loc_horiz - 1)))
                if(TRUE %in% (!adj_cells[,1] %in% 1:max_sett_row)){adj_cells <- adj_cells[-which((!adj_cells[,1] %in% 1:max_sett_row)),]}
                if(TRUE %in% (!adj_cells[,2] %in% 1:max_sett_col)){adj_cells <- adj_cells[-which((!adj_cells[,2] %in% 1:max_sett_col)),]}
                cell_nbs <- EA_mat[adj_cells]
                nb_EAs <- c(nb_EAs, cell_nbs)
              }
              nb_EAs <- unique(nb_EAs)
              nb_EAs <- nb_EAs[!is.na(nb_EAs) & !nb_EAs == small_EA]
              if(length(nb_EAs) > 0){
                if(small_EA %in% nb_EAs){nb_EAs <- nb_EAs[-which(nb_EAs == small_EA)]} ###IS THIS NEEDED, LINE ABOVE DOES THIS??? 
                # measure compactness for each of the possible combinations of blocks/EAs
                compact_vals <- NULL
                for(new_EA in nb_EAs){
                  test_EA_locs <- which(EA_mat == small_EA | EA_mat == new_EA ,arr.ind = TRUE)
                  edge_peri <- 0
                  for(df_row in 1:nrow(test_EA_locs)){
                    loc_vert <- test_EA_locs[df_row,1]
                    loc_horiz <- test_EA_locs[df_row,2]
                    adj_cells <- cbind(c((loc_vert + 1), (loc_vert - 1), loc_vert, loc_vert), c(loc_horiz, loc_horiz, (loc_horiz + 1), (loc_horiz - 1)))
                    d_test <- duplicated(rbind(test_EA_locs,adj_cells))
                    d_test <- d_test[(nrow(test_EA_locs)+1):length(d_test)] #FALSES for EA_cell_IDs that are not duplicated in cell_IDs for one side of split
                    edge_count <- length(which(d_test == FALSE))
                    edge_peri <- sum(edge_peri, edge_count)
                  }
                  compact_vals <- c(compact_vals, (edge_peri^2)/nrow(test_EA_locs)) 
                }
                #nb_EAs <- nb_EAs[order(compact_vals)] # ordered by compactness 
                #new_EA_pop <- sapply(nb_EAs, function(EA_index, small_EA_ID = small_EA){sum(subpop_mat[which(EA_mat == small_EA_ID | EA_mat == EA_index)])})
                #new_EA_size <- sapply(nb_EAs, function(EA_index, small_EA_ID = small_EA){length(which(EA_mat == small_EA_ID | EA_mat == EA_index))})
                # could consider pop and geog size but compactness is likely most important when adding on smaller chunks 
                EA_ID_to_add <- (nb_EAs[compact_vals == min(compact_vals)])[1] 
                EA_mat[EA_mat == small_EA] <- EA_ID_to_add
                #there is more code here in the small EA operation section so can't just create function and repeat 
              }
            }
          }
        }
      } 
      EA_ID_vec <- sort(unique(as.vector(EA_mat))) 
      EA_cell_IDs <- lapply(1:length(EA_ID_vec), function(EA){
        subEA_cell_IDs <- which(EA_mat == EA_ID_vec[EA], arr.ind = T)
        subEA_cell_IDs[,1] <- subEA_cell_IDs[,1] + min_clump_row - 1
        subEA_cell_IDs[,2] <- subEA_cell_IDs[,2] + min_clump_col - 1
        subEA_cell_IDs
      })
    } # end of else    
    fileID <- paste(clump_x[[2]][1], clump_x[[2]][2], clump_x[[2]][3], clump_x[[2]][4], sep="")
    save(EA_cell_IDs, file = paste(output_path,"/temp_folder", run_ID, "/EA_cell_IDs_list", fileID, ".RData",sep=""))
  } # end of function
  invisible(parLapply(cl, clump_info, EA_gen_fn, max_cells_per_EA=max_cells_per_EA, output_path = output_path, run_ID = run_ID))
  stopCluster(cl)
  print(paste(Sys.time()))
  print("EZ generation complete, compiling EZ files next")
  current_list <- list.files(path = paste(output_path,"/temp_folder", run_ID, sep=""), pattern =".RData$", full.names=TRUE)
  EA_pop_master_mat <- EA_master_mat <- matrix(NA, nrow = n_rows, ncol = n_cols)
  EA_index <- 0
  for(filenumber in 1:length(current_list)){
    load(current_list[filenumber])
    for(file_listID in 1:length(EA_cell_IDs)){   #in case we change to have multiple clumps per file
      EA_cells_IDs_temp <- EA_cell_IDs[[file_listID]]
      EA_index <- EA_index + 1
      EA_master_mat[EA_cells_IDs_temp] <- EA_index
      EA_pop_master_mat[EA_cells_IDs_temp] <- sum(population_mat[EA_cells_IDs_temp])
    }
    EA_cell_IDs <- NULL
  }  
  print(paste(Sys.time()))
  print("EZ_master_mat completed")
  EA_raster_master <- raster(EA_master_mat,crs=r_crs, xmn=r_xmn, xmx=r_xmx, ymn=r_ymn, ymx=r_ymx)
  EA_pop_raster_master <- raster(EA_pop_master_mat,crs=r_crs, xmn=r_xmn, xmx=r_xmx, ymn=r_ymn, ymx=r_ymx)
  writeRaster(EA_raster_master, paste(output_path,"/EA_raster_master", run_ID, ".tif",sep=""), datatype = 'INT4U', format = "GTiff")
  writeRaster(EA_pop_raster_master, paste(output_path,"/EA_pop_raster_master", run_ID, ".tif",sep=""), format = "GTiff")
  #return(stack(EA_raster_full, EA_pop_raster_full))
  return(NA)
  #unlink(paste(output_path,"/temp_folder", run_ID, sep=""), recursive = TRUE)
}










# a wrapper to the main function
icw_gridEA<-function(population_raster, settlement_raster, strata_raster, 
					output_path,EA_target_size = "medium", 
                     exclude_unsettled = FALSE, unsettled_ID = 0, predefined_EA_size = TRUE, 
                     EA_by_hh = FALSE, target_hh_per_EA = 100, pop_per_hh = 5,
                     EA_by_pop = FALSE, target_pop_per_EA = 500, max_cells_per_EA = 900,
                      run_number = 1) {
  

  
  # do our own preparation here and most important catch any errors			   
  #options(warn=-1)
  
  # remove the output dir and contents (if it exists)
  if(file.exists(output_path)){
    unlink(output_path,recursive=TRUE)
  }
  
  # re/create
  dir.create(output_path)
  
  
  # just catch errors : don't worry about warnings
  ret<-tryCatch(
    {
	
		population_r <- raster(population_raster)
		settlement_r<- raster(settlement_raster)
		strata_r <- raster(strata_raster)

		crs(population_r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
		crs(settlement_r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
		crs(strata_r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
		    
		using_ghs_smod_pop2015 <- TRUE
		EA_raster_stack <- gridEA(population_r,settlement_r,strata_r,exclude_unsettled,unsettled_ID,using_ghs_smod_pop2015 ,predefined_EA_size, EA_target_size,
                                EA_by_hh,target_hh_per_EA , pop_per_hh ,
                                EA_by_pop , target_pop_per_EA , max_cells_per_EA ,output_path, run_number)
      
		return(1)
      
    },
    error=function(e) {
     
      print(e)
      # Choose a return value in case of error
      return(0)
    }
    
  )
  
  #options(warn=0)
  
  return(ret)
  #return(EA_raster_stack)	
}

