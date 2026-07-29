xfun::pkg_attach2(c("bipartite", "plyr", "dplyr", "gawdis", "ggplot2", "iNEXT.3D", 
                    "megatrees", "parallel", "rtrees", "psych",
                    "sf", "sp", "terra", "tidyr", "tidyverse", "traitdata"))

#' get which cell a point belongs to (from Daijiang Li)
#' @param coord_df [data.frame] with coordinates for points, with long and lat as columns
#' @param ras_file [raster]to query with
#' @param crs_proj [character] coordinate system of projection (since obtaining it from rater gives problems to sp)
#' @note requires packages plyr, dplyr, knitr, raster, sp, terra, tibble, tidyverse
get_nth_cell_xy = function(coord_df, ras_file, crs_proj = NULL){
  require(sp)
  require(raster)
  require(plyr)
  require(dplyr)
  require(knitr)
  require(terra)
  require(tibble)
  require(tidyverse)
  
  # convert long/lat to SpatialPoint
  mypoints_sp = SpatialPoints(coord_df, CRS("+init=EPSG:4326"))
  mypoints_sp_prj = spTransform(mypoints_sp, crs(ras_file)) # convert to eck4
  coord_df$nth_cell = cellFromXY(ras_file, coordinates(mypoints_sp_prj))
  
  # another way to get the number of points within a cell
  # plot(rasterize(mypoints_sp_eck4, ref_grid_usa2, fun='count')) 
  coord_df = coord_df %>% 
    bind_cols(as_tibble(coordinates(mypoints_sp_prj)) %>% 
                setNames(c("long_prj", "lat_prj"))) %>% 
    as_tibble()
  
  coord_df = bind_cols(coord_df, 
                       as_tibble(crds(ras_file, na.rm = F)[coord_df$nth_cell, ]) %>% 
                         setNames(c("cell_long", "cell_lat")))
  
  if(!is.null(crs_proj)){
    centercells_proj = SpatialPoints(coord_df[!is.na(coord_df$nth_cell),c("cell_long", "cell_lat")], 
                                     CRS(crs_proj))
    centercells_sp = spTransform(centercells_proj, CRS("+init=EPSG:4326"))
    
    centercells_df <- bind_cols(coord_df[!is.na(coord_df$nth_cell),c("nth_cell")],
                                as_tibble(coordinates(centercells_sp))|>
                                  setNames(c("cell_decimalLongitude", "cell_decimalLatitude")))
    
    centercells_df <-  centercells_df[!duplicated(centercells_df),]
    
    coord_df = coord_df|>
      left_join(centercells_df, by = "nth_cell")
  }
  
  # # get cells with NAs in raster file
  # na_cells = which(is.na(values(ras_file)))
  # # set locations from na_cells to NA
  # coord_df$nth_cell[coord_df$nth_cell %in% na_cells] = NA
  # # remove NAs
  # coord_df = filter(coord_df, !is.na(nth_cell))
  # # add population data to df
  # coord_df$pop_km2 = values(ras_file)[coord_df$nth_cell]
  coord_df
}

####Functions to obtain WorldClim 2.1 data####
#' read WorldClim data in R (modified to get data from v.2.1 from Robert J. Hijmans)
#' @param var [character] variable to download: valid variables names are 'tavg', 'tmean', 'tmin', 'tmax', 'prec', 'bio', 'alt', 'elev', 'srad', 'wind' and 'vapr'.
#' @param res [numeric] resolution of the data: valid resolutions are 0.5, 30, 2.5, 5, 10 (note that all are in min arc but 30 that refers to secons and is equal to 0.5)
#' @param path [character] path to store the data
#' @param download [logical] download the data if it's not available?
#' @note requires packages utils and raster, and function download_wc
get_wc2.1 <- function (var, res, path, download = TRUE) {
  require(raster)
  
  #Cleaning arguments
  if(res == 30) {
    res = 0.5
  } #Substitute 30s to 0.5m for the rest of the function
  if (!res %in% c(0.5, 2.5, 5, 10)) {
    stop("resolution should be one of: 0.5, 30, 2.5, 5, 10")
  }
  
  if(var == "tmean"){
    var = "tavg"
  } #Substitute tmean to tavg for the rest of the function
  if(var == "alt"){
    var = "elev"
  } #Substitute alt to elev for the rest of the function
  stopifnot(var %in% c("tavg", "tmin", "tmax", "prec", "bio", 
                       "elev", "srad", "wind", "vapr"))
  
  path <- paste(path, "wc", res, "/", sep = "")
  dir.create(path, showWarnings = FALSE)
  
  #Path Zipfile
  if(res != 0.5){
    zip <- paste("wc2.1_", res, "m_", var, ".zip", sep = "")
  }else{
    zip <- paste("wc2.1_30s_", var, ".zip", sep = "")
  }
  
  zipfile <- paste(path, zip, sep = "")
  
  #Path individual files within zipfile
  if (var == "elev") {
    if(res != 0.5){
      files <- paste("wc2.1_", res, "m_", var, ".tif", sep = "")
    }else{
      files <- paste("wc2.1_30s_", var, ".tif", sep = "")
    }
  }
  else if (var != "bio") {
    if(res != 0.5){
      files <- paste("wc2.1_", res, "m_", var, "_", sprintf("%02d", 1:12), ".tif", sep = "")
    }else{
      files <- paste("wc2.1_30s_", var, "_", sprintf("%02d", 1:12), ".tif", sep = "")
    }
  }
  else {
    if(res != 0.5){
      files <- paste("wc2.1_", res, "m_", var, "_", 1:19, ".tif", sep = "")
    }else{
      files <- paste("wc2.1_30s_", res, "m_", var, "_", 1:19, ".tif", sep = "")
    }
  }
  theurl <- paste("https://geodata.ucdavis.edu/climate/worldclim/2_1/base/", 
                  zip, sep = "")
  
  #Get file
  #Check whether file is downloaded in the right path, otherwise download if allowed to
  fc <- sum(file.exists(files))
  
  if (fc < length(files)) {
    if (!file.exists(zipfile)) {
      if (download) {
        download_wc(theurl, zipfile)
        if (!file.exists(zipfile)) {
          message("\n Could not download file -- perhaps it does not exist")
        }
      }
      else {
        message("File not available locally. Use 'download = TRUE'")
      }
    }
    utils::unzip(zipfile, exdir = dirname(zipfile))
  }
  
  #Read file in raster format
  st <- rast(paste(path, files, sep = ""))
  return(st)
}

#' Support function to download WordClim data (modified to get data from v.2.1 from Robert J. Hijmans)
#' @param aurl [character] URL to download from
#' @param filename [character] Cfile name to store to
#' @note requires packages utils
download_wc<- function (aurl, filename) {
  fn <- paste(tempfile(), ".download", sep = "")
  res <- utils::download.file(url = aurl, destfile = fn, quiet = FALSE, 
                              mode = "wb", cacheOK = TRUE)
  if (res == 0) {
    w <- getOption("warn")
    on.exit(options(warn = w))
    options(warn = -1)
    if (!file.rename(fn, filename)) {
      file.copy(fn, filename)
      file.remove(fn)
    }
  }
  else {
    stop("could not download the file")
  }
}

