#################################
# Data Preperation functions    #
# For point vs polygon vs joint #
# Tim Lucas                     #
# 2018-05-30                    #
#################################


#' Function to collect together all data needed for the analysis. Split by country comes later.
#'
#' The paths required are:
#'   PR_path
#'   API_path
#'   pop_path
#'
#'
#'@param paths named character vector of data paths
#'@param cov_raster_paths A charac vector of paths to covariate tifs.
#'@return An object of class 'ppj' with some elements including separate training and testing data. 

load_data <- function(PR_path, API_path, pop_path, cov_raster_paths){

  check_inputs_load_data(PR_path, API_path, pop_path, cov_raster_paths)

  # Read PR data
  pr <- readr::read_csv(PR_path)

  
  # Read API data
  api <- readr::read_csv(API_path)


  # Read pop raster
  pop <- raster::raster(pop_path)


  # Read covariate rasters
  covs_list <- lapply(cov_raster_paths, raster::raster)
  crop_to <- find_smallest_extent(c(covs_list, pop))
  covs_cropped <- lapply(covs_list, function(x) crop(x, crop_to))
  covs <- do.call(stack, CombineRasters(covs_cropped))

  pop <- crop(pop, crop_to)

  data <- list(pr = pr, api = api, pop = pop, covs = covs)
  class(data) <- c('ppj_full_data', 'list)

  return(data)

}




# Factor out all args checks
check_inputs_load_data <- function(PR_path, API_path, pop_path, cov_raster_paths){

  stopifnot(inherits(PR_path, 'character'))
  stopifnot(inherits(API_path, 'character'))
  stopifnot(inherits(pop_path, 'character'))

  stopifnot(all(file.exists(PR_path, API_path, pop_path)))

  stopifnot(grepl('.tif$', pop_path))

  stopifnot(inherits(cov_raster_paths, 'character'))
  stopifnot(all(grepl('.tif$', cov_raster_paths)))
  stopifnot(all(do.call(file.exists, as.list(cov_raster_paths))))

  return(NULL)

}



find_smallest_extent <- function(raster_list){
  

  # Extract all extents
  extent_mat <- do.call(rbind, lapply(raster_list, function(x) as.vector(extent(x))))

  # Find smallest.
  extent_vec <- c(max(extent_mat[, 1]), min(extent_mat[, 2]), max(extent_mat[, 3]), min(extent_mat[, 4]))
  
  crop_to <- extent(extent_vec)

  return(crop_to)
}



