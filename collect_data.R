#################################
# Data Preperation functions    #
# For point vs polygon vs joint #
# Tim Lucas                     #
# 2018-05-30                    #
#################################


#' Function to collect together all data needed for the analysis. 
#'
#'
#'@param PR_path Path to CSV containing parasite rate data.
#'@param API_path Path to CSV containing polygon incidence data.
#'@param  pop_path Path to population raster.
#'@param cov_raster_paths vector of paths to covariate rasters.
#'@param  shapefile_path
#'@param shapefile_pattern
#'@param  api_year
#'@param  pr_year
#'@param useiso3
#'@param  standardisePR
#'@param  roundPR
#'@param   standardisePars
#'@param  api_column 
#'@param  pr_pos_column
#'@param   pr_n_column
#'@param  pr_latlon
#'@param  pr_country
#'@param pr_age_low
#'@param  pr_age_high
#'@param    shapefile_column
#'@param     shps_id_name 
#'@return An object of class 'ppj' with some elements including separate training and testing data. 

load_data <- function(PR_path, 
                      API_path, 
                      pop_path, 
                      cov_raster_paths, 
                      shapefile_path, 
                      shapefile_pattern,
                      api_year, 
                      pr_year,
                      useiso3,
                      standardisePR = c(2, 10),
                      roundPR = FALSE,
                      standardisePars = 'Pf_Smith2007',
                      api_column = 'api_mean_pf', 
                      pop_column = 'pop_at_risk_pf',
                      pr_pos_column = 'positive',
                      pr_n_column = 'examined',
                      pr_latlon = c('latitude', 'longitude'),
                      pr_country = 'country',
                      pr_age_low = 'lower_age',
                      pr_age_high = 'upper_age',
                      shapefile_column = 'shapefile_id',
                      admin_unit_level = 'ADMIN2',
                      shps_id_name = 'area_id'){
  
  check_inputs_load_data(PR_path, API_path, pop_path, cov_raster_paths)
  
  
  
  # Read API data
  admin_unit_level_char <- admin_unit_level 
  api_full <- readr::read_csv(API_path, guess_max  = 1e5)
  api <- api_full %>% filter(year %in% api_year, 
                             iso3 %in% useiso3,
                             admin_unit_level == admin_unit_level_char)
  api <- cbind(api_mean = pull(api, api_column),
               population = pull(api, pop_column),
               shapefile_id = pull(api, shapefile_column)) %>% as.data.frame
  
  api <- api %>% filter(population > 0)
  
  # Read PR data
  pr <- readr::read_csv(PR_path, guess_max  = 1e5)
  
  usecountries <- find_country_from_iso3(useiso3, api_full$iso3, api_full$country_name)
  pr <- pr %>% filter(country %in% usecountries, year_start %in% pr_year)
  
  pr_clean <- data_frame(
    prevalence = pull(pr, pr_pos_column) / pull(pr, pr_n_column),
    positive = pr %>% pull(pr_pos_column),
    examined = pr %>% pull(pr_n_column),
    latitude = pr %>% pull(pr_latlon[1]),
    longitude =  pr %>% pull(pr_latlon[2])
  )
  
  if(!is.null(standardisePR)){
    prev_stand <- convertPrevalence(pr_clean$prevalence, 
                                    pr %>% pull(pr_age_low),
                                    pr %>% pull(pr_age_high),
                                    standardisePR[1],
                                    standardisePR[2],
                                    parameters = standardisePars
    )
    if(roundPR == TRUE){
      pr_clean$positive <- round(pr_clean$examined * prev_stand)
    } else {
      pr_clean$positive <- pr_clean$examined * prev_stand
    }
    
    
  }  
  
  # Read pop raster
  pop <- raster::raster(pop_path)
  
  
  # Read covariate rasters
  covs_list <- lapply(cov_raster_paths, raster::raster)
  crop_to <- find_smallest_extent(c(covs_list, pop))
  covs_cropped <- lapply(covs_list, function(x) crop(x, crop_to))
  covs <- do.call(stack, CombineRasters(covs_cropped))
  
  pop <- crop(pop, crop_to)
  
  
  # Read in shapefiles
  shp_filepaths <- list.files(shapefile_path, pattern = shapefile_pattern, full.names = TRUE)
  stopifnot(length(shp_filepaths) > 0)
  
  shps <- lapply(shp_filepaths, raster::shapefile)
  shp_names <- shps %>% lapply(names) %>% do.call(c, .) %>% unique
  common_names <- shp_names[sapply(shp_names, function(n) all(sapply(shps, function(x) n %in% names(x))))]
  
  shps <- lapply(shps, function(x) x[, common_names]) %>% do.call(rbind, .)
  
  # Combine data.
  
  data <- list(pr = pr_clean, api = api, pop = pop, covs = covs, shapefiles = shps)
  class(data) <- c('ppj_full_data', 'list')
  
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




find_country_from_iso3 <- function(useiso3, iso3, country){
  
  
  # Need to find a country name because PR data doesn't have iso3 codes... 
  usecountries <- country[match(useiso3, iso3)]
  
  if(is.na(usecountries)) stop('No matching iso3')
  # if(!all(usecountries %in% pr$country)) stop('At least one country not in PR data. Probably a name mismatch')
  return(usecountries)
}





