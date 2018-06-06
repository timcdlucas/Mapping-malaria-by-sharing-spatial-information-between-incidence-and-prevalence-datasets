
#'@param data list of class ppj_full_data created with load_data.
#'@param Character vector of iso3s to use.
#'@param year vector of years to uuse.

process_data <- function(data, useiso3, year, api_column = 'api_mean_pf'){

  stopifnot(inherits(data, 'ppj_full_data'))
  stopifnot(inherits(useiso3, 'character'))

  # Subset polygon data
  api <- data$api
  api <- api %>% 
           filter(year %in% year, 
                  iso3 %in% useiso3)

  api <- api[!is.na(api[, api_column]), ]

  message('API data rows: ', nrow(api))

  # Test polygon data
  
  test_api(api)

  # Subset PR  

  pr <- data$pr

  # Need to find a country name because PR data doesn't have iso3 codes... 
  usecountries <- data$api$country_name[match(useiso3, data$api$iso3)]

  if(is.na(usecountries)) stop('No matching iso3')
  if(!all(usecountries %in% pr$country)) stop('At least one country not in PR data. Probably a name mismatch')

  pr <- pr %>% 
          filter(year %in% year, 
                 country %in% usecountries)
  
  message('PR data rows: ', nrow(pr))

  # Test PR data.






  
  # Subset shapefiles


  # Extract covariates

  
  # Test shape data


  # Make polygon rasters


  # Crop rasters


  data <- list(pr = pr, 
               api = api, 
               pop = pop, 
               covs = covs, 
               rasterised_shapes = rasters_shapes
               cov_rasters = cov_rasters, 
               pop_raster = pop_raster)
  class(data) <- c('ppj_data', 'list)


  return(processed_data)
}





test_api <- function(api, api_column){

  stopifnot(all(api[, api_column] >= 0))
  

  return(NULL)
}













