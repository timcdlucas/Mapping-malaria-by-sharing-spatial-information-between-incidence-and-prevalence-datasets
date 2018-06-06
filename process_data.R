


process_data <- function(data, useiso3){

  stopifnot(inherits(data, 'ppj_full_data'))
  stopifnot(inherits(useiso3, 'character'))

  # Subset PR  

  

  # Test PR data.


  # Subset polygon data


  # Test polygon data
  
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
