

process_data <- function(binomial_positive,
                         binomial_n,
                         coords,
                         response,
                         response_sd,
                         shapefile_id,
                         shapefiles,
                         shps_id_column = 'area_id',
                         pop_raster,
                         cov_rasters,
                         na.rm = FALSE){

  
  stopifnot(inherits(shapefiles, 'SpatialPolygonsDataFrame'))
  stopifnot(inherits(pop_raster, 'RasterLayer'))
  stopifnot(grepl('Raster', class(cov_rasters)))
  # Subset polygon data

  polygon <- cbind(response = response,
                   shapefile_id = shapefile_id) %>% as.data.frame

  if(na.rm == TRUE){ 
    rows <- nrow(polygon)
    polygon <- na.omit(polygon)
    newrows <- nrow(polygon)
    if(rows > newrows) message(rows - newrows, 'rows with NAs omitted from polygon data')
  }




  shapefiles@data <- shapefiles@data[, shps_id_column, drop = FALSE]


  # Subset shapefiles
  shapefiles <- shapefiles[shapefiles@data[, shps_id_column] %in% polygon$shapefile_id, ]

  if(na.rm == TRUE){
    polygon <- polygon %>% filter(shapefile_id %in% shapefiles@data[, shps_id_column])
    finalrows <- nrow(polygon)
    if(finalrows < newrows) message(newrows - finalrows, 'rows with missing polygons omitted from polygon data')
  }

  # Test polygon data
  message('API data rows: ', nrow(polygon))  
  test_api(polygon, shapefiles)

  # Subset PR  
    
  names(coords) <- c('latitude', 'longitude')
  pr <- cbind(positive = binomial_positive,
              examined = binomial_n,
              coords)

   if(na.rm == TRUE){ 
      rows <- nrow(pr)
      pr <- na.omit(pr)
      newrows <- nrow(pr)
      if(rows > newrows) message(rows-newrows, ' rows with NAs omitted from pr data')
   }
                                  
  message('PR data rows: ', nrow(pr))

  # Test PR data.

  test_pr(pr)




  # Extract covariates
  extracted <- parallelExtract(stack(pop_raster, cov_rasters), shapefiles, fun = NULL, id = 'area_id')

  raster_pts <- rasterToPoints(pop_raster %>% inset(is.na(.), value = -9999), spatial = TRUE)
  extracted <- cbind(extracted, raster_pts@coords[extracted$cellid, ])

  covs <- extracted[, -3]
  pop <- extracted[, 3]
  pop[is.na(pop)] <- 0



  data <- list(pr = pr, 
               polygon = polygon, 
               pop = pop, 
               covs = covs, 
               cov_rasters = cov_rasters, 
               shapefiles = shapefiles,
               pop_raster = pop_raster)
  class(data) <- c('ppj_data', 'list')


  return(data)
}





test_api <- function(polygon, shapefiles){

  stopifnot(all(polygon$response >= 0))
  stopifnot(!any(is.na(polygon$response)))
  stopifnot(!any(is.na(polygon$shapefile_id)))
  stopifnot(!any(!polygon$shapefile_id %in% shapefiles$area_id))

  return(NULL)
}


test_pr <- function(pr){

stopifnot(all(pr$examined >= 0))
  stopifnot(!any(is.na(pr$examined)))
  stopifnot(!any(pr$positive > pr$examined))

  return(NULL)
}




