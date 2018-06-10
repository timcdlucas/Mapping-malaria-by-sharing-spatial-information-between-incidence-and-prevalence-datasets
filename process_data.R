

process_data <- function(
                         binomial_positive,
                         binomial_n,
                         coords,
                         response,
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
  extracted <- parallelExtract(stack(pop_raster, cov_rasters), shape files, fun = NULL)

  covs <- extracted[, -3]
  pop <- extracted[, 3]


  # Make polygon rasters



  data <- list(pr = pr, 
               api = api, 
               pop = pop, 
               covs = covs, 
               cov_rasters = cov_rasters, 
               pop_raster = pop_raster)
  class(data) <- c('ppj_data', 'list)


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





find_iso3_from_country <- function(useiso3, iso3, country){


  # Need to find a country name because PR data doesn't have iso3 codes... 
  usecountries <- country[match(useiso3, iso3)]

  if(is.na(usecountries)) stop('No matching iso3')
  if(!all(usecountries %in% pr$country)) stop('At least one country not in PR data. Probably a name mismatch')
return(usecountries)
}



