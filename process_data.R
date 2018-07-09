#'@param transform Vector of indices for which covariates to log + min transform

process_data <- function(binomial_positive,
                         binomial_n,
                         coords,
                         polygon_response,
                         polygon_population,
                         shapefile_id,
                         shapefiles,
                         shps_id_column = 'area_id',
                         pop_raster,
                         cov_rasters,
                         na.rm = FALSE,
                         transform = NULL,
                         skip_extract = FALSE){
  
  
  stopifnot(inherits(shapefiles, 'SpatialPolygonsDataFrame'))
  stopifnot(inherits(pop_raster, 'RasterLayer'))
  stopifnot(grepl('Raster', class(cov_rasters)))
  # Subset polygon data
  
  polygon <- cbind(response = polygon_response,
                   population = polygon_population,
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
  
  shapefiles@data <- shapefiles@data %>% 
                       left_join(polygon, by = c('area_id' = 'shapefile_id'))
  
  
  # Test polygon data
  message('API data rows: ', nrow(polygon))  
  test_api(polygon, shapefiles)
  
  # Subset PR  
  
  names(coords) <- c('longitude', 'latitude')
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
  
  
  # Crop rasters
  
  # Make pop raster that will be raster template.
  pop_raster <- crop(pop_raster, extent(shapefiles))
  cov_rasters <- crop(cov_rasters, extent(pop_raster))
  
  
  
  # Scale rasters
  cov_rasters <- transform_rasters(cov_rasters, transform)
  cov_rasters <- scale(cov_rasters)
  
  plot_raster_histograms(cov_rasters)
  
  
  # Extract covariates
  extracted <- NULL
  if(!skip_extract){
    extracted <- parallelExtract(stack(pop_raster, cov_rasters), shapefiles, fun = NULL, id = 'area_id')
  }
  
  pr_extracted <- raster::extract(cov_rasters, SpatialPoints(coords))
  
  # Handle covariate NAs.
  
  raster_pts <- rasterToPoints(pop_raster %>% inset(is.na(.), value = -9999), spatial = TRUE)
  
  coords <- raster_pts@coords[extracted$cellid, ]
  covs <- extracted[, -3]
  pop <- extracted[, 3]
  pop[is.na(pop)] <- 0
  
  # Deal with NAs in covariates! TODO
  
  data <- list(pr = pr, 
               polygon = polygon, 
               pop = pop, 
               covs = covs, 
               pr_covs = pr_extracted,
               coords = coords,
               cov_rasters = cov_rasters, 
               shapefiles = shapefiles,
               pop_raster = pop_raster)
  class(data) <- c('ppj_data', 'list')
  
  
  return(data)
}




plot_raster_histograms <- function(covs){
  
  
  sq_sides <- ceiling(sqrt(nlayers(covs)))
  oldpar <- par(mfrow = c(sq_sides, sq_sides))
  on.exit(par(oldpar))
  
  for(i in seq_len(nlayers(covs))){
    hist(na.omit(getValues(covs[[i]])),
         main = names(covs)[i])
  }
  return(NULL)
}

transform_rasters <- function(cov_rasters, transform){
  if(is.null(transform)) return(cov_rasters)
  
  names <- names(cov_rasters)
  
  stopifnot(max(transform) <= nlayers(cov_rasters))
  
  for(i in transform){
    # If all values are above 0, min and min_nonzero get zero and we just log
    # If there are negative values, shift to min of zero then add a small amount...
    min <- minValue(cov_rasters[[i]])
    min <- ifelse(min < 0, min, 0)
    

    cov_rasters[[i]] <- log(cov_rasters[[i]] - min + 1)
  }
  names(cov_rasters) <- names
  return(cov_rasters)
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




