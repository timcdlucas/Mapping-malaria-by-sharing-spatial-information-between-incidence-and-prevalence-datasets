

#'@return a 'ppj_cv' object. Which is a list of paired train/test 'ppf_data' objects
#'

cv_spatial_folds <- function(data, k = 5, keep_pr = FALSE){
  
  # Polygons to points
  # Sort shapefiles to match order of polygon dataframe
  
  sorted_shapes <- data$shapefiles[match(data$polygon$shapefile_id, data$shapefiles$area_id), ]
  
  centroids <- sapply(seq_len(nrow(sorted_shapes)), 
                      function(x) rgeos::gCentroid(sorted_shapes[x, ])@coords)
  
  centroids <- t(centroids)
  colnames(centroids) <- c('longitude', 'latitude')
  
  all_coords <- centroids
  
  # Do k means
  folds <- kmeans(all_coords[, 1:2], k, algorithm = 'MacQueen', iter.max = 1000)
  folds <- folds$cluster
  
  
  # fold is length nrow(polygons), foldspr is length nrow(points)
  #   Should be values 1 to k.
  polygon_folds <- folds[seq_len(nrow(data$polygon))]


  # Find which pr are in which polygons.
  pr_spatial_points <- SpatialPoints(data$pr[, c('longitude', 'latitude')], CRS(projection(sorted_shapes)))

  d <- over(pr_spatial_points, sorted_shapes)

  foldspr <- polygon_folds[match(d$area_id, data$polygon$shapefile_id)]
  foldspr[is.na(foldspr)] <- 0
  
  data_cv <- list()
  class(data_cv) <- c('ppj_cv', 'list')
  
  for(i in 1:k){
    data_cv[[i]] <- subset_data_cv_spat(data, polygon_folds, foldspr, i, keep_pr)
  }
  
  return(data_cv)
  
}


subset_data_cv_spat <- function(data, polygon_folds, pr_folds, k, keep_pr){
  
  data_fold <- list(train = list(), test = list())
  
  if(keep_pr){
    data_fold$train$pr <- data$pr
    data_fold$train$pr_covs <- data$pr_covs
  } else { 
    data_fold$train$pr <- data$pr[pr_folds != k, ]
    data_fold$train$pr_covs <- data$pr_covs[pr_folds != k, ]
  }
    
  data_fold$train$polygon <- data$polygon[polygon_folds != k, ]
  cov_rows <- data$covs$area_id %in% data_fold$train$polygon$shapefile_id
  data_fold$train$covs <- data$covs[cov_rows, ]
  data_fold$train$pop <- data$pop[cov_rows]
  data_fold$train$coords <- data$coords[cov_rows, ]
  data_fold$train$cov_rasters <- data$cov_rasters
  data_fold$train$pop_raster <- data$pop_raster
  data_fold$train$shapefiles <- data$shapefiles[data$shapefiles$area_id %in% data_fold$train$covs$area_id, ]
  data_fold$train$iid_raster <- data$iid_raster  
  data_fold$train$shapefile_raster <- data$shapefile_raster  
  
  data_fold$test$pr <- data$pr
  data_fold$test$pr_covs <- data$pr_covs
  
  data_fold$test$polygon <- data$polygon[polygon_folds == k, ]
  data_fold$test$covs <- data$covs[!cov_rows, ]
  data_fold$test$pop <- data$pop[!cov_rows]
  data_fold$test$coords <- data$coords[!cov_rows, ]
  data_fold$test$cov_rasters <- data$cov_rasters
  data_fold$test$pop_raster <- data$pop_raster
  data_fold$test$shapefiles <- data$shapefiles[!(data$shapefiles$area_id %in% data_fold$train$covs$area_id), ]
  data_fold$test$iid_raster <- data$iid_raster  
  data_fold$test$shapefile_raster <- data$shapefile_raster  

  class(data_fold$test) <- 'ppj_data'
  class(data_fold$train) <- 'ppj_data'
  
  
  return(data_fold)
  
}



