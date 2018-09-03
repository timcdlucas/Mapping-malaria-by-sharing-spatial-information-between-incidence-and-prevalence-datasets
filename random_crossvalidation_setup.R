

#'@return a 'ppj_cv' object. Which is a list of paired train/test 'ppf_data' objects
#'

cv_random_folds <- function(data, k = 5){
  
  # To make sure theres's nearly exactly n / k in each fold
  f <- rep(1:k, ceiling(nrow(data$polygon) / k))
  s <- sample(1:nrow(data$polygon))
  folds <- f[s]
  

  # Same for PR
  fpr <- rep(1:k, ceiling(nrow(data$pr) / k))
  spr <- sample(1:nrow(data$pr))
  foldspr <- fpr[spr]
  
  data_cv <- list()
  class(data_cv) <- c('ppj_cv', 'list')
  
  for(i in 1:k){
    data_cv[[i]] <- subset_data_cv(data, folds, foldspr, i)
  }
  
  any_empty <- sapply(seq_along(data_cv), function(x) data_cv[[x]]$test$polygon %>% nrow) %>% 
                 equals(0) %>% any
  any_empty_pr <- sapply(seq_along(data_cv), function(x) data_cv[[x]]$test$pr %>% nrow) %>% 
    equals(0) %>% any
  

  if(any_empty) stop('Some test sets have zero polygons')
  if(any_empty_pr) stop('Some test sets have zero pr points')
  
  return(data_cv)
  
}



subset_data_cv <- function(data, polygon_folds, pr_folds, k){
  
  data_fold <- list(train = list(), test = list())
  
  data_fold$train$pr <- data$pr[pr_folds != k, ]
  data_fold$train$pr_covs <- data$pr_covs[pr_folds != k, ]
    
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
  
  data_fold$test$pr <- data$pr[pr_folds == k, ]
  data_fold$test$pr_covs <- data$pr_covs[pr_folds == k, ] # this shouldn't ever be used as far as I can see.
  
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







