# plotting functions.





# List of classes and useful plotting functions.

# ppj_data
#  Plot polygon data
#  Plot combined data.

autoplot.ppj_data <- function(object, type = 'both', trans = 'identity', limits = c(NA, NA), ...){
  
  df <- ggplot2::fortify(x$shapefiles, region = 'area_id')
  
  df <- x$shapefiles@data %>% 
          mutate(area_id = as.character(area_id)) %>% 
          left_join(df, ., by = c('id' = 'area_id'))
  
  
  p <- ggplot(df, aes(long, lat, group = group, fill = response)) + 
         geom_polygon() +
         coord_equal() +
         scale_fill_viridis_c(trans = trans, limits = limits, oob = scales::squish)
  
  
  
  if(type == 'both'){
    p <- p + geom_point(data = x$pr, 
                        aes(longitude, latitude, fill = positive / examined, group = NULL), 
                        shape = 21, colour = 'lightgrey', alpha = 0.8) 
  }
  print(p)
  p
}

# ppj_cv
# Plot polygon and pr groups

autoplot.ppj_cv <- function(object, jitter = 0.5, ...){
  
  test_poly <- lapply(seq_along(object), 
                        function(x) cbind(object[[x]]$test$shapefiles@data, fold = x)) %>% 
                   do.call(rbind, .)
  
  test_shapes <- lapply(seq_along(object), 
                        function(x) object[[x]]$test$shapefiles) %>% 
    do.call(rbind, .)
  
  test_pr <- lapply(seq_along(object), 
                        function(x) cbind(object[[x]]$test$pr, fold = x)) %>% 
    do.call(rbind, .)
  
  
  #test_shapes@data <- test_poly
  
  df <- ggplot2::fortify(test_shapes, region = 'area_id')
  
  df <- test_poly %>% 
    mutate(area_id = as.character(area_id)) %>% 
    left_join(df, ., by = c('id' = 'area_id'))
  
  
  p <- ggplot(df, aes(long, lat, group = group, fill = factor(fold))) + 
    geom_polygon() +
    geom_jitter(data = test_pr, aes(longitude, latitude, group = NULL), 
                pch = 21,
                height = jitter, width = jitter) +
    coord_equal() + 
    scale_fill_brewer(palette = 'Set3')
  
  
  
  print(p)
}




# ppj_model 
#   plot parameters
#   plot predictions - chose which.


# ppj_cv_performance
#   Plot obs vs pred for both data
#   Plot performance metrics against each other?


# ppf_cv_results
#   Plot small multiples of predictions - choose which.
#   Plot summary of predictions - choose which, choose functions.
#   Plot obs vs pred for both data
#   Plot map of predictions (assuming no spatial overlap?) and input data.


##################################################
# List of classes and summary of useful print?  ##
##################################################

# ppj_model 
#   print parameters


# ppj_cv_performance
#  Print useful summaries.
#  Summary function might be useful.

