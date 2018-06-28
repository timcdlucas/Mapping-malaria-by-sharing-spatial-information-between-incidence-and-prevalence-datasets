# plotting functions.





# List of classes and useful plotting functions.

# ppj_data
#  Plot polygon data
#  Plot combined data.

plot.ppj_data <- function(x, y = NULL, type = 'both', trans = 'identity', ...){
  
  df <- ggplot2::fortify(x$shapefiles, region = 'area_id')
  
  df <- x$shapefiles@data %>% 
          mutate(area_id = as.character(area_id)) %>% 
          left_join(df, ., by = c('id' = 'area_id'))
  
  
  p <- ggplot(df, aes(long, lat, group = group, fill = response)) + 
         geom_polygon() +
         coord_equal() +
         scale_fill_viridis_c(trans = trans)
  
  
  
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

