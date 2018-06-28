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
                height = jitter, width = jitter, show.legend = FALSE) +
    coord_equal() + 
    scale_fill_brewer(palette = 'Set3')
  
  print(p)
}




# ppj_model 
#   plot parameters
#   plot predictions - chose which.

plot.ppj_model <- function(x, layers = 'all', ...){
  
  if(layers == 'all') layers <- names(x$predictions)
  plot(stack(x$predictions[layers]), ...)
  invisible(NULL)
}


autoplot.ppj_model <- function(object, skip_node_mean = TRUE, ...){
  pars <- object$model[[1]]$par
  if(skip_node_mean)  pars <- pars[names(pars) != 'nodemean']
  
  pars.df <- data_frame(value = pars, parameter_group = names(pars))
  pars.df$parameter <- make.unique(pars.df$parameter_group, sep = '_')
  
  if(!skip_node_mean) pars.df$parameter[pars.df$parameter_group == 'nodemean'] <- ''
  
  p <- ggplot() + 
         geom_point(data = pars.df %>% dplyr::filter(parameter_group != 'nodemean'), 
                    aes(x = value, y = parameter, colour = parameter_group)) +
         geom_jitter(data = pars.df %>% dplyr::filter(parameter_group == 'nodemean'), 
               aes(x = value, y = parameter, colour = parameter_group),
               alpha = 0.6, width = 0) +
         facet_wrap(~ parameter_group, scale = 'free') + 
         scale_fill_brewer(palette = 'Set3')
  print(p)
}


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

