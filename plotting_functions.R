# plotting functions.





# List of classes and useful plotting functions.

# ppj_data
#  Plot polygon data
#  Plot combined data.

autoplot.ppj_data <- function(object, type = 'both', trans = 'identity', limits = c(NA, NA), 
                              pr_limits = c(NA, NA), ...){
  
  df <- ggplot2::fortify(object$shapefiles, region = 'area_id')
  
  df <- object$shapefiles@data %>% 
    mutate(area_id = as.character(area_id)) %>% 
    left_join(df, ., by = c('id' = 'area_id'))
  
  
  p <- ggplot(df, aes(long, lat, group = group, fill = response)) + 
    geom_polygon() +
    coord_equal() +
    scale_fill_viridis_c(trans = trans, limits = limits, oob = scales::squish)
  
  
  
  if(type == 'both'){
    p <- p + 
      geom_point(data = object$pr, 
                 aes(longitude, latitude, colour = positive / examined, 
                     fill = NULL, group = NULL), 
                 alpha = 0.8) +
      scale_colour_viridis_c(trans = trans, limits = pr_limits, oob = scales::squish)
    
  }
  print(p)
  return(p)
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
                height = jitter, width = jitter, 
                show.legend = FALSE,
                ...) +
    coord_equal() + 
    scale_fill_brewer(palette = 'Set3')
  
  print(p)
  return(p)
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
  return(p)
}


# ppj_cv_performance
#   Plot obs vs pred for both data
#   Plot performance metrics against each other?

autoplot.ppj_cv_performance <- function(object, trans = 'identity', CI = FALSE, ...){
  
  pr <- object$pr_pred_obs
  polygon <- object$polygon_pred_obs
  
  polygon_clean <- data_frame(data_type = 'polygon', 
                              predicted = polygon$pred_api,
                              observed = polygon$response,
                              lower = polygon$pred_api_lower,
                              upper = polygon$pred_api_upper)
  
  pr_clean = data_frame(data_type = 'point',
                        predicted = pr$pred_prev, 
                        observed = pr$prevalence,
                        lower = pr$prevalence_lower,
                        upper = pr$prevalence_upper)
  
  d <- rbind(polygon_clean, pr_clean)
  
  p <- ggplot(d, aes(observed, predicted)) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_point() + 
    geom_smooth(method = 'lm', alpha = 0.15, size = 0.8, colour = 'steelblue') +
    facet_wrap(~ data_type, scale = 'free') + 
    scale_x_continuous(trans = trans) +
    scale_y_continuous(trans = trans)
  
  if(CI) p <- p + geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.3)
  
  print(p)
  return(p)
}

# ppf_cv_results
#   Plot small multiples of predictions - choose which.
#   Plot summary of predictions - choose which, choose functions.
#   Plot obs vs pred for both data
#   Plot map of predictions (assuming no spatial overlap?) and input data.


autoplot.ppf_cv_results <- function(object, 
                                    type = 'layers', 
                                    layer = 'api', 
                                    fun = range,
                                    trans = 'identity',
                                    CI = TRUE,
                                    ...){
  
  if(type == 'layers'){
    p <- ppf_cv_results_layers(object, layer)
  } else if (type == 'raster_summary') {
    p <- ppf_cv_results_raster_summary(object, layer, fun)
  } else if (type == 'obs_preds') {
    p <- ppf_cv_results_obspreds(object, trans, CI)
  } else {
    stop('Availale plot types are layers, raster_summary and obs_preds')
  }
  return(p)
}

ppf_cv_results_obspreds <- function(object, trans, CI = TRUE){
  
  pr <- object$summary$combined_pr
  polygon <- object$summary$combined_aggregated
  
  polygon_clean <- data_frame(data_type = 'polygon', 
                              predicted = polygon$pred_api,
                              observed = polygon$response,
                              upper = polygon$pred_api_upper,
                              lower = polygon$pred_api_lower)
  
  pr_clean = data_frame(data_type = 'point',
                        predicted = pr$pred_prev, 
                        observed = pr$prevalence,
                        upper = pr$prevalence_upper,
                        lower = pr$prevalence_lower)
  
  d <- rbind(polygon_clean, pr_clean)
  
  p <- ggplot(d, aes(observed, predicted)) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_point() + 
    geom_smooth(method = 'lm', alpha = 0.15, size = 0.8, colour = 'steelblue') +
    facet_wrap(~ data_type, scale = 'free') + 
    scale_x_continuous(trans = trans) +
    scale_y_continuous(trans = trans)
  
  if(CI) p <- p + geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.3)
  
  print(p)
  return(p)
}

ppf_cv_results_layers <- function(object, layer, fun){
  stopifnot(layer %in% names(object$models[[1]]$predictions))
  
  r <- lapply(object$models, function(x) x$predictions[[layer]]) %>% do.call(stack, .)
  plot(r)
}

ppf_cv_results_raster_summary <- function(object, layer, fun){
  
  stopifnot(layer %in% names(object$models[[1]]$predictions))
  
  r <- lapply(object$models, function(x) x$predictions[[layer]]) %>% do.call(stack, .)
  summary <- calc(r, fun)
  plot(summary)
  return(summary)
}





###############################################
# Plot pred obs map                           #
###############################################


obspred_map <- function(cv_data, 
                        cv_preds, 
                        layer = 'api', 
                        trans = 'identity', 
                        lims = NULL, 
                        column = TRUE,
                        legend_title = 'response',
                        mask = FALSE,
                        ...){
  
  stopifnot(inherits(cv_data, 'ppj_cv'))
  stopifnot(inherits(cv_preds, 'ppf_cv_results'))
  
  
  # Create polygon data map
  
  test_poly <- lapply(seq_along(cv_data), 
                      function(x) cbind(cv_data[[x]]$test$shapefiles@data, fold = x)) %>% 
    do.call(rbind, .)
  
  test_shapes <- lapply(seq_along(cv_data), 
                        function(x) cv_data[[x]]$test$shapefiles) %>% 
    do.call(rbind, .)
  
  # test_pr <- lapply(seq_along(cv_data), 
  #                   function(x) cv_data(object[[x]]$test$pr, fold = x)) %>% 
  #   do.call(rbind, .)
  
  
  #test_shapes@data <- test_poly
  
  df <- ggplot2::fortify(test_shapes, region = 'area_id')
  
  df <- test_poly %>% 
    mutate(area_id = as.character(area_id)) %>% 
    left_join(df, ., by = c('id' = 'area_id'))
  
  # Create mosaiced predicted map
  r <- lapply(cv_preds$models, function(x) x$predictions[[layer]]) %>% 
    do.call(stack, .) %>% 
    mean
  
  if(mask){
    r <- mask(r, test_shapes)
  }
  
  test_shapes@data <- test_poly
  
  cv_raster <- rasterize(test_shapes, r, 'fold')
  
  for(i in sort(unique(test_poly$fold))){
    r[cv_raster == i] <- cv_preds$models[[i]]$predictions[[layer]][cv_raster == i]
  }
  
  r_df <- as.MAPraster(r)
  
  # Find limits
  if(is.null(lims)){
    if(trans %in% c('log', 'log10')){
      min <- min(min(r_df$z[r_df$z != 0]), min(df$response[df$response != 0]))
    } else {
      min <- min(min(r_df$z), min(df$response))
    }
    max <- max(max(r_df$z), max(df$response))
    lims <- c(min, max) 
  }
  
  
  # Create plots
  
  
  p <- ggplot(df, aes(long, lat, group = group, fill = response)) + 
    geom_polygon() +
    coord_equal() + 
    scale_fill_viridis_c(trans = trans, 
                         limits = lims, 
                         oob = scales::squish, 
                         name = legend_title,
                         ...)
  
  p2 <- ggplot(r_df, aes(x, y, fill = z)) + 
    geom_raster() +
    coord_equal() +
    scale_fill_viridis_c(trans = trans, 
                         limits = lims, 
                         oob = scales::squish, 
                         name = legend_title,
                         ...)
  
  
  # Combine in cowplot?
  ncols <- ifelse(column, 1, 2)
  grid <- plot_grid(p, p2, ncol = ncols)
  print(grid)
  invisible(list(p, p2))
}



##################################################
# List of classes and summary of useful print?  ##
##################################################

# ppj_model 
#   print parameters


# ppj_cv_performance
#  Print useful summaries.
#  Summary function might be useful.

