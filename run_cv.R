

run_cv <- function(cv_data, mesh, its = 10, model.args = NULL, parallel_delay = 0){

  stopifnot(inherits(cv_data, 'ppj_cv'))

  models <- list()
  results <- list()
  
  for(i in seq_along(cv_data)){
    message('Fitting model: ', i)
    Sys.sleep(runif(1, 0, parallel_delay))
    models[[i]] <- fit_model(cv_data[[i]]$train, mesh, its, model.args)
    #full_model <- fit_model(data_idn, mesh_idn, its = 200, model.args = arg_list)

    results[[i]] <- cv_performance(models[[i]]$predictions, cv_data[[i]]$test)
  }
  
  summary <- summarise_cv_results(results)

  return(list(summary = summary, models = models, results = results))
}



summarise_cv_results <- function(results){
  
  combined_aggregated <- 
    lapply(seq_along(results), function(x) cbind(results[[x]]$polygon_pred_obs, fold = x)) %>% 
      do.call(rbind, .)
  
  polygon_metrics <- combined_aggregated %>% 
    summarise(RMSE = sqrt(mean((pred_api - response) ^ 2)),
              MAE = mean(abs(pred_api - response)),
              pearson = cor(pred_api, response, method = 'pearson'),
              spearman = cor(pred_api, response, method = 'spearman'),
              RMSE_cases = sqrt(mean((incidence_count - pred_incidence_count) ^ 2)),
              MAE_cases = mean(abs(incidence_count - pred_incidence_count)),
              total_cases = sum(incidence_count - pred_incidence_count))
  
  combined_summaries <- 
    lapply(seq_along(results), function(x) cbind(results[[x]]$polygon_metrics, fold = x)) %>% 
    do.call(rbind, .)
  
  return(list(combined_aggregated = combined_aggregated,
              combined_summaries = combined_summaries,
              polygon_metrics = polygon_metrics))
  
}




