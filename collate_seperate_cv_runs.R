
# Functions for collating sepeartely fitted folds.



collate <- function(model_list){
  
  models <- lapply(model_list, function(x) x$models)
  results <- lapply(model_list, function(x) x$results)
  
  summary <- summarise_cv_results(results)
  out <- list(summary = summary, models = models, results = results)
  class(out) <- c('ppf_cv_results', 'list')
  return(out)
}