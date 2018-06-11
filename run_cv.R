

run_cv <- function(cv_data, mesh, model.args){

  stopifnot(inherits(cv_data, 'ppj_cv'))

  for(i in seq_along(cv_data)){
    model[[i]] <- fit_model(cv_data[[i]]$train, mesh, model.args)
    results[[i]] <- cv_performance(moodel[[i]], cv_data[[i]]$test)
  }

  return(list(model, results)
}
