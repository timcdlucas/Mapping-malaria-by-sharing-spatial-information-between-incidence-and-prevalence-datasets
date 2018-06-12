

run_cv <- function(cv_data, mesh, model.args, parallel_delay = 60){

  stopifnot(inherits(cv_data, 'ppj_cv'))

  for(i in seq_along(cv_data)){

    Sys.sleep(runif(0, 60, 1))
    models[[i]] <- fit_model(cv_data[[i]]$train, mesh, model.args)
    results[[i]] <- cv_performance(moodel[[i]], cv_data[[i]]$test)
  }

  return(list(summary = summary, models = models, results = results))
}
