


fit_model <- function(data, mesh, model.args = NULL){

  
  startendindex <- make_startend_index(data)

  # Sort out mesh bits
  spde <- (inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
  Apix <- inla.mesh.project(mesh, loc = as.matrix(data$coords))$A
  
  n_s <- nrow(spde$M0)
  

  # scale and transform variables.
  cov_matrix <- scale(data$covs[, -c(1:2)])


  
  priormean_log_kappa = -3
  priorsd_log_kappa = 0.5
  priormean_log_tau = -0.5
  priorsd_log_tau = 2
  
  priormean_intercept = -2
  priorsd_intercept = 1
  priormean_slope = 0
  priorsd_slope = 1

  

  # Compile and load the model

  
  dyn.load(dynlib("joint_model"))
  
  parameters <- list(intercept = -4,
                     slope = rep(0, nlayers(data$cov_rasters)),
                     log_tau = priormean_log_tau,
                     log_kappa = priormean_log_kappa,
                     nodemean = rep(0, n_s))

  input_data <- list(x = cov_matrix, 
                     xpop = data$pop,
                     Apixel = Apix,
                     spde = spde,
                     startendindex = startendindex,
                     polygon_mean_API = data$polygon$response,
                     prev_inc_par = c(2.616, -3.596, 1.594),
                     priormean_log_kappa = priormean_log_kappa,
                     priorsd_log_kappa = priorsd_log_kappa,
                     priormean_log_tau = priormean_log_tau,
                     priorsd_log_tau = priorsd_log_tau)
  

  obj <- MakeADFun(
    data = data, 
    parameters = parameters,
    #random = 'nodemean',
    DLL = "joint_model")
  
  opt <- nlminb(obj$par, obj$fn, obj$gr, 
                control = list(iter.max = its, eval.max = 2*its, trace = 10))



  predictions <- predict_model()
  return(list(model = modelfit,
              predictions = preds))
}




predict_model <- function(model, predict_model){


  return(predictions)
}



cv_performance <- function(predictions, holdout){



  return(list(polygon_pred_obs = polygon_pred_obs,
              pr_pred_obs = pr_pred_obs,
              polygon_metrics,
              pr_metrics))  
}





make_startend_index <- function(data){
  # Make start end index that determines whichs rows go with which polygon data.
  startendindex <- lapply(unique(data$covs$area_id), function(x) range(which(data$covs$area_id == x))) %>% 
    do.call(rbind, .)
  
  whichindices <- match(data$polygon$shapefile_id, unique(data$covs$area_id))
  
  # c++ is zero indexed.
  startendindex <- startendindex[whichindices, ] - 1L
}






