


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
  
  
  # 
  # // Pixel data. 
  # // All in long format so that if a pixel is in multiple polygons or multiple years, it will be represented by multiple rows.
  # // Environmental/other covariates data matrix
  # DATA_MATRIX(x);
  # 
  # // Population fraction of a polygon in a pixel. i.e. pixel i makes contains 0.01 of the population of polygon j
  # DATA_VECTOR(xpop);
  # 
  # // two col matrix with start end indices for each shape case. 
  # DATA_IARRAY(startendindex); 
  # 
  # // Shape data. Cases and region id.
  # DATA_VECTOR(polygon_mean_API);
  # DATA_VECTOR(polygon_sd_API);
  # 
  # // Weight polygon likelihood by this much
  # DATA_SCALAR(polyweight);
  # 
  # 
  # 
  # // ------------------------------------------------------------------------ //
  #   // Parameters
  # // ------------------------------------------------------------------------ //
  #   
  #   // regression slopes
  # // (log of) empirical mean incidence to guide intercept
  # PARAMETER(intercept); // intercept
  # PARAMETER_VECTOR(slope); 
  # 
  # 
  # 
  # DATA_SCALAR(priormean_intercept); // = -4.0; 
  # DATA_SCALAR(priorsd_intercept);// = 2.0
  # DATA_SCALAR(priormean_slope); // = 0.0;
  # DATA_SCALAR(priorsd_slope); // = 1.0;
  # 
  # // 2016 spde hyperparameters
  # // tau defines strength of random field. 
  # // kappa defines distance within which points in field affect each other. 
  # PARAMETER(log_tau);
  # PARAMETER(log_kappa);
  # 
  # // Priors on spde hyperparameters
  # //   kappa -- i.e. exp(priormean_log_kappa) -- set as approximately the width of the region being studied.
  # //   This implies prior belief in a fairly flat field.
  # //   tau -- exp(priormean_log_tau) -- set to close to zero. Betas on regression coefficients have priors of 0 so this is reasonable.
  # //Type priormean_log_kappa = -3;
  # //Type priorsd_log_kappa   = 0.5;
  # //Type priormean_log_tau   = -0.50;
  # //Type priorsd_log_tau     = 2.0;
  # DATA_SCALAR(priormean_log_kappa);
  # DATA_SCALAR(priorsd_log_kappa);
  # DATA_SCALAR(priormean_log_tau);
  # DATA_SCALAR(priorsd_log_tau);
  # 
  # // Convert hyperparameters to natural scale
  # Type tau = exp(log_tau);
  # Type kappa = exp(log_kappa);
  # 
  # 
  # // Space-time random effect parameters
  # // matrix logit_pr_offset [nrows = n_mesh, col=n_years].
  # PARAMETER_VECTOR(nodemean);
  # 
  # // Prevalence to incidence conversion parameters 
  # DATA_VECTOR(prev_inc_par); // length: 3 
  # DATA_SCALAR(max_prev);
  # 
  # 
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






