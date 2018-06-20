


fit_model <- function(data, mesh, model.args = NULL){

  
  startendindex <- make_startend_index(data)

  # Sort out mesh bits
  spde <- (inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
  Apix <- inla.mesh.project(mesh, loc = as.matrix(data$coords))$A
  
  n_s <- nrow(spde$M0)
  

  # TODO think about this more.
  data$covs[is.na(data$covs)] <- 0
  cov_matrix <- as.matrix(data$covs[, -c(1:2)])
  
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
  
  
  parameters <- list(intercept = -5,
                     slope = rep(0, nlayers(data$cov_rasters)),
                     log_tau = priormean_log_tau,
                     log_kappa = priormean_log_kappa,
                     nodemean = rep(0, n_s))

  input_data <- list(x = cov_matrix, 
                     xpop = data$pop,
                     Apixel = Apix,
                     spde = spde,
                     startendindex = startendindex,
                     polygon_cases = data$polygon$response * data$polygon$population / 1000,
                     prev_inc_par = c(2.616, -3.596, 1.594),
                     priormean_intercept = priormean_intercept,
                     priorsd_intercept = priorsd_intercept,
                     priormean_slope = priormean_slope,
                     priorsd_slope = priorsd_slope,
                     priormean_log_kappa = priormean_log_kappa,
                     priorsd_log_kappa = priorsd_log_kappa,
                     priormean_log_tau = priormean_log_tau,
                     priorsd_log_tau = priorsd_log_tau)
  

  obj <- MakeADFun(
    data = input_data, 
    parameters = parameters,
    #random = 'nodemean',
    DLL = "joint_model")
  
  its = 100
  opt <- nlminb(obj$par, obj$fn, obj$gr, 
                control = list(iter.max = its, eval.max = 2*its, trace = 0))

  # sd_out <- sdreport(obj)

  predictions <- predict_model(pars = opt$par, data, mesh)
  return(list(model = modelfit,
              predictions = preds))
}




predict_model <- function(pars, data, mesh){
  
  # Get raster coords
  raster_pts <- rasterToPoints(data$pop_raster %>% inset(is.na(.), value = -9999), spatial = TRUE)
  coords <- raster_pts@coords
  
  # Split up parameters
  pars <- split(pars, names(pars))
  
  
  # Get random field predicted
  spde <- (inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
  n_s <- nrow(spde$M0)						
  
  Amatrix <- inla.mesh.project(mesh, loc = as.matrix(coords))$A
  
  field <- MakeField(Amatrix, pars)
  field_ras <- rasterFromXYZ(cbind(coords, field))

  # Create linear predictor.
  covs_by_betas <- list()
  for(i in seq_len(nlayers(data$cov_rasters))){
    covs_by_betas[[i]] <- pars$slope[i] * data$cov_rasters[[i]]
  }
  
  cov_by_betas <- stack(covs_by_betas)
  cov_contribution <- sum(cov_by_betas) + pars$intercept
  
  linear_pred <- cov_contribution + field_ras
  
  prevalence <- 1 / (1 + exp(-1 * linear_pred))
  
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


MakeField <- function(Amatrix, pars){
  field <- (Amatrix %*% pars$nodemean)[, 1]
  return(field)
  
}





