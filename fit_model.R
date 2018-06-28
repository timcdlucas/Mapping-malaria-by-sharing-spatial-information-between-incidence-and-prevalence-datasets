


fit_model <- function(data, mesh, its = 10, model.args = NULL){

  
  startendindex <- make_startend_index(data)

  # Sort out mesh bits
  spde <- (inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
  Apix <- inla.mesh.project(mesh, loc = as.matrix(data$coords))$A
  Apoint <- inla.mesh.project(mesh, loc = as.matrix(data$pr[, c('longitude', 'latitude')]))$A
  n_s <- nrow(spde$M0)
  

  # TODO think about this more.
  data$covs[is.na(data$covs)] <- 0
  cov_matrix <- as.matrix(data$covs[, -c(1:2)])
  
  
  # Define prior defaults
  priormean_log_kappa = -3
  priorsd_log_kappa = 0.3
  priormean_log_tau = 8
  priorsd_log_tau = 0.2
  
  priormean_intercept = -2
  priorsd_intercept = 3
  priormean_slope = 0
  priorsd_slope = 0.5
  
  # Replace defaults with anything given in model.args
  if(!is.null(model.args)){
    here <- environment()
    list2env(model.args, here)
  }

  

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
                     Apoint = Apoint,
                     pointx = data$pr_covs,
                     pointcases = data$pr$positive,
                     pointtested = data$pr$examined,
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
  
  opt <- nlminb(obj$par, obj$fn, obj$gr, 
                control = list(iter.max = its, eval.max = 2*its, trace = 0))

  # sd_out <- sdreport(obj)

  predictions <- predict_model(pars = opt$par, data, mesh)
  
  out <- list(model = list(opt, obj),
           predictions = predictions)
  class(out) <- c('ppj_model', 'list')
  return(out)
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
  api <- 1000 * PrevIncConversion(prevalence)
  
  incidence_count <- api * data$pop_raster / 1000
  
  predictions <- list(api = api, 
                      prevalence = prevalence, 
                      incidence_count = incidence_count,
                      pop = data$pop_raster)
  class(predictions) <- c('ppj_preds', 'list')
  return(predictions)
}



cv_performance <- function(predictions, holdout){

  # Extract raster data
  rasters <- stack(predictions$pop, predictions$incidence_count)
  names(rasters) <- c('population', 'incidence_count')
  
  extracted <- parallelExtract(rasters, holdout$shapefiles, fun = NULL, id = 'area_id')
  
  # Calc pred incidence and API
  
  aggregated <- extracted %>% 
                  na.omit %>% 
                  group_by(area_id) %>% 
                  summarise(pred_incidence_count = sum(incidence_count),
                            pred_pop = sum(population),
                            pred_api = 1000 * sum(incidence_count) / sum(population)) %>% 
                  left_join(holdout$polygon, by = c('area_id' = 'shapefile_id')) %>% 
                  mutate(incidence_count = response * population / 1000)
  # Calc API metrics
  polygon_metrics <- aggregated %>% 
                       summarise(RMSE = sqrt(mean((pred_api - response) ^ 2)),
                                 MAE = mean(abs(pred_api - response)),
                                 pearson = cor(pred_api, response, method = 'pearson'),
                                 spearman = cor(pred_api, response, method = 'spearman'),
                                 log_pearson = cor(log1p(pred_api), log1p(response), method = 'pearson'),
                                 log_spearman = cor(log1p(pred_api), log1p(response), method = 'spearman'),
                                 RMSE_cases = sqrt(mean((incidence_count - pred_incidence_count) ^ 2)),
                                 MAE_cases = mean(abs(incidence_count - pred_incidence_count)),
                                 total_cases = sum(incidence_count - pred_incidence_count))
  
  
  # Extract PR
  pr_coords <- SpatialPoints(holdout$pr[, c('longitude', 'latitude')])
  pr_preds <- extract(predictions$prevalence, pr_coords)
  pr_pred_obs <- cbind(holdout$pr, pred_prev= pr_preds) %>% 
                   mutate(prevalence = positive / examined)
  
  pr_pred_obs <- na.omit(pr_pred_obs)
  
  pr_metrics <- pr_pred_obs %>% 
    summarise(RMSE = sqrt(mean((pred_prev - prevalence) ^ 2)),
              MAE = mean(abs(pred_prev - prevalence)),
              pearson = cor(pred_prev, prevalence, method = 'pearson'),
              spearman = cor(pred_prev, prevalence, method = 'spearman'),
              log_pearson = cor(log1p(pred_prev), log1p(prevalence), method = 'pearson'),
              log_spearman = cor(log1p(pred_prev), log1p(prevalence), method = 'spearman'))
  
  
  out <- list(polygon_pred_obs = aggregated,
              pr_pred_obs = pr_pred_obs,
              polygon_metrics = polygon_metrics,
              pr_metrics = pr_metrics)
  class(out) <- c('ppj_cv_performance', 'list')
  return(out)  
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







#' Calc kappa and tau from rho and sd.,
#' 
#' 
#' 
#'  

rsd2kt <- function(rho, sigma){
  k <- sqrt(8) / rho
  t <- 1 / (sqrt(4 * pi) * k * sigma)
  
  return(list(tau = t, kappa = k))
}

k2r <- function(k) sqrt(8) / k
r2k <- function(r) sqrt(8) / r
r2logk <- function(r) log(sqrt(8) / r)


kt2rsd <- function(kappa, tau){
  rho <- sqrt(8) / kappa
  sigma <- 1 / (sqrt(4 * pi) * kappa * tau)
  
  return(list(rho = rho, sigma = sigma))
}



