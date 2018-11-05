
#'@param data An object of class ppj_data.
#'@param mesh An INLA mesh
#'@param its Maximum iterations for the model optimisation.
#'@param model.args Names list of non-default model arguments. See below.
#'@param CI What confidence interval to return for raster realisations
#'@param N How many realisations to make.

#'
#' Mode defaults are this. 
#'   Any of these can be set manually with model.args
#' priormean_log_kappa = -3
#' priorsd_log_kappa = 0.3
#' priormean_log_tau = 8
#' priorsd_log_tau = 0.2
#'
#' priormean_intercept = -2
#' priorsd_intercept = 3
#' priormean_slope = 0
#' priorsd_slope = 0.5
#' 
#' use_polygons = 1
#' use_points = 1


fit_model <- function(data, mesh, its = 10, model.args = NULL, CI = 0.95, N = 100){

  
  startendindex <- make_startend_index(data)

  # Sort out mesh bits
  spde <- (inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
  Apix <- inla.mesh.project(mesh, loc = as.matrix(data$coords))$A
  Apoint <- inla.mesh.project(mesh, loc = as.matrix(data$pr[, c('longitude', 'latitude')]))$A
  n_s <- nrow(spde$M0)
  

  # TODO think about this more.
  data$covs[is.na(data$covs)] <- 0
  data$pr_covs[is.na(data$pr_covs)] <- 0
  
  cov_matrix <- as.matrix(data$covs[, -c(1:2)])
  
  

  prior_rho_min = 3
  prior_rho_prob = 0.00001
  prior_sigma_max = 1
  prior_sigma_prob = 0.00001
  prior_iideffect_sd_max = 0.05
  prior_iideffect_sd_prob = 0.00001
  prior_iideffect_pr_sd_max = 0.05
  prior_iideffect_pr_sd_prob = 0.00001
  priormean_intercept = -2
  priorsd_intercept = 2
  priormean_slope = 0 
  priorsd_slope = 0.4
  use_polygons = 0
  use_points = 1


  # Replace defaults with anything given in model.args
  if(!is.null(model.args)){
    here <- environment()
    list2env(model.args, here)
  }

  # Construct vector of length PR data where each value is the element of polygon enclosing that point
  overlap <- unique(c(data$polygon$shapefile_id, data$pr$shapefile_id))
  if(length(overlap) > length(data$polygon$shapefile_id)) {
    extra <- length(overlap) - length(data$polygon$shapefile_id)
    message(paste("There are", extra, "shapefiles that contain point data but not polygon data"))
  }
  # For each PR point (ordered), the map gives the corresponding index in overlap
  # map is length of PR data. Ordering of polygon data shapefile_id
  pointtopolygonmap <- match(data$pr$shapefile_id, overlap)
  
  # Compile and load the model

  #dyn.unload(dynlib("joint_model"))
  dyn.load(dynlib("joint_model"))
  
  
  parameters <- list(intercept = -5,
                     slope = rep(0, nlayers(data$cov_rasters)),
                     iideffect = rep(0, length(overlap)),
                     iideffect_log_tau = 1,
                     iideffect_pr = rep(0, nrow(data$pr)),
                     iideffect_pr_log_tau = 1,
                     log_sigma = 0,
                     log_rho = 4,
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
                     pointtopolygonmap = pointtopolygonmap,
                     prev_inc_par = c(2.616, -3.596, 1.594),
                     prior_rho_min = prior_rho_min,
                     prior_rho_prob = prior_rho_prob,
                     prior_sigma_max = prior_sigma_max,
                     prior_sigma_prob = prior_sigma_prob,
                     prior_iideffect_sd_max = prior_iideffect_sd_max,
                     prior_iideffect_sd_prob = prior_iideffect_sd_prob,
                     prior_iideffect_pr_sd_max = prior_iideffect_pr_sd_max,
                     prior_iideffect_pr_sd_prob = prior_iideffect_pr_sd_prob,
                     priormean_intercept = priormean_intercept,
                     priorsd_intercept = priorsd_intercept,
                     priormean_slope = priormean_slope, 
                     priorsd_slope = priorsd_slope,
                     use_polygons = use_polygons,
                     use_points = use_points)
  
  if(!use_points){
    fix <- list(iideffect_pr_log_tau = factor(NA), 
                iideffect_pr = factor(rep(NA, nrow(data$pr))))
    parameters$iideffect_pr_log_tau <- -100
  } else {
    fix <- list()
  }

  obj <- MakeADFun(
    data = input_data, 
    parameters = parameters,
    random = c('nodemean', 'iideffect', 'iideffect_pr'),
    map = fix,
    DLL = "joint_model")
  
  
  opt <- tryCatch({
    nlminb(obj$par, obj$fn, obj$gr, 
           control = list(iter.max = its, eval.max = 2*its, trace = 0))
    },
    error = function(e) {
      cat('Error in model. Stopping')
      return('error')
    }
  )
  
  # Rerun if needed.
  if(opt$convergence != 0){ 
    warning('Model did not converge on first try.')
    
    opt <- tryCatch({
      nlminb(obj$par + rnorm(length(obj$par)), obj$fn, obj$gr, 
             control = list(iter.max = its, eval.max = 2*its, trace = 0))
      },
      error = function(e) {
        cat('Error in model. Stopping')
        return('error')
      }
    )
    
    # Soft check
    if(opt$convergence != 0){ 
      warning('Model did not converge on second try.')
    }  
    
  } 

  cat("Optimisation has finished. Now moving onto sdreport.\n")
  
  sd_out <- sdreport(obj, getJointPrecision = TRUE)
  
  # Check sdreport worked.
  if(anyNA(sd_out$cov.fixed) | anyNA(sd_out$jointPrecision)) stop('sdreport failed. NAs in fixed SDs or joinPrecision')
  
  cat("Model fitting has finished. Now moving onto prediction.\n")
  
  
  predictions <- predict_model(pars = obj$env$last.par.best, data, mesh, overlap)
  uncertainty <- predict_uncertainty(pars = obj$env$last.par.best, 
                                     joint_pred = sd_out$jointPrecision,
                                     data, 
                                     mesh,
                                     overlap,
                                     N = N,
                                     CI)
  
  # todo to get CIs for aggregated polygons... probably need to do it on every layer.
  
  predictions <- c(predictions, uncertainty)
  
  out <- list(model = list(opt = opt, obj = obj, sd_report = sd_out),
              predictions = predictions)
  class(out) <- c('ppj_model', 'list')
  return(out)
}


predict_uncertainty <- function(pars, joint_pred, data, mesh, shapefile_ids, N, CI = 0.95){
  
  # Extract the Amatrix and coords of the field
  field_properties <- extractFieldProperties(data, mesh)

  # Initialise iid raster with values of polygon indices - rasterize is very slow, cannot be done each realisation
  # iid_ras <- raster(ncols = ncol(data$cov_rasters), nrows = nrow(data$cov_rasters), ext = extent(data$cov_rasters))
  # iid_ras <- rasterize(data$shapefiles, iid_ras, seq(1:nrow(data$polygon)))

  # Get par draws
  ch <- Cholesky(joint_pred)
  par_draws <- sparseMVN::rmvn.sparse(N, pars, ch, prec = TRUE)
                    
  prevalence <- list()
  api <- list()
  
  for(r in seq_len(N)){
    
    # Split up parameters
    p <- split(par_draws[r, ], names(pars))
    
    # Extract field values
    field <- MakeField(field_properties$Amatrix, p)
    field_ras <- rasterFromXYZ(cbind(field_properties$coords, field))
    
    linear_pred_result <- makeLinearPredictor(p, data, field_ras, data$shapefile_raster, shapefile_ids)
    linear_pred <- linear_pred_result$linear_pred

    linear_pred <- iidDraw(linear_pred, p, data$shapefile_raster, shapefile_ids)

    prevalence[[r]] <- 1 / (1 + exp(-1 * linear_pred))
    
  }
  
  prevalence <- do.call(stack, prevalence)
  
  probs <- c((1 - CI) / 2, 1 - (1 - CI) / 2)
  
  quant <- function(x) quantile(x, probs = probs, na.rm = TRUE)
  
  prevalence_ci <- calc(prevalence, fun = quant)
  api_ci <- 1000 * PrevIncConversion(prevalence_ci)
  
  api <- 1000 * PrevIncConversion(prevalence)
  incidence_count <- api * data$pop_raster / 1000
  
  predictions <- list(api_ci = api_ci, 
                      prevalence_ci = prevalence_ci,
                      api_realisations = api,
                      prevalence_realisations = prevalence,
                      incidence_count_realisations = incidence_count,
                      par_draws = par_draws)
  class(predictions) <- c('ppj_preds_ci', 'list')
  return(predictions)
}


iidDraw <- function(linear_pred, pars, shapefile_raster, shapefile_ids){

  oos_shapefile_ids <- unique(shapefile_raster)
  #oos_shapefile_ids <- unique_shapefile_ids[!unique_shapefile_ids %in% shapefile_ids]

  iid_samples <- rnorm(length(oos_shapefile_ids), 0, 1 / sqrt(exp(pars$iideffect_log_tau)))
  iid_ras <- shapefile_raster
  for(i in seq_along(iid_samples)) {
    iid_ras@data@values[which(shapefile_raster@data@values == oos_shapefile_ids[i])] <- iid_samples[i]
  }
  iid_ras[!shapefile_raster %in% oos_shapefile_ids & !is.na(iid_ras)] <- 0

  linear_pred <- linear_pred + iid_ras

  return(linear_pred)
}

predict_model <- function(pars, data, mesh, shapefile_ids){
  
  # Extract the Amatrix and coords of the field
  field_properties <- extractFieldProperties(data, mesh)
  
  # Initialise iid raster with values of polygon indices
  # We do the rasterise thing at least twice. Here and in unc. Maybe pull out.
  # iid_ras <- raster(ncols = ncol(data$cov_rasters), nrows = nrow(data$cov_rasters), ext = extent(data$cov_rasters))
  # iid_ras <- rasterize(data$shapefiles, iid_ras, seq_len(nrow(data$polygon)))
  
  # Split up parameters
  pars <- split(pars, names(pars))
  
  # Extract field values
  field <- MakeField(field_properties$Amatrix, pars)
  field_ras <- rasterFromXYZ(cbind(field_properties$coords, field))
 

  # Create linear predictor
  linear_pred_result <- makeLinearPredictor(pars, data, field_ras, data$shapefile_raster, shapefile_ids)
  linear_pred <- linear_pred_result$linear_pred
  
  prevalence <- 1 / (1 + exp(-1 * linear_pred))
  api <- 1000 * PrevIncConversion(prevalence)
  
  incidence_count <- api * data$pop_raster / 1000
  
  predictions <- list(api = api, 
                      prevalence = prevalence, 
                      incidence_count = incidence_count,
                      pop = data$pop_raster,
                      field = field_ras,
                      covariates = linear_pred_result$covariates,
                      iid = linear_pred_result$iid
                      )
  class(predictions) <- c('ppj_preds', 'list')
  return(predictions)
}


extractFieldProperties <- function(data, mesh) {
  
  # Get raster coords
  raster_pts <- rasterToPoints(data$pop_raster %>% inset(is.na(.), value = -9999), spatial = TRUE)
  coords <- raster_pts@coords
  
  # Get random field predicted
  spde <- (inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
  n_s <- nrow(spde$M0)						
  
  Amatrix <- inla.mesh.project(mesh, loc = as.matrix(coords))$A
  
  return(list(Amatrix = Amatrix, 
              coords = coords))
}

makeLinearPredictor <- function(pars, data, field_ras, shapefile_ras, shapefile_ids) {
  
  covs_by_betas <- list()
  for(i in seq_len(nlayers(data$cov_rasters))){
    covs_by_betas[[i]] <- pars$slope[i] * data$cov_rasters[[i]]
  }
  
  cov_by_betas <- stack(covs_by_betas)
  cov_contribution <- sum(cov_by_betas) + pars$intercept
  
  # Replace raster values with correct values from pars$iideffect
  # Need to make sure the indexing is correct!!!

  iid_ras <- shapefile_ras
  for(i in seq_along(pars$iideffect)) {
    #iid_ras[iid_ras == i] <- pars$iideffect[i] # perhaps neater but below is faster. 
    iid_ras@data@values[which(shapefile_ras@data@values == shapefile_ids[i])] <- pars$iideffect[i]
  }
  iid_ras[!shapefile_ras %in% shapefile_ids & !is.na(iid_ras)] <- 0
  

  linear_pred <- cov_contribution + field_ras #+ iid_ras
  
  return(list(linear_pred = linear_pred, 
              covariates = cov_contribution,
              iid = iid_ras))
}


cv_performance <- function(predictions, holdout, model_params, CI = 0.95, use_points, serial_extract = FALSE){

  # Extract raster data
  rasters <- stack(predictions$pop, predictions$incidence_count)
  names(rasters) <- c('population', 'incidence_count')
  
  if(serial_extract){
    extracted <- raster::extract(rasters, holdout$shapefiles, cellnumbers = TRUE, df = TRUE)
    extracted[, 1] <- shapefiles$area_id[extracted[, 1]]
    names(extracted)[1] <- 'area_id'
  } else {
    cl <- makeCluster(min(detectCores() - 1, 20))
    registerDoParallel(cl)
    extracted <- parallelExtract(rasters, holdout$shapefiles, fun = NULL, id = 'area_id')
    stopCluster(cl)
    registerDoSEQ()
  }
  
  # extracted <- parallelExtract(rasters, holdout$shapefiles, fun = NULL, id = 'area_id')
  #extracted <- extract(rasters, holdout$shapefiles, fun = NULL, id = 'area_id')
  
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
                                 RMSE_cases = sqrt(mean((incidence_count - pred_incidence_count) ^ 2)),
                                 MAE_cases = mean(abs(incidence_count - pred_incidence_count)),
                                 total_cases = sum(incidence_count - pred_incidence_count))
  
  
  
  # Uncertainty metrics
  
  # Todo. Think about MAP estimate for metrics vs mean of realisations.
  #   For polygons this is particularly important.
  #   Currently MAP estimate is NOT the mean of the realisaitons.
  
  rasters_real <- stack(predictions$pop, predictions$incidence_count_realisations)
  names(rasters_real)[1] <- 'population'
  
  
  if(serial_extract){
    extracted_reals <- raster::extract(rasters_real, holdout$shapefiles, cellnumbers = TRUE, df = TRUE)
    extracted_reals[, 1] <- shapefiles$area_id[extracted_reals[, 1]]
    names(extracted_reals)[1] <- 'area_id'
  } else {
    cl <- makeCluster(min(detectCores() - 1, 20))
    registerDoParallel(cl)
    extracted_reals <- parallelExtract(rasters_real, 
                                       holdout$shapefiles, fun = NULL, id = 'area_id')
    stopCluster(cl)
    registerDoSEQ()
  }
  
  
  #extracted_reals[, -c(1:3)] <- as.matrix(extracted_reals[, -c(1:3)]) * extracted_reals$population
  
  extract_reals_tidy <- gather(extracted_reals, key = layer, value = incidence_count, -cellid, -area_id, -population)
  
  # Calc pred incidence and API
  probs <- c((1 - CI) / 2, 1 - (1 - CI) / 2)
  
  aggregated_reals <- extract_reals_tidy %>%
    na.omit %>%
    group_by(area_id, layer) %>%
    summarise(pred_incidence_count = sum(incidence_count),
              pred_pop = sum(population),
              pred_api = 1000 * sum(pred_incidence_count) / sum(population))
  
  #pred_incidence_count_noise <- sapply(aggregated_reals$pred_incidence_count, function(x) rpois(1, x))
  pred_incidence_count_noise <- rpois(rep(1, nrow(aggregated_reals)), aggregated_reals$pred_incidence_count)
    #rpois(length(aggregated_reals$pred_incidence_count), aggregated_reals$pred_incidence_count)
  pred_api_noise <- 1000 * pred_incidence_count_noise / aggregated_reals$pred_pop
  
  aggregated_reals <- cbind(aggregated_reals, 
                            pred_incidence_count_noise = pred_incidence_count_noise,
                            pred_api_noise = pred_api_noise)
  
  aggregated_reals <- aggregated_reals %>%
    group_by(area_id) %>%
    summarise(pred_incidence_count_lower = quantile(pred_incidence_count_noise, probs[1]),
              pred_incidence_count_upper = quantile(pred_incidence_count_noise, probs[2]),
              pred_api_lower = quantile(pred_api_noise, probs[1]),
              pred_api_upper = quantile(pred_api_noise, probs[2])) %>%
    left_join(aggregated, by = c('area_id'))
  
  polygon_metrics_unc <- aggregated_reals %>% 
                           mutate(inCI = response > pred_api_lower & response < pred_api_upper) %>% 
                           summarise(coverage = mean(inCI))
  
  
  polygon_metrics <- cbind(polygon_metrics, polygon_metrics_unc)
  
  
  
  
  # Extract PR
  pr_coords <- SpatialPoints(holdout$pr[, c('longitude', 'latitude')])
  pr_preds <- raster::extract(predictions$prevalence, pr_coords)
  pr_pred_obs <- cbind(holdout$pr, pred_prev= pr_preds) %>% 
    mutate(prevalence = positive / examined)
  
  
  pr_preds_reals <- raster::extract(predictions$prevalence_realisations, pr_coords)
  
  # Convert prev to linear predictor, sum with iid and then convert back
  pr_preds_reals_lp <- log(pr_preds_reals/(1 - pr_preds_reals))
  pr_reals_lp <- pr_preds_reals_lp
  
  # Calculate additional uncertainty from point iid effect
  if(use_points) { 
    iid_par_index <- which(names(model_params$obj$env$last.par.best) == 'iideffect_pr_log_tau')
    pr_iid_log_tau_reals <- predictions$par_draws[, iid_par_index] # vector of N realisations
    pr_iid_sd_reals <- 1/sqrt(exp(pr_iid_log_tau_reals)) # vector of N realisations
    
    pr_iid_reals <- sapply(pr_iid_sd_reals, function(x) rnorm(length(pr_coords), 0, x)) # matrix of n(pr) by N
    
    # Total uncertainty in linear predictor
    pr_reals_lp <- pr_preds_reals_lp + pr_iid_reals
  }
  
  # Transform back to prevelance
  pr_reals_prev <- 1 / (1 + exp(-1 * pr_reals_lp))
  
  positive_reals_noise <- sapply(seq_len(nrow(pr_pred_obs)), 
                                 function(i) 
                                   rbinom(ncol(pr_reals_prev), pr_pred_obs$examined[i], pr_reals_prev[i, ])
                                 )
  pr_reals_prev_noise <- t(positive_reals_noise) / pr_pred_obs$examined
  
  pr_conf <- apply(pr_reals_prev_noise, 1, function(x) quantile(x, probs = probs, na.rm = TRUE))
  
  pr_pred_obs <- cbind(pr_pred_obs, t(pr_conf))
  names(pr_pred_obs)[c(8, 9)] <- c('prevalence_lower',  'prevalence_upper')

  pr_pred_obs <- na.omit(pr_pred_obs)
  
  pr_metrics <- pr_pred_obs %>% 
    mutate(inCI = prevalence > prevalence_lower & prevalence < prevalence_upper) %>% 
    summarise(RMSE = sqrt(mean((pred_prev - prevalence) ^ 2)),
              MAE = mean(abs(pred_prev - prevalence)),
              pearson = cor(pred_prev, prevalence, method = 'pearson'),
              spearman = cor(pred_prev, prevalence, method = 'spearman'),
              log_pearson = cor(log1p(pred_prev), log1p(prevalence), method = 'pearson'),
              coverage = mean(inCI))
  
  
  
  
  out <- list(polygon_pred_obs = aggregated_reals,
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


find_max_rho <- function(raster){
  xrange <- raster@extent[2] - raster@extent[1]
  yrange <- raster@extent[4] - raster@extent[3]
  
  rho <- sqrt(yrange^2 + xrange^2)
  rho * 2/3
}

find_max_kappa <- function(raster){
  rho <- find_max_rho(raster)
  kappa <- r2k(rho)
}

find_max_logkappa <- function(raster){
  kappa <- find_max_kappa(raster)
  log_kappa <- log(kappa)
}


