



# load all data
# Perhaps read all raster years and have a step in process data to choose the right one? Or something. Kinda annoying.
data <- load_data(PR_path, 
                  API_path, 
                  pop_path, 
                  cov_raster_paths, 
                  shapefile_path, 
                  shapefile_pattern = '.shp$', 
                  useiso3 = 'MDG', 
                  admin_unit_level = 'ADMIN1',
                  pr_year = 2013,
                  api_year = 2013)


# madagascar

# pre analysis

data_mdg <- process_data(
  binomial_positive = data$pr$positive,
  binomial_n = data$pr$examined,
  coords = data$pr[, c('longitude', 'latitude')],
  polygon_response = data$api$api_mean,
  polygon_population = data$api$population,
  shapefile_id = data$api$shapefile_id,
  shps_id_column = 'area_id',
  shapefiles = data$shapefiles,
  pop_raster = data$pop,
  cov_rasters = data$covs,
  transform = c(4:7))
closeCluster(cl)

autoplot(data_mdg)

mesh_mdg <- build_mesh(data_mdg, mesh.args = list(max.edge = c(0.3, 5), cut = 0.3))

data_cv1_mdg <- cv_folds(data_mdg, k = 3)
autoplot(data_cv1_mdg, jitter = 0.2)


# run models
# Run full model to get a handle on things.

arg_list <- list(prior_rho_min = 3, # Mean of two thirds the spatial range. rho = 27, log_kappa = -2.446
                 prior_rho_prob = 0.00001, # Want p(rho < 3) = 0.0001 -> p(log_kappa < -0.058) = 0.0001
                 prior_sigma_max = 1, # Want p(sd > 1) = 0.0001 (would explain most of prev).  Wnat mean(sd) = 0.001. Do at large rho (50).
                 prior_sigma_prob = 0.00001,
                 prior_iideffect_sd_max = 0.05, 
                 # The difference between m_low_pf and LCI(pois(m_mean_pf)), then converted to inc rate, then to prev ranges around 0-0.025. 
                 # The 0.975 quantile of that (two sided) distribution is 0.005 prevalence. 
                 # To explain 0.005 prevalence, we need a norm of 0.05. Fin.
                 prior_iideffect_sd_prob = 0.000001, # Made this stronger because too much iid.
                 prior_iideffect_pr_sd_max = 0.05,
                 prior_iideffect_pr_sd_prob = 0.000001,
                 priormean_intercept = -2,
                 priorsd_intercept = 2,  # Indonesia has prev lowish. But want intercept to take whatever value it likes.
                 priormean_slope = 0, 
                 priorsd_slope = 0.4, # Explains between 0.004 and 0.27 prevalence. 1 covariate shouldn't explain between 0 and 0.6 (range of prev).
                 use_polygons = 0,
                 # use_polygons = 1,
                 use_points = 1)

full_model <- fit_model(data_mdg, mesh_mdg, its = 600, model.args = arg_list)
autoplot(full_model)
plot(full_model, layer = 'api')

in_sample <- cv_performance(predictions = full_model$predictions, 
                            holdout = data_mdg,
                            model_params = full_model$model)
autoplot(in_sample, CI = TRUE)
autoplot(in_sample, trans = 'log1p')

# Run 3 x models with 3 x hyperpars on cv1.
arg_list[c('use_polygons', 'use_points')] <- c(0, 1)
cv1_output1 <- run_cv(data_cv1_mdg, mesh_mdg, its = 200, model.args = arg_list)
obspred_map(data_cv1_mdg, cv1_output1, column = FALSE)
ggsave('figs/mdg_points_only.png')
obspred_map(data_cv1_mdg, cv1_output1, trans = 'log10', column = FALSE)

arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv1_output2 <- run_cv(data_cv1_mdg, mesh_mdg, its = 600, model.args = arg_list)
obspred_map(data_cv1_mdg, cv1_output2, column = FALSE)
ggsave('figs/mdg_polygons_only.png')
obspred_map(data_cv1_mdg, cv1_output2, trans = 'log10', column = FALSE)


arg_list[c('use_polygons', 'use_points')] <- c(1, 1)
cv1_output3 <- run_cv(data_cv1_mdg, mesh_mdg, its = 200, model.args = arg_list)
obspred_map(data_cv1_mdg, cv1_output3, column = FALSE)
obspred_map(data_cv1_mdg, cv1_output3, trans = 'log10', column = FALSE)

save(cv1_output1, file = 'model_outputs/mdg_points_cv_1.RData')
save(cv1_output2, file = 'model_outputs/mdg_polygon_cv_1.RData')
save(cv1_output3, file = 'model_outputs/mdg_join_cv_1.RData')

cv1_output1$summary$polygon_metrics
cv1_output2$summary$polygon_metrics
cv1_output3$summary$polygon_metrics

cv1_output1$summary$pr_metrics
cv1_output2$summary$pr_metrics
cv1_output3$summary$pr_metrics








# create temp figures




# Write out data needed for final figures.