
########
# master script for point Vs polygon Vs joint analysis
# Tim Lucas
# 2018-05-30
###########

if(Sys.info()["user"] != 'anita'){
  setwd('~/timz/timothy/point_polygon_joint_comparison')
} else {
  setwd('~/Z/users/anita/point_polygon_join_comparison_analysis')
}

source("setUserInfo.R")

# define paths

PR_path <- Z('GBD2017/Processing/Stages/04b_PR_DB_Import_Export/Verified_Outputs/2018_02_15/pfpr_dump.csv')
API_path <- Z('GBD2017/Processing/Stages/04c_API_Data_Export/Checkpoint_Outputs/subnational.csv')
pop_path <- Z('GBD2017/Processing/Stages/03_Muster_Population_Figures/Verified_Outputs/Output_Pop_At_Risk_Pf_5K/ihme_corrected_frankenpop_All_Ages_3_2015_at_risk_pf.tif')
shapefile_path <- Z('master_geometries/Admin_Units/Global/GBD/GBD2017_MAP/GBD2017_MAP_MG_5K/')

cov_raster_paths <- c(
  Z('mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/5km/Synoptic/LST_Day_v6.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/EVI_v6/5km/Synoptic/EVI_v6.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/TSI-Martens2-Pf.Synoptic.Overall.Mean.5km.Data.tif'),
  Z('GBD2017/Processing/Static_Covariates/MAP/other_rasters/accessibility/accessibility.5k.MEAN.tif'),
  Z('mastergrids/Other_Global_Covariates/Elevation/SRTM-Elevation/5km/Synoptic/SRTM_elevation.Synoptic.Overall.Data.5km.mean.tif'),
  Z('mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/5km/Synoptic/LST_Day_v6.Synoptic.Overall.SD.5km.mean.tif'),
  #Z('mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCB/5km/Synoptic/TCB.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/Other_Global_Covariates/NightTimeLights/VIIRS_DNB_Composites/5km/Annual/VIIRS-SLC.2016.Annual.5km.MEDIAN.tif'),
  #Z('mastergrids/Other_Global_Covariates/UrbanAreas/Global_Urban_Footprint/From_86m/5km/Global_Urban_Footprint_5km_PropUrban.tif'),
  Z('mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/TCW_v6/5km/Synoptic/TCW_v6.Synoptic.Overall.mean.5km.mean.tif')
)

#Z('mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCB/5km/Synoptic/TCB.Synoptic.Overall.mean.5km.mean.tif'),
#Z('mastergrids/Other_Global_Covariates/UrbanAreas/Global_Urban_Footprint/From_86m/5km/Global_Urban_Footprint_5km_PropUrban.tif'),

# load packages

## Spatial packages
library(raster)
library(maptools)
library(rgeos)

## dataframe packages
library(dplyr)
library(readr)
library(magrittr)
library(tidyr)

## For standardising prevalence
library(malariaAtlas)

## For inc prevalence conversions
library(GBDutils)

## Plotting packages
library(ggplot2)
library(cowplot)
theme_set(theme_minimal())

##  Modelling packages
library(TMB)
#library(stantmb)
library(INLA)
if(Sys.info()["sysname"] != 'Windows'){
  message('using INLA unix workaround')
  INLA:::inla.dynload.workaround()
} else {
  message('Not using INLA unix workaround. Expect you are using winows.')
}
library(INLAutils)
library(sparseMVN)


# Parallel processing
library(foreach)
library(doParallel)

# load functions

source('collect_data.R')
source('process_data.R')
source('CombineRasters.R')
source('parallel-raster-extract.R')
source('build_inla_meshes.R')
source('fit_model.R')
source('run_cv.R')
source('random_crossvalidation_setup.R')
source('spatial_crossvalidation_setup.R')
source('plotting_functions.R')

# Compile the model
compile("joint_model.cpp")

set.seed(180530)

# load all data
# Perhaps read all raster years and have a step in process data to choose the right one? Or something. Kinda annoying.
data <- load_data(PR_path, 
                  API_path, 
                  pop_path, 
                  cov_raster_paths, 
                  shapefile_path, 
                  shapefile_pattern = '.shp$', 
                  useiso3 = 'MDG', 
                  pr_country = 'Madagascar',
                  admin_unit_level = 'ADMIN3',
                  pr_year = 2011:2016,
                  api_year = 2013)

#test years
#pr_year = 2016,
#api_year = 2015)

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
  useiso3 = 'MDG',
  add_pr_gp = TRUE,
  transform = 4:7)
save(data_mdg, file = 'model_outputs/mdg_full_data.RData')

autoplot(data_mdg)

mesh_mdg <- build_mesh(data_mdg, mesh.args = list(max.edge = c(0.2, 3), cut = 0.2, offset = c(2, 5)))

data_cv1_mdg <- cv_random_folds(data_mdg, k = 10)
autoplot(data_cv1_mdg, jitter = 0)
save(data_cv1_mdg, file = 'model_outputs/mdg_cv_1.RData')



# run models
# Run full model to get a handle on things.

arg_list <- list(prior_rho_min = 1, # 
                 prior_rho_prob = 0.00001, # Want p(rho < 3) = 0.0001 -> p(log_kappa < -0.058) = 0.0001
                 prior_sigma_max = 1, # Want p(sd > 1) = 0.0001 (would explain most of prev).
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
                 use_polygons = 1,
                 use_points = 1)

if(FALSE){
  full_model <- fit_model(data_mdg, mesh_mdg, its = 1000, model.args = arg_list, drop_covs = 9)
  autoplot(full_model)
  
  png('figs/mdg_full_model_in_sample_map.png')
  plot(full_model, layer = 'api')
  dev.off()
  
  in_sample <- cv_performance(predictions = full_model$predictions, 
                              holdout = data_mdg,
                              model_params = full_model$model, 
                              CI = 0.8,
                              use_points = arg_list$use_points)
  autoplot(in_sample, CI = TRUE)
  autoplot(in_sample, trans = 'log1p', CI = TRUE)
  ggsave('figs/mdg_full_model_in_sample.png')
  
  save(full_model, file = 'model_outputs/full_model_mdg.RData')
  
  
  arg_list$use_points <- 0
  prevgp_model <- fit_model(data_mdg, mesh_mdg, its = 1000, model.args = arg_list)
  autoplot(prevgp_model)
  
  png('figs/mdg_prevgp_model_in_sample_map.png')
  plot(prevgp_model, layer = 'api')
  dev.off()
  
  in_sample <- cv_performance(predictions = prevgp_model$predictions, 
                              holdout = data_mdg,
                              model_params = prevgp_model$model, 
                              CI = 0.8,
                              use_points = arg_list$use_points)
  autoplot(in_sample, CI = TRUE)
  autoplot(in_sample, trans = 'log1p', CI = TRUE)
  ggsave('figs/mdg_prevgp_model_in_sample.png')
  
  save(prevgp_model, file = 'model_outputs/full_prevgp_mdg.RData')
  
  
  
  
  
  arg_list$use_points <- 0
  baseline_model <- fit_model(data_mdg, mesh_mdg, its = 1000, model.args = arg_list, drop_covs = 9)
  autoplot(baseline_model)
  
  png('figs/mdg_baseline_model_in_sample_map.png')
  plot(baseline_model, layer = 'api')
  dev.off()
  
  in_sample <- cv_performance(predictions = baseline_model$predictions, 
                              holdout = data_mdg,
                              model_params = baseline_model$model, 
                              CI = 0.8,
                              use_points = arg_list$use_points)
  autoplot(in_sample, CI = TRUE)
  autoplot(in_sample, trans = 'log1p', CI = TRUE)
  ggsave('figs/mdg_baseline_model_in_sample.png')
  
  save(baseline_model, file = 'model_outputs/full_baseline_mdg.RData')
  
  
  
  
}

delay <- 2

cat('Start cv1 model 1')
# Run 3 x models with 3 x hyperpars on cv1.
arg_list[c('use_polygons', 'use_points')] <- c(0, 1)
cv1_output1 <- run_cv(data_cv1_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      parallel_delay = delay, cores = 10, 
                      drop_covs = 9)
obspred_map(data_cv1_mdg, cv1_output1, column = FALSE)
ggsave('figs/mdg_points_only_obspred_map.png')
obspred_map(data_cv1_mdg, cv1_output1, trans = 'log10', column = FALSE)
ggsave('figs/mdg_points_only_obspred_map_log.png')
autoplot(cv1_output1, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_points_only_obspred.png')
save(cv1_output1, file = 'model_outputs/mdg_points_cv_1.RData')


cat('Start cv1 model 2')
arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv1_output2 <- run_cv(data_cv1_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      parallel_delay = delay, cores = 10, 
                      drop_covs = 9)
obspred_map(data_cv1_mdg, cv1_output2, column = FALSE)
ggsave('figs/mdg_polygons_only_obspred_map.png')
obspred_map(data_cv1_mdg, cv1_output2, trans = 'log10', column = FALSE)
ggsave('figs/mdg_polygons_only_obspred_map_log.png')
autoplot(cv1_output2, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_polygons_only_obspred.png')
save(cv1_output2, file = 'model_outputs/mdg_polygon_cv_1.RData')

cat('Start cv1 model3')
arg_list[c('use_polygons', 'use_points')] <- c(1, 1)
cv1_output3 <- run_cv(data_cv1_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      parallel_delay = delay, cores = 10, 
                      drop_covs = 9)
obspred_map(data_cv1_mdg, cv1_output3, column = FALSE)
ggsave('figs/mdg_both_obspred_map.png')
obspred_map(data_cv1_mdg, cv1_output3, trans = 'log10', column = FALSE)
ggsave('figs/mdg_both_obspred_map_log.png')
autoplot(cv1_output3, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_both_obspred.png')

save(cv1_output3, file = 'model_outputs/mdg_joint_cv_1.RData')




cat('Start cv1 model 4')
arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv1_output4 <- run_cv(data_cv1_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      parallel_delay = delay, cores = 10)
obspred_map(data_cv1_mdg, cv1_output4, column = FALSE)
ggsave('figs/mdg_pr_gp_only_obspred_map.png')
obspred_map(data_cv1_mdg, cv1_output4, trans = 'log10', column = FALSE)
ggsave('figs/mdg_pr_gp__only_obspred_map_log.png')
autoplot(cv1_output4, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_pr_gp__only_obspred.png')
save(cv1_output4, file = 'model_outputs/mdg_pr_gp__cv_1.RData')

cv1_output1$summary$polygon_metrics
cv1_output2$summary$polygon_metrics
cv1_output3$summary$polygon_metrics
cv1_output4$summary$polygon_metrics

cv1_output1$summary$pr_metrics
cv1_output2$summary$pr_metrics
cv1_output3$summary$pr_metrics








# Run 3 x models with 3 x hyperpars on cv3 Spatial.

cat('Start cv3')
data_cv3_mdg <- cv_spatial_folds(data_mdg, k = 3, keep_pr = TRUE)
save(data_cv3_mdg, file = 'model_outputs/mdg_cv_3.RData')
autoplot(data_cv3_mdg, jitter = 0.0)
ggsave('figs/mdg_cv_spatial2.png')


cat('Start cv3 model1')
# Run 3 x models with 3 x hyperpars on cv1.
arg_list[c('use_polygons', 'use_points')] <- c(0, 1)
cv3_output1 <- run_cv(data_cv3_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      parallel_delay = delay,
                      drop_covs = 9)
obspred_map(data_cv3_mdg, cv3_output1, column = FALSE)
ggsave('figs/mdg_points_only_obspred_map3.png')
obspred_map(data_cv3_mdg, cv3_output1, trans = 'log10', column = FALSE)
ggsave('figs/mdg_points_only_obspred_map_log3.png')
autoplot(cv3_output1, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_points_only_obspred3.png')
save(cv3_output1, file = 'model_outputs/mdg_points_cv_3.RData')

cat('Start cv3 model2')
arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv3_output2 <- run_cv(data_cv3_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      parallel_delay = delay,
                      drop_covs = 9)
obspred_map(data_cv3_mdg, cv3_output2, column = FALSE)
ggsave('figs/mdg_polygons_only_obspred_map3.png')
obspred_map(data_cv3_mdg, cv3_output2, trans = 'log10', column = FALSE)
ggsave('figs/mdg_polygons_only_obspred_map_log3.png')
autoplot(cv3_output2, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_polygons_only_obspred3.png')
save(cv3_output2, file = 'model_outputs/mdg_polygon_cv_3.RData')

cat('Start cv3 model3')
arg_list[c('use_polygons', 'use_points')] <- c(1, 1)
cv3_output3 <- run_cv(data_cv3_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      parallel_delay = delay,
                      drop_covs = 9)
obspred_map(data_cv3_mdg, cv3_output3, column = FALSE)
ggsave('figs/mdg_both_obspred_map3.png')
obspred_map(data_cv3_mdg, cv3_output3, trans = 'log10', column = FALSE)
ggsave('figs/mdg_both_obspred_map_log3.png')
autoplot(cv3_output3, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_both_obspred3.png')

save(cv3_output3, file = 'model_outputs/mdg_joint_cv_3.RData')


cat('Start cv3 model4')
arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv3_output4 <- run_cv(data_cv3_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      parallel_delay = delay)
obspred_map(data_cv3_mdg, cv3_output4, column = FALSE)
ggsave('figs/mdg_pr_gps_only_obspred_map3.png')
obspred_map(data_cv3_mdg, cv3_output4, trans = 'log10', column = FALSE)
ggsave('figs/mdg_pr_gp_only_obspred_map_log3.png')
autoplot(cv3_output4, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_pr_gp_only_obspred3.png')
save(cv3_output4, file = 'model_outputs/mdg_pr_gp_cv_3.RData')


# Not sure these are run correctly.
# No drop_cov args.
# Wrong number of parameters in outputs.


cv3_output1$summary$polygon_metrics
cv3_output2$summary$polygon_metrics
cv3_output3$summary$polygon_metrics
cv3_output4$summary$polygon_metrics

cv3_output1$summary$pr_metrics
cv3_output2$summary$pr_metrics
cv3_output3$summary$pr_metrics



print('Finished')
