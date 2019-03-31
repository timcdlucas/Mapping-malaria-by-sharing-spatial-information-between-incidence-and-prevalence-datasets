
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
# API_path <- Z('GBD2017/Processing/Stages/04c_API_Data_Export/Checkpoint_Outputs/subnational.csv')
API_path <- Z('GBD2019/Processing/Stages/04c_API_Data_Export/Checkpoint_Outputs/Decomp 4/subnational.csv')
#API_path <- Z('GBD2017/Processing/Stages/04c_API_Data_Export/Checkpoint_Outputs/2018-10-17_Senegal/subnational.csv')
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
  Z('mastergrids/Other_Global_Covariates/NightTimeLights/VIIRS_DNB_Monthly/5km/Annual/VIIRS-SLC.2016.Annual.5km.MEDIAN.tif'),
  #Z('mastergrids/Other_Global_Covariates/UrbanAreas/Global_Urban_Footprint/From_86m/5km/Global_Urban_Footprint_5km_PropUrban.tif'),
  Z('mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/TCW_v6/5km/Synoptic/TCW_v6.Synoptic.Overall.mean.5km.mean.tif')
)


#which_vars <- c(1:20)[-c(19, 12, 15, 6, 14, 7, 5)] 
#which_vars <- c(2, 5, 6, 8, 10, 20)
#cov_raster_paths <- paste0('~/timz/mastergrids/Other_Global_Covariates/Pf_Covariates/PF_V', 1:20, '/5km/Monthly/PF_V', 1:20, '.2010.05.Data.5km.Data.tif')#[which_vars]


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
                  useiso3 = 'ZMB', 
                  admin_unit_level = 'ADMIN2', # todo
                  shapefile_column = 'area_id',
                  pr_year = c(2008:2012),
                  api_year = 2010)

#pr_year = 2013,
#api_year = 2013)

# madagascar

# pre analysis

data_zam <- process_data(
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
  useiso3 = 'ZMB',
  transform = 4:7)
save(data_zam, file = 'model_outputs/zam_full_data.RData')

autoplot(data_zam)

mesh_zam <- build_mesh(data_zam, mesh.args = list(max.edge = c(0.1, 0.8), cut = 0.1, offset = c(2, 2)))

data_cv1_zam <- cv_random_folds(data_zam, k = 10) # todo
autoplot(data_cv1_zam, jitter = 0)
save(data_cv1_zam, file = 'model_outputs/zam_cv_1.RData')



# run models
# Run full model to get a handle on things.

arg_list <- list(prior_rho_min = 1, # 
                 prior_rho_prob = 0.00001, # Want p(rho < 3) = 0.0001 -> p(log_kappa < -0.058) = 0.0001
                 prior_sigma_max = 0.5, # Want p(sd > 1) = 0.0001 (would explain most of prev).
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
                 # use_polygons = 1,
                 use_points = 1)

if(FALSE){
  full_model <- fit_model(data_zam, mesh_zam, its = 1000, model.args = arg_list)
  autoplot(full_model)
  
  png('figs/zam_full_model_in_sample_map.png')
  plot(full_model, layer = 'api')
  dev.off()
  
  in_sample <- cv_performance(predictions = full_model$predictions, 
                              holdout = data_zam,
                              model_params = full_model$model, 
                              CI = 0.8,
                              use_points = use_points)
  autoplot(in_sample, CI = TRUE)
  autoplot(in_sample, trans = 'log1p', CI = TRUE)
  ggsave('figs/zam_full_model_in_sample.png')
  
  save(full_model, file = 'model_outputs/full_model_zam.RData')
  
  
  
  arg_list[c('use_polygons', 'use_points')] <- c(0, 1)
  points_model <- fit_model(data_zam, mesh_zam, its = 1000, model.args = arg_list)
  autoplot(points_model)
  png('figs/zam_points_model_in_sample_map.png')
  plot(points_model, layer = 'api')
  dev.off()
  
  points_in_sample <- cv_performance(predictions = points_model$predictions, 
                                     holdout = data_zam,
                                     model_params = points_model$model, 
                                     CI = 0.8,
                                     use_points = TRUE)
  autoplot(points_in_sample, CI = TRUE)
  autoplot(points_in_sample, trans = 'log1p', CI = TRUE)
  ggsave('figs/zam_points_model_in_sample.png')
  
  
  save(points_model, file = 'model_outputs/points_model_zam.RData')
  
  
  
  
  
  arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
  polygons_model <- fit_model(data_zam, mesh_zam, its = 1000, model.args = arg_list)
  autoplot(polygons_model)
  png('figs/zam_polygons_model_in_sample_map.png')
  plot(polygons_model, layer = 'api')
  dev.off()
  
  polygons_in_sample <- cv_performance(predictions = polygons_model$predictions, 
                                       holdout = data_zam,
                                       model_params = polygons_model$model, 
                                       CI = 0.8,
                                       use_points = FALSE)
  autoplot(polygons_in_sample, CI = TRUE)
  autoplot(polygons_in_sample, trans = 'log1p', CI = TRUE)
  ggsave('figs/zam_polygon_model_in_sample.png')
  
  save(polygons_model, file = 'model_outputs/polygons_model_zam.RData')
  
  
  
}

delay <- 10

cat('Start cv1 model 1')
# Run 3 x models with 3 x hyperpars on cv1.
arg_list[c('use_polygons', 'use_points')] <- c(0, 1)
cv1_output1 <- run_cv(data_cv1_zam, mesh_zam, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      cores = 10, parallel_delay = delay)
obspred_map(data_cv1_zam, cv1_output1, column = FALSE)
ggsave('figs/zam_points_only_obspred_map.png')
obspred_map(data_cv1_zam, cv1_output1, trans = 'log10', column = FALSE)
ggsave('figs/zam_points_only_obspred_map_log.png')
autoplot(cv1_output1, type = 'obs_preds', CI = TRUE)
ggsave('figs/zam_points_only_obspred.png')

save(cv1_output1, file = 'model_outputs/zam_points_cv_1.RData')

cat('Start cv1 model 2')
arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv1_output2 <- run_cv(data_cv1_zam, mesh_zam, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      cores = 10, parallel_delay = delay)
obspred_map(data_cv1_zam, cv1_output2, column = FALSE, mask = TRUE)
ggsave('figs/zam_polygons_only_obspred_map.png')
obspred_map(data_cv1_zam, cv1_output2, trans = 'log10', column = FALSE)
ggsave('figs/zam_polygons_only_obspred_map_log.png')
autoplot(cv1_output2, type = 'obs_preds', CI = FALSE)
ggsave('figs/zam_polygons_only_obspred.png')
save(cv1_output2, file = 'model_outputs/zam_polygon_cv_1.RData')

cat('Start cv1 model3')
arg_list[c('use_polygons', 'use_points')] <- c(1, 1)
cv1_output3 <- run_cv(data_cv1_zam, mesh_zam, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      cores = 10, parallel_delay = delay)
obspred_map(data_cv1_zam, cv1_output3, column = FALSE)
ggsave('figs/zam_both_obspred_map.png')
obspred_map(data_cv1_zam, cv1_output3, trans = 'log10', column = FALSE)
ggsave('figs/zam_both_obspred_map_log.png')
autoplot(cv1_output3, type = 'obs_preds', CI = TRUE)
ggsave('figs/zam_both_obspred.png')
save(cv1_output3, file = 'model_outputs/zam_joint_cv_1.RData')






cv1_output1$summary$polygon_metrics
cv1_output2$summary$polygon_metrics
cv1_output3$summary$polygon_metrics

cv1_output1$summary$pr_metrics
cv1_output2$summary$pr_metrics
cv1_output3$summary$pr_metrics




# Run 3 x models with 3 x hyperpars on cv2 Spatial.

cat('Start cv2')
data_cv2_zam <- cv_spatial_folds(data_zam, k = 5)
save(data_cv2_zam, file = 'model_outputs/zam_cv_2.RData')
autoplot(data_cv2_zam, jitter = 0.0)
ggsave('figs/zam_cv_spatial.png')


cat('Start cv2 model1')
# Run 3 x models with 3 x hyperpars on cv1.
arg_list[c('use_polygons', 'use_points')] <- c(0, 1)
cv2_output1 <- run_cv(data_cv2_zam, mesh_zam, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = delay, cores = 5)
obspred_map(data_cv2_zam, cv2_output1, column = FALSE)
ggsave('figs/zam_points_only_obspred_map2.png')
obspred_map(data_cv2_zam, cv2_output1, trans = 'log10', column = FALSE)
ggsave('figs/zam_points_only_obspred_map_log2.png')
autoplot(cv2_output1, type = 'obs_preds', CI = TRUE)
ggsave('figs/zam_points_only_obspred2.png')
save(cv2_output1, file = 'model_outputs/zam_points_cv_2.RData')

cat('Start cv2 model2')
arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv2_output2 <- run_cv(data_cv2_zam, mesh_zam, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = delay, cores = 5)
obspred_map(data_cv2_zam, cv2_output2, column = FALSE)
ggsave('figs/zam_polygons_only_obspred_map2.png')
obspred_map(data_cv2_zam, cv2_output2, trans = 'log10', column = FALSE)
ggsave('figs/zam_polygons_only_obspred_map_log2.png')
autoplot(cv2_output2, type = 'obs_preds', CI = FALSE)
ggsave('figs/zam_polygons_only_obspred2.png')
save(cv2_output2, file = 'model_outputs/zam_polygon_cv_2.RData')

cat('Start cv2 model3')
arg_list[c('use_polygons', 'use_points')] <- c(1, 1)
cv2_output3 <- run_cv(data_cv2_zam, mesh_zam, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = delay, cores = 5)
obspred_map(data_cv2_zam, cv2_output3, column = FALSE)
ggsave('figs/zam_both_obspred_map2.png')
obspred_map(data_cv2_zam, cv2_output3, trans = 'log10', column = FALSE)
ggsave('figs/zam_both_obspred_map_log2.png')
autoplot(cv2_output3, type = 'obs_preds', CI = FALSE)
ggsave('figs/zam_both_obspred2.png')

save(cv2_output3, file = 'model_outputs/zam_joint_cv_2.RData')



cv2_output1$summary$polygon_metrics
cv2_output2$summary$polygon_metrics
cv2_output3$summary$polygon_metrics

cv2_output1$summary$pr_metrics
cv2_output2$summary$pr_metrics
cv2_output3$summary$pr_metrics






# Run 3 x models with 3 x hyperpars on cv3 Spatial.

cat('Start cv3')
data_cv3_zam <- cv_spatial_folds(data_zam, k = 5, keep_pr = TRUE)
save(data_cv3_zam, file = 'model_outputs/zam_cv_3.RData')
autoplot(data_cv3_zam, jitter = 0.0)
ggsave('figs/zam_cv_spatial.png')


cat('Start cv3 model1')
# Run 3 x models with 3 x hyperpars on cv1.
arg_list[c('use_polygons', 'use_points')] <- c(0, 1)
cv3_output1 <- run_cv(data_cv3_zam, mesh_zam, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = delay, cores = 5)
obspred_map(data_cv3_zam, cv3_output1, column = FALSE)
ggsave('figs/zam_points_only_obspred_map3,png')
obspred_map(data_cv3_zam, cv3_output1, trans = 'log10', column = FALSE)
ggsave('figs/zam_points_only_obspred_map_log3,png')
autoplot(cv3_output1, type = 'obs_preds', CI = TRUE)
ggsave('figs/zam_points_only_obspred3,png')
save(cv3_output1, file = 'model_outputs/zam_points_cv_3.RData')

cat('Start cv3 model2')
arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv3_output2 <- run_cv(data_cv3_zam, mesh_zam, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = delay, cores = 5)
obspred_map(data_cv3_zam, cv3_output2, column = FALSE)
ggsave('figs/zam_polygons_only_obspred_map3,png')
obspred_map(data_cv3_zam, cv3_output2, trans = 'log10', column = FALSE)
ggsave('figs/zam_polygons_only_obspred_map_log3,png')
autoplot(cv3_output2, type = 'obs_preds', CI = FALSE)
ggsave('figs/zam_polygons_only_obspred3,png')
save(cv3_output2, file = 'model_outputs/zam_polygon_cv_3.RData')

cat('Start cv3 model3')
arg_list[c('use_polygons', 'use_points')] <- c(1, 1)
cv3_output3 <- run_cv(data_cv3_zam, mesh_zam, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = delay, cores = 5)
obspred_map(data_cv3_zam, cv3_output3, column = FALSE)
ggsave('figs/zam_both_obspred_map3,png')
obspred_map(data_cv3_zam, cv3_output3, trans = 'log10', column = FALSE)
ggsave('figs/zam_both_obspred_map_log3,png')
autoplot(cv3_output3, type = 'obs_preds', CI = FALSE)
ggsave('figs/zam_both_obspred3,png')

save(cv3_output3, file = 'model_outputs/zam_joint_cv_3.RData')



cv3_output1$summary$polygon_metrics
cv3_output2$summary$polygon_metrics
cv3_output3$summary$polygon_metrics

cv3_output1$summary$pr_metrics
cv3_output2$summary$pr_metrics
cv3_output3$summary$pr_metrics


print('Finished')
