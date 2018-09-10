
########
# master script for point Vs polygon Vs joint analysis
# Tim Lucas
# 2018-05-30
###########


# setwd('~/Dropbox/work/analysis/point_polygon_joint_comparison')
setwd('~/timz/timothy/point_polygon_joint_comparison')

# define paths

PR_path <- '~/timz/GBD2017/Processing/Stages/04b_PR_DB_Import_Export/Verified_Outputs/2018_02_15/pfpr_dump.csv'
API_path <- '~/timz/GBD2017/Processing/Stages/04c_API_Data_Export/Checkpoint_Outputs/subnational.csv'
pop_path <- '~/timz/GBD2017/Processing/Stages/03_Muster_Population_Figures/Verified_Outputs/Output_Pop_At_Risk_Pf_5K/ihme_corrected_frankenpop_All_Ages_3_2015_at_risk_pf.tif'
shapefile_path = '~/timz/master_geometries/Admin_Units/Global/GBD/GBD2017_MAP/GBD2017_MAP_MG_5K/'

cov_raster_paths <- c(
  '~/timz/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/5km/Synoptic/LST_Day.Synoptic.Overall.mean.5km.mean.tif',
  '~/timz/mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/EVI/5km/Synoptic/EVI.Synoptic.Overall.mean.5km.mean.tif',
  '~/timz/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/TSI-Martens2-Pf.Synoptic.Overall.Mean.5km.Data.tif',
  '~/timz/GBD2017/Processing/Static_Covariates/MAP/other_rasters/accessibility/accessibility.5k.MEAN.tif',
  '~/timz/mastergrids/Other_Global_Covariates/Elevation/SRTM-Elevation/5km/Synoptic/SRTM_elevation.Synoptic.Overall.Data.5km.mean.tif',
  '~/timz/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/5km/Synoptic/LST_Day.Synoptic.Overall.SD.5km.mean.tif',
  #'~/timz/mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCB/5km/Synoptic/TCB.Synoptic.Overall.mean.5km.mean.tif',
  '~/timz/mastergrids/Other_Global_Covariates/NightTimeLights/VIIRS_DNB_Monthly/5km/Annual/VIIRS-SLC.2016.Annual.5km.MEDIAN.tif',
  #'~/timz/mastergrids/Other_Global_Covariates/UrbanAreas/Global_Urban_Footprint/From_86m/5km/Global_Urban_Footprint_5km_PropUrban.tif',
  '~/timz/mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCW/5km/Synoptic/TCW.Synoptic.Overall.mean.5km.mean.tif'
)

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
                  admin_unit_level = 'ADMIN3',
                  pr_year = 2016,
                  api_year = 2015)

#pr_year = 2013,
#api_year = 2013)

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
  transform = c(4:7))
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
                 # use_polygons = 1,
                 use_points = 1)
if(false){
full_model <- fit_model(data_mdg, mesh_mdg, its = 400, model.args = arg_list)
autoplot(full_model)
plot(full_model, layer = 'api')

in_sample <- cv_performance(predictions = full_model$predictions, 
                            holdout = data_mdg,
                            model_params = full_model$model,
                            CI = 0.8,
                            use_points = arg_list$use_points)
autoplot(in_sample, CI = TRUE)
autoplot(in_sample, trans = 'log1p')
}

cat('Start cv1 model 1')
# Run 3 x models with 3 x hyperpars on cv1.
arg_list[c('use_polygons', 'use_points')] <- c(0, 1)
cv1_output1 <- run_cv(data_cv1_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      cores = 10, parallel_delay = 1500)
obspred_map(data_cv1_mdg, cv1_output1, column = FALSE)
ggsave('figs/mdg_points_only_obspred_map.png')
obspred_map(data_cv1_mdg, cv1_output1, trans = 'log10', column = FALSE)
ggsave('figs/mdg_points_only_obspred_map_log.png')
autoplot(cv1_output1, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_points_only_obspred.png')



cat('Start cv1 model 2')
arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv1_output2 <- run_cv(data_cv1_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      cores = 10, parallel_delay = 1500)
obspred_map(data_cv1_mdg, cv1_output2, column = FALSE)
ggsave('figs/mdg_polygons_only_obspred_map.png')
obspred_map(data_cv1_mdg, cv1_output2, trans = 'log10', column = FALSE)
ggsave('figs/mdg_polygons_only_obspred_map_log.png')
autoplot(cv1_output2, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_polygons_only_obspred.png')


cat('Start cv1 model3')
arg_list[c('use_polygons', 'use_points')] <- c(1, 1)
cv1_output3 <- run_cv(data_cv1_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      cores = 10, parallel_delay = 1500)
obspred_map(data_cv1_mdg, cv1_output3, column = FALSE)
ggsave('figs/mdg_both_obspred_map.png')
obspred_map(data_cv1_mdg, cv1_output3, trans = 'log10', column = FALSE)
ggsave('figs/mdg_both_obspred_map_log.png')
autoplot(cv1_output3, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_both_obspred.png')

save(cv1_output1, file = 'model_outputs/mdg_points_cv_1.RData')
save(cv1_output2, file = 'model_outputs/mdg_polygon_cv_1.RData')
save(cv1_output3, file = 'model_outputs/mdg_join_cv_1.RData')

cv1_output1$summary$polygon_metrics
cv1_output2$summary$polygon_metrics
cv1_output3$summary$polygon_metrics

cv1_output1$summary$pr_metrics
cv1_output2$summary$pr_metrics
cv1_output3$summary$pr_metrics




# Run 3 x models with 3 x hyperpars on cv2 Spatial.

cat('Start cv2')
data_cv2_mdg <- cv_spatial_folds(data_mdg, k = 3)
save(data_cv2_mdg, file = 'model_outputs/mdg_cv_2.RData')
autoplot(data_cv2_mdg, jitter = 0.0)
ggsave('figs/idn_cv_spatial.png')


cat('Start cv2 model1')
# Run 3 x models with 3 x hyperpars on cv1.
arg_list[c('use_polygons', 'use_points')] <- c(0, 1)
cv2_output1 <- run_cv(data_cv2_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = 500)
obspred_map(data_cv2_mdg, cv2_output1, column = FALSE)
ggsave('figs/mdg_points_only_obspred_map2.png')
obspred_map(data_cv2_mdg, cv2_output1, trans = 'log10', column = FALSE)
ggsave('figs/mdg_points_only_obspred_map_log2.png')
autoplot(cv2_output1, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_points_only_obspred2.png')


cat('Start cv2 model2')
arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv2_output2 <- run_cv(data_cv2_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = 500)
obspred_map(data_cv2_mdg, cv2_output2, column = FALSE)
ggsave('figs/mdg_polygons_only_obspred_map2.png')
obspred_map(data_cv2_mdg, cv2_output2, trans = 'log10', column = FALSE)
ggsave('figs/mdg_polygons_only_obspred_map_log2.png')
autoplot(cv2_output2, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_polygons_only_obspred2.png')


cat('Start cv2 model3')
arg_list[c('use_polygons', 'use_points')] <- c(1, 1)
cv2_output3 <- run_cv(data_cv2_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, parallel_delay = 500)
obspred_map(data_cv2_mdg, cv2_output3, column = FALSE)
ggsave('figs/mdg_both_obspred_map2.png')
obspred_map(data_cv2_mdg, cv2_output3, trans = 'log10', column = FALSE)
ggsave('figs/mdg_both_obspred_map_log2.png')
autoplot(cv2_output3, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_both_obspred2.png')

save(cv2_output1, file = 'model_outputs/mdg_points_cv_2.RData')
save(cv2_output2, file = 'model_outputs/mdg_polygon_cv_2.RData')
save(cv2_output3, file = 'model_outputs/mdg_joint_cv_2.RData')

cv2_output1$summary$polygon_metrics
cv2_output2$summary$polygon_metrics
cv2_output3$summary$polygon_metrics

cv2_output1$summary$pr_metrics
cv2_output2$summary$pr_metrics
cv2_output3$summary$pr_metrics

print('Finished')
