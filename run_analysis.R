########
# master script for point Vs polygon Vs joint analysis
# Tim Lucas
# 2018-05-30
###########


# setwd('~/Dropbox/work/analysis/point_polygon_joint_comparison')
# setwd('~/timz/timothy/point_polygon_joint_comparison')

# define paths

PR_path <- '~/timz/GBD2017/Processing/Stages/04b_PR_DB_Import_Export/Verified_Outputs/2018_02_15/pfpr_dump.csv'
API_path <- '~/timz/GBD2017/Processing/Stages/04c_API_Data_Export/Checkpoint_Outputs/subnational.csv'
pop_path <- '~/timz/GBD2017/Processing/Stages/03_Muster_Population_Figures/Verified_Outputs/Output_Pop_At_Risk_Pf_5K/ihme_corrected_frankenpop_All_Ages_3_2012_at_risk_pf.tif'
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
cl <- makeCluster(min(detectCores() - 1, 20))
registerDoParallel(cl)


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
                  useiso3 = 'IDN', 
                  pr_year = 2008,
                  api_year = 2014)
                  # pr_year = 2010,
                  # api_year = 2012)


# indonesia

# pre analysis

data_idn <- process_data(
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

autoplot(data_idn, pr_limits = c(0, 0.3))
ggsave('figs/idn_input_data.png')

mesh_idn <- build_mesh(data_idn, mesh.args = list(max.edge = c(0.5, 5), cut = 0.5))
autoplot(mesh_idn)

data_cv1_idn <- cv_random_folds(data_idn, k = 10)
autoplot(data_cv1_idn, jitter = 0.7)
ggsave('figs/idn_cv_random.png')

#autoplot(data_cv1_idn[[1]]$train, pr_limits = c(0, 0.3))


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
                 prior_iideffect_sd_prob = 0.0000001, # Made this stronger because too much iid.
                 prior_iideffect_pr_sd_max = 0.3, # Max difference between PR points within a cell (with n > 500)
                 prior_iideffect_pr_sd_prob = 0.0000001,
                 priormean_intercept = -2,
                 priorsd_intercept = 2,  # Indonesia has prev lowish. But want intercept to take whatever value it likes.
                 priormean_slope = 0, 
                 priorsd_slope = 0.4, # Explains between 0.004 and 0.27 prevalence. 1 covariate shouldn't explain between 0 and 0.6 (range of prev).
                 use_polygons = 0,
                 # use_polygons = 1,
                 use_points = 1)

full_model <- fit_model(data_idn, mesh_idn, its = 500, model.args = arg_list)
autoplot(full_model)
plot(full_model, layer = 'api')

in_sample <- cv_performance(predictions = full_model$predictions, 
                            holdout = data_idn,
                            model_params = full_model$model, 
                            CI = 0.8)
autoplot(in_sample, CI = TRUE)
autoplot(in_sample, trans = 'log1p', CI = TRUE)


# Run 3 x models with 3 x hyperpars on cv1.
arg_list[c('use_polygons', 'use_points')] <- c(0, 1)
cv1_output1 <- run_cv(data_cv1_idn, mesh_idn, its = 200, model.args = arg_list)

autoplot(cv1_output1, type = 'obs_preds', trans = 'log1p')
obspred_map(data_cv1_idn, cv1_output1)
obspred_map(data_cv1_idn, cv1_output1, trans = 'log10', lims = c(1e-5, 620))

arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv1_output2 <- run_cv(data_cv1_idn, mesh_idn, its = 700, model.args = arg_list)
obspred_map(data_cv1_idn, cv1_output2)
obspred_map(data_cv1_idn, cv1_output2, trans = 'log10', lims = c(1e-5, 620))


arg_list[c('use_polygons', 'use_points')] <- c(1, 1)
cv1_output3 <- run_cv(data_cv1_idn, mesh_idn, its = 200, model.args = arg_list)
obspred_map(data_cv1_idn, cv1_output3)
obspred_map(data_cv1_idn, cv1_output3, trans = 'log10', lims = c(1e-5, 620))

save(cv1_output1, file = 'model_outputs/points_cv_1.RData')
save(cv1_output2, file = 'model_outputs/polygon_cv_1.RData')
save(cv1_output3, file = 'model_outputs/join_cv_1.RData')

cv1_output1$summary$polygon_metrics
cv1_output2$summary$polygon_metrics
cv1_output3$summary$polygon_metrics

cv1_output1$summary$pr_metrics
cv1_output2$summary$pr_metrics
cv1_output3$summary$pr_metrics



data_cv1_idn <- cv_spatial_folds(data_idn, k = 7)
autoplot(data_cv1_idn, jitter = 0.7)
ggsave('figs/idn_cv_spatial.png')





# Run 3 x models with 3 x hyperpars on cv2


#data_cv2_idn <- cv_spat_folds(data_idn)

# cv1_model <- fit_model(data_cv1_idn[[1]]$train, mesh_idn, model.args = arg_list)
# cv1_test <- cv_performance(predictions = cv1_model$predictions, 
#                             holdout = data_cv1_idn[[1]]$test)




#cv2_output <- run_cv(data_cv2_idn, mesh, model.args = arg_list)





# Choose best hyperpar for each model, for each CV and collate.









# create temp figures




# Write out data needed for final figures.





















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



# Run 3 x models with 3 x hyperpars on cv2


#data_cv2_mdg <- cv_spat_folds(data_mdg)

# cv1_model <- fit_model(data_cv1_mdg[[1]]$train, mesh_mdg, model.args = arg_list)
# cv1_test <- cv_performance(predictions = cv1_model$predictions, 
#                             holdout = data_cv1_mdg[[1]]$test)




#cv2_output <- run_cv(data_cv2_mdg, mesh, model.args = arg_list)





# Choose best hyperpar for each model, for each CV and collate.









# create temp figures




# Write out data needed for final figures.

