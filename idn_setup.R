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
  Z('mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/5km/Synoptic/LST_Day.Synoptic.Overall.mean.5km.mean.tif'),
  #Z('mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/EVI/5km/Synoptic/EVI.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/TSI-Martens2-Pf.Synoptic.Overall.Mean.5km.Data.tif'),
  Z('GBD2017/Processing/Static_Covariates/MAP/other_rasters/accessibility/accessibility.5k.MEAN.tif'),
  Z('mastergrids/Other_Global_Covariates/Elevation/SRTM-Elevation/5km/Synoptic/SRTM_elevation.Synoptic.Overall.Data.5km.mean.tif'),
  Z('mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/5km/Synoptic/LST_Day.Synoptic.Overall.SD.5km.mean.tif'),
  #Z('mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCB/5km/Synoptic/TCB.Synoptic.Overall.mean.5km.mean.tif'),
  Z('mastergrids/Other_Global_Covariates/NightTimeLights/VIIRS_DNB_Monthly/5km/Annual/VIIRS-SLC.2016.Annual.5km.MEDIAN.tif'),
  #Z('mastergrids/Other_Global_Covariates/UrbanAreas/Global_Urban_Footprint/From_86m/5km/Global_Urban_Footprint_5km_PropUrban.tif'),
  Z('mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/TCW/5km/Synoptic/TCW.Synoptic.Overall.mean.5km.mean.tif')
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
                  useiso3 = 'IDN', 
                  # pr_year = 2008,
                  # api_year = 2014)
                  pr_year = 2010,
                  api_year = 2012)


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
  useiso3 = 'IDN',
  transform = c(4:7))
save(data_idn, file = 'model_outputs/idn_full_data.RData')

autoplot(data_idn, pr_limits = c(0, 0.3))
ggsave('figs/idn_input_data.png')

mesh_idn <- build_mesh(data_idn, mesh.args = list(max.edge = c(0.5, 5), cut = 0.5))
autoplot(mesh_idn)
save(mesh_idn, file = 'model_outputs/idn_mesh.RData')



# Define cross validation strategies
data_cv1_idn <- cv_random_folds(data_idn, k = 10)
autoplot(data_cv1_idn, jitter = 0.7)
ggsave('figs/idn_cv_random.png')
save(data_cv1_idn, file = 'model_outputs/idn_cv_1.RData')


# Spatial
data_cv2_idn <- cv_spatial_folds(data_idn, k = 7)
autoplot(data_cv2_idn, jitter = 0.7)
ggsave('figs/idn_cv_spatial2.png')
save(data_cv2_idn, file = 'model_outputs/idn_cv_2.RData')


#autoplot(data_cv1_idn[[1]]$train, pr_limits = c(0, 0.3))

use_points <- 1
use_polygons <- 1
# run models
# Run full model to get a handle on things.

arg_list <- list(prior_rho_min = 3, # 
                 prior_rho_prob = 0.00001, # Want p(rho < 3) = 0.0001
                 prior_sigma_max = 1, # Want p(sd > 1) = 0.0001 (would explain most of prev). 
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
                 use_polygons = use_polygons,
                 use_points = use_points)



if(FALSE){
  full_model <- fit_model(data_idn, mesh_idn, its = 1000, model.args = arg_list)
  autoplot(full_model)
  
  png('figs/full_model_in_sample_map.png')
  plot(full_model, layer = 'api')
  dev.off()
  
  in_sample <- cv_performance(predictions = full_model$predictions, 
                              holdout = data_idn,
                              model_params = full_model$model, 
                              CI = 0.8,
                              use_points = use_points)
  autoplot(in_sample, CI = TRUE)
  autoplot(in_sample, trans = 'log1p', CI = TRUE)
  ggsave('figs/idn_full_model_in_sample.png')
  
  save(full_model, file = 'model_outputs/full_model_idn.RData')
  
  
  
  arg_list[c('use_polygons', 'use_points')] <- c(0, 1)
  points_model <- fit_model(data_idn, mesh_idn, its = 1000, model.args = arg_list)
  autoplot(points_model)
  png('figs/points_model_in_sample_map.png')
  plot(points_model, layer = 'api')
  dev.off()
  
  points_in_sample <- cv_performance(predictions = points_model$predictions, 
                                     holdout = data_idn,
                                     model_params = points_model$model, 
                                     CI = 0.8,
                                     use_points = use_points)
  autoplot(points_in_sample, CI = TRUE)
  autoplot(points_in_sample, trans = 'log1p', CI = TRUE)
  ggsave('figs/idn_points_model_in_sample.png')
  
  
  save(points_model, file = 'model_outputs/points_model_idn.RData')
  
  
  
  
  
  arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
  polygons_model <- fit_model(data_idn, mesh_idn, its = 1000, model.args = arg_list)
  autoplot(polygons_model)
  png('figs/polygons_model_in_sample_map.png')
  plot(polygons_model, layer = 'api')
  dev.off()
  
  polygons_in_sample <- cv_performance(predictions = polygons_model$predictions, 
                                       holdout = data_idn,
                                       model_params = polygons_model$model, 
                                       CI = 0.8,
                                       use_points = FALSE)
  autoplot(polygons_in_sample, CI = TRUE)
  autoplot(polygons_in_sample, trans = 'log1p', CI = TRUE)
  ggsave('figs/idn_polygon_model_in_sample.png')
  
  save(polygons_model, file = 'model_outputs/polygons_model_idn.RData')
  
  
  
}





