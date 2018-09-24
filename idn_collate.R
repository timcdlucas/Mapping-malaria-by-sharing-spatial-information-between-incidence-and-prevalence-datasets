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
source('collate_seperate_cv_runs.R')

# Compile the model
compile("joint_model.cpp")

set.seed(180530)

arg_list <- list()

# Load cv data

load('model_outputs/idn_cv_1.RData')
load('model_outputs/idn_cv_2.RData')





# Run 3 x models on cv1.
cat('Start cv1 model 1')


arg_list[c('use_polygons', 'use_points')] <- c(0, 1)


files1 <- paste0('model_outputs/idn-random-points-', 1:10, '.RData')
cv1_output1_list <- lapply(files1, function(x) get(load(x)))
cv1_output1 <- collate(cv1_output1_list)

# cv1_output1 <- run_cv(data_cv1_idn, mesh_idn, its = 1000, 
#                       model.args = arg_list, CI = 0.8, parallel_delay = 0, cores = 1)
obspred_map(data_cv1_idn, cv1_output1, column = TRUE)
ggsave('figs/idn_points_only_obspred_map.png')
obspred_map(data_cv1_idn, cv1_output1, trans = 'log10', column = TRUE)
ggsave('figs/idn_points_only_obspred_map_log.png')
autoplot(cv1_output1, type = 'obs_preds', CI = TRUE)
ggsave('figs/idn_points_only_obspred.png')

cat('Start cv1 model 2')

files2 <- paste0('model_outputs/idn-random-polygons-', 1:10, '.RData')
cv1_output2_list <- lapply(files2, function(x) get(load(x)))
cv1_output2 <- collate(cv1_output2_list)


obspred_map(data_cv1_idn, cv1_output2, column = TRUE)
ggsave('figs/idn_polygons_only_obspred_map.png')
obspred_map(data_cv1_idn, cv1_output2, trans = 'log10', column = TRUE)
ggsave('figs/idn_polygons_only_obspred_map_log.png')
autoplot(cv1_output2, type = 'obs_preds', CI = TRUE)
ggsave('figs/idn_polygons_only_obspred.png')

cat('Start cv1 model 3')

files3 <- paste0('model_outputs/idn-random-joint-', 1:10, '.RData')
cv1_output3_list <- lapply(files3, function(x) get(load(x)))
cv1_output3 <- collate(cv1_output3_list)

obspred_map(data_cv1_idn, cv1_output3, column = TRUE)
ggsave('figs/idn_joint_obspred_map.png')
obspred_map(data_cv1_idn, cv1_output3, trans = 'log10', column = TRUE)
ggsave('figs/idn_joint_obspred_map_log.png')
autoplot(cv1_output3, type = 'obs_preds', CI = TRUE)
ggsave('figs/idn_joint_obspred.png')


save(cv1_output1, file = 'model_outputs/idn_points_cv_1.RData')
save(cv1_output2, file = 'model_outputs/idn_polygon_cv_1.RData')
save(cv1_output3, file = 'model_outputs/idn_joint_cv_1.RData')

cv1_output1$summary$polygon_metrics
cv1_output2$summary$polygon_metrics
cv1_output3$summary$polygon_metrics

cv1_output1$summary$pr_metrics
cv1_output2$summary$pr_metrics
cv1_output3$summary$pr_metrics




# Run 3 x models on cv2.
cat('Start cv2 model 1')

files4 <- paste0('model_outputs/idn-spatial-points-', 1:7, '.RData')
cv2_output1_list <- lapply(files4, function(x) get(load(x)))
cv2_output1 <- collate(cv2_output1_list)

obspred_map(data_cv2_idn, cv2_output1, column = FALSE)
ggsave('figs/idn_points_only_obspred_map2.png')
obspred_map(data_cv2_idn, cv2_output1, trans = 'log10', column = FALSE)
ggsave('figs/idn_points_only_obspred_map_log2.png')
autoplot(cv2_output1, type = 'obs_preds', CI = TRUE)
ggsave('figs/idn_points_only_obspred2.png')

cat('Start cv2 model 2')

files5 <- paste0('model_outputs/idn-spatial-polygons-', 1:7, '.RData')
cv2_output2_list <- lapply(files5, function(x) get(load(x)))
cv2_output2 <- collate(cv2_output2_list)

obspred_map(data_cv2_idn, cv2_output2, column = FALSE)
ggsave('figs/idn_polygons_only_obspred_map2.png')
obspred_map(data_cv2_idn, cv2_output2, trans = 'log10', column = FALSE)
ggsave('figs/idn_polygons_only_obspred_map_log2.png')
autoplot(cv2_output2, type = 'obs_preds', CI = TRUE)
ggsave('figs/idn_polygons_only_obspred2.png')

cat('Start cv2 model 3')

files6 <- paste0('model_outputs/idn-spatial-joint-', 1:7, '.RData')
cv2_output3_list <- lapply(files6, function(x) get(load(x)))
cv2_output3 <- collate(cv2_output3_list)

obspred_map(data_cv2_idn, cv2_output3, column = FALSE)
ggsave('figs/idn_joint_obspred_map2.png')
obspred_map(data_cv2_idn, cv2_output3, trans = 'log10', column = FALSE)
ggsave('figs/idn_joint_obspred_map_log2.png')
autoplot(cv2_output3, type = 'obs_preds', CI = TRUE)
ggsave('figs/idn_joint_obspred2.png')


save(cv2_output1, file = 'model_outputs/idn_points_cv_2.RData')
save(cv2_output2, file = 'model_outputs/idn_polygon_cv_2.RData')
save(cv2_output3, file = 'model_outputs/idn_joint_cv_2.RData')

cv2_output1$summary$polygon_metrics
cv2_output2$summary$polygon_metrics
cv2_output3$summary$polygon_metrics

cv2_output1$summary$pr_metrics
cv2_output2$summary$pr_metrics
cv2_output3$summary$pr_metrics







