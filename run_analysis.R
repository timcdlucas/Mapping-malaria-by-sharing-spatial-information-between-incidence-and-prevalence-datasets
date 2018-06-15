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
pop_path <- '~/timz/GBD2017/Processing/Stages/03_Muster_Population_Figures/Verified_Outputs/Output_Pop_At_Risk_Pf_5K/ihme_corrected_frankenpop_All_Ages_3_2015_at_risk_pf.tif'
shapefile_path = '~/timz/master_geometries/Admin_Units/Global/GBD/GBD2017_MAP/GBD2017_MAP_MG_5K/'

cov_raster_paths <- c(
  '~/timz/mastergrids/MODIS_Global/MOD11A2_LST/LST_Day/5km/Annual/LST_Day.2015.Annual.mean.5km.Mean.tif',
  '~/timz/mastergrids/MODIS_Global/MCD43B4_BRDF_Reflectance/EVI/5km/Annual/EVI_Overall_Mean_0.tif',
  '~/timz/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/TSI-Martens2-Pf.Synoptic.01.Min.5km.Data.tif'
  )

# load packages

## Spatial packages
library(raster)
library(maptools)

## dataframe packages
library(dplyr)
library(readr)
library(magrittr)

## For standardising prevalence
library(malariaAtlas)

## For inc prevalence conversions
library(GBDutils)

## Plotting packages
library(ggplot2)


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
                  api_year = 2012)


# indonesia

# pre analysis

data_idn <- process_data(
  binomial_positive = data$pr$positive,
  binomial_n = data$pr$examined,
  coords = data$pr[, c('latitude', 'longitude')],
  polygon_response = data$api$api_mean,
  polygon_population = data$api$population,
  shapefile_id = data$api$shapefile_id,
  shps_id_column = 'area_id',
  shapefiles = data$shapefiles,
  pop_raster = data$pop,
  cov_rasters = data$covs)


mesh_idn <- build_mesh(data_idn, mesh.args = list(max.edge = c(0.8, 5), cut = 0.8))


data_cv1_idn <- cv_folds(data_idn)
data_cv2_idn <- cv_spat_folds(data_idn)


# run models

arg_list <- NULL
full_model <- fit_models(data_idn, mesh_idn, model.args = arg_list)

cv1_output <- run_cv(data_cv1_idn, mesh, model.args = arg_list)

cv2_output <- run_cv(data_cv2_idn, mesh, model.args = arg_list)



# create figures
