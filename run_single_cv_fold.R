args <- commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])
cv_type <- args[2]
model_type <- args[3]

if(Sys.info()["user"] != 'anita'){
  setwd('~/timz/timothy/point_polygon_joint_comparison')
} else {
  setwd('~/Z/users/anita/point_polygon_join_comparison_analysis')
}

source("setUserInfo.R")

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

set.seed(i + 20)

# load all data

if(cv_type == 'random'){
  cv_data <- get(load('model_outputs/idn_cv_1.RData'))
} else {
  cv_data <- get(load('model_outputs/idn_cv_2.RData'))
}

load('model_outputs/idn_mesh.RData')


if(model_type == 'points'){
  use_points <- 1
  use_polygons <- 0
} else if(model_type == 'polygons') {
  use_points <- 0
  use_polygons <- 1
} else if(model_type == 'joint'){
  use_points <- 1
  use_polygons <- 1
} else {
  stop('Wrong type of model requested')
}


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


data <- cv_data[[i]]$train


model <- fit_model(data, mesh_idn, its = 400, model.args = arg_list)
#save(full_model, file = 'model_outputs/full_model_idn.RData')


results <- cv_performance(model$predictions, 
                          cv_data[[i]]$test, 
                          model$model,
                          use_points = arg_list$use_points,
                          CI = 0.8)


#summary <- summarise_cv_results(list(results))
out <- list(models = model, results = results)

save(out, file = paste0('model_outputs/idn-', cv_type, '-', model_type, '-', i, '.RData'))











