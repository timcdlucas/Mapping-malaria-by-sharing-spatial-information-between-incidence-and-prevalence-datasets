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

source('plotting_functions.R')
source('run_cv.R')
source('collate_seperate_cv_runs.R')

set.seed(180530)

arg_list <- list()

# Load cv data

load('model_outputs/idn_cv_1.RData')
load('model_outputs/idn_cv_2.RData')
load('model_outputs/idn_cv_3.RData')

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









# Run 3 x models on cv3.

cat('Start cv3 model 2')

files5 <- paste0('model_outputs/idn-spatialkeeppr-polygons-', 1:7, '.RData')
cv3_output2_list <- lapply(files5, function(x) get(load(x)))
cv3_output2 <- collate(cv3_output2_list)

obspred_map(data_cv3_idn, cv3_output2, column = FALSE)
ggsave('figs/idn_polygons_only_obspred_map3.png')
obspred_map(data_cv3_idn, cv3_output2, trans = 'log10', column = FALSE)
ggsave('figs/idn_polygons_only_obspred_map_log3.png')
autoplot(cv3_output2, type = 'obs_preds', CI = TRUE)
ggsave('figs/idn_polygons_only_obspred3.png')

cat('Start cv3 model 3')

files6 <- paste0('model_outputs/idn-spatialkeeppr-joint-', 1:7, '.RData')
cv3_output3_list <- lapply(files6, function(x) get(load(x)))
cv3_output3 <- collate(cv3_output3_list)

obspred_map(data_cv3_idn, cv3_output3, column = FALSE)
ggsave('figs/idn_joint_obspred_map3.png')
obspred_map(data_cv3_idn, cv3_output3, trans = 'log10', column = FALSE)
ggsave('figs/idn_joint_obspred_map_log3.png')
autoplot(cv3_output3, type = 'obs_preds', CI = TRUE)
ggsave('figs/idn_joint_obspred3.png')



save(cv3_output2, file = 'model_outputs/idn_polygon_cv_3.RData')
save(cv3_output3, file = 'model_outputs/idn_joint_cv_3.RData')

cv3_output2$summary$polygon_metrics
cv3_output3$summary$polygon_metrics






