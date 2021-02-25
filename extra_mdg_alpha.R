
########
# master script for point Vs polygon Vs joint analysis
# Tim Lucas
# 2018-05-30
###########

# define paths


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


load('model_outputs/mdg_full_data.RData')

mesh_mdg <- build_mesh(data_mdg, mesh.args = list(max.edge = c(0.2, 3), cut = 0.2, offset = c(2, 5)))


load('model_outputs/mdg_cv_1.RData')



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
                 priorsd_log_prev_inc_extra_slope = 1,
                 use_points = 1)


delay <- 2

cat('Start cv1 model3')
arg_list[c('use_polygons', 'use_points')] <- c(1, 1)
cv1_output3 <- run_cv(data_cv1_mdg, mesh_mdg, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      parallel_delay = delay, cores = 5, drop_covs = 9)
obspred_map(data_cv1_mdg, cv1_output3, column = FALSE)
ggsave('figs/mdg_both_obspred_map_new_priors.png')
obspred_map(data_cv1_mdg, cv1_output3, trans = 'log10', column = FALSE)
ggsave('figs/mdg_both_obspred_map_log_new_priors.png')
autoplot(cv1_output3, type = 'obs_preds', CI = TRUE)
ggsave('figs/mdg_both_obspred_new_priors.png')

save(cv1_output3, file = 'model_outputs/mdg_joint_cv_1_testalphaprior.RData')

cv1_output3_new <- cv1_output3




load('model_outputs/mdg_joint_cv_1_oldpriors.RData')


cv1_output3_old <- cv1_output3


cv1_output3_old$summary$polygon_metrics
# MAE 39.1
cv1_output3_new$summary$polygon_metrics
# MAE 39.1



alphaold <- sapply(cv1_output3_old$models, function(x) x$model$opt$par['log_prev_inc_extra_slope'])
alphanew <- sapply(cv1_output3_new$models, function(x) x$model$opt$par['log_prev_inc_extra_slope'])
png('figs/robust_prev_inc_extra_mdg.pdf')
plot(alphaold, alphanew, main = 'spatial range')
abline(0,1)
dev.off()





