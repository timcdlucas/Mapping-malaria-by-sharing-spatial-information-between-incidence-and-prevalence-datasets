
########
# master script for point Vs polygon Vs joint analysis
# Tim Lucas
# 2018-05-30
###########

# define paths


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



load('model_outputs/sen_full_data.RData')

autoplot(data_sen)


# Force small n prev.

pr <- data_sen$pr$positive / data_sen$pr$examined
data_sen$pr$examined <- 5

pos <- (5 * pr)
round_ud <- 
  (pos - floor(pos)) > runif(length(pr))

pos[round_ud] <- ceiling(pos[round_ud])
pos[!round_ud] <- floor(pos[!round_ud])
data_sen$pr$positive <- pos



mesh_sen <- build_mesh(data_sen, mesh.args = list(max.edge = c(0.1, 0.8), cut = 0.1, offset = c(2, 2), concave = -0.03))

data_cv1_sen <- cv_random_folds(data_sen, k = 10) # todo
autoplot(data_cv1_sen, jitter = 0)
save(data_cv1_sen, file = 'model_outputs/senpsn_cv_1.RData')



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



delay <- 10

cat('Start cv1 model 1')
# Run 3 x models with 3 x hyperpars on cv1.

cat('Start cv1 model 2')
arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv1_output2 <- run_cv(data_cv1_sen, mesh_sen, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      cores = 10, parallel_delay = delay, 
                      drop_covs = 9)
obspred_map(data_cv1_sen, cv1_output2, column = FALSE)
ggsave('figs/senpsn_polygons_only_obspred_map.png')
obspred_map(data_cv1_sen, cv1_output2, trans = 'log10', column = FALSE)
ggsave('figs/senpsn_polygons_only_obspred_map_log.png')
autoplot(cv1_output2, type = 'obs_preds', CI = FALSE)
ggsave('figs/senpsn_polygons_only_obspred.png')
save(cv1_output2, file = 'model_outputs/senpsn_polygon_cv_1.RData')

cat('Start cv1 model3')
arg_list[c('use_polygons', 'use_points')] <- c(1, 1)
cv1_output3 <- run_cv(data_cv1_sen, mesh_sen, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      cores = 10, parallel_delay = delay, drop_covs = 9)
obspred_map(data_cv1_sen, cv1_output3, column = FALSE)
ggsave('figs/senpsn_both_obspred_map.png')
obspred_map(data_cv1_sen, cv1_output3, trans = 'log10', column = FALSE)
ggsave('figs/senpsn_both_obspred_map_log.png')
autoplot(cv1_output3, type = 'obs_preds', CI = TRUE)
ggsave('figs/senpsn_both_obspred.png')
save(cv1_output3, file = 'model_outputs/senpsn_joint_cv_1.RData')




cat('Start cv1 model 4')
arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv1_output4 <- run_cv(data_cv1_sen, mesh_sen, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      cores = 10, parallel_delay = delay)
obspred_map(data_cv1_sen, cv1_output4, column = FALSE)
ggsave('figs/senpsn_pr_gp_obspred_map.png')
obspred_map(data_cv1_sen, cv1_output4, trans = 'log10', column = FALSE)
ggsave('figs/senpsn_pr_gp_only_obspred_map_log.png')
autoplot(cv1_output4, type = 'obs_preds', CI = FALSE)
ggsave('figs/senpsn_pr_gp_only_obspred.png')
save(cv1_output4, file = 'model_outputs/senpsn_pr_gp_cv_1.RData')


cv1_output2$summary$polygon_metrics
cv1_output4$summary$polygon_metrics
cv1_output3$summary$polygon_metrics








# Run 3 x models with 3 x hyperpars on cv3 Spatial.

cat('Start cv3')
data_cv3_sen <- cv_spatial_folds(data_sen, k = 5, keep_pr = TRUE)
save(data_cv3_sen, file = 'model_outputs/senpsn_cv_3.RData')
autoplot(data_cv3_sen, jitter = 0.0)
ggsave('figs/senpsn_cv_spatial.png')


cat('Start cv3 model1')
# Run 3 x models with 3 x hyperpars on cv1.

cat('Start cv3 model2')
arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv3_output2 <- run_cv(data_cv3_sen, mesh_sen, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      parallel_delay = delay, cores = 5, drop_covs = 9)
obspred_map(data_cv3_sen, cv3_output2, column = FALSE)
ggsave('figs/senpsn_polygons_only_obspred_map3.png')
obspred_map(data_cv3_sen, cv3_output2, trans = 'log10', column = FALSE)
ggsave('figs/senpsn_polygons_only_obspred_map_log3.png')
autoplot(cv3_output2, type = 'obs_preds', CI = FALSE)
ggsave('figs/senpsn_polygons_only_obspred3.png')
save(cv3_output2, file = 'model_outputs/senpsn_polygon_cv_3.RData')

cat('Start cv3 model3')
arg_list[c('use_polygons', 'use_points')] <- c(1, 1)
cv3_output3 <- run_cv(data_cv3_sen, mesh_sen, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      parallel_delay = delay, cores = 5, drop_covs = 9)
obspred_map(data_cv3_sen, cv3_output3, column = FALSE)
ggsave('figs/senpsn_both_obspred_map3.png')
obspred_map(data_cv3_sen, cv3_output3, trans = 'log10', column = FALSE)
ggsave('figs/senpsn_both_obspred_map_log3.png')
autoplot(cv3_output3, type = 'obs_preds', CI = FALSE)
ggsave('figs/senpsn_both_obspred3.png')

save(cv3_output3, file = 'model_outputs/senpsn_joint_cv_3.RData')




cat('Start cv3 model4')
arg_list[c('use_polygons', 'use_points')] <- c(1, 0)
cv3_output4 <- run_cv(data_cv3_sen, mesh_sen, its = 1000, 
                      model.args = arg_list, CI = 0.8, 
                      parallel_delay = delay, cores = 5, drop_covs = NULL)
obspred_map(data_cv3_sen, cv3_output4, column = FALSE)
ggsave('figs/senpsn_pr_gp_obspred_map3.png')
obspred_map(data_cv3_sen, cv3_output4, trans = 'log10', column = FALSE)
ggsave('figs/senpsn_pr_gp_obspred_map_log3.png')
autoplot(cv3_output4, type = 'obs_preds', CI = FALSE)
ggsave('figs/senpsn_pr_gp_obspred3.png')

save(cv3_output4, file = 'model_outputs/senpsn_pr_gp_cv_3.RData')



cv3_output2$summary$polygon_metrics
cv3_output4$summary$polygon_metrics
cv3_output3$summary$polygon_metrics




senpsn1df <- rbind(cv1_output2$summary$combined_aggregated %>% mutate(model = 'polys'),
                  cv1_output3$summary$combined_aggregated %>% mutate(model = 'both'),
                  cv1_output4$summary$combined_aggregated %>% mutate(model = 'prgp'))

write.csv(senpsn1df, 'model_outputs/senpsn1df.csv')





senpsn3df <- rbind(cv3_output2$summary$combined_aggregated %>% mutate(model = 'polys'),
                  cv3_output3$summary$combined_aggregated %>% mutate(model = 'both'),
                  cv3_output4$summary$combined_aggregated %>% mutate(model = 'prgp'))

write.csv(senpsn3df, 'model_outputs/senpsn3df.csv')


get_sigs <- function(df){
  summary <- 
    df %>% 
      mutate(error = abs(pred_api - response)) %>% 
      group_by(fold, model) %>% 
      summarise(mae = mean(error))
  
  summary_prgp <- 
    summary %>% 
      filter(model != 'both') %>% 
      pivot_wider(names_from = model, values_from = mae)
  prgp <- t.test(summary_prgp$polys, summary_prgp$prgp, paired = TRUE)$p.value
  
  summary_both <- 
    summary %>% 
    filter(model != 'prgp') %>% 
    pivot_wider(names_from = model, values_from = mae)
  both <- t.test(summary_both$polys, summary_both$both, paired = TRUE)$p.value
  c(prgp = prgp, both = both)
}


get_sigs(sen1k2df)
get_sigs(sen3k2df)




print('Finished')
