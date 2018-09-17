
# figure 1.cross validation. %% Do fig 1 and 2, random and spatial cv. IDN on top, MDG and SEN below in each.
# figure 3 and 4. data and predicted incidence maps. Indonesia and Senegal only. Fig 3 ind, fig 4 sen Data, Rand, Spatial for best model? Joint model?
# figure 5, 6. Spat and random cv. PR vs Poly columns, countries as rows, model as colour?


if(Sys.info()["user"] != 'anita'){
  setwd('~/timz/timothy/point_polygon_joint_comparison')
} else {
  setwd('~/Z/users/anita/point_polygon_join_comparison_analysis')
}

# Libs

## Spatial packages
library(raster)
library(maptools)
library(rgeos)

## dataframe packages
library(dplyr)
library(readr)
library(magrittr)
library(tidyr)

library(malariaAtlas)

library(ggplot2)
library(cowplot)
theme_set(theme_minimal())

source('plotting_functions.R')

# Paths

## IDN

### Full data object

full_data_idn_path <- 'model_outputs/idn_full_data.RData'

### Cross validation object

data_cv1_idn_path <- 'model_outputs/idn_cv_1.RData'
data_cv2_idn_path <- 'model_outputs/idn_cv_2.RData'

### CV 1 output

cv1_points_idn_path <- 'model_outputs/idn_points_cv_1.RData'
cv1_polys_idn_path <- 'model_outputs/idn_polygon_cv_1.RData'
cv1_both_idn_path <- 'model_outputs/idn_joint_cv_1.RData'

### CV 2 output

cv2_points_idn_path <- 'model_outputs/idn_points_cv_2.RData'
cv2_polys_idn_path <- 'model_outputs/idn_polygon_cv_2.RData'
cv2_both_idn_path <- 'model_outputs/idn_joint_cv_2.RData'


## SEN

### Full data object

full_data_sen_path <- 'model_outputs/sen_full_data.RData'


### Cross validation object

data_cv1_sen_path <- 'model_outputs/sen_cv_1.RData'
data_cv2_sen_path <- 'model_outputs/sen_cv_2.RData'


### CV 1 output

cv1_points_sen_path <- 'model_outputs/sen_points_cv_1.RData'
cv1_polys_sen_path <- 'model_outputs/sen_polygon_cv_1.RData'
cv1_both_sen_path <- 'model_outputs/sen_joint_cv_1.RData'


### CV 2 output

cv2_points_sen_path <- 'model_outputs/sen_points_cv_2.RData'
cv2_polys_sen_path <- 'model_outputs/sen_polygon_cv_2.RData'
cv2_both_sen_path <- 'model_outputs/sen_joint_cv_2.RData' # todo



## MDG


### Full data object

full_data_mdg_path <- 'model_outputs/mdg_full_data.RData'


### Cross validation object

data_cv1_mdg_path <- 'model_outputs/mdg_cv_1.RData'
data_cv2_mdg_path <- 'model_outputs/mdg_cv_2.RData'


### CV 1 output

cv1_points_mdg_path <- 'model_outputs/mdg_points_cv_1.RData'
cv1_polys_mdg_path <- 'model_outputs/mdg_polygon_cv_1.RData'
cv1_both_mdg_path <- 'model_outputs/mdg_joint_cv_1.RData'


### CV 2 output

cv2_points_mdg_path <- 'model_outputs/mdg_points_cv_2.RData'
cv2_polys_mdg_path <- 'model_outputs/mdg_polygon_cv_2.RData'
cv2_both_mdg_path <- 'model_outputs/mdg_joint_cv_2.RData'


# figure 1.cross validation. %% Do fig 1 and 2, random and spatial cv. IDN on top, MDG and SEN below in each.

# Fig 1 - random cross validation for each country
# Load data
data_cv1_idn <- get(load(data_cv1_idn_path))
data_cv1_sen <- get(load(data_cv1_sen_path))
data_cv1_mdg <- get(load(data_cv1_mdg_path))

p1 <- autoplot(data_cv1_idn) + guides(fill = FALSE)
p2 <- autoplot(data_cv1_sen) + guides(fill = FALSE)
p3 <- autoplot(data_cv1_mdg) + guides(fill = FALSE)

bottom_row <- plot_grid(p2, p3, labels = c('B', 'C'))

full_plot <- plot_grid(p1, bottom_row, ncol = 1, labels = c('A', ''))

png('figs/random_crossvalidation_full.png', height = 1000, width = 1000)
print(full_plot)
dev.off()

rm(data_cv1_mdg)
gc()


# Fig 2 - spatial cross validation for each country
data_cv2_idn <- get(load(data_cv2_idn_path))
data_cv2_sen <- get(load(data_cv2_sen_path))
data_cv2_mdg <- get(load(data_cv2_mdg_path))


p1 <- autoplot(data_cv2_idn) + guides(fill = FALSE)
p2 <- autoplot(data_cv2_sen) + guides(fill = FALSE)
p3 <- autoplot(data_cv2_mdg) + guides(fill = FALSE)

bottom_row <- plot_grid(p2, p3, labels = c('B', 'C'))

full_plot <- plot_grid(p1, bottom_row, ncol = 1, labels = c('A', ''))

png('figs/spatial_crossvalidation_full.png', height = 1000, width = 1000)
print(full_plot)
dev.off()

rm(data_cv2_mdg)
gc()



# figure 3 data and predicted incidence maps. Indonesia only. Data, Rand, Spatial for best model? Joint model?
# todo add prevalence points

# Fig 3 - IDN: a) Data, predicted incidence from joint model for b) random cv, and c) spatial cv
#full_data_idn <- get(load(full_data_idn_path))
cv1_both_idn <- get(load(cv1_both_idn_path))
cv2_both_idn <- get(load(cv2_both_idn_path))


p1 <- obspred_map(data_cv1_idn, cv1_both_idn, trans = 'log1p')
p2 <- obspred_map(data_cv2_idn, cv2_both_idn, trans = 'log1p')


idn_preds_plot <- plot_grid(p1[[1]], p1[[2]], p2[[2]], labels = LETTERS[1:3], ncol = 1)

png('figs/idn_both_cv12_preds.png', height = 1500, width = 1200)
print(idn_preds_plot)
dev.off()


#rm(full_data_idn)
rm(cv1_both_idn)
rm(cv2_both_idn)
gc()


# figure 4 data and predicted incidence maps. Senegal only. Data, Rand, Spatial for best model? Joint model?

# Fig 4 - SEN: a) Data, predicted incidence from joint model for b) random cv, and c) spatial cv
#full_data_sen <- get(load(full_data_sen_path))
cv1_both_sen <- get(load(cv1_both_sen_path))
cv2_both_sen <- get(load(cv2_both_sen_path))


p1 <- obspred_map(data_cv1_sen, cv1_both_sen, trans = 'log1p')
p2 <- obspred_map(data_cv2_sen, cv2_both_sen, trans = 'log1p')

# Todo! switch back to p1
sen_preds_plot <- plot_grid(p2[[1]], p2[[2]], p2[[2]], labels = LETTERS[1:3], ncol = 3)

png('figs/sen_both_cv12_preds.png', height = 700, width = 1200)
print(sen_preds_plot)
dev.off()


#rm(full_data_sen)
rm(cv1_both_sen)
rm(cv2_both_sen)
gc()


# figure 5, random cv. PR vs Poly columns, countries as rows, model as colour?


cv1_points_idn <- get(load(cv1_points_idn_path))
cv1_polys_idn <- get(load(cv1_polys_idn_path))
cv1_both_idn <- get(load(cv1_both_idn_path))

idn_cv1_poly_df <- rbind(cv1_points_idn$summary$combined_aggregated %>% cbind(model = 'points'),
                         cv1_polys_idn$summary$combined_aggregated %>% cbind(model = 'polygons'),
                         cv1_both_idn$summary$combined_aggregated %>% cbind(model = 'both'))
idn_cv1_pr_df <-  rbind(cv1_points_idn$summary$combined_pr %>% cbind(model = 'points'),
                        cv1_polys_idn$summary$combined_pr %>% cbind(model = 'polygons'),
                        cv1_both_idn$summary$combined_pr %>% cbind(model = 'both'))


idn_cv1_metrics <- list(rbind(cv1_points_idn$summary$polygon_metrics %>% cbind(model = 'points'),
                              cv1_polys_idn$summary$polygon_metrics %>% cbind(model = 'polygons'),
                              cv1_both_idn$summary$polygon_metrics %>% cbind(model = 'both')),
                        rbind(cv1_points_idn$summary$pr_metrics %>% cbind(model = 'points'),
                              cv1_polys_idn$summary$pr_metrics %>% cbind(model = 'polygons'),
                              cv1_both_idn$summary$pr_metrics %>% cbind(model = 'both')))

rm(cv1_points_idn)
rm(cv1_polys_idn)
rm(cv1_both_idn)
gc()


cv1_points_sen <- get(load(cv1_points_sen_path))
cv1_polys_sen <- get(load(cv1_polys_sen_path))
cv1_both_sen <- get(load(cv1_both_sen_path))


sen_cv1_poly_df <- rbind(cv1_points_sen$summary$combined_aggregated %>% cbind(model = 'points'),
                         cv1_polys_sen$summary$combined_aggregated %>% cbind(model = 'polygons'),
                         cv1_both_sen$summary$combined_aggregated %>% cbind(model = 'both'))
sen_cv1_pr_df <-  rbind(cv1_points_sen$summary$combined_pr %>% cbind(model = 'points'),
                        cv1_polys_sen$summary$combined_pr %>% cbind(model = 'polygons'),
                        cv1_both_sen$summary$combined_pr %>% cbind(model = 'both'))

sen_cv1_metrics <- list(rbind(cv1_points_sen$summary$polygon_metrics %>% cbind(model = 'points'),
                              cv1_polys_sen$summary$polygon_metrics %>% cbind(model = 'polygons'),
                              cv1_both_sen$summary$polygon_metrics %>% cbind(model = 'both')),
                        rbind(cv1_points_sen$summary$pr_metrics %>% cbind(model = 'points'),
                              cv1_polys_sen$summary$pr_metrics %>% cbind(model = 'polygons'),
                              cv1_both_sen$summary$pr_metrics %>% cbind(model = 'both')))



rm(cv1_points_sen)
rm(cv1_polys_sen)
rm(cv1_both_sen)
gc()


cv1_points_mdg <- get(load(cv1_points_mdg_path))
cv1_polys_mdg <- get(load(cv1_polys_mdg_path))
cv1_both_mdg <- get(load(cv1_both_mdg_path))



mdg_cv1_poly_df <- rbind(cv1_points_mdg$summary$combined_aggregated %>% cbind(model = 'points'),
                         cv1_polys_mdg$summary$combined_aggregated %>% cbind(model = 'polygons'),
                         cv1_both_mdg$summary$combined_aggregated %>% cbind(model = 'both'))
mdg_cv1_pr_df <-  rbind(cv1_points_mdg$summary$combined_pr %>% cbind(model = 'points'),
                        cv1_polys_mdg$summary$combined_pr %>% cbind(model = 'polygons'),
                        cv1_both_mdg$summary$combined_pr %>% cbind(model = 'both'))



mdg_cv1_metrics <- list(rbind(cv1_points_mdg$summary$polygon_metrics %>% cbind(model = 'points'),
                              cv1_polys_mdg$summary$polygon_metrics %>% cbind(model = 'polygons'),
                              cv1_both_mdg$summary$polygon_metrics %>% cbind(model = 'both')),
                        rbind(cv1_points_mdg$summary$pr_metrics %>% cbind(model = 'points'),
                              cv1_polys_mdg$summary$pr_metrics %>% cbind(model = 'polygons'),
                              cv1_both_mdg$summary$pr_metrics %>% cbind(model = 'both')))

rm(cv1_points_mdg)
rm(cv1_polys_mdg)
rm(cv1_both_mdg)
gc()






idn_poly <- ggplot(idn_cv1_poly_df, aes(response, pred_api, colour = model)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  scale_y_log10() + 
  scale_x_log10() +
  guides(colour = FALSE)
idn_points <- ggplot(idn_cv1_pr_df, aes(prevalence, pred_prev, colour = model)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)  +
  guides(colour = FALSE)

sen_poly <- ggplot(sen_cv1_poly_df, aes(response, pred_api, colour = model)) + 
  geom_jitter(width = 0.5, height = 0.5) + #todo
  geom_abline(slope = 1, intercept = 0) +
  scale_y_log10() + 
  scale_x_log10() +
  guides(colour = FALSE)
sen_points <- ggplot(sen_cv1_pr_df, aes(prevalence, pred_prev, colour = model)) + 
  geom_jitter(width = 0.02, height = 0.02) + #todo
  geom_abline(slope = 1, intercept = 0)  +
  guides(colour = FALSE)


mdg_poly <- ggplot(mdg_cv1_poly_df, aes(response, pred_api, colour = model)) + 
  geom_jitter(width = 0.5, height = 0.5) + #todo
  geom_abline(slope = 1, intercept = 0) +
  scale_y_log10() + 
  scale_x_log10() +
  guides(colour = FALSE)
mdg_points <- ggplot(mdg_cv1_pr_df, aes(prevalence, pred_prev, colour = model)) + 
  geom_jitter(width = 0.02, height = 0.02) + #todo
  geom_abline(slope = 1, intercept = 0)  +
  guides(colour = FALSE)


scatter_list <- list(idn_poly, idn_points, sen_poly, sen_points, mdg_poly, mdg_points)
full_obs_pred_scatter <- plot_grid(plotlist = scatter_list,
                                   labels = LETTERS[1:6], ncol = 2)

png('figs/random_cv_scatter.png', height = 1500, width = 1100)
print(full_obs_pred_scatter)
dev.off()



pdf('figs/random_cv_scatter.pdf', height = 15, width = 11)
print(full_obs_pred_scatter)
dev.off()


# Fig 4 and 6


cv2_points_idn <- get(load(cv2_points_idn_path))
cv2_polys_idn <- get(load(cv2_polys_idn_path))
cv2_both_idn <- get(load(cv2_both_idn_path))

idn_cv2_poly_df <- rbind(cv2_points_idn$summary$combined_aggregated %>% cbind(model = 'points'),
                         cv2_polys_idn$summary$combined_aggregated %>% cbind(model = 'polygons'),
                         cv2_both_idn$summary$combined_aggregated %>% cbind(model = 'both'))
idn_cv2_pr_df <-  rbind(cv2_points_idn$summary$combined_pr %>% cbind(model = 'points'),
                        cv2_polys_idn$summary$combined_pr %>% cbind(model = 'polygons'),
                        cv2_both_idn$summary$combined_pr %>% cbind(model = 'both'))


idn_cv2_metrics <- list(rbind(cv2_points_idn$summary$polygon_metrics %>% cbind(model = 'points'),
                              cv2_polys_idn$summary$polygon_metrics %>% cbind(model = 'polygons'),
                              cv2_both_idn$summary$polygon_metrics %>% cbind(model = 'both')),
                        rbind(cv2_points_idn$summary$pr_metrics %>% cbind(model = 'points'),
                              cv2_polys_idn$summary$pr_metrics %>% cbind(model = 'polygons'),
                              cv2_both_idn$summary$pr_metrics %>% cbind(model = 'both')))

rm(cv2_points_idn)
rm(cv2_polys_idn)
rm(cv2_both_idn)
gc()


cv2_points_sen <- get(load(cv2_points_sen_path))
cv2_polys_sen <- get(load(cv2_polys_sen_path))
cv2_both_sen <- get(load(cv2_both_sen_path))


sen_cv2_poly_df <- rbind(cv2_points_sen$summary$combined_aggregated %>% cbind(model = 'points'),
                         cv2_polys_sen$summary$combined_aggregated %>% cbind(model = 'polygons'),
                         cv2_both_sen$summary$combined_aggregated %>% cbind(model = 'both'))
sen_cv2_pr_df <-  rbind(cv2_points_sen$summary$combined_pr %>% cbind(model = 'points'),
                        cv2_polys_sen$summary$combined_pr %>% cbind(model = 'polygons'),
                        cv2_both_sen$summary$combined_pr %>% cbind(model = 'both'))

sen_cv2_metrics <- list(rbind(cv2_points_sen$summary$polygon_metrics %>% cbind(model = 'points'),
                              cv2_polys_sen$summary$polygon_metrics %>% cbind(model = 'polygons'),
                              cv2_both_sen$summary$polygon_metrics %>% cbind(model = 'both')),
                        rbind(cv2_points_sen$summary$pr_metrics %>% cbind(model = 'points'),
                              cv2_polys_sen$summary$pr_metrics %>% cbind(model = 'polygons'),
                              cv2_both_sen$summary$pr_metrics %>% cbind(model = 'both')))



rm(cv2_points_sen)
rm(cv2_polys_sen)
rm(cv2_both_sen)
gc()


cv2_points_mdg <- get(load(cv2_points_mdg_path))
cv2_polys_mdg <- get(load(cv2_polys_mdg_path))
cv2_both_mdg <- get(load(cv2_both_mdg_path))



mdg_cv2_poly_df <- rbind(cv2_points_mdg$summary$combined_aggregated %>% cbind(model = 'points'),
                         cv2_polys_mdg$summary$combined_aggregated %>% cbind(model = 'polygons'),
                         cv2_both_mdg$summary$combined_aggregated %>% cbind(model = 'both'))
mdg_cv2_pr_df <-  rbind(cv2_points_mdg$summary$combined_pr %>% cbind(model = 'points'),
                        cv2_polys_mdg$summary$combined_pr %>% cbind(model = 'polygons'),
                        cv2_both_mdg$summary$combined_pr %>% cbind(model = 'both'))



mdg_cv2_metrics <- list(rbind(cv2_points_mdg$summary$polygon_metrics %>% cbind(model = 'points'),
                              cv2_polys_mdg$summary$polygon_metrics %>% cbind(model = 'polygons'),
                              cv2_both_mdg$summary$polygon_metrics %>% cbind(model = 'both')),
                        rbind(cv2_points_mdg$summary$pr_metrics %>% cbind(model = 'points'),
                              cv2_polys_mdg$summary$pr_metrics %>% cbind(model = 'polygons'),
                              cv2_both_mdg$summary$pr_metrics %>% cbind(model = 'both')))

rm(cv2_points_mdg)
rm(cv2_polys_mdg)
rm(cv2_both_mdg)
gc()






idn_poly <- ggplot(idn_cv2_poly_df, aes(response, pred_api, colour = model)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  scale_y_log10() + 
  scale_x_log10() +
  guides(colour = FALSE)
idn_points <- ggplot(idn_cv2_pr_df, aes(prevalence, pred_prev, colour = model)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)  +
  guides(colour = FALSE)

sen_poly <- ggplot(sen_cv2_poly_df, aes(response, pred_api, colour = model)) + 
  geom_jitter(width = 0.5, height = 0.5) + #todo
  geom_abline(slope = 1, intercept = 0) +
  scale_y_log10() + 
  scale_x_log10() +
  guides(colour = FALSE)
sen_points <- ggplot(sen_cv2_pr_df, aes(prevalence, pred_prev, colour = model)) + 
  geom_jitter(width = 0.02, height = 0.02) + #todo
  geom_abline(slope = 1, intercept = 0)  +
  guides(colour = FALSE)


mdg_poly <- ggplot(mdg_cv2_poly_df, aes(response, pred_api, colour = model)) + 
  geom_jitter(width = 0.5, height = 0.5) + #todo
  geom_abline(slope = 1, intercept = 0) +
  scale_y_log10() + 
  scale_x_log10() +
  guides(colour = FALSE)
mdg_points <- ggplot(mdg_cv2_pr_df, aes(prevalence, pred_prev, colour = model)) + 
  geom_jitter(width = 0.02, height = 0.02) + #todo
  geom_abline(slope = 1, intercept = 0)  +
  guides(colour = FALSE)


scatter_list <- list(idn_poly, idn_points, sen_poly, sen_points, mdg_poly, mdg_points)
full_obs_pred_scatter <- plot_grid(plotlist = scatter_list,
                                   labels = LETTERS[1:6], ncol = 2)

png('figs/spatial_cv_scatter.png', height = 1500, width = 700)
print(full_obs_pred_scatter)
dev.off()


pdf('figs/spatial_cv_scatter.pdf', height = 15, width = 11)
print(full_obs_pred_scatter)
dev.off()



# Useful summary tables
table1_skeleton <- 
  "Incidence & Pearson & Indonesia & %s & %s &  %s\\\\
&& Senegal & %s & %s &  %s\\\\
&& Madagascar & %s & %s &  %s\\vspace{1mm}\\\\
& Spearman & Indonesia & %s & %s &  %s\\\\
&& Senegal & %s & %s &  %s\\\\
&& Madagascar & %s & %s &  %s\\vspace{3mm} \\\\
Prevalence & Pearson & Indonesia & %s & %s &  %s\\\\
&& Senegal & %s & %s &  %s\\\\
&& Madagascar & %s & %s &  %s\\vspace{1mm}\\\\
& Spearman & Indonesia & %s & %s &  %s\\\\\
&& Senegal & %s & %s &  %s\\\\
&& Madagascar & %s & %s &  %s\\\\"

r <- c(idn_cv1_metrics[[1]]$pearson, sen_cv1_metrics[[1]]$pearson, mdg_cv1_metrics[[1]]$pearson,
       idn_cv1_metrics[[1]]$spearman, sen_cv1_metrics[[1]]$spearman, mdg_cv1_metrics[[1]]$spearman,
       idn_cv1_metrics[[2]]$pearson, sen_cv1_metrics[[2]]$pearson, mdg_cv1_metrics[[2]]$pearson,
       idn_cv1_metrics[[2]]$spearman, sen_cv1_metrics[[2]]$spearman, mdg_cv1_metrics[[2]]$spearman)

r <- format(round(r, 2), nsmall = 2)

table1 <- do.call(sprintf, c(table1_skeleton, as.list(r)))

write(table1, 'figs/table1.txt')




r_spat <- c(idn_cv1_metrics[[1]]$pearson, sen_cv1_metrics[[1]]$pearson, mdg_cv1_metrics[[1]]$pearson,
            idn_cv1_metrics[[1]]$spearman, sen_cv1_metrics[[1]]$spearman, mdg_cv1_metrics[[1]]$spearman,
            idn_cv1_metrics[[2]]$pearson, sen_cv1_metrics[[2]]$pearson, mdg_cv1_metrics[[2]]$pearson,
            idn_cv1_metrics[[2]]$spearman, sen_cv1_metrics[[2]]$spearman, mdg_cv1_metrics[[2]]$spearman)

r_spat <- format(round(r_spat, 2), nsmall = 2)

table2 <- do.call(sprintf, c(table1_skeleton, as.list(r_spat)))

write(table2, 'figs/table2.txt')



coverage_skeleton <-
  "Incidence & Random & Indonesia & %s & %s &  %s\\\\
&& Senegal & %s & %s &  %s\\\\
&& Madagascar & %s & %s &  %s\\vspace{1mm}\\\\
& Spatial & Indonesia & %s & %s &  %s\\\\
&& Senegal & %s & %s &  %s\\\\
&& Madagascar & %s & %s &  %s\\vspace{3mm} \\\\
Prevalence & Random & Indonesia & %s & %s &  %s\\\\
&& Senegal & %s & %s &  %s\\\\
&& Madagascar & %s & %s &  %s\\vspace{1mm}\\\\
& Spatial & Indonesia & %s & %s &  %s\\\\\
&& Senegal & %s & %s &  %s\\\\
&& Madagascar & %s & %s &  %s\\\\"

cov <- c(idn_cv1_metrics[[1]]$coverage, sen_cv1_metrics[[1]]$coverage, mdg_cv1_metrics[[1]]$coverage,
         idn_cv1_metrics[[1]]$coverage, sen_cv1_metrics[[1]]$coverage, mdg_cv1_metrics[[1]]$coverage,
         idn_cv1_metrics[[2]]$coverage, sen_cv1_metrics[[2]]$coverage, mdg_cv1_metrics[[2]]$coverage,
         idn_cv1_metrics[[2]]$spearman, sen_cv1_metrics[[2]]$coverage, mdg_cv1_metrics[[2]]$coverage)

cov <- format(round(cov, 2), nsmall = 2)

table3 <- do.call(sprintf, c(coverage_skeleton, as.list(cov)))

write(table3, 'figs/table3.txt')


# Further SI figures.