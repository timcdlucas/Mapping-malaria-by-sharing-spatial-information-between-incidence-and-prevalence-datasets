
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
theme_update(text = element_text(size = 10))


# Metrics 



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

p1 <- autoplot(data_cv1_idn, jitter = 0, size = 0.7) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')
p2 <- autoplot(data_cv1_sen, jitter = 0, size = 0.7) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')
p3 <- autoplot(data_cv1_mdg, jitter = 0, size = 0.7) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')

bottom_row <- plot_grid(p2, p3, labels = c('B', 'C'), rel_widths = c(0.6, 0.4))

full_plot <- plot_grid(p1, bottom_row, ncol = 1, labels = c('A', ''), rel_heights = c(0.45, 0.55))

png('figs/summaries/random_crossvalidation_full.png', height = 100, width = 100, unit = 'mm', res = 720)
print(full_plot)
dev.off()

gc()


# Fig 2 - spatial cross validation for each country
data_cv2_idn <- get(load(data_cv2_idn_path))
data_cv2_sen <- get(load(data_cv2_sen_path))
data_cv2_mdg <- get(load(data_cv2_mdg_path))


p1 <- autoplot(data_cv2_idn, jitter = 0, size = 0.7) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')
p2 <- autoplot(data_cv2_sen, jitter = 0, size = 0.7) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')
p3 <- autoplot(data_cv2_mdg, jitter = 0, size = 0.7) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')

bottom_row <- plot_grid(p2, p3, labels = c('B', 'C'), rel_widths = c(0.6, 0.4))

full_plot <- plot_grid(p1, bottom_row, ncol = 1, labels = c('A', ''), rel_heights = c(0.45, 0.55))

png('figs/summaries/spatial_crossvalidation_full.png', height = 100, width = 100, unit = 'mm', res = 720)
print(full_plot)
dev.off()

gc()



# figure 3 data and predicted incidence maps. Indonesia only. Data, Rand, Spatial for best model? Joint model?
# todo add prevalence points

# Fig 3 - IDN: a) Data, predicted incidence from joint model for b) random cv, and c) spatial cv
#full_data_idn <- get(load(full_data_idn_path))
cv1_both_idn <- get(load(cv1_both_idn_path))
cv2_both_idn <- get(load(cv2_both_idn_path))

xx= data_cv2_idn[1:3]
class(xx) = 'ppj_cv' #todo
p1 <- obspred_map(data_cv1_idn, cv1_both_idn, trans = 'log1p', 
                  legend_title = 'API', 
                  breaks = c(1, 10, 100, 300, 500))
p2 <- obspred_map(xx, cv2_both_idn, trans = 'log1p', legend_title = 'API')


panel1 <- p1[[1]] +
  guides(fill = FALSE) + 
  labs(x = '', y = '')
panel2 <- p1[[2]] + 
  guides(fill = FALSE) + 
  labs(x = '', y = 'Latitude')
panel3 <- p2[[2]] + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = '')

legend <- get_legend(p1[[1]])

idn_preds_plot <- plot_grid(panel1, panel2, panel3, labels = LETTERS[1:3], ncol = 1)
full_plot <- plot_grid(idn_preds_plot, legend, ncol = 2, rel_widths = c(4, 1))

png('figs/summaries/idn_both_cv12_preds.png', height = 130, width = 100, unit = 'mm', res = 720)
print(full_plot)
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


p1 <- obspred_map(data_cv1_sen, cv1_both_sen, trans = 'log1p', 
                  legend_title = 'API', 
                  mask = TRUE,
                  breaks = c(1, 10, 100, 300))
p2 <- obspred_map(data_cv2_sen, cv2_both_sen, 
                  trans = 'log1p', 
                  legend_title = 'API',
                  mask = TRUE)


panel1 <- p1[[1]] +
  guides(fill = FALSE) + 
  labs(x = '', y = '') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
panel2 <- p1[[2]] + 
  guides(fill = FALSE) + 
  labs(x = '', y = 'Latitude') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
panel3 <- p2[[2]] + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = '')+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

legend <- get_legend(p1[[1]])

sen_preds_plot <- plot_grid(panel1, panel2, panel3, labels = LETTERS[1:3], ncol = 1)
full_plot <- plot_grid(sen_preds_plot, legend, ncol = 2, rel_widths = c(4, 1))




png('figs/summaries/sen_both_cv12_preds.png', height = 150, width = 100, unit = 'mm', res = 720)
print(full_plot)
dev.off()


#rm(full_data_sen)
rm(cv1_both_sen)
rm(cv2_both_sen)
gc()

# 
# # Fig 5 - MDG: a) Data, predicted incidence from joint model for b) random cv, and c) spatial cv
# #full_data_mdg <- get(load(full_data_mdg_path))
# cv1_both_mdg <- get(load(cv1_both_mdg_path))
# cv2_both_mdg <- get(load(cv2_both_mdg_path))
# 
# 
# p1 <- obspred_map(data_cv1_mdg, cv1_both_mdg, trans = 'log1p')
# p2 <- obspred_map(data_cv2_mdg, cv2_both_mdg, trans = 'log1p')
# 
# 
# mdg_preds_plot <- plot_grid(p1[[1]], p1[[2]], p2[[2]], labels = LETTERS[1:3], ncol = 1)
# 
# png('figs/summaries/mdg_both_cv12_preds.png', height = 1500, width = 1200)
# print(mdg_preds_plot)
# dev.off()
# 
# 
# #rm(full_data_mdg)
# rm(cv1_both_mdg)
# rm(cv2_both_mdg)
# gc()


# figure 6, random cv. PR vs Poly columns, countries as rows, model as colour?

#################################################################################
## Make model comparisons. Random CV.                                          ##
#################################################################################


cv1_points_idn <- get(load(cv1_points_idn_path))
cv1_polys_idn <- get(load(cv1_polys_idn_path))
cv1_both_idn <- get(load(cv1_both_idn_path))

idn_cv1_poly_df <- rbind(cv1_points_idn$summary$combined_aggregated %>% cbind(model = 'points'),
                         cv1_polys_idn$summary$combined_aggregated %>% cbind(model = 'polygons'),
                         cv1_both_idn$summary$combined_aggregated %>% cbind(model = 'both'))
idn_cv1_pr_df <-  rbind(cv1_points_idn$summary$combined_pr %>% cbind(model = 'points'),
                        cv1_polys_idn$summary$combined_pr %>% cbind(model = 'polygons'),
                        cv1_both_idn$summary$combined_pr %>% cbind(model = 'both'))


idn_cv1_new_metrics <- 
  rbind(cv1_points_idn$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'points'),
        cv1_polys_idn$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'polygons'),
        cv1_both_idn$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'both')
  )




idn_cv1_metrics <- list(rbind(cv1_points_idn$summary$polygon_metrics %>% cbind(model = 'points'),
                              cv1_polys_idn$summary$polygon_metrics %>% cbind(model = 'polygons'),
                              cv1_both_idn$summary$polygon_metrics %>% cbind(model = 'both')),
                        rbind(cv1_points_idn$summary$pr_metrics %>% cbind(model = 'points'),
                              cv1_polys_idn$summary$pr_metrics %>% cbind(model = 'polygons'),
                              cv1_both_idn$summary$pr_metrics %>% cbind(model = 'both')))
idn_cv1_metrics[[2]] <- left_join(idn_cv1_metrics[[2]], idn_cv1_new_metrics)

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

sen_cv1_new_metrics <- 
  rbind(cv1_points_sen$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'points'),
        cv1_polys_sen$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'polygons'),
        cv1_both_sen$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'both')
  )




sen_cv1_metrics <- list(rbind(cv1_points_sen$summary$polygon_metrics %>% cbind(model = 'points'),
                              cv1_polys_sen$summary$polygon_metrics %>% cbind(model = 'polygons'),
                              cv1_both_sen$summary$polygon_metrics %>% cbind(model = 'both')),
                        rbind(cv1_points_sen$summary$pr_metrics %>% cbind(model = 'points'),
                              cv1_polys_sen$summary$pr_metrics %>% cbind(model = 'polygons'),
                              cv1_both_sen$summary$pr_metrics %>% cbind(model = 'both')))
sen_cv1_metrics[[2]] <- left_join(sen_cv1_metrics[[2]], sen_cv1_new_metrics)



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


mdg_cv1_new_metrics <- 
  rbind(cv1_points_mdg$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'points'),
        cv1_polys_mdg$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'polygons'),
        cv1_both_mdg$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'both')
        )




mdg_cv1_metrics <- list(rbind(cv1_points_mdg$summary$polygon_metrics %>% cbind(model = 'points'),
                              cv1_polys_mdg$summary$polygon_metrics %>% cbind(model = 'polygons'),
                              cv1_both_mdg$summary$polygon_metrics %>% cbind(model = 'both')),
                        rbind(cv1_points_mdg$summary$pr_metrics %>% cbind(model = 'points'),
                              cv1_polys_mdg$summary$pr_metrics %>% cbind(model = 'polygons'),
                              cv1_both_mdg$summary$pr_metrics %>% cbind(model = 'both')))
mdg_cv1_metrics[[2]] <- left_join(mdg_cv1_metrics[[2]], mdg_cv1_new_metrics)


rm(cv1_points_mdg)
rm(cv1_polys_mdg)
rm(cv1_both_mdg)
gc()




a <- 0.2
s <- 16

idn_poly <- ggplot(idn_cv1_poly_df, aes(response, pred_api, colour = model)) + 
  geom_point(alpha = a, size = 3, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_y_sqrt() + 
  scale_x_sqrt() +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)


idn_points <- ggplot(idn_cv1_pr_df, aes(prevalence, pred_prev, colour = model, size = examined)) + 
  geom_point(alpha = a, shape = s) + 
  geom_abline(slope = 1, intercept = 0)  +
  scale_y_sqrt() + 
  scale_x_sqrt() +
  labs(x = 'Observed Prevalence', y = 'Predicted Prevalence', size = 'Sample size') +
  guides(colour = FALSE)


sen_poly <- ggplot(sen_cv1_poly_df, aes(response, pred_api, colour = model)) + 
  geom_point(alpha = a + 0.2, size = 3, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_y_sqrt() + 
  scale_x_sqrt() +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)

sen_points <- sen_cv1_pr_df %>% 
  mutate(noise = abs(rnorm(nrow(.), 0, 2e-5))) %>% 
  mutate(prevalence = ifelse(prevalence == min(prevalence), prevalence + noise, prevalence)) %>% 
  ggplot(aes(prevalence, pred_prev, colour = model, size = examined)) + 
    geom_point(alpha = a, shape = s) + 
    geom_abline(slope = 1, intercept = 0)  +
    scale_y_sqrt() + 
    scale_x_sqrt() +
    labs(x = 'Observed Prevalence', y = 'Predicted Prevalence', size = 'Sample size') +
    guides(colour = FALSE)


mdg_poly <- ggplot(mdg_cv1_poly_df, aes(response, pred_api, colour = model)) + 
  geom_point(alpha = a + 0.2, size = 3, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_y_sqrt() + 
  scale_x_sqrt() +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)
mdg_points <- mdg_cv1_pr_df %>% 
  mutate(noise = abs(rnorm(nrow(.), 0, 1e-4))) %>% 
  mutate(prevalence = ifelse(prevalence == min(prevalence), prevalence + noise, prevalence)) %>% 
  ggplot(aes(prevalence, pred_prev, colour = model, size = examined)) + 
    geom_point(alpha = a, shape = s) + 
    geom_abline(slope = 1, intercept = 0)  +
    scale_y_sqrt() + 
    scale_x_sqrt() +
    labs(x = 'Observed Prevalence', y = 'Predicted Prevalence', size = 'Sample size') +
    guides(colour = FALSE)


scatter_list <- list(idn_poly, idn_points, sen_poly, sen_points, mdg_poly, mdg_points)
full_obs_pred_scatter <- plot_grid(plotlist = scatter_list,
                                   labels = LETTERS[1:6], ncol = 2)

png('figs/summaries/random_cv_scatter.png', height = 1500, width = 1100)
print(full_obs_pred_scatter)
dev.off()



pdf('figs/summaries/random_cv_scatter.pdf', height = 15, width = 11)
print(full_obs_pred_scatter)
dev.off()


# Fig 4 and 6
#################################################################################
## Make model comparisons. Spatial CV.                                         ##
#################################################################################

cv2_points_idn <- get(load(cv2_points_idn_path))
cv2_polys_idn <- get(load(cv2_polys_idn_path))
cv2_both_idn <- get(load(cv2_both_idn_path))

idn_cv2_poly_df <- rbind(cv2_points_idn$summary$combined_aggregated %>% cbind(model = 'points'),
                         cv2_polys_idn$summary$combined_aggregated %>% cbind(model = 'polygons'),
                         cv2_both_idn$summary$combined_aggregated %>% cbind(model = 'both'))
idn_cv2_pr_df <-  rbind(cv2_points_idn$summary$combined_pr %>% cbind(model = 'points'),
                        cv2_polys_idn$summary$combined_pr %>% cbind(model = 'polygons'),
                        cv2_both_idn$summary$combined_pr %>% cbind(model = 'both'))


idn_cv2_new_metrics <- 
  rbind(cv2_points_idn$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'points'),
        cv2_polys_idn$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'polygons'),
        cv2_both_idn$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'both')
  )




idn_cv2_metrics <- list(rbind(cv2_points_idn$summary$polygon_metrics %>% cbind(model = 'points'),
                              cv2_polys_idn$summary$polygon_metrics %>% cbind(model = 'polygons'),
                              cv2_both_idn$summary$polygon_metrics %>% cbind(model = 'both')),
                        rbind(cv2_points_idn$summary$pr_metrics %>% cbind(model = 'points'),
                              cv2_polys_idn$summary$pr_metrics %>% cbind(model = 'polygons'),
                              cv2_both_idn$summary$pr_metrics %>% cbind(model = 'both')))
idn_cv2_metrics[[2]] <- left_join(idn_cv2_metrics[[2]], idn_cv2_new_metrics)

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


sen_cv2_new_metrics <- 
  rbind(cv2_points_sen$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'points'),
        cv2_polys_sen$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'polygons'),
        cv2_both_sen$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'both')
  )




sen_cv2_metrics <- list(rbind(cv2_points_sen$summary$polygon_metrics %>% cbind(model = 'points'),
                              cv2_polys_sen$summary$polygon_metrics %>% cbind(model = 'polygons'),
                              cv2_both_sen$summary$polygon_metrics %>% cbind(model = 'both')),
                        rbind(cv2_points_sen$summary$pr_metrics %>% cbind(model = 'points'),
                              cv2_polys_sen$summary$pr_metrics %>% cbind(model = 'polygons'),
                              cv2_both_sen$summary$pr_metrics %>% cbind(model = 'both')))
sen_cv2_metrics[[2]] <- left_join(sen_cv2_metrics[[2]], sen_cv2_new_metrics)



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


mdg_cv2_new_metrics <- 
  rbind(cv2_points_mdg$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'points'),
        cv2_polys_mdg$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'polygons'),
        cv2_both_mdg$summary$combined_pr %>%
          mutate(abs_error = abs(prevalence - pred_prev)) %>% 
          summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
          cbind(model = 'both')
  )




mdg_cv2_metrics <- list(rbind(cv2_points_mdg$summary$polygon_metrics %>% cbind(model = 'points'),
                              cv2_polys_mdg$summary$polygon_metrics %>% cbind(model = 'polygons'),
                              cv2_both_mdg$summary$polygon_metrics %>% cbind(model = 'both')),
                        rbind(cv2_points_mdg$summary$pr_metrics %>% cbind(model = 'points'),
                              cv2_polys_mdg$summary$pr_metrics %>% cbind(model = 'polygons'),
                              cv2_both_mdg$summary$pr_metrics %>% cbind(model = 'both')))
mdg_cv2_metrics[[2]] <- left_join(mdg_cv2_metrics[[2]], mdg_cv2_new_metrics)

rm(cv2_points_mdg)
rm(cv2_polys_mdg)
rm(cv2_both_mdg)
gc()







a <- 0.2
s <- 16

idn_poly <- ggplot(idn_cv2_poly_df, aes(response, pred_api, colour = model)) + 
  geom_point(alpha = a, size = 3, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_y_sqrt() + 
  scale_x_sqrt() +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)


idn_points <- ggplot(idn_cv2_pr_df, aes(prevalence, pred_prev, colour = model, size = examined)) + 
  geom_point(alpha = a, shape = s) + 
  geom_abline(slope = 1, intercept = 0)  +
  scale_y_sqrt() + 
  scale_x_sqrt() +
  labs(x = 'Observed Prevalence', y = 'Predicted Prevalence', size = 'Sample size') +
  guides(colour = FALSE)


sen_poly <- ggplot(sen_cv2_poly_df, aes(response, pred_api, colour = model)) + 
  geom_point(alpha = a + 0.2, size = 3, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_y_sqrt() + 
  scale_x_sqrt() +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)

sen_points <- sen_cv2_pr_df %>% 
  mutate(noise = abs(rnorm(nrow(.), 0, 2e-5))) %>% 
  mutate(prevalence = ifelse(prevalence == min(prevalence), prevalence + noise, prevalence)) %>% 
  ggplot(aes(prevalence, pred_prev, colour = model, size = examined)) + 
  geom_point(alpha = a, shape = s) + 
  geom_abline(slope = 1, intercept = 0)  +
  scale_y_sqrt() + 
  scale_x_sqrt() +
  labs(x = 'Observed Prevalence', y = 'Predicted Prevalence', size = 'Sample size') +
  guides(colour = FALSE)


mdg_poly <- ggplot(mdg_cv2_poly_df, aes(response, pred_api, colour = model)) + 
  geom_point(alpha = a + 0.2, size = 3, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_y_sqrt() + 
  scale_x_sqrt() +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)
mdg_points <- mdg_cv2_pr_df %>% 
  mutate(noise = abs(rnorm(nrow(.), 0, 1e-4))) %>% 
  mutate(prevalence = ifelse(prevalence == min(prevalence), prevalence + noise, prevalence)) %>% 
  ggplot(aes(prevalence, pred_prev, colour = model, size = examined)) + 
  geom_point(alpha = a, shape = s) + 
  geom_abline(slope = 1, intercept = 0)  +
  scale_y_sqrt() + 
  scale_x_sqrt() +
  labs(x = 'Observed Prevalence', y = 'Predicted Prevalence', size = 'Sample size') +
  guides(colour = FALSE)

scatter_list <- list(idn_poly, idn_points, sen_poly, sen_points, mdg_poly, mdg_points)
full_obs_pred_scatter <- plot_grid(plotlist = scatter_list,
                                   labels = LETTERS[1:6], ncol = 2)

png('figs/summaries/spatial_cv_scatter.png', height = 1500, width = 700)
print(full_obs_pred_scatter)
dev.off()


pdf('figs/summaries/spatial_cv_scatter.pdf', height = 15, width = 11)
print(full_obs_pred_scatter)
dev.off()



# Useful summary tables
table1_skeleton <- 
  "Incidence & Indonesia & %s & %s &  %s\\\\
& Senegal & %s & %s &  %s\\\\
& Madagascar & %s & %s &  %s\\vspace{3mm}\\\\
Prevalence & Indonesia & %s & %s &  %s\\\\
& Senegal & %s & %s &  %s\\\\
& Madagascar & %s & %s &  %s\\\\"

r <- c(idn_cv1_metrics[[1]]$MAE, sen_cv1_metrics[[1]]$MAE, mdg_cv1_metrics[[1]]$MAE,
       idn_cv1_metrics[[2]]$weightMAE, sen_cv1_metrics[[2]]$weightMAE, mdg_cv1_metrics[[2]]$weightMAE)

r <- c(format(round(r[1:9], 2), nsmall = 2), format(round(r[10:18], 3), nsmall = 3))

table1 <- do.call(sprintf, c(table1_skeleton, as.list(r)))

write(table1, 'figs/summaries/table1.txt')




r_spat <- c(idn_cv2_metrics[[1]]$MAE, sen_cv2_metrics[[1]]$MAE, mdg_cv2_metrics[[1]]$MAE,
            idn_cv2_metrics[[2]]$weightMAE, sen_cv2_metrics[[2]]$weightMAE, mdg_cv2_metrics[[2]]$weightMAE)

r_spat <- c(format(round(r_spat[1:9], 2), nsmall = 2), format(round(r_spat[10:18], 3), nsmall = 3))

table2 <- do.call(sprintf, c(table1_skeleton, as.list(r_spat)))

write(table2, 'figs/summaries/table2.txt')



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
         idn_cv2_metrics[[2]]$coverage, sen_cv2_metrics[[2]]$coverage, mdg_cv2_metrics[[2]]$coverage,
         idn_cv2_metrics[[2]]$coverage, sen_cv2_metrics[[2]]$coverage, mdg_cv2_metrics[[2]]$coverage)

cov <- format(round(cov, 2), nsmall = 2)

table3 <- do.call(sprintf, c(coverage_skeleton, as.list(cov)))

write(table3, 'figs/summaries/table3.txt')


# Further SI figures.