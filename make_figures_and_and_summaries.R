
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
library(Matrix)

library(malariaAtlas)

library(xtable)
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
#data_cv2_idn_path <- 'model_outputs/idn_cv_2.RData'
data_cv3_idn_path <- 'model_outputs/idn_cv_3.RData'

### CV 1 output

#cv1_points_idn_path <- 'model_outputs/idn_points_cv_1.RData'
cv1_polys_idn_path <- 'model_outputs/idn_polygon_cv_1.RData'
cv1_both_idn_path <- 'model_outputs/idn_joint_cv_1.RData'
cv1_prgp_idn_path <- 'model_outputs/idn_prgp_cv_1.RData'


### CV 3 output

#cv3_points_idn_path <- 'model_outputs/idn_points_cv_2.RData'
cv3_polys_idn_path <- 'model_outputs/idn_polygon_cv_3.RData'
cv3_both_idn_path <- 'model_outputs/idn_joint_cv_3.RData'
cv3_prgp_idn_path <- 'model_outputs/idn_prgp_cv_3.RData'


### Full models
idn_baseline <- 'model_outputs/full_baseline_idn.RData'
idn_prgp <- 'model_outputs/full_prevgp_idn.RData'
idn_joint <- 'model_outputs/full_model_idn.RData'



## SEN

### Full data object

full_data_sen_path <- 'model_outputs/sen_full_data.RData'


### Cross validation object

data_cv1_sen_path <- 'model_outputs/sen_cv_1.RData'
#data_cv2_sen_path <- 'model_outputs/sen_cv_2.RData'
data_cv3_sen_path <- 'model_outputs/sen_cv_3.RData'


### CV 1 output

cv1_polys_sen_path <- 'model_outputs/sen_polygon_cv_1.RData'
cv1_both_sen_path <- 'model_outputs/sen_joint_cv_1.RData' 
cv1_prgp_sen_path <- 'model_outputs/sen_pr_gp_cv_1.RData' 

### CV 3 output

#cv3_points_sen_path <- 'model_outputs/sen_points_cv_3.RData'
cv3_polys_sen_path <- 'model_outputs/sen_polygon_cv_3.RData'
cv3_both_sen_path <- 'model_outputs/sen_joint_cv_3.RData'
cv3_prgp_sen_path <- 'model_outputs/sen_pr_gp_cv_3.RData'

### Full models
sen_baseline <- 'model_outputs/full_baseline_sen.RData'
sen_prgp <- 'model_outputs/full_prevgp_sen.RData'
sen_joint <- 'model_outputs/full_model_sen.RData'

## MDG


### Full data object

full_data_mdg_path <- 'model_outputs/mdg_full_data.RData'


### Cross validation object

data_cv1_mdg_path <- 'model_outputs/mdg_cv_1.RData'
#data_cv2_mdg_path <- 'model_outputs/mdg_cv_2.RData'
data_cv3_mdg_path <- 'model_outputs/mdg_cv_3.RData'

### CV 1 output

#cv1_points_mdg_path <- 'model_outputs/mdg_points_cv_1.RData'
cv1_polys_mdg_path <- 'model_outputs/mdg_polygon_cv_1.RData'
cv1_both_mdg_path <- 'model_outputs/mdg_joint_cv_1.RData'
cv1_prgp_mdg_path <- 'model_outputs/mdg_pr_gp__cv_1.RData'


### CV 3 output

#cv3_points_mdg_path <- 'model_outputs/mdg_points_cv_3.RData'
cv3_polys_mdg_path <- 'model_outputs/mdg_polygon_cv_3.RData'
cv3_both_mdg_path <- 'model_outputs/mdg_joint_cv_3.RData'
cv3_prgp_mdg_path <- 'model_outputs/mdg_pr_gp_cv_3.RData'

### Full models
mdg_baseline <- 'model_outputs/full_baseline_mdg.RData'
mdg_prgp <- 'model_outputs/full_prevgp_mdg.RData'
mdg_joint <- 'model_outputs/full_model_mdg.RData'


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

png('figs/summaries/random_crossvalidation_full.png', height = 200, width = 200, unit = 'mm', res = 720)
print(full_plot)
dev.off()





p1 <- autoplot(data_cv1_idn, jitter = 0, size = 0.3) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')
p2 <- autoplot(data_cv1_sen, jitter = 0, size = 0.3) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')
p3 <- autoplot(data_cv1_mdg, jitter = 0, size = 0.3) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')

bottom_row <- plot_grid(p2, p3, labels = c('B', 'C'), rel_widths = c(0.6, 0.4))

full_plot <- plot_grid(p1, bottom_row, ncol = 1, labels = c('A', ''), rel_heights = c(0.45, 0.55))

png('figs/summaries/random_crossvalidation_full_big.png', height = 100, width = 100, unit = 'mm', res = 720)
print(full_plot)
dev.off()



# Fig 2 - spatial cross validation for each country
data_cv3_idn <- get(load(data_cv3_idn_path))
data_cv3_sen <- get(load(data_cv3_sen_path))
data_cv3_mdg <- get(load(data_cv3_mdg_path))


p1 <- autoplot(data_cv3_idn, jitter = 0, size = 0.7) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')
p2 <- autoplot(data_cv3_sen, jitter = 0, size = 0.7) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')
p3 <- autoplot(data_cv3_mdg, jitter = 0, size = 0.7) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')

bottom_row <- plot_grid(p2, p3, labels = c('B', 'C'), rel_widths = c(0.6, 0.4))

full_plot <- plot_grid(p1, bottom_row, ncol = 1, labels = c('A', ''), rel_heights = c(0.45, 0.55))

png('figs/summaries/spatial_crossvalidation_full.png', height = 200, width = 200, unit = 'mm', res = 720)
print(full_plot)
dev.off()

gc()




p1 <- autoplot(data_cv3_idn, jitter = 0, size = 0.4) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')
p2 <- autoplot(data_cv3_sen, jitter = 0, size = 0.4) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')
p3 <- autoplot(data_cv3_mdg, jitter = 0, size = 0.4) + 
  guides(fill = FALSE) + 
  labs(x = 'Longitude', y = 'Latitude')

bottom_row <- plot_grid(p2, p3, labels = c('B', 'C'), rel_widths = c(0.6, 0.4))

full_plot <- plot_grid(p1, bottom_row, ncol = 1, labels = c('A', ''), rel_heights = c(0.45, 0.55))

png('figs/summaries/spatial_crossvalidation_full_big.png', height = 100, width = 100, unit = 'mm', res = 720)
print(full_plot)
dev.off()

# figure 3 data and predicted incidence maps. Indonesia only. Data, Rand, Spatial for best model? Joint model?
# todo add prevalence points

# Fig 3 - IDN: a) Data, predicted incidence from joint model for b) M2 , and c) M3
full_data_idn <- get(load(full_data_idn_path))
cv3_both_idn <- get(load(cv3_both_idn_path))
cv3_prgp_idn <- get(load(cv3_prgp_idn_path))

p1 <- obspred_map(data_cv3_idn, cv3_prgp_idn, trans = 'log1p',
                  legend_title = 'Cases per 1000',
                  breaks = c(1, 10, 100, 300, 500), mask = TRUE)
p2 <- obspred_map(data_cv3_idn, cv3_both_idn, trans = 'log1p', 
                  legend_title = 'API',
                  mask = TRUE)

WorldData <- map_data('world') %>% filter(region != "Antarctica") %>% fortify
idn_extent <- extent(cv3_both_idn$models[[1]]$predictions$api)
WorldData <- WorldData %>% 
               filter(region == 'Indonesia')

panel1 <- p1[[1]] +
  guides(fill = FALSE) +
  labs(x = '', y = '') + 
  geom_map(data = WorldData, map = WorldData,
           aes(x = long, y = lat, group = group, map_id=region),
           fill = NA, colour = 'darkgrey', size=0.5) +
  lims(x = idn_extent[1:2], y = idn_extent[3:4]) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))




panel2 <- p1[[2]] +
  guides(fill = FALSE) +
  geom_map(data = WorldData, map = WorldData,
           aes(x = long, y = lat, group = group, map_id=region),
           fill = NA, colour = 'darkgrey', size=0.5) +
  labs(x = '', y = 'Latitude')+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

panel3 <- p2[[2]] +
  geom_map(data = WorldData, map = WorldData,
           aes(x = long, y = lat, group = group, map_id=region),
           fill = NA, colour = 'darkgrey', size=0.5) +
  guides(fill = FALSE) +
  labs(x = 'Longitude', y = '')+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

legend <- get_legend(p1[[1]])

idn_preds_plot <- plot_grid(panel1, panel2, panel3, labels = LETTERS[1:3], ncol = 1)
full_plot <- plot_grid(idn_preds_plot, legend, ncol = 2, rel_widths = c(4, 1.2))

png('figs/summaries/idn_both_cv12_preds.png', height = 130, width = 110, unit = 'mm', res = 720)
print(full_plot)
dev.off()


#rm(full_data_idn)
rm(cv3_both_idn)
rm(cv3_prgp_idn)
gc()




# figure 4 data and predicted incidence maps. Senegal only. Data, Rand, Spatial for best model? Joint model?

# Fig 4 - SEN: a) Data, predicted incidence from joint model for b) random cv, and c) spatial cv
full_data_sen <- get(load(full_data_sen_path))
cv3_both_sen <- get(load(cv3_both_sen_path))
cv3_prgp_sen <- get(load(cv3_prgp_sen_path))


sen_y_max <- 
  max(
   c(full_data_sen$polygon$response),
     sapply(seq_along(cv3_both_sen$models), 
       function(i) maxValue(cv3_both_sen$models[[i]]$predictions$api)),
     sapply(seq_along(cv3_both_sen$models), 
       function(i) maxValue(cv3_prgp_sen$models[[i]]$predictions$api))
  )



p1 <- obspred_map(data_cv3_sen, cv3_prgp_sen, trans = 'log1p',
                  legend_title = 'API',
                  mask = TRUE,
                  lims = c(0, sen_y_max),
                  breaks = c(1, 10, 100, 300))
p2 <- obspred_map(data_cv3_sen, cv3_both_sen,
                  trans = 'log1p',
                  legend_title = 'API',
                  lims = c(0, sen_y_max),
                  mask = TRUE)

WorldData <- map_data('world') %>% filter(region != "Antarctica") %>% fortify
sen_extent <- extent(cv3_both_sen$models[[1]]$predictions$api)
WorldData <- WorldData %>% 
               filter(region == 'Senegal')

panel1 <- p1[[1]] +
  scale_fill_viridis_c(trans = 'log1p', 
                       limits = c(0, sen_y_max), 
                       oob = scales::squish, 
                       name = 'Cases per 1000') +
  geom_map(data = WorldData, map = WorldData,
           aes(x = long, y = lat, group = group, map_id=region),
           fill = NA, colour = 'darkgrey', size=0.5) +
  guides(fill = FALSE) +
  labs(x = '', y = 'Latitude') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

panel2 <- p1[[2]] +
  scale_fill_viridis_c(trans = 'log1p', 
                       limits = c(0, sen_y_max), 
                       oob = scales::squish, 
                       name = 'Cases per 1000') +
  geom_map(data = WorldData, map = WorldData,
           aes(x = long, y = lat, group = group, map_id=region),
           fill = NA, colour = 'darkgrey', size=0.5) +  
  guides(fill = FALSE) +
  labs(x = 'Longitude', y = '') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
panel3 <- p2[[2]] +
  scale_fill_viridis_c(trans = 'log1p', 
                       limits = c(0, sen_y_max), 
                       oob = scales::squish, 
                       name = 'Cases per 1000') +
  geom_map(data = WorldData, map = WorldData,
           aes(x = long, y = lat, group = group, map_id=region),
           fill = NA, colour = 'darkgrey', size=0.5) +  
  guides(fill = FALSE) +
  labs(x = '', y = '')+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

legend <- get_legend(p1[[1]] +
                       scale_fill_viridis_c(trans = 'log1p', 
                                            limits = c(0, sen_y_max), 
                                            oob = scales::squish, 
                                            name = 'Cases per 1000',
                                            breaks = c(0, 10, 50, 100, 200)))

sen_preds_plot <- plot_grid(panel1, panel2, panel3, labels = LETTERS[1:3], ncol = 3)
full_plot <- plot_grid(sen_preds_plot, legend, ncol = 2, rel_widths = c(5, 1))



png('figs/summaries/sen_both_cv12_preds.png', height = 90, width = 300, unit = 'mm', res = 720)
print(full_plot)
dev.off()


#rm(full_data_sen)
rm(cv3_both_sen)
rm(cv3_prgp_sen)
gc()


# Fig 5 - MDG: a) Data, predicted incidence from joint model for b) random cv, and c) spatial cv
full_data_mdg <- get(load(full_data_mdg_path))
cv3_both_mdg <- get(load(cv3_both_mdg_path))
cv3_prgp_mdg <- get(load(cv3_prgp_mdg_path))




mdg_y_max <- 
  max(
    c(full_data_mdg$polygon$response),
    sapply(seq_along(cv3_both_mdg$models), 
           function(i) maxValue(cv3_both_mdg$models[[i]]$predictions$api)),
    sapply(seq_along(cv3_both_mdg$models), 
           function(i) maxValue(cv3_prgp_mdg$models[[i]]$predictions$api))
  )

p1 <- obspred_map(data_cv3_mdg, cv3_prgp_mdg, trans = 'log1p',
                  legend_title = 'API',
                  mask = TRUE,
                  lims = c(0, sen_y_max),
                  breaks = c(1, 10, 100, 300))
p2 <- obspred_map(data_cv3_mdg, cv3_both_mdg,
                  trans = 'log1p',
                  legend_title = 'API',
                  lims = c(0, sen_y_max),
                  mask = TRUE)


panel1 <- p1[[1]] +
  scale_fill_viridis_c(trans = 'log1p', 
                       limits = c(0, mdg_y_max), 
                       oob = scales::squish, 
                       name = 'Cases per 1000') +
  guides(fill = FALSE) +
  labs(x = '', y = 'Latitude') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
panel2 <- p1[[2]] +
  scale_fill_viridis_c(trans = 'log1p', 
                       limits = c(0, mdg_y_max), 
                       oob = scales::squish, 
                       name = 'Cases per 1000') +
  
  guides(fill = FALSE) +
  labs(x = 'Longitude', y = '') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
panel3 <- p2[[2]] +
  scale_fill_viridis_c(trans = 'log1p', 
                       limits = c(0, mdg_y_max), 
                       oob = scales::squish, 
                       name = 'Cases per 1000') +
  
  guides(fill = FALSE) +
  labs(x = '', y = '')+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

legend <- get_legend(p1[[1]] +
                       scale_fill_viridis_c(trans = 'log1p', 
                                            limits = c(0, mdg_y_max), 
                                            oob = scales::squish, 
                                            name = 'Cases per 1000',
                                            breaks = c(0, 10, 25, 100, 250, 500)))

mdg_preds_plot <- plot_grid(panel1, panel2, panel3, labels = LETTERS[1:3], ncol = 3)
full_plot <- plot_grid(mdg_preds_plot, legend, ncol = 2, rel_widths = c(5, 1))




png('figs/summaries/mdg_both_cv12_preds.png', height = 150, width = 300, unit = 'mm', res = 720)
print(full_plot)
dev.off()


#rm(full_data_mdg)
rm(cv1_both_mdg)
rm(cv2_both_mdg)
gc()




# ------------------------------------------------------------------- #


# figure 3 data and predicted incidence maps. Indonesia only. Data, Rand, Spatial for best model? Joint model?
# todo add prevalence points

# Fig 3 - IDN: a) Data, predicted incidence from joint model for b) M2 , and c) M3
cv1_both_idn <- get(load(cv1_both_idn_path))
cv1_prgp_idn <- get(load(cv1_prgp_idn_path))

p1 <- obspred_map(data_cv1_idn, cv1_prgp_idn, trans = 'log1p',
                  legend_title = 'API',
                  breaks = c(1, 10, 100, 300, 500))
p2 <- obspred_map(data_cv1_idn, cv1_both_idn, trans = 'log1p', legend_title = 'API')


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

png('figs/summaries/idn_both_cv12_random_preds.png', height = 130, width = 100, unit = 'mm', res = 720)
print(full_plot)
dev.off()


#rm(full_data_idn)
rm(cv1_both_idn)
rm(cv1_prgp_idn)
gc()




# figure 4 data and predicted incidence maps. Senegal only. Data, Rand, Spatial for best model? Joint model?

# Fig 4 - SEN: a) Data, predicted incidence from joint model for b) random cv, and c) spatial cv
cv1_both_sen <- get(load(cv1_both_sen_path))
cv1_prgp_sen <- get(load(cv1_prgp_sen_path))


sen_y_max <- 
  max(
    c(full_data_sen$polygon$response),
    sapply(seq_along(cv1_both_sen$models), 
           function(i) maxValue(cv1_both_sen$models[[i]]$predictions$api)),
    sapply(seq_along(cv3_both_sen$models), 
           function(i) maxValue(cv1_prgp_sen$models[[i]]$predictions$api))
  )

p1 <- obspred_map(data_cv1_sen, cv1_prgp_sen, trans = 'log1p',
                  legend_title = 'API',
                  mask = TRUE,
                  lims = c(0, sen_y_max),
                  breaks = c(1, 10, 100, 300))
p2 <- obspred_map(data_cv1_sen, cv1_both_sen,
                  trans = 'log1p',
                  legend_title = 'API',
                  lims = c(0, sen_y_max),
                  mask = TRUE)


panel1 <- p1[[1]] +
  scale_fill_viridis_c(trans = 'log1p', 
                       limits = c(0, sen_y_max), 
                       oob = scales::squish, 
                       name = 'Cases per 1000') +
  guides(fill = FALSE) +
  labs(x = '', y = 'Latitude') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
panel2 <- p1[[2]] +
  scale_fill_viridis_c(trans = 'log1p', 
                       limits = c(0, sen_y_max), 
                       oob = scales::squish, 
                       name = 'Cases per 1000') +
  
  guides(fill = FALSE) +
  labs(x = 'Longitude', y = '') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
panel3 <- p2[[2]] +
  scale_fill_viridis_c(trans = 'log1p', 
                       limits = c(0, sen_y_max), 
                       oob = scales::squish, 
                       name = 'Cases per 1000') +
  
  guides(fill = FALSE) +
  labs(x = '', y = '')+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

legend <- get_legend(p1[[1]] +
                       scale_fill_viridis_c(trans = 'log1p', 
                                            limits = c(0, sen_y_max), 
                                            oob = scales::squish, 
                                            name = 'Cases per 1000',
                                            breaks = c(0, 10, 50, 100, 200)))

sen_preds_plot <- plot_grid(panel1, panel2, panel3, labels = LETTERS[1:3], ncol = 3)
full_plot <- plot_grid(sen_preds_plot, legend, ncol = 2, rel_widths = c(5, 1))



png('figs/summaries/sen_both_cv12_random_preds.png', height = 90, width = 300, unit = 'mm', res = 720)
print(full_plot)
dev.off()


#rm(full_data_sen)
rm(cv1_both_sen)
rm(cv1_prgp_sen)
gc()


# Fig 5 - MDG: a) Data, predicted incidence from joint model for b) random cv, and c) spatial cv
cv1_both_mdg <- get(load(cv1_both_mdg_path))
cv1_prgp_mdg <- get(load(cv1_prgp_mdg_path))




mdg_y_max <- 
  max(
    c(full_data_mdg$polygon$response),
    sapply(seq_along(cv1_both_mdg$models), 
           function(i) maxValue(cv1_both_mdg$models[[i]]$predictions$api)),
    sapply(seq_along(cv1_both_mdg$models), 
           function(i) maxValue(cv1_prgp_mdg$models[[i]]$predictions$api))
  )

p1 <- obspred_map(data_cv1_mdg, cv1_prgp_mdg, trans = 'log1p',
                  legend_title = 'API',
                  mask = TRUE,
                  lims = c(0, sen_y_max),
                  breaks = c(1, 10, 100, 300))
p2 <- obspred_map(data_cv1_mdg, cv1_both_mdg,
                  trans = 'log1p',
                  legend_title = 'API',
                  lims = c(0, sen_y_max),
                  mask = TRUE)


panel1 <- p1[[1]] +
  scale_fill_viridis_c(trans = 'log1p', 
                       limits = c(0, mdg_y_max), 
                       oob = scales::squish, 
                       name = 'Cases per 1000') +
  guides(fill = FALSE) +
  labs(x = '', y = 'Latitude') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
panel2 <- p1[[2]] +
  scale_fill_viridis_c(trans = 'log1p', 
                       limits = c(0, mdg_y_max), 
                       oob = scales::squish, 
                       name = 'Cases per 1000') +
  
  guides(fill = FALSE) +
  labs(x = 'Longitude', y = '') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
panel3 <- p2[[2]] +
  scale_fill_viridis_c(trans = 'log1p', 
                       limits = c(0, mdg_y_max), 
                       oob = scales::squish, 
                       name = 'Cases per 1000') +
  
  guides(fill = FALSE) +
  labs(x = '', y = '')+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

legend <- get_legend(p1[[1]] +
                       scale_fill_viridis_c(trans = 'log1p', 
                                            limits = c(0, mdg_y_max), 
                                            oob = scales::squish, 
                                            name = 'Cases per 1000',
                                            breaks = c(0, 10, 25, 100, 250, 500)))

mdg_preds_plot <- plot_grid(panel1, panel2, panel3, labels = LETTERS[1:3], ncol = 3)
full_plot <- plot_grid(mdg_preds_plot, legend, ncol = 2, rel_widths = c(5, 1))




png('figs/summaries/mdg_both_cv12_random_preds.png', height = 150, width = 300, unit = 'mm', res = 720)
print(full_plot)
dev.off()


#rm(full_data_mdg)
rm(cv1_both_mdg)
rm(cv1_prgp_mdg)
gc()






# figure 6, random cv. PR vs Poly columns, countries as rows, model as colour?

#################################################################################
## Make model comparisons. Random CV.                                          ##
#################################################################################


#cv1_points_idn <- get(load(cv1_points_idn_path))
cv1_polys_idn <- get(load(cv1_polys_idn_path))
cv1_prgp_idn <- get(load(cv1_prgp_idn_path))
cv1_both_idn <- get(load(cv1_both_idn_path))


idn_cv1_poly_df <- rbind(cv1_polys_idn$summary$combined_aggregated %>% cbind(model = 'Baseline'),
                         cv1_prgp_idn$summary$combined_aggregated %>% cbind(model = 'Prev GP'),
                         cv1_both_idn$summary$combined_aggregated %>% cbind(model = 'Joint'))
idn_cv1_pr_df <-  rbind(cv1_polys_idn$summary$combined_pr %>% cbind(model = 'Baseline'),
                        cv1_prgp_idn$summary$combined_pr %>% cbind(model = 'Prev GP'),                        
                        cv1_both_idn$summary$combined_pr %>% cbind(model = 'Joint'))
# 
# 
# idn_cv1_new_metrics <- 
#   rbind(cv1_points_idn$summary$combined_pr %>%
#           mutate(abs_error = abs(prevalence - pred_prev)) %>% 
#           summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
#           cbind(model = 'points'),
#         cv1_polys_idn$summary$combined_pr %>%
#           mutate(abs_error = abs(prevalence - pred_prev)) %>% 
#           summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
#           cbind(model = 'polygons'),
#         cv1_both_idn$summary$combined_pr %>%
#           mutate(abs_error = abs(prevalence - pred_prev)) %>% 
#           summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
#           cbind(model = 'both')
#   )
# 



idn_cv1_metrics <- list(rbind(cv1_polys_idn$summary$polygon_metrics %>% cbind(model = 'Baseline'),
                              cv1_prgp_idn$summary$polygon_metrics %>% cbind(model = 'Prev GP'),
                              cv1_both_idn$summary$polygon_metrics %>% cbind(model = 'Joint')),
                        rbind(cv1_polys_idn$summary$pr_metrics %>% cbind(model = 'Baseline'),
                              cv1_prgp_idn$summary$pr_metrics %>% cbind(model = 'Prev GP'),                              
                              cv1_both_idn$summary$pr_metrics %>% cbind(model = 'Joint')))


rm(cv1_prgp_idn)
rm(cv1_polys_idn)
rm(cv1_both_idn)
gc()


cv1_prgp_sen <- get(load(cv1_prgp_sen_path))
cv1_polys_sen <- get(load(cv1_polys_sen_path))
cv1_both_sen <- get(load(cv1_both_sen_path))


sen_cv1_poly_df <- rbind(cv1_polys_sen$summary$combined_aggregated %>% cbind(model = 'Baseline'),
                         cv1_prgp_sen$summary$combined_aggregated %>% cbind(model = 'Prev GP'),
                         cv1_both_sen$summary$combined_aggregated %>% cbind(model = 'Joint'))
sen_cv1_pr_df <-  rbind(cv1_polys_sen$summary$combined_pr %>% cbind(model = 'Baseline'),
                        cv1_prgp_sen$summary$combined_pr %>% cbind(model = 'Prev GP'),
                        cv1_both_sen$summary$combined_pr %>% cbind(model = 'Joint'))
# 
# sen_cv1_new_metrics <- 
#   rbind(cv1_points_sen$summary$combined_pr %>%
#           mutate(abs_error = abs(prevalence - pred_prev)) %>% 
#           summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
#           cbind(model = 'points'),
#         cv1_polys_sen$summary$combined_pr %>%
#           mutate(abs_error = abs(prevalence - pred_prev)) %>% 
#           summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
#           cbind(model = 'polygons'),
#         cv1_both_sen$summary$combined_pr %>%
#           mutate(abs_error = abs(prevalence - pred_prev)) %>% 
#           summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
#           cbind(model = 'both')
#   )




sen_cv1_metrics <- list(rbind(cv1_polys_sen$summary$polygon_metrics %>% cbind(model = 'Baseline'),
                              cv1_prgp_sen$summary$polygon_metrics %>% cbind(model = 'Prev GP'),
                              cv1_both_sen$summary$polygon_metrics %>% cbind(model = 'Joint')),
                        rbind(cv1_polys_sen$summary$pr_metrics %>% cbind(model = 'Baseline'),
                              cv1_prgp_sen$summary$pr_metrics %>% cbind(model = 'Prev GP'),
                              cv1_both_sen$summary$pr_metrics %>% cbind(model = 'Joint')))



rm(cv1_points_sen)
rm(cv1_polys_sen)
rm(cv1_both_sen)
gc()



cv1_prgp_mdg <- get(load(cv1_prgp_mdg_path))
cv1_polys_mdg <- get(load(cv1_polys_mdg_path))
cv1_both_mdg <- get(load(cv1_both_mdg_path))



mdg_cv1_poly_df <- rbind(cv1_polys_mdg$summary$combined_aggregated %>% cbind(model = 'Baseline'),
                         cv1_prgp_mdg$summary$combined_aggregated %>% cbind(model = 'Prev GP'),
                         cv1_both_mdg$summary$combined_aggregated %>% cbind(model = 'Joint'))
mdg_cv1_pr_df <-  rbind(cv1_polys_mdg$summary$combined_pr %>% cbind(model = 'Baseline'),
                        cv1_prgp_mdg$summary$combined_pr %>% cbind(model = 'Prev GP'),                        
                        cv1_both_mdg$summary$combined_pr %>% cbind(model = 'Joint'))



mdg_cv1_metrics <- list(rbind(cv1_polys_mdg$summary$polygon_metrics %>% cbind(model = 'Baseline'),
                              cv1_prgp_mdg$summary$polygon_metrics %>% cbind(model = 'Prev GP'),
                              cv1_both_mdg$summary$polygon_metrics %>% cbind(model = 'Joint')),
                        rbind(cv1_polys_mdg$summary$pr_metrics %>% cbind(model = 'Baseline'),
                              cv1_prgp_mdg$summary$pr_metrics %>% cbind(model = 'Prev GP'),
                              cv1_both_mdg$summary$pr_metrics %>% cbind(model = 'Joint')))


rm(cv1_points_mdg)
rm(cv1_polys_mdg)
rm(cv1_both_mdg)
gc()




a <- 0.5
s <- 16

idn_poly <- ggplot(idn_cv1_poly_df, aes(pred_api, response, colour = model)) + 
  geom_point(alpha = a, size = 3, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  #scale_y_sqrt() + 
  #scale_x_sqrt() +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)

idn_poly_facet <- ggplot(idn_cv1_poly_df, aes(pred_api, response, colour = model)) + 
  geom_point(alpha = a, size = 2, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  #scale_y_sqrt() + 
  #scale_x_sqrt() +
  facet_wrap(~model) + 
  geom_smooth(method = 'lm', aes(colour = model), se = FALSE) +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)



sen_poly <- ggplot(sen_cv1_poly_df, aes(pred_api, response, colour = model)) + 
  geom_point(alpha = a + 0.2, size = 2, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  #scale_y_sqrt() + 
  #scale_x_sqrt() +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)

sen_poly_facet <- ggplot(sen_cv1_poly_df, aes(pred_api, response, colour = model)) + 
  geom_point(alpha = a + 0.2, size = 2, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  #scale_y_sqrt() + 
  #scale_x_sqrt() +
  facet_wrap(~ model) +
  geom_smooth(method = 'lm', aes(colour = model), se = FALSE) +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)


mdg_poly <- ggplot(mdg_cv1_poly_df, aes(pred_api, response, colour = model)) + 
  geom_point(alpha = a + 0.2, size = 2, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  #scale_y_sqrt() + 
  #scale_x_sqrt() +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)


mdg_poly_facet <- ggplot(mdg_cv1_poly_df, aes(pred_api, response, colour = model)) + 
  geom_point(alpha = a + 0.2, size = 2, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  #scale_y_sqrt() + 
  #scale_x_sqrt() +  
  facet_wrap(~ model) + 
  geom_smooth(method = 'lm', aes(colour = model), se = FALSE) +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)



poly_facet_list <- list(idn_poly_facet, sen_poly_facet, mdg_poly_facet)
full_obs_pred_poly_facet <- plot_grid(plotlist = poly_facet_list,
                                       labels = c('A) Indonesia', 'B) Senegal', 'C) Madagascar'), ncol = 1, hjust = -0.1, label_size = 13)

png('figs/summaries/random_cv_poly_facet.png', height = 1500, width = 1100)
print(full_obs_pred_poly_facet)
dev.off()



pdf('figs/summaries/random_cv_poly_facet.pdf', height = 10, width = 7)
print(full_obs_pred_poly_facet)
dev.off()


pdf('figs/summaries/random_cv_poly_facet_big.pdf', height = 7, width = 6)
print(full_obs_pred_poly_facet)
dev.off()




# Fig 4 and 6
#################################################################################
## Make model comparisons. Spatial CV keeppr                                   ##
#################################################################################


#cv3_points_idn <- get(load(cv3_points_idn_path))
cv3_polys_idn <- get(load(cv3_polys_idn_path))
cv3_prgp_idn <- get(load(cv3_prgp_idn_path))
cv3_both_idn <- get(load(cv3_both_idn_path))


idn_cv3_poly_df <- rbind(cv3_polys_idn$summary$combined_aggregated %>% cbind(model = 'Baseline'),
                         cv3_prgp_idn$summary$combined_aggregated %>% cbind(model = 'Prev GP'),
                         cv3_both_idn$summary$combined_aggregated %>% cbind(model = 'Joint'))
idn_cv3_pr_df <-  rbind(cv3_polys_idn$summary$combined_pr %>% cbind(model = 'Baseline'),
                        cv3_prgp_idn$summary$combined_pr %>% cbind(model = 'Prev GP'),                        
                        cv3_both_idn$summary$combined_pr %>% cbind(model = 'Joint'))
# 
# 
# idn_cv3_new_metrics <- 
#   rbind(cv3_points_idn$summary$combined_pr %>%
#           mutate(abs_error = abs(prevalence - pred_prev)) %>% 
#           summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
#           cbind(model = 'points'),
#         cv3_polys_idn$summary$combined_pr %>%
#           mutate(abs_error = abs(prevalence - pred_prev)) %>% 
#           summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
#           cbind(model = 'polygons'),
#         cv3_both_idn$summary$combined_pr %>%
#           mutate(abs_error = abs(prevalence - pred_prev)) %>% 
#           summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
#           cbind(model = 'both')
#   )
# 



idn_cv3_metrics <- list(rbind(cv3_polys_idn$summary$polygon_metrics %>% cbind(model = 'Baseline'),
                              cv3_prgp_idn$summary$polygon_metrics %>% cbind(model = 'Prev GP'),
                              cv3_both_idn$summary$polygon_metrics %>% cbind(model = 'Joint')),
                        rbind(cv3_polys_idn$summary$pr_metrics %>% cbind(model = 'Baseline'),
                              cv3_prgp_idn$summary$pr_metrics %>% cbind(model = 'Prev GP'),                              
                              cv3_both_idn$summary$pr_metrics %>% cbind(model = 'Joint')))


rm(cv3_prgp_idn)
rm(cv3_polys_idn)
rm(cv3_both_idn)
gc()


cv3_prgp_sen <- get(load(cv3_prgp_sen_path))
cv3_polys_sen <- get(load(cv3_polys_sen_path))
cv3_both_sen <- get(load(cv3_both_sen_path))


sen_cv3_poly_df <- rbind(cv3_polys_sen$summary$combined_aggregated %>% cbind(model = 'Baseline'),
                         cv3_prgp_sen$summary$combined_aggregated %>% cbind(model = 'Prev GP'),
                         cv3_both_sen$summary$combined_aggregated %>% cbind(model = 'Joint'))
sen_cv3_pr_df <-  rbind(cv3_polys_sen$summary$combined_pr %>% cbind(model = 'Baseline'),
                        cv3_prgp_sen$summary$combined_pr %>% cbind(model = 'Prev GP'),
                        cv3_both_sen$summary$combined_pr %>% cbind(model = 'Joint'))
# 
# sen_cv3_new_metrics <- 
#   rbind(cv3_points_sen$summary$combined_pr %>%
#           mutate(abs_error = abs(prevalence - pred_prev)) %>% 
#           summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
#           cbind(model = 'points'),
#         cv3_polys_sen$summary$combined_pr %>%
#           mutate(abs_error = abs(prevalence - pred_prev)) %>% 
#           summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
#           cbind(model = 'polygons'),
#         cv3_both_sen$summary$combined_pr %>%
#           mutate(abs_error = abs(prevalence - pred_prev)) %>% 
#           summarise(weightMAE = weighted.mean(abs_error, w = examined)) %>% 
#           cbind(model = 'both')
#   )




sen_cv3_metrics <- list(rbind(cv3_polys_sen$summary$polygon_metrics %>% cbind(model = 'Baseline'),
                              cv3_prgp_sen$summary$polygon_metrics %>% cbind(model = 'Prev GP'),
                              cv3_both_sen$summary$polygon_metrics %>% cbind(model = 'Joint')),
                        rbind(cv3_polys_sen$summary$pr_metrics %>% cbind(model = 'Baseline'),
                              cv3_prgp_sen$summary$pr_metrics %>% cbind(model = 'Prev GP'),
                              cv3_both_sen$summary$pr_metrics %>% cbind(model = 'Joint')))



rm(cv3_points_sen)
rm(cv3_polys_sen)
rm(cv3_both_sen)
gc()



cv3_prgp_mdg <- get(load(cv3_prgp_mdg_path))
cv3_polys_mdg <- get(load(cv3_polys_mdg_path))
cv3_both_mdg <- get(load(cv3_both_mdg_path))



mdg_cv3_poly_df <- rbind(cv3_polys_mdg$summary$combined_aggregated %>% cbind(model = 'Baseline'),
                         cv3_prgp_mdg$summary$combined_aggregated %>% cbind(model = 'Prev GP'),
                         cv3_both_mdg$summary$combined_aggregated %>% cbind(model = 'Joint'))
mdg_cv3_pr_df <-  rbind(cv3_polys_mdg$summary$combined_pr %>% cbind(model = 'Baseline'),
                        cv3_prgp_mdg$summary$combined_pr %>% cbind(model = 'Prev GP'),                        
                        cv3_both_mdg$summary$combined_pr %>% cbind(model = 'Joint'))



mdg_cv3_metrics <- list(rbind(cv3_polys_mdg$summary$polygon_metrics %>% cbind(model = 'Baseline'),
                              cv3_prgp_mdg$summary$polygon_metrics %>% cbind(model = 'Prev GP'),
                              cv3_both_mdg$summary$polygon_metrics %>% cbind(model = 'Joint')),
                        rbind(cv3_polys_mdg$summary$pr_metrics %>% cbind(model = 'Baseline'),
                              cv3_prgp_mdg$summary$pr_metrics %>% cbind(model = 'Prev GP'),
                              cv3_both_mdg$summary$pr_metrics %>% cbind(model = 'Joint')))


rm(cv3_points_mdg)
rm(cv3_polys_mdg)
rm(cv3_both_mdg)
gc()




a <- 0.5
s <- 16

idn_poly <- ggplot(idn_cv3_poly_df, aes(pred_api, response, colour = model)) + 
  geom_point(alpha = a, size = 3, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  #scale_y_sqrt() + 
  #scale_x_sqrt() +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)

idn_poly_facet <- ggplot(idn_cv3_poly_df, aes(pred_api, response, colour = model)) + 
  geom_point(alpha = a, size = 2, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  #scale_y_sqrt() + 
  #scale_x_sqrt() +
  facet_wrap(~model) + 
  geom_smooth(method = 'lm', aes(colour = model), se = FALSE) +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)



sen_poly <- ggplot(sen_cv3_poly_df, aes(pred_api, response, colour = model)) + 
  geom_point(alpha = a + 0.2, size = 2, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  #scale_y_sqrt() + 
  #scale_x_sqrt() +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)

sen_poly_facet <- ggplot(sen_cv3_poly_df, aes(pred_api, response, colour = model)) + 
  geom_point(alpha = a + 0.2, size = 2, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  #scale_y_sqrt() + 
  #scale_x_sqrt() +
  facet_wrap(~ model) +
  geom_smooth(method = 'lm', aes(colour = model), se = FALSE) +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)


mdg_poly <- ggplot(mdg_cv3_poly_df, aes(pred_api, response, colour = model)) + 
  geom_point(alpha = a + 0.2, size = 2, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  #scale_y_sqrt() + 
  #scale_x_sqrt() +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)


mdg_poly_facet <- ggplot(mdg_cv3_poly_df, aes(pred_api, response, colour = model)) + 
  geom_point(alpha = a + 0.2, size = 2, shape = s) + 
  geom_abline(slope = 1, intercept = 0) +
  #scale_y_sqrt() + 
  #scale_x_sqrt() +  
  facet_wrap(~ model) + 
  geom_smooth(method = 'lm', aes(colour = model), se = FALSE) +
  labs(x = 'Observed API', y = 'Predicted API') +
  guides(colour = FALSE)



poly_facet_list <- list(idn_poly_facet, sen_poly_facet, mdg_poly_facet)
full_obs_pred_poly_facet <- plot_grid(plotlist = poly_facet_list,
                                       labels = c('A) Indonesia', 'B) Senegal', 'C) Madagascar'), ncol = 1, hjust = -0.1, label_size = 13)

png('figs/summaries/spatialkeeppr_cv_poly_facet.png', height = 1500, width = 1100)
print(full_obs_pred_poly_facet)
dev.off()



pdf('figs/summaries/spatialkeeppr_cv_poly_facet.pdf', height = 10, width = 8)
print(full_obs_pred_poly_facet)
dev.off()


pdf('figs/summaries/spatialkeeppr_cv_poly_facet_big.pdf', height = 7, width = 6)
print(full_obs_pred_poly_facet)
dev.off()







# Useful summary tables
table1_skeleton <- 
"Random & Indonesia  & %s &  %s &  %s\\\\
& Senegal  & %s &  %s &  %s\\\\
& Madagascar  & %s &  %s &  %s\\vspace{3mm}\\\\
Spatial & Indonesia & %s &  %s &  %s\\\\
& Senegal  & %s &  %s &  %s\\\\
& Madagascar & %s &  %s &  %s\\\\"



r <- c(idn_cv1_metrics[[1]]$MAE,
       sen_cv1_metrics[[1]]$MAE,
       mdg_cv1_metrics[[1]]$MAE,
       idn_cv3_metrics[[1]]$MAE, 
       sen_cv3_metrics[[1]]$MAE, 
       mdg_cv3_metrics[[1]]$MAE)

r <- format(round(r, 2), nsmall = 2)

table1 <- do.call(sprintf, c(table1_skeleton, as.list(r)))

write(table1, 'figs/summaries/table1.txt')




coverage_skeleton <-
"Random & Indonesia  & %s &  %s &  %s\\\\
& Senegal  & %s &  %s &  %s\\\\
& Madagascar  & %s &  %s &  %s\\vspace{1mm}\\\\
 Spatial  & Indonesia & %s &  %s &  %s\\\\
& Senegal  & %s &  %s &  %s\\\\
& Madagascar  & %s &  %s &  %s\\\\
"

cov <- c(idn_cv1_metrics[[1]]$coverage, sen_cv1_metrics[[1]]$coverage, mdg_cv1_metrics[[1]]$coverage,
         idn_cv3_metrics[[1]]$coverage, sen_cv3_metrics[[1]]$coverage, mdg_cv3_metrics[[1]]$coverage)

cov <- format(round(cov, 2), nsmall = 2)

table3 <- do.call(sprintf, c(coverage_skeleton, as.list(cov)))

write(table3, 'figs/summaries/table3.txt')


# Further SI figures.

## prevgp fields.

sen_full <- get(load(full_data_sen_path))
idn_full <- get(load(full_data_idn_path))
mdg_full <- get(load(full_data_mdg_path))

png('figs/summaries/sen_prevgp.png', height = 1000, width = 1100)
plot(sen_full$cov_rasters[[9]])
dev.off()
png('figs/summaries/idn_prevgp.png', height = 800, width = 1100)
plot(idn_full$cov_rasters[[9]])
dev.off()
png('figs/summaries/mdg_prevgp.png', height = 1500, width = 1100)
plot(mdg_full$cov_rasters[[9]])
dev.off()



## Disaggregation field

mdg_joint_model <- get(load(mdg_joint))
mdg_baseline_model <- get(load(mdg_baseline))
mdg_prgp_model <- get(load(mdg_prgp))


png('figs/summaries/mdg_baseline_field.png', height = 1500, width = 1100)
plot(mask(mdg_baseline_model$predictions$field, mdg_baseline_model$predictions$api))
dev.off()
png('figs/summaries/mdg_prgp_field.png', height = 1500, width = 1100)
plot(mask(mdg_prgp_model$predictions$field, mdg_baseline_model$predictions$api))
dev.off()
png('figs/summaries/mdg_joint_field.png', height = 1500, width = 1100)
plot(mask(mdg_joint_model$predictions$field, mdg_baseline_model$predictions$api))
dev.off()



sen_joint_model <- get(load(sen_joint))
sen_baseline_model <- get(load(sen_baseline))
sen_prgp_model <- get(load(sen_prgp))


png('figs/summaries/sen_baseline_field.png', height = 1100, width = 1100)
plot(mask(sen_baseline_model$predictions$field, sen_baseline_model$predictions$api))
dev.off()
png('figs/summaries/sen_prgp_field.png', height = 1100, width = 1100)
plot(mask(sen_prgp_model$predictions$field, sen_baseline_model$predictions$api))
dev.off()
png('figs/summaries/sen_joint_field.png', height = 1100, width = 1100)
plot(mask(sen_joint_model$predictions$field, sen_baseline_model$predictions$api))
dev.off()



idn_joint_model <- get(load(idn_joint))
idn_baseline_model <- get(load(idn_baseline))
idn_prgp_model <- get(load(idn_prgp))


png('figs/summaries/idn_baseline_field.png', height = 800, width = 1100)
plot(mask(idn_baseline_model$predictions$field, idn_baseline_model$predictions$api))
dev.off()
png('figs/summaries/idn_prgp_field.png', height = 800, width = 1100)
plot(mask(idn_prgp_model$predictions$field, idn_baseline_model$predictions$api))
dev.off()
png('figs/summaries/idn_joint_field.png', height = 800, width = 1100)
plot(mask(idn_joint_model$predictions$field, idn_joint_model$predictions$api))
dev.off()


## Full model predictions

png('figs/summaries/mdg_baseline_api.png', height = 1500, width = 1100)
plot(mask(mdg_baseline_model$predictions$api, mdg_baseline_model$predictions$api))
dev.off()
png('figs/summaries/mdg_prgp_api.png', height = 1500, width = 1100)
plot(mask(mdg_prgp_model$predictions$api, mdg_baseline_model$predictions$api))
dev.off()
png('figs/summaries/mdg_joint_api.png', height = 1500, width = 1100)
plot(mask(mdg_joint_model$predictions$api, mdg_baseline_model$predictions$api))
dev.off()



png('figs/summaries/sen_baseline_api.png', height = 1100, width = 1100)
plot(mask(sen_baseline_model$predictions$api, sen_baseline_model$predictions$api))
dev.off()
png('figs/summaries/sen_prgp_api.png', height = 1100, width = 1100)
plot(mask(sen_prgp_model$predictions$api, sen_baseline_model$predictions$api))
dev.off()
png('figs/summaries/sen_joint_api.png', height = 1100, width = 1100)
plot(mask(sen_joint_model$predictions$api, sen_baseline_model$predictions$api))
dev.off()



png('figs/summaries/idn_baseline_api.png', height = 800, width = 1100)
plot(mask(idn_baseline_model$predictions$api, idn_baseline_model$predictions$api))
dev.off()
png('figs/summaries/idn_prgp_api.png', height = 800, width = 1100)
plot(mask(idn_prgp_model$predictions$api, idn_baseline_model$predictions$api))
dev.off()
png('figs/summaries/idn_joint_api.png', height = 800, width = 1100)
plot(mask(idn_joint_model$predictions$api, idn_joint_model$predictions$api))
dev.off()



# parameters


sen_sd_baseline <- summary(sen_baseline_model$model$sd_report, select= 'fixed')[, 2]
sen_sd_prgp <- summary(sen_prgp_model$model$sd_report, select= 'fixed')[, 2]
sen_sd_joint <- summary(sen_joint_model$model$sd_report, select= 'fixed')[, 2]

sen_mean_baseline <- summary(sen_baseline_model$model$sd_report, select= 'fixed')[, 1]
sen_mean_prgp <- summary(sen_prgp_model$model$sd_report, select= 'fixed')[, 1]
sen_mean_joint <- summary(sen_joint_model$model$sd_report, select= 'fixed')[, 1]

sen_df <- data.frame(Parameters = c('$\\beta_0$', '$\\beta_p$', 
                                    'LST day mean', 'EVI', 'TSI', 'accessability', 'elevation', 'LST day SD', 'Night lights', 'TCW', 'Prev GP',
                                    '$\\omega_v$', '$\\omega_w$', '$\\log(\\sigma_u)$', '$\\log(\\rho)$', '$\\log(\\alpha)$'),
                     Baseline = c(sen_mean_baseline[1:10], NA, sen_mean_baseline[11], NA, sen_mean_baseline[12:14]),
                     BaselineSD = c(sen_sd_baseline[1:10], NA, sen_sd_baseline[11], NA, sen_sd_baseline[12:14]),
                     PrevalenceGP = c(sen_mean_prgp[1:12], NA, sen_mean_prgp[13:15]),
                     PrevalenceGPSD = c(sen_sd_prgp[1:12], NA, sen_sd_prgp[13:15]),
                     Joint = c(sen_mean_joint[1:10], NA, sen_mean_joint[11:15]),
                     JointSD = c(sen_sd_joint[1:10], NA, sen_sd_joint[11:15]))

sen_df[16, 2:5] <- NA
sen_df[2, 2:5] <- NA




tableS1 <- xtable(sen_df,
                  caption = 'Parameter estimates for Senegal. $\\omega_v = \\log\\left(\\frac{1}{{\\sigma_v}^2}\\right)$ and $\\omega_w = \\log\\left(\\frac{1}{{\\sigma_w}^2}\\right)$')

print(tableS1, file = 'figs/summaries/tableS1.txt', include.rownames = FALSE, sanitize.text.function = function(x)x)



mdg_sd_baseline <- summary(mdg_baseline_model$model$sd_report, select= 'fixed')[, 2]
mdg_sd_prgp <- summary(mdg_prgp_model$model$sd_report, select= 'fixed')[, 2]
mdg_sd_joint <- summary(mdg_joint_model$model$sd_report, select= 'fixed')[, 2]

mdg_mean_baseline <- summary(mdg_baseline_model$model$sd_report, select= 'fixed')[, 1]
mdg_mean_prgp <- summary(mdg_prgp_model$model$sd_report, select= 'fixed')[, 1]
mdg_mean_joint <- summary(mdg_joint_model$model$sd_report, select= 'fixed')[, 1]

mdg_df <- data.frame(Parameters = c('$\\beta_0$', '$\\beta_p$', 
                                    'LST day mean', 'EVI', 'TSI', 'accessability', 'elevation', 'LST day SD', 'Night lights', 'TCW', 'Prev GP',
                                    '$\\omega_v$', '$\\omega_w$', '$\\log(\\sigma_u)$', '$\\log(\\rho)$', '$\\log(\\alpha)$'),
                     Baseline = c(mdg_mean_baseline[1:10], NA, mdg_mean_baseline[11], NA, mdg_mean_baseline[12:14]),
                     BaselineSD = c(mdg_sd_baseline[1:10], NA, mdg_sd_baseline[11], NA, mdg_sd_baseline[12:14]),
                     PrevalenceGP = c(mdg_mean_prgp[1:12], NA, mdg_mean_prgp[13:15]),
                     PrevalenceGPSD = c(mdg_sd_prgp[1:12], NA, mdg_sd_prgp[13:15]),
                     Joint = c(mdg_mean_joint[1:10], NA, mdg_mean_joint[11:15]),
                     JointSD = c(mdg_sd_joint[1:10], NA, mdg_sd_joint[11:15]),
                     stringsAsFactors = FALSE)

mdg_df[16, 2:5] <- NA
mdg_df[2, 2:5] <- NA


tableS2 <- xtable(mdg_df,
                  caption = 'Parameter estimates for Madagascar. $\\omega_v = \\log\\left(\\frac{1}{{\\sigma_v}^2}\\right)$ and $\\omega_w = \\log\\left(\\frac{1}{{\\sigma_w}^2}\\right)$')

print(tableS2, file = 'figs/summaries/tableS2.txt', include.rownames = FALSE, sanitize.text.function = function(x)x)





idn_sd_baseline <- summary(idn_baseline_model$model$sd_report, select= 'fixed')[, 2]
idn_sd_prgp <- summary(idn_prgp_model$model$sd_report, select= 'fixed')[, 2]
idn_sd_joint <- summary(idn_joint_model$model$sd_report, select= 'fixed')[, 2]

idn_mean_baseline <- summary(idn_baseline_model$model$sd_report, select= 'fixed')[, 1]
idn_mean_prgp <- summary(idn_prgp_model$model$sd_report, select= 'fixed')[, 1]
idn_mean_joint <- summary(idn_joint_model$model$sd_report, select= 'fixed')[, 1]

idn_df <- data.frame(Parameters = c('$\\beta_0$', '$\\beta_p$', 
                                    'LST day mean', 'EVI', 'TSI', 'accessability', 'elevation', 'LST day SD', 'Night lights', 'TCW', 'Prev GP',
                                    '$\\omega_v$', '$\\omega_w$', '$\\log(\\sigma_u)$', '$\\log(\\rho)$', '$\\log(\\alpha)$'),
                     Baseline = c(idn_mean_baseline[1:10], NA, idn_mean_baseline[11], NA, idn_mean_baseline[12:14]),
                     BaselineSD = c(idn_sd_baseline[1:10], NA, idn_sd_baseline[11], NA, idn_sd_baseline[12:14]),
                     PrevalenceGP = c(idn_mean_prgp[1:12], NA, idn_mean_prgp[13:15]),
                     PrevalenceGPSD = c(idn_sd_prgp[1:12], NA, idn_sd_prgp[13:15]),
                     Joint = c(idn_mean_joint[1:10], NA, idn_mean_joint[11:15]),
                     JointSD = c(idn_sd_joint[1:10], NA, idn_sd_joint[11:15]),
                     stringsAsFactors = FALSE)

idn_df[16, 2:5] <- NA
idn_df[2, 2:5] <- NA


tableS3 <- xtable(idn_df,
                  caption = 'Parameter estimates for Indonesia. $\\omega_v = \\log\\left(\\frac{1}{{\\sigma_v}^2}\\right)$ and $\\omega_w = \\log\\left(\\frac{1}{{\\sigma_w}^2}\\right)$')

print(tableS3, file = 'figs/summaries/tableS3.txt', include.rownames = FALSE, sanitize.text.function = function(x)x)






# Data summaries.

full_data_mdg$pr %>% dim
full_data_mdg$pr$examined %>% sum
full_data_mdg$polygon %>% dim

full_data_idn$pr %>% dim
full_data_idn$pr$examined %>% sum
full_data_idn$polygon %>% dim

full_data_sen$pr %>% dim
full_data_sen$pr$examined %>% sum
full_data_sen$polygon %>% dim

