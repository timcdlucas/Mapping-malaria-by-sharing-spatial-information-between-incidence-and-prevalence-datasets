
# figure 1.cross validation. %% Do fig 1 and 2, random and spatial cv. IDN on top, MDG and SEN below in each.
# figure 3 and 4. data and predicted incidence maps. Indonesia and Senegal only. Fig 3 ind, fig 4 sen Data, Rand, Spatial for best model? Joint model?
# figure 5, 6. Spat and random cv. PR vs Poly columns, countries as rows, model as colour?



setwd('~/timz/timothy/point_polygon_joint_comparison')

# Libs

library(ggplot2)
library(cowplot)
source('plotting_functions.R')

# Paths

## IDN

### Full data object

full_data_idn_path <- 'model_outputs/idn_full_data.RData'

### Cross validation object

data_cv1_idn_path <- 'model_outputs/cv_1.RData'
data_cv2_idn_path <- 'model_outputs/cv_2.RData'

### CV 1 output

cv1_points_idn_path <- 'model_outputs/points_cv_1.RData'
cv1_polys_idn_path <- 'model_outputs/polygons_cv_1.RData'
cv1_both_idn_path <- 'model_outputs/join_cv_1.RData'

### CV 2 output

cv1_points_idn_path <- 'model_outputs/points_cv_1.RData'
cv1_polys_idn_path <- 'model_outputs/polygons_cv_1.RData'
cv1_both_idn_path <- 'model_outputs/join_cv_1.RData'


## SEN

### Full data object

full_data_sen_path <- 'model_outputs/mdg_full_data.RData'


### Cross validation object

data_cv1_sen_path <- 'model_outputs/mdg_cv_1.RData'
data_cv2_sen_path <- 'model_outputs/mdg_cv_2.RData'


### CV 1 output

cv1_points_sen_path <- 'model_outputs/mdg_points_cv_1.RData'
cv1_polys_sen_path <- 'model_outputs/mdg_polygon_cv_1.RData'
cv1_both_sen_path <- 'model_outputs/mdg_join_cv_1.RData'


### CV 2 output

cv1_points_sen_path <- 'model_outputs/mdg_points_cv_1.RData'
cv1_polys_sen_path <- 'model_outputs/mdg_polygon_cv_1.RData'
cv1_both_sen_path <- 'model_outputs/mdg_join_cv_1.RData'




## MDG


### Full data object

# full_data_mdg_path <- 'model_outputs/mdg_full_data.RData'


### Cross validation object

data_cv1_mdg_path <- 'model_outputs/mdg_cv_1.RData'
data_cv2_mdg_path <- 'model_outputs/mdg_cv_2.RData'


### CV 1 output

cv1_points_mdg_path <- 'model_outputs/mdg_points_cv_1.RData'
cv1_polys_mdg_path <- 'model_outputs/mdg_polygon_cv_1.RData'
cv1_both_mdg_path <- 'model_outputs/mdg_join_cv_1.RData'


### CV 2 output

cv1_points_mdg_path <- 'model_outputs/mdg_points_cv_1.RData'
cv1_polys_mdg_path <- 'model_outputs/mdg_polygon_cv_1.RData'
cv1_both_mdg_path <- 'model_outputs/mdg_join_cv_1.RData'


# figure 1.cross validation. %% Do fig 1 and 2, random and spatial cv. IDN on top, MDG and SEN below in each.

# Fig1 
# Load data
data_cv1_idn <- get(load(data_cv1_idn_path))
data_cv1_sen <- get(load(data_cv1_sen_path))
data_cv1_mdg <- get(load(data_cv1_mdg_path))



rm(data_cv1_idn)
rm(data_cv1_sen)
rm(data_cv1_mdg)
gc()

# Fig2
data_cv2_idn <- get(load(data_cv2_idn_path))
data_cv2_sen <- get(load(data_cv2_sen_path))
data_cv2_mdg <- get(load(data_cv2_mdg_path))


rm(data_cv2_idn)
rm(data_cv2_sen)
rm(data_cv2_mdg)
gc()



# figure 3 data and predicted incidence maps. Indonesia only. Data, Rand, Spatial for best model? Joint model?

full_data_idn <- get(load(full_data_idn_path))
cv1_both_idn <- get(load(cv1_both_idn_path))
cv2_both_idn <- get(load(cv2_both_idn_path))






rm(full_data_idn)
rm(cv1_both_idn)
rm(cv2_both_idn)
gc()


# figure 4 data and predicted incidence maps. Senegal only. Data, Rand, Spatial for best model? Joint model?

full_data_sen <- get(load(full_data_sen_path))
cv1_both_sen <- get(load(cv1_both_sen_path))
cv2_both_sen<- get(load(cv2_both_sen_path))


rm(full_data_sen)
rm(cv1_both_sen)
rm(cv2_both_sen)
gc()

# figure 5, Spat and random cv. PR vs Poly columns, countries as rows, model as colour?


cv1_points_idn <- get(load(cv1_points_idn_path))
cv1_polys_idn <- get(load(cv1_polys_idn_path))
cv1_both_idn <- get(load(cv1_both_idn_path))

cv1_points_sen <- get(load(cv1_points_sen_path))
cv1_polys_sen <- get(load(cv1_polys_sen_path))
cv1_both_sen <- get(load(cv1_both_sen_path))

cv1_points_mdg <- get(load(cv1_points_mdg_path))
cv1_polys_mdg <- get(load(cv1_polys_mdg_path))
cv1_both_mdg <- get(load(cv1_both_mdg_path))





rm(cv1_points_idn)
rm(cv1_polys_idn)
rm(cv1_both_idn)

rm(cv1_points_sen)
rm(cv1_polys_sen)
rm(cv1_both_sen)

rm(cv1_points_mdg)
rm(cv1_polys_mdg)
rm(cv1_both_mdg)
gc()

# Fig 4 and 6


cv2_points_idn <- get(load(cv2_points_idn_path))
cv2_polys_idn <- get(load(cv2_polys_idn_path))
cv2_both_idn <- get(load(cv2_both_idn_path))

cv2_points_sen <- get(load(cv2_points_sen_path))
cv2_polys_sen <- get(load(cv2_polys_sen_path))
cv2_both_sen <- get(load(cv2_both_sen_path))

cv2_points_mdg <- get(load(cv2_points_mdg_path))
cv2_polys_mdg <- get(load(cv2_polys_mdg_path))
cv2_both_mdg <- get(load(cv2_both_mdg_path))





rm(cv2_points_idn)
rm(cv2_polys_idn)
rm(cv2_both_idn)

rm(cv2_points_sen)
rm(cv2_polys_sen)
rm(cv2_both_sen)

rm(cv2_points_mdg)
rm(cv2_polys_mdg)
rm(cv2_both_mdg)
gc()

# Useful summary tables