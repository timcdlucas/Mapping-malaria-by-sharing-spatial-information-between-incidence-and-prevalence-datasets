
setwd('~/timz/timothy/point_polygon_joint_comparison')

# Libs

library(ggplot2)
library(cowplot)
source('plotting_functions.R')

# Paths

## IDN

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








# Load data

## IDN

### Cross validation object

data_cv1_idn <- get(load(data_cv1_idn_path))
data_cv2_idn <- get(load(data_cv2_idn_path))


### CV 1 output

cv1_points_idn <- get(load(cv1_points_idn_path))
cv1_polys_idn <- get(load(cv1_polys_idn_path))
cv1_both_idn <- get(load(cv1_both_idn_path))


### CV 2 output

cv1_points_idn <- get(load(cv1_points_idn_path))
cv1_polys_idn <- get(load(cv1_polys_idn_path))
cv1_both_idn <- get(load(cv1_both_idn_path))


## SEN

### Cross validation object

data_cv1_sen <- get(load(data_cv1_sen_path))
data_cv2_sen <- get(load(data_cv2_sen_path))


### CV 1 output

cv1_points_sen <- get(load(cv1_points_sen_path))
cv1_polys_sen <- get(load(cv1_polys_sen_path))
cv1_both_sen <- get(load(cv1_both_sen_path))


### CV 2 output

cv1_points_sen <- get(load(cv1_points_sen_path))
cv1_polys_sen <- get(load(cv1_polys_sen_path))
cv1_both_sen <- get(load(cv1_both_sen_path))




## MDG

### Cross validation object

data_cv1_mdg <- get(load(data_cv1_mdg_path))
data_cv2_mdg <- get(load(data_cv2_mdg_path))


### CV 1 output

cv1_points_mdg <- get(load(cv1_points_mdg_path))
cv1_polys_mdg <- get(load(cv1_polys_mdg_path))
cv1_both_mdg <- get(load(cv1_both_mdg_path))


### CV 2 output

cv1_points_mdg <- get(load(cv1_points_mdg_path))
cv1_polys_mdg <- get(load(cv1_polys_mdg_path))
cv1_both_mdg <- get(load(cv1_both_mdg_path))



# figure 1.cross validation. %% Do fig 1 and 2, random and spatial cv. IDN on top, MDG and SEN below in each.






# figure 3 and 4. data and predicted incidence maps. Indonesia and Senegal only. Fig 3 ind, fig 4 rand. Data, Rand, Spatial for best model? Joint model?





# figure 5, 6. Spat and random cv. PR vs Poly columns, countries as rows, model as colour?



# Useful summary tables