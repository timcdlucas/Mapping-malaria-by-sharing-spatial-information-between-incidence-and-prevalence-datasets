

library(raster)
library(dplyr)
library(sp)

load('~/timz/timothy/point_polygon_joint_comparison/model_outputs/sen_full_data.RData')


# Detect correct user

Z <- function(path){
  if(Sys.info()["user"] != 'anita'){
    fullpath <- paste0('~/timz/', path)
  } else {
    fullpath <- paste0( '~/Z/', path)
  }
  if(Sys.info()["sysname"] == 'Windows'){
    fullpath <- paste0( 'Z://', path)
  }
  return(fullpath)
}

covs <- data_sen$cov_rasters
crs(covs) <- projection(covs)

writeRaster(covs, 
            filename = 'data/senegal_covariate',
            bylayer = TRUE,                 
            format="GTiff", overwrite = TRUE, 
            options = c('COMPRESS' = 'LZW'))

pop <- data_sen$pop_raster
crs(pop) <- projection(pop)

writeRaster(data_sen$pop_raster, 
            filename = 'data/senegal_population.tif',
            format="GTiff", overwrite = TRUE, 
            options = c('COMPRESS' = 'LZW'))




shapefile(data_sen$shapefiles, filename = 'data/senegal-shapefiles.shp')

PR_path <- Z('GBD2017/Processing/Stages/04b_PR_DB_Import_Export/Verified_Outputs/2018_02_15/pfpr_dump.csv')
API_path <- Z('GBD2017/Processing/Stages/04c_API_Data_Export/Checkpoint_Outputs/2018-10-17_Senegal/subnational.csv')

pr <- read.csv(PR_path)

pr_sim <- pr %>% 
            filter(country == 'Senegal', year_start %in% c(2007:2011)) %>% 
            select(country, year_start, lower_age, upper_age, examined, positive, latitude, longitude)


pr_sim$examined <- rpois(nrow(pr_sim), 300)
pr_sim$positive <- round(pr_sim$examined * runif(nrow(pr_sim)))

write.csv(pr_sim, 'data/senegal_simulated_pr.csv')

api <- read.csv(API_path)
api_sim <- api %>% 
             filter(iso3 == 'SEN', year == 2009) %>% 
             select(iso3, country_name, admin_unit_level, admin_unit_name, shapefile_id, year, pop_at_risk_pf, api_mean_pf)

api_sim$api_mean_pf <- runif(nrow(api_sim), 1, 40)

write.csv(api_sim, 'data/senegal_simulated_incidence.csv')








