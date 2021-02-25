
# figure 1.cross validation. %% Do fig 1 and 2, random and spatial cv. IDN on top, MDG and SEN below in each.
# figure 3 and 4. data and predicted incidence maps. Indonesia and Senegal only. Fig 3 ind, fig 4 sen Data, Rand, Spatial for best model? Joint model?
# figure 5, 6. Spat and random cv. PR vs Poly columns, countries as rows, model as colour?


# Libs

## dataframe packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(xtable)

### CV 3 output

idn3dfr <- read.csv('model_outputs/idn3df.csv')

sen3dfr <- read.csv('model_outputs/sen3df.csv')

mdg3dfr <- read.csv('model_outputs/mdg3df.csv')




idn3df %>% head

idn3df <- 
  idn3dfr %>% 
    group_by(model) %>% 
    summarise(RMSE = sqrt(mean((pred_api - response) ^ 2)),
              MAE = mean(abs(pred_api - response)),
              AE = mean(pred_api - response)) %>% 
    mutate(country = 'IDN')
    
    



idn3dfr %>% 
  group_by(model) %>% 
  filter(response > 30) %>% 
  summarise(RMSE = sqrt(mean((pred_api - response) ^ 2)),
            MAE = mean(abs(pred_api - response)),
            AE = mean(pred_api - response)) %>% 
  mutate(country = 'IDN')





    
mdg3df <- 
  mdg3dfr %>% 
    group_by(model) %>% 
    summarise(RMSE = sqrt(mean((pred_api - response) ^ 2)),
              MAE = mean(abs(pred_api - response)),
              AE = mean(pred_api - response)) %>% 
    mutate(country = 'MDG')


sen3df <-
  sen3dfr %>% 
    group_by(model) %>% 
    summarise(RMSE = sqrt(mean((pred_api - response) ^ 2)),
              MAE = mean(abs(pred_api - response)),
              AE = mean(pred_api - response)) %>% 
    mutate(country = 'SEN')
          
    

sen3dfr %>% 
  filter(model != 'prgp') %>% 
  mutate(error = pred_api - response) %>%
  wilcox.test(error ~ model, data = .)


sen3dfr %>% 
  filter(model != 'both') %>% 
  mutate(error = pred_api - response) %>%
  wilcox.test(error ~ model, data = .)


sen3dfr %>% 
  filter(model == 'both') %>% 
  mutate(error = pred_api - response) %>%
  wilcox.test(.$error, data = .)

sen3dfr %>% 
  filter(model == 'polys') %>% 
  mutate(error = pred_api - response) %>%
  wilcox.test(.$error, data = .)





idn3dfr %>% 
  filter(model != 'prgp') %>% 
  mutate(error = pred_api - response) %>%
  wilcox.test(error ~ model, data = .)

idn3dfr %>% 
  filter(model == 'both') %>% 
  mutate(error = pred_api - response) %>%
  wilcox.test(.$error, data = .)

idn3dfr %>% 
  filter(model == 'polys') %>% 
  mutate(error = pred_api - response) %>%
  wilcox.test(.$error, data = .)


idn3dfr %>% 
  sample_frac(1) %>% 
  mutate(error = pred_api - response) %>%
  mutate(model = recode(model, 
                    both = 'Joint', 
                    prgp = 'PrevGP', 
                    polys = 'Base')) %>% 
  ggplot(aes(x = response, y = error, colour = model)) + 
    geom_point(alpha = 0.6, size = 2) + 
    geom_abline(slope = 0, intercept = 0, alpha = 0.4) +
    xlim(00, 100) +
    ylim(-100, 50) +
    ylab('Model error') + 
    xlab('Observed API')
ggsave('figs/idn_model_error_vs_observed.pdf')


idn3dfr %>% 
  sample_frac(1) %>% 
  mutate(error = pred_api - response) %>%
  mutate(model = recode(model, 
                    both = 'Joint', 
                    prgp = 'PrevGP', 
                    polys = 'Base')) %>% 
  ggplot(aes(x = response, y = error, colour = model)) + 
    geom_point(alpha = 0.6, size = 2) + 
    geom_abline(slope = 0, intercept = 0, alpha = 0.4) +
    xlim(00, 10) +
    ylim(-10, 40) +
    ylab('Model error') + 
    xlab('Observed API')
ggsave('figs/idn_model_error_vs_observed_small.pdf')




sen3dfr %>% 
  sample_frac(1) %>% 
  mutate(model = recode(model, 
                    both = 'Joint', 
                    prgp = 'PrevGP', 
                    polys = 'Base')) %>% 
  mutate(error = pred_api - response) %>%
  ggplot(aes(x = response, y = error, colour = model)) + 
    geom_point(alpha = 0.6, size = 2) + 
    geom_abline(slope = 0, intercept = 0, alpha = 0.4) +
    ylab('Model error') + 
    xlab('Observed API')
ggsave('figs/sen_model_error_vs_observed.pdf')



mdg3dfr %>% 
  sample_frac(1) %>% 
  mutate(model = recode(model, 
                    both = 'Joint', 
                    prgp = 'PrevGP', 
                    polys = 'Base')) %>% 
  mutate(error = pred_api - response) %>%
  ggplot(aes(x = response, y = error, colour = model)) + 
    geom_point(alpha = 0.6, size = 2) + 
    geom_abline(slope = 0, intercept = 0, alpha = 0.4) +
    ylab('Model error') + 
    xlab('Observed API')
ggsave('figs/mdg_model_error_vs_observed.pdf')



idn3dfr %>% 
  filter(model != 'prgp') %>% 
  mutate(error = pred_api - response) %>%
  dplyr::select(model, error, response) %>% 
  group_by(model) %>% 
  mutate(id = sequence(n())) %>% 
  pivot_wider(id_cols = id, 
              names_from = model, 
              values_from = error) %>% 
  mutate(error_diff = polys - both) %>% 
  cbind(idn3dfr[seq(nrow(.)), ]) %>% 
  ggplot(aes(x = response, y = error_diff)) + 
    geom_point(alpha = 0.6, size = 2) + 
    geom_abline(slope = 0, intercept = 0, alpha = 0.4) +
    xlim(0, 100) 



rbind(
  idn3dfr %>% mutate(country = 'IDN'),
  sen3dfr %>% mutate(country = 'SEN'),
  mdg3dfr %>% mutate(country = 'MDG')
  ) %>% 
  mutate(model = recode(model, 
                      both = 'Joint', 
                      prgp = 'PrevGP', 
                      polys = 'Base')) %>% 
  mutate(error = pred_api - response) %>%
  filter(!(country == 'IDN' & error < -40)) %>% 
  ggplot(aes(x = error, colour = model)) + 
    geom_density(size = 1) +
    facet_wrap(~ country, scales = 'free') +
    geom_rug() +
    ggtitle('Density plots of spatial cv model error')

ggsave('figs/spatial_cv_model_error.pdf',
       height = 4, width = 7)







rbind(
  idn1dfr %>% mutate(country = 'IDN'),
  sen1dfr %>% mutate(country = 'SEN'),
  mdg1dfr %>% mutate(country = 'MDG')
  ) %>% 
  mutate(model = recode(model, 
                        both = 'Joint', 
                        prgp = 'PrevGP', 
                        polys = 'Base')) %>% 
  mutate(error = pred_api - response) %>%
  filter(!(country == 'IDN' & abs(error) > 1)) %>% 
  ggplot(aes(x = error, colour = model)) + 
    geom_density(size = 1) +
    facet_wrap(~ country, scales = 'free') 

ggsave('figs/random_cv_model_error.pdf',
       height = 4, width = 7)



sen3dfr %>% 
  mutate(error = pred_api - response) %>%
  ggplot(aes(x = error, colour = model)) + 
    geom_density() 




d <- rbind(sen3df, idn3df, mdg3df)


ggplot(d, aes(MAE, AE, colour = model, shape = country)) + 
  geom_point()




ggplot(mdg3dfr, aes(pred_api, response, colour = model)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  #scale_y_sqrt() + 
  #scale_x_sqrt() +
  facet_wrap(~model) + 
  geom_smooth(method = 'lm', aes(colour = model), se = FALSE) +
  labs(y = 'Observed API', x = 'Predicted API') +
  guides(colour = FALSE) 



ggplot(idn3dfr, aes(pred_api, response, colour = model)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() + 
  scale_y_log10() +
  facet_wrap(~model) + 
  geom_smooth(method = 'lm', aes(colour = model), se = FALSE) +
  geom_smooth(aes(colour = model), se = FALSE) +
  labs(y = 'Observed API', x = 'Predicted API') +
  guides(colour = FALSE) 


ggplot(sen3dfr, aes(pred_api, response, colour = model)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() + 
  scale_y_log10() +
  facet_wrap(~model) + 
  geom_smooth(method = 'lm', aes(colour = model), se = FALSE) +
  geom_smooth(aes(colour = model), se = FALSE) +
  labs(y = 'Observed API', x = 'Predicted API') +
  guides(colour = FALSE)


ggplot(mdg3dfr, aes(pred_api, response, colour = model)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() + 
  scale_y_log10() +
  facet_wrap(~model) + 
  geom_smooth(method = 'lm', aes(colour = model), se = FALSE) +
  geom_smooth(aes(colour = model), se = FALSE) +
  labs(y = 'Observed API', x = 'Predicted API') +
  guides(colour = FALSE)








idn1dfr <- read.csv('model_outputs/idn1df.csv')

sen1dfr <- read.csv('model_outputs/sen1df.csv')

mdg1dfr <- read.csv('model_outputs/mdg1df.csv')




idn1df %>% head

idn1df <- 
  idn1dfr %>% 
    group_by(model) %>% 
    summarise(RMSE = sqrt(mean((pred_api - response) ^ 2)),
              MAE = mean(abs(pred_api - response)),
              AE = mean(pred_api - response)) %>% 
    mutate(country = 'IDN')
    
    
    
mdg1df <- 
  mdg1dfr %>% 
    group_by(model) %>% 
    summarise(RMSE = sqrt(mean((pred_api - response) ^ 2)),
              MAE = mean(abs(pred_api - response)),
              AE = mean(pred_api - response)) %>% 
    mutate(country = 'MDG')


sen1df <-
  sen1dfr %>% 
    group_by(model) %>% 
    summarise(RMSE = sqrt(mean((pred_api - response) ^ 2)),
              MAE = mean(abs(pred_api - response)),
              AE = mean(pred_api - response)) %>% 
    mutate(country = 'SEN')
          
    


d1 <- rbind(sen1df, idn1df, mdg1df)


d1 %>% 
  select(country, model, RMSE, MAE, AE) %>% 
  mutate(model = recode(model, 
                  both = 'Joint', 
                  prgp = 'PrevGP', 
                  polys = 'Base')) %>% 
  xtable(type = 'latex') %>% 
  print(include.rownames = FALSE)


d %>% 
  select(country, model, RMSE, MAE, AE) %>% 
  mutate(model = recode(model, 
                  both = 'Joint', 
                  prgp = 'PrevGP', 
                  polys = 'Base')) %>% 
  xtable(type = 'latex') %>% 
  print(include.rownames = FALSE)








load('model_outputs/idn_joint_cv_3.RData', verbose = T)
load('model_outputs/idn_polygon_cv_3.RData')


cv3_output2$models[[1]]$model$sd_report %>% 
  summary %>% `[`(3:10, 1)



sd_joint <- 
  lapply(cv3_output3$models, 
         function(x) x$model$sd_report %>% 
           summary %>% `[`(3:10, 1)
  ) %>% do.call(c, .)

sd_base <- 
  lapply(cv3_output2$models, 
         function(x) x$model$sd_report %>% 
           summary %>% `[`(3:10, 1)
  ) %>% do.call(c, .)

plot(sd_joint, sd_base)
abline(0, 1)

data.frame(joint = sd_joint, base = sd_base) %>%
  mutate(sign_joint = sign(joint), 
         sign_base = sign(base)) %>% 
  filter(sign_joint == sign_base) %>% 
  ggplot(aes(abs(joint), abs(base))) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1)
  


data.frame(joint = sd_joint, base = sd_base) %>%
  mutate(sign_joint = sign(joint), 
         sign_base = sign(base)) %>% 
  mutate(same_sign = sign_joint == sign_base) %>% 
  ggplot(aes(joint, base, colour = same_sign)) + 
    geom_abline(intercept = 0, slope = 1) + 
    geom_hline(yintercept = 0, colour = 'darkgrey') + 
    geom_vline(xintercept = 0, colour = 'darkgrey') +
    geom_point()



data.frame(joint = sd_joint, base = sd_base) %>%
  mutate(sign_joint = sign(joint), 
         sign_base = sign(base)) %>% 
  mutate(same_sign = sign_joint == sign_base) %>% 
  ggplot(aes(abs(joint), abs(base), colour = same_sign)) + 
    geom_abline(intercept = 0, slope = 1) + 
    geom_point()

ggsave('idn_cov_regression_pars.pdf')




load('model_outputs/sen_joint_cv_3.RData', verbose = T)
load('model_outputs/sen_polygon_cv_3.RData')


cv3_output2$models[[1]]$model$sd_report %>% 
  summary %>% `[`(3:10, 1)



sd_joint <- 
  lapply(cv3_output3$models, 
         function(x) x$model$sd_report %>% 
           summary %>% `[`(3:10, 1)
  ) %>% do.call(c, .)

sd_base <- 
  lapply(cv3_output2$models, 
         function(x) x$model$sd_report %>% 
           summary %>% `[`(3:10, 1)
  ) %>% do.call(c, .)

plot(sd_joint, sd_base)
abline(0, 1)

data.frame(joint = sd_joint, base = sd_base) %>%
  mutate(sign_joint = sign(joint), 
         sign_base = sign(base)) %>% 
  filter(sign_joint == sign_base) %>% 
  ggplot(aes(abs(joint), abs(base))) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1)
  


data.frame(joint = sd_joint, base = sd_base) %>%
  mutate(sign_joint = sign(joint), 
         sign_base = sign(base)) %>% 
  mutate(same_sign = sign_joint == sign_base) %>% 
  ggplot(aes(joint, base, colour = same_sign)) + 
    geom_abline(intercept = 0, slope = 1) + 
    geom_hline(yintercept = 0, colour = 'darkgrey') + 
    geom_vline(xintercept = 0, colour = 'darkgrey') +
    geom_point()



data.frame(joint = sd_joint, base = sd_base) %>%
  mutate(sign_joint = sign(joint), 
         sign_base = sign(base)) %>% 
  mutate(same_sign = sign_joint == sign_base) %>% 
  ggplot(aes(abs(joint), abs(base), colour = same_sign)) + 
    geom_abline(intercept = 0, slope = 1) + 
    geom_point()

ggsave('sen_cov_regression_pars.pdf')
