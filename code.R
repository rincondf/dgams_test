
library(mvgam)

all(levels(blh_dat$series) %in% 
      unique(blh_dat$series))

get_mvgam_priors(counts ~ 1,
                 data = blh_dat,
                 family = Gamma()
)


plot_mvgam_series(data = blh_dat, y = "counts", series = "all")

dat_train <- blh_dat[which(blh_dat$time <= 360), ]
dat_test <- blh_dat[which(blh_dat$time > 360), ]


plot_mvgam_series(data = dat_train, 
                  newdata = dat_test,
                  y = "counts")


mod <- mvgam(counts ~
               
               # Hierarchical intercepts capture variation in average
               s(series, bs = "re") +
               # relative abundances
               s(DDs, k = 6, bs = 'fs') +
               s(year, k = 4, bs = 'fs'),
             
             # Condition on the training data
             data = dat_train,
             
             # Automatically compute forecasts for the test data
             newdata = dat_test,
             
             tred_model = AR(3),
             
             # Beta observations with independent precisions
             family = Gamma(),
             
             # cmdstanr is highly recommended over rstan as 
             # it is much more up-to-date with the latest 
             # features in the Stan universe
             backend = 'cmdstanr')
summary(mod)
plot(mod, type = 're')

gratia::draw(mod)
conditional_effects(mod)
conditional_effects(mod,
                    type = 'link')


marginaleffects::avg_predictions(mod, variable = 'series')
pp_check(mod,
         type = "ecdf_overlay_grouped",
         group = "series",
         ndraws = 50)



hcs <- hindcast(mod)
class(hcs)
?mvgam::`mvgam_forecast-class`
methods(class = "mvgam_forecast")

layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot(hcs, series = 1)
plot(hcs, series = 2)
plot(hcs, series = 3)
plot(hcs, series = 4)
plot(hcs, series = 5)
plot(hcs, series = 6)
plot(hcs, series = 7)


layout(1)



fcs <- forecast(mod)

plot(fcs, series = 1)
plot(fcs, series = 2)
plot(fcs, series = 3)
plot(fcs, series = 4)
plot(fcs, series = 5)
plot(fcs, series = 6)
plot(fcs, series = 7)




plot(mod, type = 'forecast', series = 1)
plot(mod, type = 'forecast', series = 2)
plot(mod, type = 'forecast', series = 3)
plot(mod, type = 'forecast', series = 4)
plot(mod, type = 'forecast', series = 5)
plot(mod, type = 'forecast', series = 6)
plot(mod, type = 'forecast', series = 7)

plot_predictions(mod, by = "time", points = 0.5)
plot_predictions(mod, condition = c("time", "series", "series"), points = 0.5)


year_seq <- seq(min(dat_train$year), 
                max(dat_train$year), 
                length.out = 20)

plot_predictions(
  model = mod, 
  
  # Predict over the sequence of times in 'year_seq'
  newdata = datagrid(year = year_seq,
                     series = unique),
  
  # Use by = 'year' to ensure predictions are averaged
  # over the sequence of times in year_seq
  by = 'year',
  
  # Compute predictions on the link scale
  type = 'link'
)


plot_comparisons(
  model = mod, 
  
  # Predict over the sequence of times in 'year_seq'
  newdata = datagrid(year = year_seq,
                     series = unique),
  
  # At each time in 'year_seq', calculate contrasts between
  # each pair of time series
  variables = list(series = 'pairwise'),
  
  # Plot the difference trends using year as the x-axis
  by = 'year',
  
  # Compute predictions on the link scale
  type = 'link'
)


plot_predictions(mod,
                 condition = c('year', 'series'),
                 type = 'link',
                 conf_level = 0.5)
