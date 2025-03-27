
library(mvgam)

# Checking all series are complete

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
             
             trend_model = AR(3),
             
             family = Gamma(),
             
             # cmdstanr is highly recommended over rstan as 
             # it is much more up-to-date with the latest 
             # features in the Stan universe
             noncentred = TRUE,
             backend = 'cmdstanr')
summary(mod)
plot(mod, type = 're')

gratia::draw(mod)
conditional_effects(mod, points = 0.5)
conditional_effects(mod,
                    type = 'link')


marginaleffects::avg_predictions(mod, variable = 'series')
pp_check(mod,
         type = "ecdf_overlay_grouped",
         group = "series",
         ndraws = 50)



hcs <- hindcast(mod)


plot(hcs, series = 1)
plot(hcs, series = 2)
plot(hcs, series = 3)
plot(hcs, series = 4)
plot(hcs, series = 5)
plot(hcs, series = 6)
plot(hcs, series = 7)



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


DD_seq <- seq(min(dat_train$DDs), 
                max(dat_train$DDs), 
                length.out = 20)



plot_predictions(
  model = mod, 
  
  # Predict over the sequence of times in 'year_seq'
  newdata = datagrid(DDs = DD_seq,
                     series = unique),
  
  # Use by = 'year' to ensure predictions are averaged
  # over the sequence of times in year_seq
  by = 'DDs',
  
  # Compute predictions on the link scale
  type = 'link'
)


plot_comparisons(
  model = mod, 
  
  # Predict over the sequence of times in 'year_seq'
  newdata = datagrid(DDs = DD_seq,
                     series = unique),
  
  # At each time in 'year_seq', calculate contrasts between
  # each pair of time series
  variables = list(series = 'pairwise'),
  
  # Plot the difference trends using year as the x-axis
  by = 'DDs',
  
  # Compute predictions on the link scale
  type = 'link'
)


plot_predictions(mod,
                 condition = c('DDs', 'series'),
                 type = 'link',
                 conf_level = 0.5)








mod2 <- mvgam(counts ~ 0,
              
              trend_formula = ~
                s(DDs, k = 6, bs = 'fs') +
                s(year, k = 4, bs = 'fs') - 1,
              
              # Condition on the training data
              data = dat_train,
              
              # Automatically compute forecasts for the test data
              newdata = dat_test,
              
              trend_model = AR(3),
              
              family = Gamma(),
              
              # cmdstanr is highly recommended over rstan as 
              # it is much more up-to-date with the latest 
              # features in the Stan universe
              noncentred = TRUE,
              backend = 'cmdstanr')
summary(mod2)

gratia::draw(mod2)
conditional_effects(mod2, points = 0.5)
conditional_effects(mod2,
                    type = 'link')


marginaleffects::avg_predictions(mod2, variable = 'series')
pp_check(mod2,
         type = "ecdf_overlay_grouped",
         group = "series",
         ndraws = 50)



hcs2 <- hindcast(mod2)

plot(hcs2, series = 1)
plot(hcs2, series = 2)
plot(hcs2, series = 3)
plot(hcs2, series = 4)
plot(hcs2, series = 5)
plot(hcs2, series = 6)
plot(hcs2, series = 7)


plot(mod2, type = 'forecast', series = 1)
plot(mod2, type = 'forecast', series = 2)
plot(mod2, type = 'forecast', series = 3)
plot(mod2, type = 'forecast', series = 4)
plot(mod2, type = 'forecast', series = 5)
plot(mod2, type = 'forecast', series = 6)
plot(mod2, type = 'forecast', series = 7)

plot_predictions(mod2, by = "time", points = 0.5)
plot_predictions(mod2, condition = c("time", "series", "series"), points = 0.5)



plot_predictions(
  model = mod2, 
  
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
  model = mod2, 
  
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


plot_predictions(mod2,
                 condition = c('year', 'series'),
                 type = 'link',
                 conf_level = 0.5)




plot_predictions(
  model = mod2, 
  
  # Predict over the sequence of times in 'year_seq'
  newdata = datagrid(DDs = DD_seq,
                     series = unique),
  
  # Use by = 'year' to ensure predictions are averaged
  # over the sequence of times in year_seq
  by = 'DDs',
  
  # Compute predictions on the link scale
  type = 'link'
)


plot_comparisons(
  model = mod2, 
  
  # Predict over the sequence of times in 'year_seq'
  newdata = datagrid(DDs = DD_seq,
                     series = unique),
  
  # At each time in 'year_seq', calculate contrasts between
  # each pair of time series
  variables = list(series = 'pairwise'),
  
  # Plot the difference trends using year as the x-axis
  by = 'DDs',
  
  # Compute predictions on the link scale
  type = 'link'
)


plot_predictions(mod2,
                 condition = c('DDs', 'series'),
                 type = 'link',
                 conf_level = 0.5)













mod3 <- mvgam(counts ~
                
                # Hierarchical intercepts capture variation in average
                s(series, bs = "re") +
                # relative abundances
                s(DDs, k = 6, bs = 'fs') +
                s(year, k = 4, bs = 'fs'),
              
              # Condition on the training data
              data = dat_train,
              
              # Automatically compute forecasts for the test data
              newdata = dat_test,
              
              trend_model = AR(1, cor = TRUE),
              
              
              family = Gamma(),
              
              # cmdstanr is highly recommended over rstan as 
              # it is much more up-to-date with the latest 
              # features in the Stan universe
              noncentred = TRUE,
              backend = 'cmdstanr')
summary(mod3)
plot(mod3, type = 're')

gratia::draw(mod3)
conditional_effects(mod3, points = 0.5)
conditional_effects(mod3,
                    type = 'link')


marginaleffects::avg_predictions(mod3, variable = 'series')
pp_check(mod3,
         type = "ecdf_overlay_grouped",
         group = "series",
         ndraws = 50)



hcs3<- hindcast(mod3)


plot(hcs3, series = 1)
plot(hcs3, series = 2)
plot(hcs3, series = 3)
plot(hcs3, series = 4)
plot(hcs3, series = 5)
plot(hcs3, series = 6)
plot(hcs3, series = 7)



plot(mod3, type = 'forecast', series = 1, ylim = c(0, 120))
plot(mod3, type = 'forecast', series = 2, ylim = c(0, 120))
plot(mod3, type = 'forecast', series = 3, ylim = c(0, 120))
plot(mod3, type = 'forecast', series = 4, ylim = c(0, 120))
plot(mod3, type = 'forecast', series = 5, ylim = c(0, 120))
plot(mod3, type = 'forecast', series = 6, ylim = c(0, 120))
plot(mod3, type = 'forecast', series = 7, ylim = c(0, 120))

plot_predictions(mod3, by = "time", points = 0.5)
plot_predictions(mod3, condition = c("time", "series", "series"), points = 0.5)


year_seq <- seq(min(dat_train$year), 
                max(dat_train$year), 
                length.out = 20)

plot_predictions(
  model = mod3, 
  
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
  model = mod3, 
  
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


plot_predictions(mod3,
                 condition = c('year', 'series'),
                 type = 'link',
                 conf_level = 0.5)


DD_seq <- seq(min(dat_train$DDs), 
              max(dat_train$DDs), 
              length.out = 20)



plot_predictions(
  model = mod3, 
  
  # Predict over the sequence of times in 'year_seq'
  newdata = datagrid(DDs = DD_seq,
                     series = unique),
  
  # Use by = 'year' to ensure predictions are averaged
  # over the sequence of times in year_seq
  by = 'DDs',
  
  # Compute predictions on the link scale
  type = 'link'
)


plot_comparisons(
  model = mod3, 
  
  # Predict over the sequence of times in 'year_seq'
  newdata = datagrid(DDs = DD_seq,
                     series = unique),
  
  # At each time in 'year_seq', calculate contrasts between
  # each pair of time series
  variables = list(series = 'pairwise'),
  
  # Plot the difference trends using year as the x-axis
  by = 'DDs',
  
  # Compute predictions on the link scale
  type = 'link'
)


plot_predictions(mod3,
                 condition = c('DDs', 'series'),
                 type = 'link',
                 conf_level = 0.5)

loo_compare(mod,
            mod3)










mod4 <- mvgam(counts ~
                
                s(DDs, bs = "fs", k = 8) +
                s(year, by = series, bs = "cr", k = 5),
                
                
              
              # Condition on the training data
              data = dat_train,
              
              # Automatically compute forecasts for the test data
              newdata = dat_test,
              
              trend_model = "None",
              
              
              family = Gamma(),
              
              # cmdstanr is highly recommended over rstan as 
              # it is much more up-to-date with the latest 
              # features in the Stan universe
              backend = 'cmdstanr')
summary(mod4)

gratia::draw(mod4)
conditional_effects(mod4, points = 0.5)
conditional_effects(mod4,
                    type = 'link')


marginaleffects::avg_predictions(mod4, variable = 'series')
pp_check(mod4,
         type = "ecdf_overlay_grouped",
         group = "series",
         ndraws = 50)



hcs4 <- hindcast(mod4)

plot(hcs4, series = 1)
plot(hcs4, series = 2)
plot(hcs4, series = 3)
plot(hcs4, series = 4)
plot(hcs4, series = 5)
plot(hcs4, series = 6)
plot(hcs4, series = 7)


plot(mod4, type = 'forecast', series = 1)
plot(mod4, type = 'forecast', series = 2)
plot(mod4, type = 'forecast', series = 3)
plot(mod4, type = 'forecast', series = 4)
plot(mod4, type = 'forecast', series = 5)
plot(mod4, type = 'forecast', series = 6)
plot(mod4, type = 'forecast', series = 7)

plot_predictions(mod4, by = "time", points = 0.5)
plot_predictions(mod4, condition = c("time", "series", "series"), points = 0.5)



plot_predictions(
  model = mod4, 
  
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
  model = mod4, 
  
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


plot_predictions(mod4,
                 condition = c('year', 'series'),
                 type = 'link',
                 conf_level = 0.5)




plot_predictions(
  model = mod4, 
  
  # Predict over the sequence of times in 'year_seq'
  newdata = datagrid(DDs = DD_seq,
                     series = unique),
  
  # Use by = 'year' to ensure predictions are averaged
  # over the sequence of times in year_seq
  by = 'DDs',
  
  # Compute predictions on the link scale
  type = 'link'
)



plot_predictions(mod4,
                 condition = c('DDs', 'series'),
                 type = 'link',
                 conf_level = 0.5)





















mod5 <- mvgam(counts ~ s(series, bs = "re"),
              
              # Hierarchical intercepts capture variation in average
              
              
              trend_formula = ~
                # relative abundances
                s(DDs, k = 6, bs = 'fs') +
                s(year, k = 4, bs = 'fs') - 1,
              # Condition on the training data
              data = dat_train,
              
              # Automatically compute forecasts for the test data
              newdata = dat_test,
              
              trend_model = AR(1, cor = TRUE),
              
              
              family = Gamma(),
              
              # cmdstanr is highly recommended over rstan as 
              # it is much more up-to-date with the latest 
              # features in the Stan universe
              noncentred = TRUE,
              backend = 'cmdstanr')
summary(mod5)

gratia::draw(mod5)
conditional_effects(mod5, points = 0.5)
conditional_effects(mod5,
                    type = 'link')


marginaleffects::avg_predictions(mod4, variable = 'series')
pp_check(mod4,
         type = "ecdf_overlay_grouped",
         group = "series",
         ndraws = 50)



hcs4 <- hindcast(mod4)

plot(hcs4, series = 1)
plot(hcs4, series = 2)
plot(hcs4, series = 3)
plot(hcs4, series = 4)
plot(hcs4, series = 5)
plot(hcs4, series = 6)
plot(hcs4, series = 7)


plot(mod4, type = 'forecast', series = 1)
plot(mod4, type = 'forecast', series = 2)
plot(mod4, type = 'forecast', series = 3)
plot(mod4, type = 'forecast', series = 4)
plot(mod4, type = 'forecast', series = 5)
plot(mod4, type = 'forecast', series = 6)
plot(mod4, type = 'forecast', series = 7)

plot_predictions(mod4, by = "time", points = 0.5)
plot_predictions(mod4, condition = c("time", "series", "series"), points = 0.5)



plot_predictions(
  model = mod4, 
  
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
  model = mod4, 
  
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


plot_predictions(mod4,
                 condition = c('year', 'series'),
                 type = 'link',
                 conf_level = 0.5)




plot_predictions(
  model = mod4, 
  
  # Predict over the sequence of times in 'year_seq'
  newdata = datagrid(DDs = DD_seq,
                     series = unique),
  
  # Use by = 'year' to ensure predictions are averaged
  # over the sequence of times in year_seq
  by = 'DDs',
  
  # Compute predictions on the link scale
  type = 'link'
)



plot_predictions(mod4,
                 condition = c('DDs', 'series'),
                 type = 'link',
                 conf_level = 0.5)



















mod6 <- mvgam(counts ~ s(series, bs = "re"),
               
               # Hierarchical intercepts capture variation in average
               
              
              trend_formula = ~
                # relative abundances
                s(DDs, k = 6, bs = 'fs') +
                s(year, k = 4, bs = 'fs') - 1,
             
             # Condition on the training data
             data = dat_train,
             
             # Automatically compute forecasts for the test data
             newdata = dat_test,
             
             trend_model = AR(3),
             
             
             family = Gamma(),
             
             # cmdstanr is highly recommended over rstan as 
             # it is much more up-to-date with the latest 
             # features in the Stan universe
             noncentred = TRUE,
             backend = 'cmdstanr')
summary(mod6)


conditional_effects(mod6, points = 0.5)
conditional_effects(mod6,
                    type = 'link')


marginaleffects::avg_predictions(mod6, variable = 'series')
pp_check(mod6,
         type = "ecdf_overlay_grouped",
         group = "series",
         ndraws = 50)



hcs6 <- hindcast(mod6)


plot(hcs6, series = 1)
plot(hcs6, series = 2)
plot(hcs6, series = 3)
plot(hcs6, series = 4)
plot(hcs6, series = 5)
plot(hcs6, series = 6)
plot(hcs6, series = 7)



plot(mod6, type = 'forecast', series = 1)
plot(mod6, type = 'forecast', series = 2)
plot(mod6, type = 'forecast', series = 3)
plot(mod6, type = 'forecast', series = 4)
plot(mod6, type = 'forecast', series = 5)
plot(mod6, type = 'forecast', series = 6)
plot(mod6, type = 'forecast', series = 7)

plot_predictions(mod6, by = "time", points = 0.5)
plot_predictions(mod6, condition = c("time", "series", "series"), points = 0.5)


plot_predictions(
  model = mod6, 
  
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


plot_predictions(mod6,
                 condition = c('year', 'series'),
                 type = 'link',
                 conf_level = 0.5)


DD_seq <- seq(min(dat_train$DDs), 
              max(dat_train$DDs), 
              length.out = 20)



plot_predictions(
  model = mod6, 
  
  # Predict over the sequence of times in 'year_seq'
  newdata = datagrid(DDs = DD_seq,
                     series = unique),
  
  # Use by = 'year' to ensure predictions are averaged
  # over the sequence of times in year_seq
  by = 'DDs',
  
  # Compute predictions on the link scale
  type = 'link'
)





plot_predictions(mod6,
                 condition = c('DDs', 'series'),
                 type = 'link',
                 conf_level = 0.5)


loo_compare(mod,
            mod2,
            mod4,
            mod5,
            mod6)






fc_mod <- forecast(mod)
fc_mod2 <- forecast(mod2)
fc_mod3 <- forecast(mod3)
fc_mod4 <- forecast(mod4)
fc_mod5 <- forecast(mod5)
fc_mod6 <- forecast(mod6)




sc_mod <- score(fc_mod)
sc_mod2 <- score(fc_mod2)
sc_mod3 <- score(fc_mod3)
sc_mod4 <- score(fc_mod4)
sc_mod5 <- score(fc_mod5)
sc_mod6 <- score(fc_mod6)

par(mar = c(5, 5, 2, 2))
plot(sc_mod$all_series$eval_horizon, sc_mod$all_series$score, 
     type = "o", lwd = 2, xlab = "Prediction horizon", ylab = "Error score", 
     cex.lab = 1.8, cex.axis = 1.8, ylim = c(0, 150))
abline(h = 0, lty = 2)



par(mar = c(5, 5, 2, 2))
plot(sc_mod6$all_series$eval_horizon, sc_mod6$all_series$score, 
     type = "o", lwd = 2, xlab = "Prediction horizon", ylab = "Error score", 
     cex.lab = 1.8, cex.axis = 1.8, ylim = c(0, 150))
abline(h = 0, lty = 2)


plot(sc_mod2$all_series$eval_horizon, sc_mod2$all_series$score)
plot(sc_mod3$all_series$eval_horizon, sc_mod3$all_series$score)
plot(sc_mod4$all_series$eval_horizon, sc_mod4$all_series$score)
plot(sc_mod5$all_series$eval_horizon, sc_mod5$all_series$score)
plot(sc_mod6$all_series$eval_horizon, sc_mod6$all_series$score)


sum(sc_mod$all_series$score, na.rm = TRUE)
sum(sc_mod2$all_series$score, na.rm = TRUE)
sum(sc_mod3$all_series$score, na.rm = TRUE)
sum(sc_mod4$all_series$score, na.rm = TRUE)
sum(sc_mod5$all_series$score, na.rm = TRUE)
sum(sc_mod6$all_series$score, na.rm = TRUE)


ss <- "Basin City"
aa <- "fc_mod"
fc_mod$forecasts$paste(ss)

paste0(aa, "$", "forecasts", "$", ss)


plot_blh <- function(model, region) {
  
  
  lowBC <- rep(NA, 31)
  
  for(i in 1: 31) {
    lowBC[i] <- as.numeric(summary(model$forecasts[region][[1]][, i])[2])
  }
  
  hiBC <- rep(NA, 31)
  
  for(i in 1: 31) {
    hiBC[i] <- as.numeric(summary(model$forecasts[region][[1]][, i])[5])
  }
  
  
  maxBC <- rep(NA, 31)
  
  for(i in 1: 31) {
    maxBC[i] <- as.numeric(quantile(model$forecasts[region][[1]][, i], probs = 0.95))
  }
  
  minBC <- rep(NA, 31)
  
  for(i in 1: 31) {
    minBC[i] <- as.numeric(quantile(model$forecasts[region][[1]][, i], probs = 0.05))
  }
  
  
  par(mar = c(5, 5, 2, 2))
  plot(seq(361, 391), maxBC, xlim = c(360, 395), type = "l", lty = 3,
       xlab = "Time (weeks)", ylab = "Counts", cex.lab = 1.8)
  lines(seq(361, 391), lowBC, lty = 2)
  points(dat_test$time[which(dat_test$series == 
                              levels(dat_test$series)[region])], 
        dat_test$counts[which(dat_test$series == 
                                levels(dat_test$series)[region])],
        lwd = 2)
  lines(seq(361, 391), hiBC, lty = 2)
  lines(seq(361, 391), minBC, lty = 3)
  lines(seq(361, 391), colMeans(model$forecasts[region][[1]]), lwd = 2)
  
}

plot_blh(fc_mod, region = 1)
plot_blh(fc_mod, region = 2)
plot_blh(fc_mod, region = 3)
plot_blh(fc_mod, region = 4)
plot_blh(fc_mod, region = 5)
plot_blh(fc_mod, region = 6)
plot_blh(fc_mod, region = 7)






plot_blh(fc_mod6, region = 1)
plot_blh(fc_mod6, region = 2)
plot_blh(fc_mod6, region = 3)
plot_blh(fc_mod6, region = 4)
plot_blh(fc_mod6, region = 5)
plot_blh(fc_mod6, region = 6)
plot_blh(fc_mod6, region = 7)

