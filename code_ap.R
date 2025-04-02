
library(mvgam)

# Checking all series are complete

get_mvgam_priors(count ~ 1,
                 data = aphids2,
                 family = Gamma()
)


plot_mvgam_series(data = aphids2, y = "count")

train_ap <- aphids2[which(aphids2$time <= 800), ]
test_ap <- aphids2[which(aphids2$time > 800), ]


plot_mvgam_series(data = train_ap, 
                  newdata = test_ap,
                  y = "count")


mod_ap <- mvgam(count ~
                  
                  # relative abundances
                  s(DDs, k = 6, bs = 'tp') +
                  s(year, k = 4, bs = 'tp'),
                
                # Condition on the training data
                data = train_ap,
                
                # Automatically compute forecasts for the test data
                newdata = test_ap,
                
                trend_model = AR(2),
                
                family = Gamma(),
                
                # cmdstanr is highly recommended over rstan as 
                # it is much more up-to-date with the latest 
                # features in the Stan universe
                noncentred = TRUE,
                backend = 'cmdstanr')
summary(mod_ap)
plot(mod_ap)

gratia::draw(mod_ap)
conditional_effects(mod_ap, points = 0.5)
conditional_effects(mod_ap,
                    type = 'link')





hcs_ap <- hindcast(mod_ap)

plot(hcs_ap)

plot(mod_ap, type = 'forecast')



