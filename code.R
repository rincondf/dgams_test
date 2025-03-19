Y2016 <- blh_regionWA[which(blh_regionWA$year == 2016), -c(9, 10)]
Y2017 <- blh_regionWA[which(blh_regionWA$year == 2017), -c(9, 10)]
Y2018 <- blh_regionWA[which(blh_regionWA$year == 2018), -c(9, 10)]
Y2019 <- blh_regionWA[which(blh_regionWA$year == 2019), -c(9, 10)]
Y2020 <- blh_regionWA[which(blh_regionWA$year == 2020), -c(9, 10)]
Y2021 <- blh_regionWA[which(blh_regionWA$year == 2021), -c(9, 10)]
Y2022 <- blh_regionWA[which(blh_regionWA$year == 2022), -c(9, 10)]
Y2023 <- blh_regionWA[which(blh_regionWA$year == 2023), -c(9, 10)]




convi <- function(x, yr) {
  sems <- seq(1, (x$JULIAN[which.min(x$JULIAN)] - 7), 7)
  
  
  xA <- data.frame("region" = rep(unique(x$region), length(sems)),
                   "LATITUDE" = rep(unique(x$LATITUDE), length(sems)),
                   "LONGITUDE" = rep(unique(x$LONGITUDE), length(sems)),
                   "SURVEY_DATE" = as.Date(rep(sems, 
                                               each = length(unique(x$region))), 
                                           origin = as.Date(paste0((yr - 1), "-12-31"))),
                   "TRAP_COUNT_RAW" = rep(NA, length(sems)),
                   "year" = rep(yr, length(sems)),
                   "JULIAN" = rep(sems, each = length(unique(x$region))),
                   "DDs" = rep(NA, length(sems))
                   )
  
  
  sems1 <- seq((x$JULIAN[which.max(x$JULIAN)] + 7), (360-7), 7)
  
  
  xB <- data.frame("region" = rep(unique(x$region), length(sems1)),
                   "LATITUDE" = rep(unique(x$LATITUDE), length(sems1)),
                   "LONGITUDE" = rep(unique(x$LONGITUDE), length(sems1)),
                   "SURVEY_DATE" = as.Date(rep(sems1, 
                                               each = length(unique(x$region))), 
                                           origin = as.Date(paste0((yr - 1), "-12-31"))),
                   "TRAP_COUNT_RAW" = rep(NA, length(sems1)),
                   "year" = rep(yr, length(sems1)),
                   "JULIAN" = rep(sems1, each = length(unique(x$region))),
                   "DDs" = rep(NA, length(sems1))
  )
  
  xF <- rbind(xA, x, xB)
  
  xF
}



Y2016A <- convi(Y2016, 2016)
Y2017A <- convi(Y2017, 2017)
Y2018A <- convi(Y2018, 2018)
Y2019A <- convi(Y2019, 2019)
Y2020A <- convi(Y2020, 2020)
Y2021A <- convi(Y2021, 2021)
Y2022A <- convi(Y2022, 2022)
Y2023A <- convi(Y2023, 2023)


blh_dat <- rbind(Y2016A,
                 Y2017A,
                 Y2018A,
                 Y2019A,
                 Y2020A,
                 Y2021A,
                 Y2022A,
                 Y2023A)
blh_dat$time <- as.numeric(blh_dat$SURVEY_DATE) - (as.numeric(as.Date("2016-01-01")) - 1)

blh_dat <- blh_dat[-seq(2719, 2802), ]
blh_dat$series <- as.factor(blh_dat$region)

semas <- seq(1, length(unique(blh_dat$time)))


blh_dat$time <- rep(NA, length(blh_dat$region))

for(i in 1: length(unique(unique(blh_dat$time1)))) {
  blh_dat$time[which(blh_dat$time1 == unique(blh_dat$time1)[i])] <-
    rep(semas[i], length(unique(blh_dat$time1)[i]))
}

library(mvgam)

all(levels(blh_dat$series) %in% 
      unique(blh_dat$series))

blh_dat %>%
  dplyr::right_join(expand.grid(
    time = seq(
      min(blh_dat$time),
      max(blh_dat$time)
    ),
    series = factor(unique(blh_dat$series),
                    levels = levels(blh_dat$series)
    )
  )) %>%
  dplyr::arrange(time) -> blh_dat1


blh_dat1$counts <- blh_dat1$TRAP_COUNT_RAW + 0.00001

get_mvgam_priors(counts ~ 1,
                 data = blh_dat1,
                 family = Gamma()
)




plot_mvgam_series(data = blh_dat1, y = "counts", series = "all")

dat_train <- blh_dat1[which(blh_dat1$time <= 313), ]
dat_test <- blh_dat1[which(blh_dat1$time > 313), ]


plot_mvgam_series(data = dat_train, 
                  newdata = dat_test,
                  y = "counts")


mod <- mvgam(counts ~ 
               
               # Hierarchical intercepts capture variation in average
               # relative abundances
               s(series, bs = 're') +
               s(DDs, k = 7, bs = 'fs'),
             
             # Condition on the training data
             data = dat_train,
             
             # Automatically compute forecasts for the test data
             newdata = dat_test,
             
             tred_model = AR(p = 3),
             
             # Beta observations with independent precisions
             family = Gamma(),
             
             # cmdstanr is highly recommended over rstan as 
             # it is much more up-to-date with the latest 
             # features in the Stan universe
             backend = 'cmdstanr')
summary(mod)
plot(mod, type = 're')

plot(mod, type = 'forecast', series = 1)
plot(mod, type = 'forecast', series = 2)
plot(mod, type = 'forecast', series = 3)
plot(mod, type = 'forecast', series = 4)
plot(mod, type = 'forecast', series = 5)
plot(mod, type = 'forecast', series = 6)
plot(mod, type = 'forecast', series = 7)


plot_predictions(mod, by = "time", points = 0.5)

