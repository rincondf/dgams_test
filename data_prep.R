#------------------------------------------------------------------------------
# Proper formatting
#------------------------------------------------------------------------------

# Import BLH dataset already pooled by region

blh_regionWA <- read.csv("blh_regionWA.csv")
blh_regionWA <- blh_regionWA[, -1]

# Load functions to calculate degree days (used to fill gaps in unsampled weeks)

source("DD_clacFool.R")

# DGAMS can take NAs in the response variable, but not in the time or covariates
# So we'll fill the times and degree days in between seasons. We only use data 
# from 2016 to 2023, since sample size is too small from previous years. 

Y2016 <- blh_regionWA[which(blh_regionWA$year == 2016), -c(9, 10)]
Y2017 <- blh_regionWA[which(blh_regionWA$year == 2017), -c(9, 10)]
Y2018 <- blh_regionWA[which(blh_regionWA$year == 2018), -c(9, 10)]
Y2019 <- blh_regionWA[which(blh_regionWA$year == 2019), -c(9, 10)]
Y2020 <- blh_regionWA[which(blh_regionWA$year == 2020), -c(9, 10)]
Y2021 <- blh_regionWA[which(blh_regionWA$year == 2021), -c(9, 10)]
Y2022 <- blh_regionWA[which(blh_regionWA$year == 2022), -c(9, 10)]
Y2023 <- blh_regionWA[which(blh_regionWA$year == 2023), -c(9, 10)]

# This function fills between-season gaps fro each year

convi <- function(x, yr) {
  sems <- seq(1, (x$JULIAN[which.min(x$JULIAN)] - 7), 7)
  
  
  xA <- data.frame("region" = rep(unique(x$region), length(sems)),
                   "LATITUDE" = rep(unique(x$LATITUDE), length(sems)),
                   "LONGITUDE" = rep(unique(x$LONGITUDE), length(sems)),
                   "SURVEY_DATE" = as.Date(rep(sems, 
                                               each = length(unique(x$region))), 
                                           origin = as.Date(paste0((yr - 1), 
                                                                   "-12-31"))),
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
                                           origin = as.Date(paste0((yr - 1), 
                                                                   "-12-31"))),
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

# Merge all years again

blh_dat <- rbind(Y2016A,
                 Y2017A,
                 Y2018A,
                 Y2019A,
                 Y2020A,
                 Y2021A,
                 Y2022A,
                 Y2023A)

# Create time variable with julian days starting jan 1, 2026

blh_dat$time1 <- as.numeric(blh_dat$SURVEY_DATE) - 
  (as.numeric(as.Date("2016-01-01")) - 1)

# Remove last rows with only NAs
blh_dat <- blh_dat[-seq(2719, 2802), ]

# Create a "series" factor variable with the regions
blh_dat$series <- as.factor(blh_dat$region)

# Create a time variable with weeks starting jan 1, 2026. This is the time
# variable we'll use

semas <- seq(1, length(unique(blh_dat$time1)))

blh_dat$time <- rep(NA, length(blh_dat$region))

for(i in 1: length(unique(unique(blh_dat$time1)))) {
  blh_dat$time[which(blh_dat$time1 == unique(blh_dat$time1)[i])] <-
    rep(semas[i], length(unique(blh_dat$time1)[i]))
}


#------------------------------------------------------------------------------
# Calculation of degree days in between seasons
#------------------------------------------------------------------------------

# Find unique location-year combinations

locations_WA <- subset(blh_dat, 
                       duplicated(blh_dat[c("LATITUDE", 
                                            "LONGITUDE", "year")]) == 
                         FALSE)[, c(1, 2, 3, 6)]


# Extract max and min temp

library(daymetr)

temp_recs <- list()

for(i in 1: length(locations_WA$region)) {
  temp_recs[[i]] <- download_daymet(
    site = "Daymet",
    lat = locations_WA$LATITUDE[i],
    lon = locations_WA$LONGITUDE[i],
    start = locations_WA$year[i],
    end = locations_WA$year[i],
    path = tempdir(),
    internal = TRUE,
    silent = FALSE,
    force = FALSE,
    simplify = FALSE
  )
}


tmaxWA <- list()

for(i in 1: length(locations_WA$region)) {
  tmaxWA[[i]] <- as.numeric(temp_recs[[i]]$data[,7])
}

tminWA <- list()

for(i in 1: length(locations_WA$region)) {
  tminWA[[i]] <- as.numeric(temp_recs[[i]]$data[,8])
}


# Calculate cumulative degree days using retrieved max and min temps 
# and the DAS function for DDs calculation

# upper and lower temp thresholds for BLH

upT_BLH <- 10 + 19.8
baseT_BLH <- 10

# Calculate degree days for all years and locations

DDs_WA <- list()

for(i in 1: length(locations_WA$region)) {
  DDs_WA[[i]] <- calc_dd_vec(tmax = tmaxWA[[i]], tmin = tminWA[[i]], 
                             lower_threshold = baseT_BLH, 
                             upper_threshold = upT_BLH, 
                             cutoff = "horizontal")
  DDs_WA[[i]] <- cumsum(DDs_WA[[i]])
}


# Matching data with days

blh_dat$DDS <- rep(NA, length(blh_dat$SURVEY_DATE))


for(i in 1: length(locations_WA$region)) {
  blh_dat$DDS[which(blh_dat$LATITUDE == locations_WA$LATITUDE[i] & 
                      blh_dat$LONGITUDE == locations_WA$LONGITUDE[i])] <-
    DDs_WA[[i]][blh_dat$JULIAN[which(blh_dat$LATITUDE == 
                                       locations_WA$LATITUDE[i] & 
                                       blh_dat$LONGITUDE == 
                                       locations_WA$LONGITUDE[i])]]
}

blh_dat$DDs[which(is.na(blh_dat$DDs))] <- blh_dat$DDS[which(is.na(blh_dat$DDs))]

# Cleaning up some unnecessary columns
blh_dat <- subset(blh_dat, select = -c(region, SURVEY_DATE, LATITUDE, LONGITUDE, DDS))


# As some of the regions have missing counts at some time points, they must be 
# filled with NAs


#------------------------------------------------------------------------------
# Making sure the time variable is complete and all times have a row
#------------------------------------------------------------------------------

library(dplyr)

blh_dat <- right_join(blh_dat, expand.grid(
  time = seq(
    min(blh_dat$time),
    max(blh_dat$time)
  ),
  series = factor(unique(blh_dat$series),
                  levels = levels(blh_dat$series)
  )
))

blh_dat <- arrange(blh_dat, blh_dat$time)

# Create a new response variable without zeros, so that a Gamma pdf can be used
blh_dat$counts <- blh_dat$TRAP_COUNT_RAW + 0.00001


# Fill NAs with data for covariates

blh_dat$year[which(is.na(blh_dat$DDs))[1:4]] <- 2017
blh_dat$JULIAN[which(is.na(blh_dat$DDs))[1:4]] <- 125
blh_dat$DDs[which(is.na(blh_dat$DDs))[1:4]] <- blh_dat$DDs[473:476] + 30
blh_dat$time1[which(is.na(blh_dat$time1))[1:4]] <- 491


blh_dat$year[which(is.na(blh_dat$DDs))[1:2]] <- 2018
blh_dat$JULIAN[which(is.na(blh_dat$DDs))[1:2]] <- 124
blh_dat$DDs[which(is.na(blh_dat$DDs))[1:2]] <- blh_dat$DDs[825:826] + 30
blh_dat$time1[which(is.na(blh_dat$time1))[1:2]] <- 855


blh_dat$year[which(is.na(blh_dat$DDs))[1:4]] <- 2019
blh_dat$JULIAN[which(is.na(blh_dat$DDs))[1:4]] <- 123
blh_dat$DDs[which(is.na(blh_dat$DDs))[1:4]] <- blh_dat$DDs[c(1173, 1176, 1174, 1175)] + 30
blh_dat$time1[which(is.na(blh_dat$time1))[1:4]] <- 1219

blh_dat$year[which(is.na(blh_dat$DDs))[1]] <- 2019
blh_dat$JULIAN[which(is.na(blh_dat$DDs))[1]] <- 130
blh_dat$DDs[which(is.na(blh_dat$DDs))[1]] <- blh_dat$DDs[1181] + 30
blh_dat$time1[which(is.na(blh_dat$time1))[1]] <- 1226



blh_dat$year[which(is.na(blh_dat$DDs))[1:4]] <- 2020
blh_dat$JULIAN[which(is.na(blh_dat$DDs))[1:4]] <- 122
blh_dat$DDs[which(is.na(blh_dat$DDs))[1:4]] <- blh_dat$DDs[1523:1526] + 30
blh_dat$time1[which(is.na(blh_dat$time1))[1:4]] <- 1583


blh_dat$year[which(is.na(blh_dat$DDs))[1:4]] <- 2021
blh_dat$JULIAN[which(is.na(blh_dat$DDs))[1:4]] <- 127
blh_dat$DDs[which(is.na(blh_dat$DDs))[1:4]] <- blh_dat$DDs[1887:1890] + 30
blh_dat$time1[which(is.na(blh_dat$time1))[1:4]] <- 1954

# Cleaning up an unnecessary column

blh_dat <- subset(blh_dat, select = -c(TRAP_COUNT_RAW))
