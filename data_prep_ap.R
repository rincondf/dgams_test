#------------------------------------------------------------------------------
# Proper formatting
#------------------------------------------------------------------------------

library(dplyr)

# Import aphid dataset

aphids <- read.csv("aphids_pan.csv")
aphids <- aphids[, -1]

# Reconstruction of sample dates

aphids$date <- as.Date((aphids$JulianDay - 1), 
                       origin = paste0(aphids$Year, "-01-01"))

# name columns according to subsequent functions

names(aphids) <- c("Site", "LATITUDE", "LONGITUDE", "year", "JULIAN", 
                   "COUNT", "DDs", "SURVEY_DATE")

# Find unique location-year combinations

locs <- subset(aphids, 
               duplicated(aphids[c("LATITUDE", 
                                   "LONGITUDE", "year")]) == 
                 FALSE)[, seq(1, 4)]

locs <- arrange(locs, locs$year)

# clean up datset

aphids <- aggregate(aphids[, -seq(1,3)], list(aphids$SURVEY_DATE), mean, na.rm = TRUE)[, -1]

# Load functions to calculate degree days (used to fill gaps in unsampled weeks)

source("DD_clacFool.R")

# DGAMs can take NAs in the response variable, but not in time or covariates.
# So, we'll fill in the times and degree days in between seasons with data. 
# We only use data from 2016 to 2023, since sample size is too small in previous 
# years. 

Y2007 <- aphids[which(aphids$year == 2007), ]
Y2008 <- aphids[which(aphids$year == 2008), ]
Y2009 <- aphids[which(aphids$year == 2009), ]
Y2010 <- aphids[which(aphids$year == 2010), ]
Y2011 <- aphids[which(aphids$year == 2011), ]
Y2012 <- aphids[which(aphids$year == 2012), ]
Y2013 <- aphids[which(aphids$year == 2013), ]
Y2014 <- aphids[which(aphids$year == 2014), ]
Y2015 <- aphids[which(aphids$year == 2015), ]
Y2016 <- aphids[which(aphids$year == 2016), ]
Y2017 <- aphids[which(aphids$year == 2017), ]
Y2018 <- aphids[which(aphids$year == 2018), ]
Y2019 <- aphids[which(aphids$year == 2019), ]
Y2020 <- aphids[which(aphids$year == 2020), ]
Y2021 <- aphids[which(aphids$year == 2021), ]
Y2022 <- aphids[which(aphids$year == 2022), ]




# This function fills between-season gaps for each year

convi2 <- function(x, yr) {
  sems <- seq(1, (x$JULIAN[which.min(x$JULIAN)] - 7), 7)
  
  
  xA <- data.frame(
    "SURVEY_DATE" = as.Date(sems, 
                            origin = as.Date(paste0((yr - 1), 
                                                    "-12-31"))),
    "COUNT" = rep(NA, length(sems)),
    "year" = rep(yr, length(sems)),
    "JULIAN" = sems,
    "DDs" = rep(NA, length(sems))
  )
  
  
  sems1 <- seq((x$JULIAN[which.max(x$JULIAN)] + 7), (360-7), 7)
  
  
  xB <- data.frame(
    "SURVEY_DATE" = as.Date(sems1, 
                            origin = as.Date(paste0((yr - 1), 
                                                    "-12-31"))),
    "COUNT" = rep(NA, length(sems1)),
    "year" = rep(yr, length(sems1)),
    "JULIAN" = sems1,
    "DDs" = rep(NA, length(sems1))
  )
  
  xF <- rbind(xA, x, xB)
  
  xF
}


Y2007A <- convi2(Y2007, 2007)
Y2008A <- convi2(Y2008, 2008)
Y2009A <- convi2(Y2009, 2009)
Y2010A <- convi2(Y2010, 2010)
Y2011A <- convi2(Y2011, 2011)
Y2012A <- convi2(Y2012, 2012)
Y2013A <- convi2(Y2013, 2013)
Y2014A <- convi2(Y2014, 2014)
Y2015A <- convi2(Y2015, 2015)
Y2016A <- convi2(Y2016, 2016)
Y2017A <- convi2(Y2017, 2017)
Y2018A <- convi2(Y2018, 2018)
Y2019A <- convi2(Y2019, 2019)
Y2020A <- convi2(Y2020, 2020)
Y2021A <- convi2(Y2021, 2021)
Y2022A <- convi2(Y2022, 2022)




# Merge all years again

aphids2 <- rbind(
  Y2007A,
  Y2008A,
  Y2009A,
  Y2010A,
  Y2011A,
  Y2012A,
  Y2013A,
  Y2014A,
  Y2015A,
  Y2016A,
  Y2017A,
  Y2018A,
  Y2019A,
  Y2020A,
  Y2021A,
  Y2022A
)

# Create time variable with julian days starting jan 1, 2016

aphids2$time1 <- as.numeric(aphids2$SURVEY_DATE) - 
  (as.numeric(as.Date("2007-01-01")) - 1)

aphids2 <- arrange(aphids2, aphids2$time1)

# Remove last rows with only NAs
aphids2 <- aphids2[-seq(827, 846), ]


# Create a time variable with weeks starting jan 1, 2007. This is the time
# variable we'll use

semas <- seq(1, length(unique(aphids2$time1)))

aphids2$time <- rep(NA, length(aphids2$SURVEY_DATE))

for(i in 1: length(unique(aphids2$time1))) {
  aphids2$time[which(aphids2$time1 == unique(aphids2$time1)[i])] <-
    rep(semas[i], length(unique(aphids2$time1)[i]))
}


#------------------------------------------------------------------------------
# Calculation of degree days in-between seasons
#------------------------------------------------------------------------------


# Extract max and min temp

library(daymetr)

temp_recs <- list()

for(i in 1: length(locs$Site)) {
  temp_recs[[i]] <- download_daymet(
    site = "Daymet",
    lat = locs$LATITUDE[i],
    lon = locs$LONGITUDE[i],
    start = locs$year[i],
    end = locs$year[i],
    path = tempdir(),
    internal = TRUE,
    silent = FALSE,
    force = FALSE,
    simplify = FALSE
  )
}


tmaxAp <- list()

for(i in 1: length(locs$Site)) {
  tmaxAp[[i]] <- as.numeric(temp_recs[[i]]$data[,7])
}

tminAp <- list()

for(i in 1: length(locs$Site)) {
  tminAp[[i]] <- as.numeric(temp_recs[[i]]$data[,8])
}


# Calculate cumulative degree days using retrieved max and min temps 
# and the DAS function for DDs calculation

# upper and lower temp thresholds for BLH

upT_aph <- 30
baseT_aph <- 5.5

# Calculate degree days for all years and locations

DDs_Aph <- list()

for(i in 1: length(locs$Site)) {
  DDs_Aph[[i]] <- calc_dd_vec(tmax = tmaxAp[[i]], tmin = tminAp[[i]], 
                             lower_threshold = baseT_aph, 
                             upper_threshold = upT_aph, 
                             cutoff = "horizontal")
  DDs_Aph[[i]] <- cumsum(DDs_Aph[[i]])
}

years_ap <- seq(2007, 2022)
dd2 <- list()

for(k in 1: length(years_ap)) {
  dd1 <- rep(NA, 365)
  
  dd <- rep(NA, length(which(locs$year == years_ap[k])))
  
  for(j in 1: 365) {
    for(i in 1: length(which(locs$year == years_ap[k]))) {
      dd[i] <- DDs_Aph[which(locs$year == years_ap[k])][[i]][j]
    }
    dd1[j] <- mean(dd)
  }
  
  dd2[[k]] <- dd1
}


leapY <- seq(2008, 2020, 4)

for(i in 1: length(leapY)) {
  dd2[[which(years_ap %in% leapY)[i]]] <- c(dd2[[which(years_ap %in% leapY)[i]]], dd2[[which(years_ap %in% leapY)[i]]][365])
}

dd2 <- unlist(dd2)


# Matching data with days

aphids2$DDS <- dd2[aphids2$time1]

aphids2$DDs[which(is.na(aphids2$DDs))] <- aphids2$DDS[which(is.na(aphids2$DDs))]



#------------------------------------------------------------------------------
# Making sure the time variable is complete and all times have a row
#------------------------------------------------------------------------------



# Create a new response variable without zeros, so that a Gamma pdf can be used
aphids2$count <- aphids2$COUNT + 0.00001

