
# 2. Find unique locations

locations_WA <- subset(blh_dat, 
                       duplicated(blh_dat[c("LATITUDE", 
                                              "LONGITUDE", "year")]) == 
                         FALSE)[, c(1, 2, 3, 6)]

write.csv(locations_WA, file = "locationsBLH.csv")

# 3. Extract max and min temp

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



# 4. Calculate cumulative degree days using retrieved max and min temps 
# and the DAS function for DDs calculation

upT_BLH <- 10 + 19.8
baseT_BLH <- 10

DDs_WA <- list()

for(i in 1: length(locations_WA$region)) {
  DDs_WA[[i]] <- calc_dd_vec(tmax = tmaxWA[[i]], tmin = tminWA[[i]], 
                             lower_threshold = baseT_BLH, 
                             upper_threshold = upT_BLH, 
                             cutoff = "horizontal")
  DDs_WA[[i]] <- cumsum(DDs_WA[[i]])
}


# Matching data collection days

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
blh_dat <- blh_dat[, -c(2, 3, 11)]
