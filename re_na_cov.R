blh_dat1$year[which(is.na(blh_dat1$DDs))[1:4]] <- 2017
blh_dat1$JULIAN[which(is.na(blh_dat1$DDs))[1:4]] <- 125
blh_dat1$DDs[which(is.na(blh_dat1$DDs))[1:4]] <- blh_dat1$DDs[473:476] + 30
blh_dat1$time1[which(is.na(blh_dat1$time1))[1:4]] <- 491


blh_dat1$year[which(is.na(blh_dat1$DDs))[1:2]] <- 2018
blh_dat1$JULIAN[which(is.na(blh_dat1$DDs))[1:2]] <- 124
blh_dat1$DDs[which(is.na(blh_dat1$DDs))[1:2]] <- blh_dat1$DDs[825:826] + 30
blh_dat1$time1[which(is.na(blh_dat1$time1))[1:2]] <- 855


blh_dat1$year[which(is.na(blh_dat1$DDs))[1:4]] <- 2019
blh_dat1$JULIAN[which(is.na(blh_dat1$DDs))[1:4]] <- 123
blh_dat1$DDs[which(is.na(blh_dat1$DDs))[1:4]] <- blh_dat1$DDs[c(1173, 1176, 1174, 1175)] + 30
blh_dat1$time1[which(is.na(blh_dat1$time1))[1:4]] <- 1219

blh_dat1$year[which(is.na(blh_dat1$DDs))[1]] <- 2019
blh_dat1$JULIAN[which(is.na(blh_dat1$DDs))[1]] <- 130
blh_dat1$DDs[which(is.na(blh_dat1$DDs))[1]] <- blh_dat1$DDs[1181] + 30
blh_dat1$time1[which(is.na(blh_dat1$time1))[1]] <- 1226



blh_dat1$year[which(is.na(blh_dat1$DDs))[1:4]] <- 2020
blh_dat1$JULIAN[which(is.na(blh_dat1$DDs))[1:4]] <- 122
blh_dat1$DDs[which(is.na(blh_dat1$DDs))[1:4]] <- blh_dat1$DDs[1523:1526] + 30
blh_dat1$time1[which(is.na(blh_dat1$time1))[1:4]] <- 1583


blh_dat1$year[which(is.na(blh_dat1$DDs))[1:4]] <- 2021
blh_dat1$JULIAN[which(is.na(blh_dat1$DDs))[1:4]] <- 127
blh_dat1$DDs[which(is.na(blh_dat1$DDs))[1:4]] <- blh_dat1$DDs[1887:1890] + 30
blh_dat1$time1[which(is.na(blh_dat1$time1))[1:4]] <- 1954

blh_dat1 <- blh_dat1[, -c(1,2)]
