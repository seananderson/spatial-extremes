library(dplyr)
library(mgcv)

rm(list=ls())

dat = readRDS("crabs/data/joined_crab.rds")
dat$cpue = dat$cpue_kgperha * dat$area_swept_ha

dat$bottom_tempc = log(dat$bottom_tempc)
dat$o2_mlperL = log(dat$o2_mlperL)

mod = gam(cpue_kgperha ~ as.factor(year) + te(lon,lat, by = year) +
    I(bottom_depth) + I(bottom_depth^2) + (o2_mlperL) + (bottom_tempc),
  family = Tweedie(1.25),data = dat)
