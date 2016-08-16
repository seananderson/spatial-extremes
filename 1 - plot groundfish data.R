library(dplyr)
library(PBSmapping)
library(ggplot2)

dat = readRDS("data/triennial_annual_joined.rds")

group_by(dat, species) %>%
  summarize(n = length(which(!is.na(age_yrs)))) %>%
  arrange(-n)

data(nepacLL)

drkb = dat[which(dat$species=="darkblotched rockfish"),]
drkb$PID = 1
drkb$POS = seq(1,dim(drkb)[1])
drkb$X = drkb$haul_longitude_dd
drkb$Y = drkb$haul_latitude_dd
attr(drkb, "projection") = "LL"
attr(drkb, "zone") = 7
xlim = range(drkb$X)

ggplot(convUL(nepacLL), aes(X, Y, group = PID)) + geom_polygon() +
  geom_point(data = convUL(drkb), aes(X, Y,
    colour = log(age_yrs), aes = 0.02), size=0.1) +
  scale_colour_gradient(low = "blue") +
  coord_cartesian(xlim = c(1000,2500), ylim=c(4000,6000))
