# data from https://github.com/USGS-CIDA/CIDA-Viz
# see http://cida.usgs.gov/ca_drought

sites = read.csv('examples/snow-water/data/sites.csv', as.is=TRUE)
all_data = data.frame()

for(i in 1:nrow(sites)){
  fname = paste0('examples/snow-water/data/data_sane/', sites$Station[i], '.csv')

  data = read.table(fname, sep=',', header=TRUE, as.is=TRUE)
  if(nrow(data)==0){
    next
  }
  data$site = sites$Station[i]
  all_data = rbind(all_data, data)
}

head(sites)
sites <- rename(sites, site = Station)
names(sites) <- tolower(names(sites))

library(tidyverse)
library(lubridate)
d <- all_data %>% mutate(date = ymd(datetime), year = year(date), month = month(date), day = day(date)) %>%
  inner_join(sites) %>%
  rename(snow = snow.water.eq.inches)

ggplot(d, aes(date, log(snow+1))) + geom_line() +
  facet_wrap(~site)

hist(d$year)
hist(d$month)
hist(d$day)

d2 <- filter(d, month == 4) %>%
  group_by(year) %>%
  filter(day == min(day)) %>% ungroup() %>%
  filter(year >= 1979) # data sparse spatially before

ggplot(d2, aes(lon, lat, colour = log(snow+1))) + geom_point() +
  facet_wrap(~year)

d2$lat <- d2$lat * 10
d2$lon <- d2$lon * 10

hist(d2$snow)
min(d2$snow)
sum(d2$snow==0)

length(unique(d2$site))

library(rrfields)
m <- rrfield(log(snow+1) ~ 1, data = d2, time = "year", lon = "lon", lat = "lat",
  nknots = 15, obs_error = "gamma", station = "site",
  estimate_df = TRUE, estimate_ar = TRUE, fixed_ar_value = 0, algorithm = "sampling",
  chains = 2, iter = 800,
  prior_intercept = student_t(3, 0, 10), prior_gp_sigma = half_t(3, 0, 1),
  control = list(adapt_delta = 0.9), cores = 4)
m
