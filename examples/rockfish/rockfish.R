
library(dplyr)
library(ggplot2)

catch = read.csv("examples/rockfish/TrawlCatch.csv")
haul_chars = read.csv("examples/rockfish/TrawlHaulChars.csv")

catch = left_join(haul_chars, filter(catch, common_name == "canary rockfish"))
catch$year = as.numeric(substr(catch$date_yyyymmdd,1,4))
# filter out: data since 2003, positive values, no missing temp
catch = filter(catch, year >= 2012 & !is.na(cpue_kg_per_ha_der) & !is.na(catch$temperature_at_gear_c_der))

# convert to UTM
catch$ID = seq(1,nrow(catch))
sp::coordinates(catch) = c("longitude_dd","latitude_dd")
sp::proj4string(catch) <- sp::CRS("+proj=longlat +datum=WGS84")  ## for example
catch = as.data.frame(sp::spTransform(catch, sp::CRS(paste("+proj=utm +zone=10"," ellps=WGS84",sep=''))))

# rescale
catch$latitude_dd = catch$latitude_dd/1000000
catch$longitude_dd = catch$longitude_dd/1000000

library(rrfields)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

holdout = sample(1:nrow(catch), size=round(nrow(catch)*0.1,0), replace=F)

ggplot(catch, aes(longitude_dd, latitude_dd, color = log(cpue_kg_per_ha_der))) + geom_point() + facet_wrap(~year)

mvt_gamma = rrfield(cpue_kg_per_ha_der ~ temperature_at_gear_c_der + depth_m + as.factor(year),
  data = catch[-holdout,],
  time = "year", lon = "longitude_dd", lat = "latitude_dd", station = "ID",
  nknots = 30,
  obs_error = "gamma",
  prior_gp_sigma = half_t(100, 0, 5),
  prior_gp_scale = half_t(100, 0, 5),
  prior_intercept = student_t(100, 0, 10),
  prior_beta = student_t(100, 0, 1),
  prior_sigma = half_t(100, 0, 3),
  estimate_df = TRUE,
  chains = 3, iter = 700)

mvn_gamma = rrfield(HgResult ~ Length,
  data = d[-holdout,],
  time = "Year", lon = "lon", lat = "lat", station = "StationID",
  nknots = 15,
  obs_error = "gamma",
  prior_gp_sigma = half_t(100, 0, 1),
  prior_gp_scale = half_t(100, 0, 1),
  prior_intercept = student_t(100, 0, 1),
  prior_beta = student_t(100, 0, 1),
  prior_sigma = half_t(100, 0, 1),
  estimate_df = FALSE,
  chains = 4, iter = 5000)
