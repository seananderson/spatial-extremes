library(dplyr)
library(ggplot2)
library(rrfields)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#catch = read.csv("examples/rockfish/TrawlCatch.csv")
catch = read.csv("examples/rockfish/rockfishCatch.csv")
haul_chars = read.csv("examples/rockfish/TrawlHaulChars.csv")
unique(catch$common_name)
# halfbanded, shortbellt: < 40 lat
catch = left_join(haul_chars, filter(catch, common_name == "widow rockfish"))
catch$year = as.numeric(substr(catch$date_yyyymmdd,1,4))

# filter out: data since 2003, positive values, no missing temp
catch = filter(catch, year >= 2003 & !is.na(cpue_kg_per_ha_der))
#catch = catch[-which((catch$latitude_dd < 40 & catch$longitude_dd > -124.5) | catch$longitude_dd > -123 | catch$longitude_dd < -125),]
ggplot(catch, aes(longitude_dd, latitude_dd, color = log(cpue_kg_per_ha_der))) +geom_point() + facet_wrap(~year)

# convert to UTM
catch$ID = seq(1,nrow(catch))
sp::coordinates(catch) = c("longitude_dd","latitude_dd")
sp::proj4string(catch) <- sp::CRS("+proj=longlat +datum=WGS84")  ## for example
catch = as.data.frame(sp::spTransform(catch, sp::CRS(paste("+proj=utm +zone=10"," ellps=WGS84",sep=''))))

# rescale
catch$latitude_dd = catch$latitude_dd/1000000
catch$longitude_dd = catch$longitude_dd/1000000

ggplot(catch, aes(longitude_dd, latitude_dd, color = log(cpue_kg_per_ha_der))) + geom_point() + facet_wrap(~year)

#catch$depth_m = scale(catch$depth_m, center=TRUE, scale=FALSE)
#catch$temperature_at_gear_c_der = scale(catch$temperature_at_gear_c_der, center=TRUE, scale=FALSE)
#temperature_at_gear_c_der + depth_m

#catch$depth = scale(log(catch$depth_m), scale=FALSE)

catch = select(catch, ID, latitude_dd, longitude_dd, temperature_at_gear_c_der, cpue_kg_per_ha_der,
  depth_m, year)
holdout = sample(1:nrow(catch), size=round(nrow(catch)*0.1,0), replace=F)

mvt_gamma = rrfield(log(cpue_kg_per_ha_der) ~ as.factor(year), data = catch[-holdout,],
  time = "year", lon = "longitude_dd", lat = "latitude_dd", station = "ID",
  nknots = 20,
  obs_error = "normal",
  prior_gp_sigma = half_t(100, 0, 3),
  prior_gp_scale = half_t(100, 0, 5),
  prior_intercept = student_t(100, 0, 50),
  prior_beta = student_t(100, 0, 10),
  prior_sigma = half_t(100, 0, 3),
  estimate_ar = FALSE,
  estimate_df = TRUE,
  chains = 2, iter = 700)

