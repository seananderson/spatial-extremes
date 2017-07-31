library(dplyr)
library(ggplot2)
library(rrfields)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

catch = read.csv("examples/rockfish/rockfishCatch.csv")
haul_chars = read.csv("examples/rockfish/TrawlHaulChars.csv")
unique(catch$common_name)
# halfbanded, shortbellt: < 40 lat
catch = left_join(haul_chars, filter(catch, common_name == "darkblotched rockfish"))
catch$year = as.numeric(substr(catch$date_yyyymmdd,1,4))

# filter out: data since 2003, positive values, no missing temp
catch = filter(catch, year >= 2003 & !is.na(cpue_kg_per_ha_der))
catch = catch[-which((catch$latitude_dd < 40 & catch$longitude_dd > -124.5) | catch$longitude_dd > -123 | catch$longitude_dd < -125),]
ggplot(catch, aes(longitude_dd, latitude_dd, color = log(cpue_kg_per_ha_der))) +
  geom_point() + facet_wrap(~year) +
  viridis::scale_color_viridis()

# convert to UTM
catch$ID = seq(1,nrow(catch))
sp::coordinates(catch) = c("longitude_dd","latitude_dd")
sp::proj4string(catch) <- sp::CRS("+proj=longlat +datum=WGS84")  ## for example
catch = as.data.frame(sp::spTransform(catch, sp::CRS(paste("+proj=utm +zone=10"," ellps=WGS84",sep=''))))

# rescale
catch$latitude_dd = catch$latitude_dd/100000
catch$longitude_dd = catch$longitude_dd/100000

holdout = sample(1:nrow(catch), size=round(nrow(catch)*0.1,0), replace=F)

ggplot(catch, aes(longitude_dd, latitude_dd, color = log(cpue_kg_per_ha_der))) +
  geom_point() + facet_wrap(~year) + viridis::scale_color_viridis()

#catch$depth_m = scale(catch$depth_m, center=TRUE, scale=FALSE)
#catch$temperature_at_gear_c_der = scale(catch$temperature_at_gear_c_der, center=TRUE, scale=FALSE)
#temperature_at_gear_c_der + depth_m

catch$depth_scaled = as.numeric(arm::rescale(log(catch$depth_m)))

ggplot(catch, aes(longitude_dd, latitude_dd, color = depth_scaled)) +
  geom_point() + facet_wrap(~year) + viridis::scale_color_viridis()

mvt_gamma = rrfield(log(cpue_kg_per_ha_der*10000) ~ as.factor(year) + poly(depth_scaled, 2), data = catch,
  time = "year", lon = "longitude_dd", lat = "latitude_dd",
  nknots = 20, fixed_ar_value = 1,
  obs_error = "normal",
  prior_gp_sigma = half_t(7, 0, 3),
  prior_gp_scale = half_t(7, 0, 3),
  prior_intercept = student_t(100, 0, 10),
  prior_beta = student_t(100, 0, 2),
  prior_sigma = half_t(7, 0, 3),
  estimate_ar = FALSE,
  estimate_df = TRUE,
  chains = 6, iter = 400)

library(nlme)
mod = gls(log(cpue_kg_per_ha_der) ~ 1, correlation = corExp(form = ~longitude_dd +
    latitude_dd, nugget = TRUE), data = catch[catch$year==2005,])


