library(dplyr)
library(ggplot2)
library(rrfields)
library(nlme)

# read in data on 31 rockfish spp
catch_dat = read.csv("examples/rockfish/rockfishCatch.csv")
haul_chars = read.csv("examples/rockfish/TrawlHaulChars.csv")

spp = unique(catch_dat$common_name)

# filter data since 2003 -- consistent survey
catch_dat$year = as.numeric(substr(catch_dat$date_dim.full_date,1,4))
catch_dat = filter(catch_dat, catch_dat$year >= 2003)

for(i in 22:length(spp)) {

# merge df for a single spp
catch = left_join(haul_chars, filter(catch_dat, common_name == spp[i])) %>%
  filter(!is.na(cpue_kg_per_ha_der)) %>%
  select(latitude_dd, longitude_dd, temperature_at_gear_c_der, cpue_kg_per_ha_der,
    depth_m, year)

# filter data to remove spatial extremes /outliers in coverage
catch = filter(catch, longitude_dd < quantile(longitude_dd,0.975) &
    longitude_dd > quantile(longitude_dd,0.025) & quantile(latitude_dd,0.975) &
    latitude_dd > quantile(latitude_dd,0.025))
ggplot(catch, aes(longitude_dd, latitude_dd, color = log(cpue_kg_per_ha_der))) + geom_point() + facet_wrap(~year)

# convert to UTM
catch$ID = seq(1,nrow(catch))
sp::coordinates(catch) = c("longitude_dd","latitude_dd")
sp::proj4string(catch) <- sp::CRS("+proj=longlat +datum=WGS84")  ## for example
catch = as.data.frame(sp::spTransform(catch, sp::CRS(paste("+proj=utm +zone=10"," ellps=WGS84",sep=''))))

# rescale lat & lon
catch$lat = catch$latitude_dd/1000000
catch$lon = catch$longitude_dd/1000000

# fit logged cpue to spatial model
mod = gls(log(cpue_kg_per_ha_der) ~ -1 + as.factor(year),
  correlation = corExp(form = ~lon + lat, nugget = TRUE),
  data = catch)

# something like RMSE by year should so anomalies where t-dist would
# improve things / fit outliers better
catch$residuals = mod$residuals

summary = group_by(catch, year) %>%
  summarize(rmse = sqrt(mean(residuals^2)))
summary$spp = spp[i]
g1 = ggplot(summary, aes(year,rmse)) +
  geom_line() + ggtitle(paste0(spp[i]))
if(i == 1) {
  out_df = summary
}
else {
  out_df = rbind(out_df, summary)
}
ggsave(filename = paste0("examples/rockfish/",spp[i],".png"), plot=g1)
}

saveRDS(out_df,"explore_out_df.rds")
group_by(out_df,spp) %>%
  mutate(rmse2 = rmse - mean(rmse)) %>%
ggplot(aes(year, rmse2, group = spp, col=spp)) + geom_line() + facet_wrap(~spp)
