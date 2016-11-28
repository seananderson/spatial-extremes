
library(sp)
library(rgdal)
library(ggplot2)
library(rrfields)

d = readRDS("examples/starling/bbs.rds")# load in 50 - stop data, summed at route level
d = d[which(d$AOU=="4930"),] # use only starling data
d = d[which(d$statenum%in%c(33,69,89) & d$year > 1997),] # data from ID, WA, OR

# Project coordinates to UTM using sp()
d$ID = seq(1,nrow(d))
coordinates(d) = c("Longitude", "Latitude")
proj4string(d) <- CRS("+proj=longlat +datum=WGS84")  ## for example
d = as.data.frame(spTransform(d, CRS(paste("+proj=utm +zone=11"," ellps=WGS84",sep=''))))

d$RouteID = as.numeric(as.factor(d$Route))
#plot of data
#ggplot(d, aes(Longitude,Latitude,color=log(sum))) + geom_point(alpha=1,size=1) + facet_wrap(~year)

# fit model with MVT random field
mvt_gamma = rrfield(sum ~ -1, data=d, time = "year", lon="Longitude", lat="Latitude",
  station = "RouteID", nknots = 30L, obs_error = "gamma", covariance="squared-exponential",
  estimate_df = TRUE, estimate_ar = TRUE, algorithm="sampling", year_re = TRUE,
  chains = 3L, iter=1000, control = list(adapt_delta = 0.99))

# fit model with MVN random field
mvn_gamma = rrfield(sum ~ -1, data=d, time = "year", lon="Longitude", lat="Latitude",
  station = "RouteID", nknots = 30L, obs_error = "gamma", covariance="squared-exponential",
  estimate_df = FALSE, fixed_df_value = 100, estimate_ar = TRUE, algorithm="sampling",
  year_re = TRUE, chains = 3L, iter=1000, control = list(adapt_delta = 0.99))

save.image("starling_rrfields.Rdata")

# plot df
df=data.frame("x"=extract(mvt_gamma)[["scaledf"]])
pdf("starling_df.pdf")
ggplot(data=df, aes(x=x)) + geom_histogram(alpha=0.5) + xlab("MVT degrees of freedom")
dev.off()

# plot year effects with ribbon plot
mvt_year=extract(mvt_gamma)[["yearEffects"]]
mvt_year = data.frame("year"=seq(min(d$year),max(d$year)), "mean"=apply(mvt_year,2,mean), "lower"=apply(mvt_year,2,quantile,0.025),
  "upper"=apply(mvt_year,2,quantile,0.975),"model"="mvt")
mvn_year=extract(mvn_gamma)[["yearEffects"]]
mvn_year = data.frame("year"=seq(min(d$year),max(d$year)), "mean"=apply(mvn_year,2,mean), "lower"=apply(mvn_year,2,quantile,0.025),
  "upper"=apply(mvn_year,2,quantile,0.975),"model"="mvn")
mvt_year = rbind(mvt_year,mvn_year)

pdf("starling_yearEffects.pdf")
ggplot(data=mvt_year) + geom_ribbon(aes(x=year,ymin=lower,ymax=upper,group=model,fill=model),alpha=0.2) + xlab("Year") + ylab("Estimated year effect")
dev.off()

