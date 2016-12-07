library(sp)
library(rgdal)
library(ggplot2)
library(rrfields)
library(rstan)

# Load data
d = readRDS("examples/starling/bbs.rds")# load in 50 - stop data, summed at route level
d = d[which(d$AOU=="4930"),] # use only starling data
d = d[which(d$statenum%in%c(33,69,89) & d$year > 1997),] # data from ID, WA, OR

# Project coordinates to UTM using sp()
d$ID = seq(1,nrow(d))
coordinates(d) = c("Longitude", "Latitude")
proj4string(d) <- CRS("+proj=longlat +datum=WGS84")  ## for example
d = as.data.frame(spTransform(d, CRS(paste("+proj=utm +zone=11"," ellps=WGS84",sep=''))))

# Hold out 10% for test data set
set.seed(123)
s = sample(seq(1,nrow(d)), size = floor(nrow(d)*0.1), replace=F)
holdout = d[s,]
d = d[-s,]

d$RouteID = as.numeric(as.factor(d$Route))
#plot of data
#ggplot(d, aes(Longitude,Latitude,color=log(sum))) + geom_point(alpha=1,size=1) + facet_wrap(~year)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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

#save.image("starling_rrfields.Rdata")

# plot df
df=data.frame("x"=extract(mvt_gamma$model)[["df"]])
pdf("figs/starling_df.pdf")
ggplot(data=df, aes(x=x)) + geom_histogram(alpha=0.5) + xlab("MVT degrees of freedom")
dev.off()

# Plot the predictions from each object to the fitted data
mvt_pred = predict(mvt_gamma, newdata = holdout, mcmc_draws = 1000, quantiles=c(0.025,0.975))$summary
mvn_pred = predict(mvn_gamma, newdata = holdout, mcmc_draws = 1000, quantiles=c(0.025,0.975))$summary

# Calculate coverage of confidence intervals (including observation model). Defaults to 95%
mvn_coverage = length(which(holdout[,"sum"] > mvn_pred$confint_lower & holdout[,"sum"] < mvn_pred$confint_upper))
mvt_coverage = length(which(holdout[,"sum"] > mvt_pred$confint_lower & holdout[,"sum"] < mvt_pred$confint_upper))

pdf("figs/predictions_vs_observed.pdf")
plot(log(mvn_pred$mean), log(holdout[,"sum"]), col="red",
  pch = 16, cex=0.6, xlab = "Prediction", ylab = "Observed")
points(log(mvt_pred$mean), log(holdout[,"sum"]), cex=0.8, col="blue")
abline(0,1)
dev.off()

# Make plot of the ratio of confidence intervals from the mvn to mvt
ratio = log((mvn_pred$confint_upper - mvn_pred$confint_lower) / (mvt_pred$confint_upper - mvt_pred$confint_lower))
pdf("histogram of ratio of CI widths.pdf")
hist(ratio, 40, col="grey", xlab = paste0("MVN CI widths : MVT CI widths"), main="")
dev.off()

# plot year effects with ribbon plot
mvt_year=extract(mvt_gamma$model)[["yearEffects"]]
mvt_year = data.frame("year"=seq(min(d$year),max(d$year)), "mean"=apply(mvt_year,2,mean), "lower"=apply(mvt_year,2,quantile,0.025),
  "upper"=apply(mvt_year,2,quantile,0.975),"model"="mvt")
mvn_year=extract(mvn_gamma$model)[["yearEffects"]]
mvn_year = data.frame("year"=seq(min(d$year),max(d$year)), "mean"=apply(mvn_year,2,mean), "lower"=apply(mvn_year,2,quantile,0.025),
  "upper"=apply(mvn_year,2,quantile,0.975),"model"="mvn")
mvt_year = rbind(mvt_year,mvn_year)

pdf("figs/starling_yearEffects.pdf")
ggplot(data=mvt_year) + geom_ribbon(aes(x=year,ymin=lower,ymax=upper,group=model,fill=model),alpha=0.2) + xlab("Year") + ylab("Estimated year effect")
dev.off()

