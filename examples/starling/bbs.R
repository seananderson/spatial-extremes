
library(sp)
library(rgdal)
library(ggplot2)

d = readRDS("examples/bbs.rds")
d = d[which(d$AOU=="4930"),] # starlings
d = d[which(d$statenum%in%c(33,69,89) & d$year > 1997),] # ID, WA, OR

d$ID = seq(1,nrow(d))
coordinates(d) = c("Longitude", "Latitude")
proj4string(d) <- CRS("+proj=longlat +datum=WGS84")  ## for example
d = as.data.frame(spTransform(d, CRS(paste("+proj=utm +zone=11"," ellps=WGS84",sep=''))))

#ggplot(d, aes(Longitude,Latitude,color=log(sum))) + geom_point(alpha=1,size=1) + facet_wrap(~year)

nKnots = 30
knots = cluster::pam(d[,c("Longitude","Latitude")],nKnots)$medoids
distKnots = as.matrix(dist(knots)/100000)
distKnotsSq = distKnots^2 # squared distances

# Calculate distance from knots to grid
distAll = as.matrix(dist(rbind(d[,c("Longitude","Latitude")], knots))/100000)^2
# this is the transpose of the lower left corner
nLocs = nrow(d)
distKnots21Sq = t(distAll[-c(1:nLocs), -c((nLocs+1):ncol(distAll))])

Y = as.numeric(d$sum)
yearID = as.numeric((d$year))
yearID = 1 + yearID - min(yearID)
stationID = seq(1,nrow(d))

# create list for STAN
spatglm_data = list("nKnots"=nKnots, "nLocs"=nLocs, "nT" = length(unique(yearID)), "N" = length(Y), "stationID" = stationID, "yearID" = yearID, "y" = Y, "distKnotsSq" = distKnotsSq, "distKnots21Sq" = distKnots21Sq, "x"=rep(0,length(Y)))
spatglm_pars = c("scaledf","yearEffects", "CV", "gp_sigmaSq", "gp_scale", "year_sigma","ar","spatialEffectsKnots")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# estimate model. This model is modified from the simulation model by (1) including indices to allow NAs in the inputted data, and (2) including estimated year effects (intercepts)
mvt_gamma = stan(file = 'stan_models/mvtGamma_estSigma_index_cov_yr_ar1.stan',data = spatglm_data,
  verbose = TRUE, chains = 3, thin = 1, warmup = 500, iter = 1000, pars = spatglm_pars)

spatglm_pars = c("yearEffects", "CV", "gp_sigmaSq", "gp_scale", "year_sigma","ar","spatialEffectsKnots")
mvn_gamma = stan(file = 'stan_models/mvnGamma_estSigma_index_cov_yr_ar1.stan',data = spatglm_data,
  verbose = TRUE, chains = 3, thin = 1, warmup = 500, iter = 1000, pars = spatglm_pars)

save.image("starling_bbs_mvtmvn.Rdata")

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

