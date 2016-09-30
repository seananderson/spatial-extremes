# devtools::install_github("ropensci/rerddap")
# install.packages("rerddap")

library(dplyr)

rerddap::ed_search(query = "temp") %>% .$info

(out <- rerddap::info("osuBloomsModisChla"))

d <- list()
m <- c("02-01", "02-15", "03-01", "03-15",
  "04-01", "04-15", "05-01", "05-15",
  "06-01", "06-15", "07-01", "07-15",
  "08-01", "08-15", "09-01", "09-15")
t1 <- paste0("2014-", m, "T00:00:00Z")
t2 <- paste0("2014-", m, "T00:00:01Z")

for (i in seq_along(m)) {
  message(m[[i]])
  d[[i]] <- rerddap::griddap(out,
    time = c(t1[[i]], t2[[i]]),
    longitude = c(-126, -124),
    latitude = c(42, 48.5))$data
  d[[i]]$month <- m[[i]]
  d[[i]] <- as_data_frame(d[[i]])
}
d <- bind_rows(d)

x <- d$prc_chla+abs(min(d$prc_chla, na.rm = T))
x[x == 0] <- NA
hist(log(x))
d$prc_chla_trans <- log(x)

d <- d %>% group_by(month) %>%
  mutate(prc_chla_trans = prc_chla_trans - mean(prc_chla_trans, na.rm = TRUE)) %>%
  ungroup()

library(ggplot2)
filter(d, prc_chla_trans < 3, prc_chla_trans > -3) %>%
  ggplot(aes(lon, lat, fill = prc_chla_trans)) +
  geom_raster() +
  scale_fill_gradient2() +
  facet_wrap(~month)


# convert to UTMs
# sample grid - at ?
# fiddle with downsampling (500?), fiddle with knots (30-50)

# try downsampling of every 20th lat/lon
dowsample=20
lons = unique(d$lon)[seq(1, length(unique(d$lon)), dowsample)]
lats = unique(d$lat)[seq(1, length(unique(d$lat)), dowsample)]
dsub = d[which(d$lat%in% lats & d$lon %in% lons),]

dsub = dsub[which(is.na(dsub$prc_chla_trans)==FALSE),]

nKnots = 75
knots = cluster::pam(dsub[,c("lon","lat")],nKnots)$medoids
filter(dsub, prc_chla_trans < 3, prc_chla_trans > -3) %>%
  ggplot() + geom_point(data=dsub, aes(lon, lat, color = prc_chla_trans)) +
  facet_wrap(~month)
distKnots = as.matrix(dist(knots))
distKnotsSq = distKnots^2 # squared distances

# Calculate distance from knots to grid
distAll = as.matrix(dist(rbind(dsub[,c("lon","lat")], knots)))^2
# this is the transpose of the lower left corner
nLocs = nrow(dsub)
distKnots21Sq = t(distAll[-c(1:nLocs), -c((nLocs+1):ncol(distAll))])

Y = as.numeric(dsub$prc_chla_trans)
yearID = as.numeric(as.factor(dsub$month))
stationID = seq(1,nrow(dsub))

# create list for STAN
spatglm_data = list("nKnots"=nKnots, "nLocs"=nLocs, "nT" = length(unique(yearID)), "N" = length(Y), "stationID" = stationID, "yearID" = yearID, "y" = Y, "distKnotsSq" = distKnotsSq, "distKnots21Sq" = distKnots21Sq, "x" = rep(0, length(Y)))
spatglm_pars = c("scaledf","yearEffects", "CV", "gp_sigmaSq", "gp_scale", "year_sigma","ar","B","spatialEffectsKnots")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# estimate model. This model is modified from the simulation model by (1) including indices to allow NAs in the inputted data, and (2) including estimated year effects (intercepts)
stanMod_gamma = stan(file = 'stan_models/mvtNorm_estSigma_index_yr_ar1.stan',data = spatglm_data,
  verbose = TRUE, chains = 3, thin = 1, warmup = 10, iter = 20, pars = spatglm_pars)
