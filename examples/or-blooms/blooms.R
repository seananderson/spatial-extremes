# devtools::install_github("ropensci/rerddap")
# install.packages("rerddap")

library(dplyr)
library(ggplot2)

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

# shortbelly
# english sole, petrole

# filter(d, prc_chla_trans < 3, prc_chla_trans > -3) %>%
#   ggplot(aes(lon, lat, fill = prc_chla_trans)) +
#   geom_raster() +
#   scale_fill_gradient2() +
#   facet_wrap(~month)


# convert to UTMs
# sample grid - at ?
# fiddle with downsampling (500?), fiddle with knots (30-50)

# try downsampling of every nth lat/lon
dowsample=20
lons = unique(d$lon)[seq(1, length(unique(d$lon)), dowsample)]
lats = unique(d$lat)[seq(1, length(unique(d$lat)), dowsample)]
dsub = d[which(d$lat%in% lats & d$lon %in% lons),]

# convert to UTM
dsub$ID = seq(1,nrow(dsub))
sp::coordinates(dsub) = c("lon","lat")
sp::proj4string(dsub) <- sp::CRS("+proj=longlat +datum=WGS84")  ## for example
dsub = as.data.frame(sp::spTransform(dsub, sp::CRS(paste("+proj=utm +zone=10"," ellps=WGS84",sep=''))))

dsub = dsub[which(is.na(dsub$prc_chla_trans)==FALSE),]

# filter(dsub, prc_chla_trans < 3, prc_chla_trans > -3) %>%
#   ggplot(aes(lon, lat, fill = prc_chla_trans)) +
#   geom_raster() +
#   scale_fill_gradient2() +
#   facet_wrap(~month)

nKnots = 20L

dsub <- mutate(dsub, lat_scaled = lat / 100000, lon_scaled = lon / 100000,
  timeID = as.integer(as.factor(month)))

dsub %>%
  ggplot() + geom_point(data=dsub, aes(lon, lat, color = prc_chla_trans)) +
  facet_wrap(~month) + scale_color_gradient2()

# ------------------------------
# new model

library(rrfields)
library(rstan)
options(mc.cores = parallel::detectCores())
# mvt_norm <- rrfield(prc_chla_trans ~ 0, data = filter(dsub, timeID > 4),
#   time = "timeID", lon="lon_scaled", lat="lat_scaled",
#   nknots = nKnots, estimate_df = TRUE,
#   algorithm = "sampling",
#   chains = 3L, iter = 1000L)

mvt_norm <- rrfield(prc_chla_trans ~ 0, data = dsub,
  time = "timeID", lon="lon_scaled", lat="lat_scaled",
  nknots = nKnots, estimate_df = TRUE,
  algorithm = "sampling",
  chains = 4L, iter = 800L)

mvn_norm <- rrfield(prc_chla_trans ~ 0, data = dsub,
  time = "timeID", lon="lon_scaled", lat="lat_scaled",
  nknots = nKnots, estimate_df = FALSE, fixed_df_value = 1e6,
  algorithm = "sampling",
  chains = 4L, iter = 800L)

saveRDS(mvt_norm, file = "examples/OR_blooms/blooms-2016-12-14.rds")
mvt_norm <- readRDS("examples/OR_blooms/mvt-blooms-2016-12-14.rds")

mvt_norm

# plot(mvt_norm)

obs <- mvt_norm$data$prc_chla_trans
p <- predict(mvt_norm)
# pp <- predict(mvt_norm2, interval = "prediction")
plot(obs, p$estimate, col = "#00000030")
plot(p$estimate, obs - p$estimate, col = "#00000030");abline(h = 0)
cor(obs, p$estimate)
# segments(obs, pp$conf_low, obs, pp$conf_high, lwd = 0.5, col = "#00000020")
# segments(obs, p$conf_low, obs, p$conf_high, lwd = 0.5, col = "#00000030")
# abline(a = 0, b = 1)
# (coverage <- mean(obs > pp$conf_low & obs < pp$conf_high))

# ------------------
# old model
knots = cluster::pam(dsub[,c("lon_scaled","lat_scaled")],nKnots)$medoids
distKnots = as.matrix(dist(knots))
distKnotsSq = distKnots^2 # squared distances
# Calculate distance from knots to grid
distAll = as.matrix(dist(rbind(dsub[,c("lon_scaled","lat_scaled")], knots)))^2
# this is the transpose of the lower left corner
nLocs = nrow(dsub)
distKnots21Sq = t(distAll[-c(1:nLocs), -c((nLocs+1):ncol(distAll))])

Y = as.numeric(dsub$prc_chla_trans)
yearID = as.numeric(as.factor(dsub$month))

# create list for STAN
spatglm_data = list("nKnots"=nKnots, "nLocs"=nLocs, "nT" = length(unique(yearID)),
  "N" = length(Y), "yearID" = yearID, "y" = Y,
  "distKnotsSq" = distKnotsSq, "distKnots21Sq" = distKnots21Sq, "x" = rep(0, length(Y)))
spatglm_pars = c("scaledf","yearEffects", "sigma", "gp_sigmaSq", "gp_scale",
  "year_sigma","ar","spatialEffectsKnots")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# estimate model. This model is modified from the simulation model by (1) including indices to allow NAs in the inputted data, and (2) including estimated year effects (intercepts)
stanMod_mvt_norm = stan(file = 'stan_models/mvtNorm_estSigma_index_yr_ar1.stan',
  data = spatglm_data, chains = 4L, warmup = 750, iter = 1500, pars = spatglm_pars,
  control = list(adapt_delta = 0.99))

saveRDS(stanMod_mvt_norm,"examples/OR_blooms/stanMod_mvt_norm.rds")

spatglm_pars = c("yearEffects", "sigma", "gp_sigmaSq", "gp_scale",
  "year_sigma","ar","spatialEffectsKnots")
stanMod_mvn_norm = stan(file = 'stan_models/mvnNorm_estSigma_index_yr_ar1.stan',
  data = spatglm_data, chains = 4L, warmup = 750, iter = 1500, pars = spatglm_pars,
  control = list(adapt_delta = 0.99))

saveRDS(stanMod_mvn_norm,"examples/OR_blooms/stanMod_mvn_norm.rds")

# --------------------------------------
# check posterior / priors
stanMod_mvt_norm <- readRDS("examples/OR_blooms/stanMod_mvt_norm.rds")

e <- extract(stanMod_mvt_norm)
names(e)
dt2 <- function(x, df, mu, a) 1/a * dt((x - mu)/a, df)

plot_prior <- function(par_name, dens_fun = dnorm, xlim = c(0, 5),
  xtrans = I, ...) {
  hist(e[[par_name]], xlim = xlim, main = par_name, xlab = "")
  xx <- seq(xlim[1], xlim[2], length.out = 200)
  par(new = TRUE)
  yy <- dens_fun(xx, ...)
  plot(xtrans(xx), yy, type = "l", ylim = c(0, max(yy)), axes = FALSE,
    ylab = "", xlab = "", xlim = xlim)
}

pdf("examples/OR_blooms/or-blooms-priors.pdf", width = 8, height = 4)
par(mfrow = c(2, 3))
par(cex = 0.8, mar = c(3, 3, 1, 1))
plot_prior("scaledf", dgamma, shape = 2, scale = 1/0.1, xlim = c(1, 40))
plot_prior("gp_sigmaSq", dnorm, xlim = c(0, 5))
plot_prior("sigma", dnorm, xlim = c(0, 5))
plot_prior("ar", dnorm, xlim = c(-2, 2))
plot_prior("gp_scale", dnorm, xlim = c(0, 5))
plot_prior("year_sigma", dnorm, xlim = c(0, 5))
dev.off()

# alternatives:
pdf("examples/OR_blooms/or-blooms-priors-proposed.pdf", width = 8, height = 4)
par(mfrow = c(2, 3))
par(cex = 0.8, mar = c(3, 3, 1, 1))
a <- 2 # scale
df <- 3
plot_prior("scaledf", dgamma, shape = 2, scale = 1/0.1, xlim = c(1, 40))
plot_prior("gp_sigmaSq", dt2, df = df, mu = 0, a = a, xlim = c(0, 4))
plot_prior("sigma", dt2, df = df, mu = 0, a = a, xlim = c(0, 4))
plot_prior("ar", dbeta, shape1 = 2, shape2 = 2, xlim = c(-1.5, 1.5),
  xtrans = function(x) 2 * (x - 1/2))
plot_prior("gp_scale", dt2, df = df, mu = 0, a = 20, xlim = c(0, 40))
plot_prior("year_sigma", dt2, df = df, mu = 0, a = a, xlim = c(0, 4))
dev.off()
