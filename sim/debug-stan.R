library(cluster)
library(mvtnorm)
library(fields)
library(ggplot2)
library(MASS)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# simulate
nDataPoints = 31
grid = cbind("lon" = runif(nDataPoints, 5, 15),
  "lat" = runif(nDataPoints, 5, 15))
nLocs = dim(grid)[1]

simulateData = function(nKnots = 10, gp_scale = 0.3, sigma_t = 0.15,
  nDraws = 1000, mvt = TRUE, df = 3) {

  # cluster analysis to determine knot locations
  knots = jitter(pam(grid,nKnots)$medoids)
  distKnots = as.matrix(dist(knots))
  distKnotsSq = distKnots^2 # squared distances

  corKnots = exp(-gp_scale*distKnotsSq)
  sigmaKnots = sigma_t * sigma_t * corKnots
  invSigmaKnots = solve(sigmaKnots)
  # Calculate distance from knots to grid
  distAll = as.matrix(dist(rbind(grid, knots)))^2
  # this is the transpose of the lower left corner
  distKnots21Sq = t(distAll[-c(1:nLocs), -c((nLocs+1):ncol(distAll))])
  Sigma21 = exp(-gp_scale*distKnots21Sq) * sigma_t * sigma_t

  # Generate vector of random effects
  # each 'draw' here is hypothetical draw from posterior
  if(mvt==TRUE) reKnots = mvtnorm::rmvt(nDraws, sigma = sigmaKnots, df = df)
  if(mvt==FALSE) reKnots = mvtnorm::rmvnorm(nDraws, sigma = sigmaKnots)

  # if(mvt==TRUE) {
  #   MVN <- MASS::mvrnorm(nDraws, rep(0, ncol(sigmaKnots)), sigmaKnots)
  #   reKnots <- t(t(MVN / sqrt(df / stats::rchisq(n, df))))
  # } else {
  #   reKnots = MASS::rmvnorm(nDraws, rep(0, ncol(sigmaKnots)), sigmaKnots)
  # }

  # Project random effects to locations of the data
  proj = t((Sigma21 %*% invSigmaKnots) %*% t(reKnots))
  return(list(knots = knots, reKnots = reKnots, proj = proj, distKnotsSq = distKnotsSq, distKnots21Sq = distKnots21Sq, sigmaKnots = sigmaKnots))
}

set.seed(1)
dataSim = simulateData(nKnots = 30, nDraws = 10, gp_scale = 0.5, mvt=TRUE)

if(dim(dataSim$reKnots)[1]==1) {
  dataSim$reKnots = t(dataSim$reKnots)
  dataSim$proj = t(dataSim$proj)
}

k <- data.frame(dataSim$knots, v = t(dataSim$reKnots)[,1], type = "knots")
p <- data.frame(grid, v = t(dataSim$proj)[,1], type = "projection")
d <- rbind(k, p)
ggplot(d, aes(lon, lat, colour = v)) + facet_wrap(~type) +
  geom_point(size=3) + scale_colour_gradient2()


# fit the model
spatglm_data = list(
  "nKnots"=nrow(dataSim$knots),
  "nLocs"=nLocs,
  "nT" = nrow(dataSim$reKnots),
  "y" = dataSim$proj,
  "distKnotsSq" = dataSim$distKnotsSq,
  "distKnots21Sq" = dataSim$distKnots21Sq)

spatglm_pars = c("spatialEffectsKnots", "gp_scale", "gp_sigmaSq")
stanMod_mvnmvt_1 = stan(file = 'simulationTesting/mvnNormal.stan',
  data = spatglm_data, iter = 200, chains = 4,
  pars = spatglm_pars)

stanMod_mvnmvt_1

broom::tidyMCMC(stanMod_mvnmvt_1, rhat = TRUE, ess = TRUE)

library(dplyr)
broom::tidyMCMC(stanMod_mvnmvt_1, rhat = TRUE, ess = TRUE) %>%
  dplyr::filter(!grepl("Knots", term))

# jags

spatglm_data = list(
  "nKnots"=nrow(dataSim$knots),
  "nLocs"=nLocs,
  "nT" = nrow(dataSim$reKnots),
  "y" = dataSim$proj,
  "distKnotsSq" = dataSim$distKnotsSq,
  "distKnots21Sq" = dataSim$distKnots21Sq,
  "muZeros" = rep(0, nrow(dataSim$knots)),
  "diagKnots" = diag(nrow(dataSim$knots)))

estModel = R2jags::jags(data = spatglm_data, inits=NULL,
  parameters.to.save=c("gp_scale"),
  model.file="simulationTesting/archive_jags/recover_rf_mvnGaussian.txt",
 n.chains=2, n.iter=2000, n.burnin=1000,
 n.thin=1)

estModel
