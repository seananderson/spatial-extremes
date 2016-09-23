# This script just does a simple speed test of matrix multiplication in STAN -- is it better
# to multiply a giant (1000 x 30 matrix) or 10 smaller 100x30 matrices?


nDataPoints = 100
grid = cbind("lon" = runif(nDataPoints, 5, 15),
  "lat" = runif(nDataPoints, 5, 15))
nLocs = dim(grid)[1]
set.seed(1)
dataSim <-simulateData(nKnots = 30, nDraws = 20, gp_scale = 0.5, mvt = TRUE, df = 1)

allknots = seq(5, 35, 5)
secs = 0
for(i in 1:length(allknots)) {

  ptm = proc.time()
  # generate data with observation error
  dataSim <-simulateData(nKnots = allknots[i], nDraws = 20, gp_scale = 0.5, mvt = TRUE, df = 1)
  Y = dataSim$proj + matrix(rnorm(ncol(dataSim$proj) * nrow(dataSim$proj), sd = 0.1), nrow(dataSim$proj), ncol(dataSim$proj))

  # create list for STAN
  spatglm_data = list("nKnots"=nrow(dataSim$knots), "nLocs"=nLocs, "nT" = nrow(dataSim$reKnots), "y" = Y, "distKnotsSq" = dataSim$distKnotsSq, "distKnots21Sq" = dataSim$distKnots21Sq)
  spatglm_pars = c("scaledf")
  # estimate model
  stanMod = stan(file = 'stan_models/mvtNormal_estSigma.stan',data = spatglm_data,
    verbose = TRUE, chains = 1, thin = 1, warmup = 300, iter = 600,
    pars = spatglm_pars)
  secs[i] = (proc.time() - ptm)[3]
}

