# This script just does a simple speed test of matrix multiplication in STAN -- is it better
# to multiply a giant (1000 x 30 matrix) or 10 smaller 100x30 matrices?

set.seed(1)

y = rnorm(10, 0, sd=3)

#y = rnorm(N, 0, 1)
nKnots = 50
distKnotsSq = matrix(runif(nKnots*nKnots), nKnots, nKnots)
nLocs = 100
nYears = 10
distKnots21Sq = matrix(runif(nLocs*nYears*nKnots), nLocs*nYears, nKnots)

spatglm_data = list("nKnots"=nKnots, "nLocs"=nLocs, "y"=y, "distKnotsSq"=distKnotsSq, "distKnots21Sq"=distKnots21Sq)
spatglm_pars = c("sigma")

# estimate model. This model is modified from the simulation model by (1) including indices to allow NAs in the inputted data, and (2) including estimated year effects (intercepts)
ptm <- proc.time()
for(i in 1:10) model1 = stan(file = 'stan_models/matrix_test_full.stan',data = spatglm_data,
  verbose = TRUE, chains = 3, thin = 1, warmup = 1000, iter = 2000, pars = spatglm_pars)
proc.time() - ptm

#user  system elapsed
#4.318   1.387 128.415

distKnots21Sq1 = matrix(runif(nLocs*nKnots), nLocs, nKnots)
distKnots21Sq2 = matrix(runif(nLocs*nKnots), nLocs, nKnots)
distKnots21Sq3 = matrix(runif(nLocs*nKnots), nLocs, nKnots)
distKnots21Sq4 = matrix(runif(nLocs*nKnots), nLocs, nKnots)
distKnots21Sq5 = matrix(runif(nLocs*nKnots), nLocs, nKnots)
distKnots21Sq6 = matrix(runif(nLocs*nKnots), nLocs, nKnots)
distKnots21Sq7 = matrix(runif(nLocs*nKnots), nLocs, nKnots)
distKnots21Sq8 = matrix(runif(nLocs*nKnots), nLocs, nKnots)
distKnots21Sq9 = matrix(runif(nLocs*nKnots), nLocs, nKnots)
distKnots21Sq10 = matrix(runif(nLocs*nKnots), nLocs, nKnots)
spatglm_data = list("nKnots"=nKnots, "nLocs"=nLocs, "y"=y, "distKnotsSq"=distKnotsSq, "distKnots21Sq1"=distKnots21Sq1,
  "distKnots21Sq2"=distKnots21Sq2, "distKnots21Sq3"=distKnots21Sq3, "distKnots21Sq4"=distKnots21Sq4, "distKnots21Sq5"=distKnots21Sq5,
  "distKnots21Sq6"=distKnots21Sq6, "distKnots21Sq7"=distKnots21Sq7, "distKnots21Sq8"=distKnots21Sq8,
  "distKnots21Sq9"=distKnots21Sq9, "distKnots21Sq10"=distKnots21Sq10)
spatglm_pars = c("sigma")

ptm2 <- proc.time()
for(i in 1:10) model2 = stan(file = 'stan_models/matrix_test_split.stan',data = spatglm_data,
  verbose = TRUE, chains = 3, thin = 1, warmup = 1000, iter = 2000, pars = spatglm_pars)
proc.time() - ptm2
# user  system elapsed
# 5.055   1.556 127.371



###########################################
# Test 2: how much does adding more data affect speed of computations?

library(rstan)

locs = seq(20, 200, 20)
sec = 0
for(i in 1:length(locs)) {
  nLocs = locs[i]
  #y = rnorm(N, 0, 1)
  nKnots = 50
  distKnotsSq = matrix(runif(nKnots*nKnots), nKnots, nKnots)
  nYears = 10
  distKnots21Sq = matrix(runif(nLocs*nYears*nKnots), nLocs*nYears, nKnots)

  spatglm_data = list("nKnots"=nKnots, "nLocs"=nLocs, "y"=y, "distKnotsSq"=distKnotsSq, "distKnots21Sq"=distKnots21Sq)
  spatglm_pars = c("sigma")

  ptm2 <- proc.time()
  model2 = stan(file = 'stan_models/matrix_test_full.stan',data = spatglm_data,
    verbose = TRUE, chains = 3, thin = 1, warmup = 1000, iter = 2000, pars = spatglm_pars)
  proc.time() - ptm2
  sec[i] = (proc.time() - ptm2)[3]
}

pdf("simulationTesting/cost_of_more_data.pdf")
plot(locs, sec, xlab="Observations", ylab="Seconds", type="b",lwd=3, main="Computational cost of adding more data")
dev.off()



