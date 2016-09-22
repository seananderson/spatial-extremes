set.seed(1)
N = 1000
y = rt(N, df=2)
#y = rnorm(N, 0, 1)
spatglm_data = list("N"=N, "y"=y)
spatglm_pars = c("scaledf", "W")

# estimate model. This model is modified from the simulation model by (1) including indices to allow NAs in the inputted data, and (2) including estimated year effects (intercepts)
univ = stan(file = 'stan_models/univariate_trans.stan',data = spatglm_data,
  verbose = TRUE, chains = 3, thin = 1, warmup = 1000, iter = 2000, pars = spatglm_pars)

extract(univ, permuted = TRUE)[["scaledf"]]

apply(extract(univ, permuted = TRUE)[["W"]], 2, median)
