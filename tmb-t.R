library(TMB)
compile("tmb_t.cpp")
dyn.load(dynlib("tmb_t"))

set.seed(1)
y <- rnorm(1000)

y <- rt(1000, 2)

dat <- list(y = y)

# parameters <- list(log_df = log(50), location = mean(y), log_scale = log(sd(y)),
  # w = rep(0, length(y)))

# parameters <- list(log_df = log(50), location = mean(y), w = rep(0, length(y)))

parameters <- list(log_df = log(50), w = rep(0, length(y)))

obj <- MakeADFun(
  data = dat,
  parameters = parameters,
  DLL = "tmb_t", random = "w")

opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr,
  lower = c(log(1), -10, -20), upper = c(log(100), 10, 20))

opt$par
exp(opt$par[1])

I(opt$par[2])
exp(opt$par[3])
mean(y)
sd(y)
