library(tidyverse)
source("simulationTesting/sim_mvt_rf.R")
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4L)

nDataPoints <- 100
g <- data.frame(lon = runif(nDataPoints, 0, 10),
  lat = runif(nDataPoints, 0, 10))
nLocs <- dim(g)[1]

#g <- expand.grid(lon = seq(1, 10, length.out = 10),
#  lat = seq(1, 10, length.out = 10))
#nLocs <- nrow(g)

draws <- 10
simulation_data <- sim_mvt_rf(df = 2, grid = g, n_pts = nrow(g), seed = NULL,
      n_draws = draws, n_knots = 30, gp_scale = 0.5, sigma_t = 0.5)

out <- reshape2::melt(simulation_data$proj)
names(out) <- c("i", "pt", "re")
out <- arrange(out, i, pt)
out$lon <- rep(g$lon, draws)
out$lat <- rep(g$lat, draws)

p1 <- out %>%
  ggplot(aes(x = lon, y = lat, z = re, colour = re)) +
  facet_wrap(~i, ncol = 5) +
  geom_point(size = 3) +
  viridis::scale_color_viridis() +
  theme_light()

CV <- 0.2
gamma.a <- 1/(CV^2)
gamma.b <- gamma.a/exp(out$re)
out$obs <- rgamma(nrow(out), shape = gamma.a, rate = c(gamma.b))

# out <- mutate(out, obs = re + rnorm(nrow(out), sd = 0.5))
p2 <- out %>%
  ggplot(aes(x = lon, y = lat, z = log(obs), colour = log(obs))) +
  facet_wrap(~i, ncol = 5) +
  geom_point(size = 3) +
  viridis::scale_color_viridis() +
  theme_light()
gridExtra::grid.arrange(p1, p2)

gamma.a <- 1/(CV^2)
gamma.b <- gamma.a/exp(simulation_data$proj)
y.gamma <- rgamma(nrow(gamma.b)*ncol(gamma.b), shape = gamma.a, rate = c(gamma.b))
y <- matrix(y.gamma, ncol = ncol(gamma.b))

# y <- simulation_data$proj + matrix(
#     rnorm(ncol(simulation_data$proj) * nrow(simulation_data$proj), sd = 0.1),
#     nrow(simulation_data$proj), ncol(simulation_data$proj))

model_data <- list(nKnots = nrow(simulation_data$knots), nLocs = nLocs,
  nT = nrow(simulation_data$re_knots), y = y,
  distKnotsSq = simulation_data$dist_knots_sq,
  distKnots21Sq = simulation_data$dist_knots21_sq)

pars <- c("scaledf", "gp_scale", "gp_sigmaSq", "CV")
m <- stan(file = 'stan_models/mvtGamma_estSigma.stan',
  data = model_data, chains = 4L, warmup = 200L, iter = 400L, pars = pars)
m
# ----------
# now try across multiple arguments


