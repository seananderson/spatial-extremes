library(tidyverse)
source("simulationTesting/sim_mvt_rf.R")
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4L)

nDataPoints <- 100
g <- data.frame(lon = runif(nDataPoints, 0, 10),
  lat = runif(nDataPoints, 0, 10))
nLocs <- dim(g)[1]

# draws <- 10
# simulation_data <- sim_mvt_rf(df = 2, grid = g, n_pts = nrow(g), seed = NULL,
#       n_draws = draws, n_knots = 30, gp_scale = 0.5, sigma_t = 0.5)
#
# out <- reshape2::melt(simulation_data$proj)
# names(out) <- c("i", "pt", "re")
# out <- arrange(out, i, pt)
# out$lon <- rep(g$lon, draws)
# out$lat <- rep(g$lat, draws)
#
# p1 <- out %>%
#   ggplot(aes(x = lon, y = lat, z = re, colour = re)) +
#   facet_wrap(~i, ncol = 5) +
#   geom_point(size = 3) +
#   viridis::scale_color_viridis() +
#   theme_light()
#
# CV <- 0.2
# gamma.a <- 1/(CV^2)
# gamma.b <- gamma.a/exp(out$re)
# out$obs <- rgamma(nrow(out), shape = gamma.a, rate = c(gamma.b))
#
# # out <- mutate(out, obs = re + rnorm(nrow(out), sd = 0.5))
# p2 <- out %>%
#   ggplot(aes(x = lon, y = lat, z = log(obs), colour = log(obs))) +
#   facet_wrap(~i, ncol = 5) +
#   geom_point(size = 3) +
#   viridis::scale_color_viridis() +
#   theme_light()
# gridExtra::grid.arrange(p1, p2)
#
# gamma.a <- 1/(CV^2)
# gamma.b <- gamma.a/exp(simulation_data$proj)
# y.gamma <- rgamma(nrow(gamma.b)*ncol(gamma.b), shape = gamma.a, rate = c(gamma.b))
# y <- matrix(y.gamma, ncol = ncol(gamma.b))
#
# # y <- simulation_data$proj + matrix(
# #     rnorm(ncol(simulation_data$proj) * nrow(simulation_data$proj), sd = 0.1),
# #     nrow(simulation_data$proj), ncol(simulation_data$proj))
#
# model_data <- list(nKnots = nrow(simulation_data$knots), nLocs = nLocs,
#   nT = nrow(simulation_data$re_knots), y = y,
#   distKnotsSq = simulation_data$dist_knots_sq,
#   distKnots21Sq = simulation_data$dist_knots21_sq)
#
# pars <- c("scaledf", "gp_scale", "gp_sigmaSq", "CV")
# m <- stan(file = 'stan_models/mvtGamma_estSigma.stan',
#   data = model_data, chains = 4L, warmup = 200L, iter = 400L, pars = pars)
# m

# ----------
# now try across multiple arguments

stan_mod <- rstan::stan_model("stan_models/mvtGamma_estSigma.stan")

sim_fit <- function(df = 2, n_draws, n_knots = 30, gp_scale = 0.5, cv_obs = 0.2,
  comment = "", sigma_t = 0.5) {

  nDataPoints <- 100
  g <- data.frame(lon = runif(nDataPoints, 0, 10),
    lat = runif(nDataPoints, 0, 10))
  n_pts <- nrow(g)

  simulation_data <- sim_mvt_rf(df = df, grid = g, n_pts = n_pts, seed = NULL,
    n_draws = n_draws, n_knots = n_knots, gp_scale = gp_scale, sigma_t = sigma_t)

  gamma.a <- 1/(cv_obs^2)
  gamma.b <- gamma.a/exp(simulation_data$proj)
  y.gamma <- rgamma(nrow(gamma.b)*ncol(gamma.b), shape = gamma.a, rate = c(gamma.b))
  y <- matrix(y.gamma, ncol = ncol(gamma.b))
  model_data <- list(nKnots = nrow(simulation_data$knots), nLocs = n_pts,
    nT = nrow(simulation_data$re_knots), y = y,
    distKnotsSq = simulation_data$dist_knots_sq,
    distKnots21Sq = simulation_data$dist_knots21_sq)

  pars <- c("scaledf", "gp_scale", "gp_sigmaSq", "CV")
  m <-  rstan::sampling(stan_mod, data = model_data, chains = 4L, warmup = 200L,
    iter = 400L, pars = pars)
  m
}

set.seed(1)
arguments <- readxl::read_excel("simulationTesting/simulation-arguments.xlsx")
# arguments$case <- c(0, rep(c(-1, 1), 5))
arguments$count <- 20L
arguments <- arguments[rep(seq_len(nrow(arguments)), arguments$count), ]
arguments_apply <- dplyr::select(arguments, -count, -case)
nrow(arguments)

out <- plyr::mlply(arguments_apply, sim_fit)

saveRDS(out, file = "simulationTesting/gamma-sim-testing.rds")

out <- readRDS("simulationTesting/gamma-sim-testing.rds")
out_print <- out %>%
  map_df(function(x) {
    broom::tidyMCMC(x, estimate.method = "median") %>%
      dplyr::select(-std.error) %>%
      tidyr::spread(term, estimate)})
names(out_print) <- paste0(names(out_print), "_est")


rhat <- out %>%
  map_df(function(x) {
    broom::tidyMCMC(x, estimate.method = "median", rhat = TRUE, ess = TRUE) %>%
      summarise(rhat = max(rhat), ess = min(ess))
    })


out_summary <- data.frame(arguments, out_print, rhat) %>%
  filter(rhat < 1.05, ess > 100) %>%
  select(-count, -rhat, -ess) %>%
  mutate(sigma_t = 0.5^2)

plot_viol <- function(term, term_true) {
  x <- tidyr::gather(out_summary, parameter, estimate, -df, -n_knots, -n_draws,
    -gp_scale, -cv_obs, -case, -comment, -sigma_t) %>%
    filter(parameter == term)

  ggplot(x, aes(paste(comment, case), estimate)) +
    geom_violin(draw_quantiles = c(0.5), trim = TRUE, fill = "grey93") +
    coord_flip() +
    geom_point(aes_string(y = term_true), colour = "red", size = 2) +
    theme_light() +
    labs(title = term_true, x = "")
}

p1 <- plot_viol("scaledf_est", "df")
p2 <- plot_viol("gp_scale_est", "gp_scale")
p3 <- plot_viol("gp_sigmaSq_est", "sigma_t")
p4 <- plot_viol("CV_est", "cv_obs")

pdf("simulationTesting/sim-gamma-pars.pdf", width = 7, height = 6)
gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()
