library(tidyverse)
library(rrfields)
library(rstan)
options(mc.cores = min(c(4L, parallel::detectCores())))
library(ggsidekick) # devtools::install_github("seananderson/ggsidekick")
library(assertthat)

# ------------------------------------------------------------
# Pick reasonable values:
if (interactive()) {
  library(manipulate)
  manipulate({
    set.seed(seed)
    simulation_data <- sim_rrfield(df = df, n_data_points = 50, seed = NULL,
      n_draws = 6, n_knots = 7, gp_scale = gp_scale, gp_sigma = gp_sigma,
      obs_error = "gamma", sd_obs = CV)
    print(simulation_data$plot)
  }, gp_scale = slider(0.05, 10, 1, step = 0.25),
    gp_sigma = slider(0.05, 10, 0.5, step = 0.25),
    df = slider(2, 50, 4, step = 1),
    CV = slider(0.01, 1, 0.05, step = 0.05),
    seed = slider(1, 300, 10, step = 1))
}

# ------------------------------------------------------------
# Now run across multiple arguments

i <<- 0

sim_fit <- function(n_draws, df = 2, n_knots = 30, gp_scale = 0.5, sd_obs = 0.2,
  comment = "", gp_sigma = 0.5) {

  i <<- i + 1
  message(i)

  s <- sim_rrfield(df = df, n_data_points = 100, seed = NULL,
    n_draws = n_draws, n_knots = n_knots, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sd_obs, obs_error = "gamma")

  fit_model <- function(iter) {
    rrfield(y ~ 0, data = s$dat, time = "time", lon = "lon", lat = "lat",
      nknots = n_knots,
      station = "station_id",
      chains = 4L, iter = iter, obs_error = "gamma",
      prior_gp_scale = half_t(3, 0, 3),
      prior_gp_sigma = half_t(3, 0, 3),
      prior_sigma = half_t(3, 0, 3),
      prior_intercept = student_t(1e6, 0, 1),
      prior_beta = student_t(1e6, 0, 1))
  }

  m <- fit_model(iter = 400L)
  b <- broom::tidyMCMC(m$model, rhat = TRUE, ess = TRUE)
  if (any(b$ess < 200) | any(b$rhat > 1.05)) {
    m <- fit_model(iter = 1200L)
  }
  m
}

set.seed(123)
# arguments <- readxl::read_excel("simulationTesting/simulation-arguments.xlsx")
# arguments$count <- 5L
# arguments <- arguments[rep(seq_len(nrow(arguments)), arguments$count), ]
# arguments_apply <- dplyr::select(arguments, -count, -case)
# nrow(arguments)

arguments <- expand.grid(
  df = c(2, 5, 20),
  n_knots = 15,
  n_draws = c(10, 20, 30),
  gp_scale = 1,
  gp_sigma = 1,
  sd_obs = c(0.1, 0.5, 0.9, 1.3)
)
nrow(arguments)
arguments$count <- 50L
arguments <-
  arguments[rep(seq_len(nrow(arguments)), arguments$count),]
arguments_apply <- dplyr::select(arguments,-count)
nrow(arguments_apply)

arguments <- expand.grid(
  df = 3,
  n_knots = 15,
  n_draws = 20,
  gp_scale = 1,
  gp_sigma = 1,
  sd_obs = 1.3)

out <- plyr::mlply(arguments_apply, sim_fit)
saveRDS(out, file = "simulationTesting/mvt-norm-sim-testing.rds")

out <- readRDS("simulationTesting/mvt-norm-sim-testing.rds")
out_print <- out %>%
  map_df(function(x) {
    broom::tidyMCMC(x$model, estimate.method = "median") %>%
      dplyr::select(-std.error) %>%
      tidyr::spread(term, estimate)})
names(out_print) <- paste0(names(out_print), "_est")
out_print <- select(out_print, -starts_with("spatialEffectsKnots"))
names(out_print) <- sub("\\[1\\]", "", names(out_print))
saveRDS(out_print, file = "simulationTesting/mvt-norm-sim-testing-summary.rds")

rhat <- out %>%
  map_df(function(x) {
    broom::tidyMCMC(x$model, rhat = TRUE, ess = TRUE) %>%
      summarise(rhat = max(rhat), ess = min(ess))
    })
assert_that(max(rhat$rhat) < 1.05)
assert_that(max(rhat$ess) > 200)

out_summary <- data.frame(arguments, out_print, rhat) %>%
  filter(rhat < 1.05, ess > 200) %>%
  select(-count, -rhat, -ess)# %>%
#  mutate(gp_sigmaSq = 0.5^2)

ggplot(out_summary, aes(sd_obs, df_est, group = sd_obs, fill = as.factor(df))) +
  facet_grid(df~n_draws) + theme_sleek() +
  geom_violin(colour = NA) +
  geom_jitter(colour = "#00000020", height = 0, width = 0.02) +
  geom_hline(aes(yintercept = df), colour = "grey50") +
  scale_fill_brewer(palette = "YlOrRd", direction = -1)
ggsave("figs/sim-recapture.pdf", width = 10, height = 7)

# plot_viol <- function(term, term_true) {
#   x <- tidyr::gather(out_summary, parameter, estimate, -df, -n_knots, -n_draws,
#     -gp_scale, -sd_obs, -gp_sigma) %>%
#     filter(parameter == term)
#
#   ggplot(x, aes(paste(comment, case), estimate)) +
#     geom_violin(draw_quantiles = c(0.5), trim = TRUE, fill = "grey93") +
#     coord_flip() +
#     geom_point(aes_string(y = term_true), colour = "red", size = 2) +
#     theme_sleek() +
#     labs(title = term_true, x = "")
# }
#
# p1 <- plot_viol("df_est", "df")
# p2 <- plot_viol("gp_scale_est", "gp_scale")
# p3 <- plot_viol("gp_sigma_est", "gp_sigma")
# p4 <- plot_viol("sigma_est", "sd_obs")
#
# pdf("simulationTesting/sim-mvt-norm-pars.pdf", width = 7, height = 6)
# gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
# dev.off()
