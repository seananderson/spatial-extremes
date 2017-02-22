library(tidyverse)
library(rrfields)
library(rstan)
options(mc.cores = min(c(4L, parallel::detectCores())))
library(ggsidekick) # devtools::install_github("seananderson/ggsidekick")

# ------------------------------------------------------------
# Pick reasonable values:
if (interactive()) {
  library(manipulate)
  manipulate({
    set.seed(seed)
    simulation_data <- sim_rrfield(df = df, n_data_points = 100, seed = NULL,
      n_draws = 6, n_knots = 7, gp_scale = gp_scale, gp_sigma = gp_sigma,
      obs_error = "normal", sd_obs = CV)
    print(simulation_data$plot)
  }, gp_scale = slider(0.05, 10, 1, step = 0.25),
    gp_sigma = slider(0.05, 10, 0.5, step = 0.25),
    df = slider(2, 50, 4, step = 1),
    CV = slider(0.01, 1, 0.05, step = 0.05),
    seed = slider(1, 300, 10, step = 1))
}

# ------------------------------------------------------------
# Now run across multiple arguments

sim_fit <- function(n_draws, df = 2, n_knots = 30, gp_scale = 0.5, sd_obs = 0.2,
  comment = "", gp_sigma = 0.5) {

  s <- sim_rrfield(df = df, n_data_points = 100, seed = NULL,
    n_draws = n_draws, n_knots = n_knots, gp_scale = gp_scale, gp_sigma = gp_sigma,
    sd_obs = sd_obs, obs_error = "normal")

  fit_model <- function(iter) {
    rrfield(y ~ 1, data = s$dat, time = "time", lon = "lon", lat = "lat",
      nknots = n_knots, chains = 4L, iter = iter, obs_error = "normal",
      prior_gp_scale = half_t(3, 0, 5),
      prior_gp_sigma = half_t(3, 0, 3),
      prior_sigma = half_t(3, 0, 3),
      prior_intercept = student_t(1e6, 0, 1),
      prior_beta = student_t(1e6, 0, 1))
  }

  m <- fit_model(iter = 500L)
  b <- broom::tidyMCMC(m$model, rhat = TRUE, ess = TRUE)
  if (any(b$ess < 100) | any(b$rhat > 1.05)) {
    m <- fit_model(iter = 1000L)
  }
  m
}

set.seed(123)
arguments <- readxl::read_excel("simulationTesting/simulation-arguments.xlsx")
arguments$count <- 5L
arguments <- arguments[rep(seq_len(nrow(arguments)), arguments$count), ]
arguments_apply <- dplyr::select(arguments, -count, -case)
nrow(arguments)

out <- plyr::mlply(arguments_apply[1,], sim_fit)
saveRDS(out, file = "simulationTesting/mvt-norm-sim-testing.rds")

# TODO pick only 4 params:
out <- readRDS("simulationTesting/mvt-norm-sim-testing.rds")
out_print <- out %>%
  map_df(function(x) {
    broom::tidyMCMC(x, estimate.method = "median") %>%
      dplyr::select(-std.error) %>%
      tidyr::spread(term, estimate)})
names(out_print) <- paste0(names(out_print), "_est")

rhat <- out %>%
  map_df(function(x) {
    broom::tidyMCMC(x, rhat = TRUE, ess = TRUE) %>%
      summarise(rhat = max(rhat), ess = min(ess))
    })

out_summary <- data.frame(arguments, out_print, rhat) %>%
  filter(rhat < 1.05, ess > 100) %>%
  select(-count, -rhat, -ess) %>%
  mutate(gp_sigmaSq = 0.5^2)

plot_viol <- function(term, term_true) {
  x <- tidyr::gather(out_summary, parameter, estimate, -df, -n_knots, -n_draws,
    -gp_scale, -cv_obs, -case, -comment, -gp_sigmaSq) %>%
    filter(parameter == term)

  ggplot(x, aes(paste(comment, case), estimate)) +
    geom_violin(draw_quantiles = c(0.5), trim = TRUE, fill = "grey93") +
    coord_flip() +
    geom_point(aes_string(y = term_true), colour = "red", size = 2) +
    theme_sleek() +
    labs(title = term_true, x = "")
}

p1 <- plot_viol("df_est", "df")
p2 <- plot_viol("gp_scale_est", "gp_scale")
p3 <- plot_viol("gp_sigma_est", "gp_sigmaSq")
p4 <- plot_viol("sigma_est", "sd_obs")

pdf("simulationTesting/sim-mvt-norm-pars.pdf", width = 7, height = 6)
gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()
