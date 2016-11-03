library(dplyr)
library(ggplot2)
source("simulationTesting/sim_mvt_rf.R")
source("simulationTesting/theme_gg.R")
g <- expand.grid(lon = seq(1, 10, 0.5),
  lat = seq(1, 10, 0.5))

exp_cor <- function(delta_ij, phi) {
  exp(-phi * delta_ij)
}

library(manipulate)
manipulate({
  d <- seq(0, max, length.out = 200)
  r <- diff(range(d))
  scaled_phi <- phi_raw / r
  par(mfrow = c(2, 1), mar = c(3, 3, 1, 1))
  plot(d, exp_cor(d, phi = scaled_phi), type = "l",
    main = paste("scaled phi =", round(scaled_phi, 2)))
  x <- seq(0, (sd_prior_raw / r) * 3, length.out = 200)
  plot(x, dnorm(x, sd = sd_prior_raw / r), type = "l")
  abline(v = scaled_phi)
}, phi_raw = slider(1, 20, 5),
  max = slider(1, 20, 10),
  sd_prior_raw = slider(1, 50, 10))

gauss_cor <- function(delta_ij, phi) {
  exp(-phi * (delta_ij^2))
}

manipulate({
  d <- seq(0, max, length.out = 200)
  r <- diff(range(d))
  scaled_phi <- phi_raw / (r^2)
  par(mfrow = c(2, 1), mar = c(3, 3, 1, 1))
  plot(d, gauss_cor(d, phi = scaled_phi), type = "l",
    main = paste("scaled phi =", round(scaled_phi, 2)))
  x <- seq(0, (sd_prior_raw / (r^2)) * 3, length.out = 200)
  # plot(x, dnorm(x, sd = sd_prior_raw / (r^2)), type = "l")
  plot(x, metRology::dt.scaled(x, df = df_prior, sd = sd_prior_raw / (r^2)),
    type = "l")
  abline(v = scaled_phi)
}, phi_raw = slider(0.1, 400, 5, step = 0.1),
  max = slider(1, 20, 10),
  sd_prior_raw = slider(1, 600, 10),
  df_prior = slider(1, 50, 3))

manipulate({
  d <- seq(0, max, length.out = 200)
  r <- diff(range(d))
  par(mfrow = c(2, 1), mar = c(3, 3, 1, 1))
  plot(d, gauss_cor(d, phi = scaled_phi / (r^2)), type = "l",
    main = paste("scaled phi =", round(scaled_phi, 2)))
  x <- seq(0, (sd_prior * 3), length.out = 200)
  plot(x, metRology::dt.scaled(x, df = df_prior, sd = sd_prior, ncp = 0),
    type = "l")
  abline(v = scaled_phi)
}, scaled_phi = slider(0.1, 4, 0.2, step = 0.05),
  max = slider(1, 20, 10),
  sd_prior = slider(1, 30, 10),
  df_prior = slider(1, 50, 3))

# ds <- sim_mvt_rf(n_pts = nrow(g), grid = g)

# n_pts = 100, n_knots = 30, n_draws = 20, gp_scale = 0.3,
# sigma_t = 0.2, mvt = TRUE, df = 3,
# grid = cbind(lon = runif(n_pts, 5, 15), lat = runif(n_pts, 5, 15)))

draws <- lapply(c(2, 1e9),
  function(x) {
    draws <- 5
    s <- sim_mvt_rf(df = x, grid = g, n_pts = nrow(g), seed = 29,
      n_draws = draws, gp_scale = 2, sigma_t = 0.3, n_knots = 30)
    out <- reshape2::melt(s$proj)
    names(out) <- c("i", "pt", "re")
    out <- arrange(out, i, pt)
    out$nu <- x
    out$lon <- rep(g$lon, draws)
    out$lat <- rep(g$lat, draws)
    out
}) %>% bind_rows()

draws <- mutate(draws, i = paste("Draw", i))

labels <- tibble::tibble(
  nu = c(2, 1e9),
  nu_lab = c("MVT, df = 2", "MVN")
)

draws <- inner_join(draws, labels, by = "nu")

p <- draws %>%
  ggplot(aes(x = lon, y = lat, z = re, fill = re)) +
  facet_grid(nu_lab~i) +
  geom_raster() +
  # scale_fill_gradient2(low = "#01665e", mid = "#f5f5f5", high = "#8c510a") +
  viridis::scale_fill_viridis() +
  theme_gg() +
  theme(axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank()) +
  theme(panel.spacing = unit(0.001, "lines"))
print(p)
ggsave("figs/nu-rf-illustration-high-gpscale.pdf", width = 7, height = 3)

# for irregular spacing interpolation:
# ds <- akima::interp(x = lon, y = lat, z = v))
# with(fld, persp(x, y, z))

# df <- reshape2::melt(fld$z, na.rm = TRUE)
# df$lon <- fld$x[df$x]
# df$lat <- fld$y[df$y]
