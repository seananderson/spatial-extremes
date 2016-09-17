library(dplyr)
library(ggplot2)
source("simulationTesting/sim_mvt_rf.R")
source("simulationTesting/theme_gg.R")
g <- expand.grid(lon = seq(0, 1, length.out = 25),
  lat = seq(0, 1, length.out = 25))

# ds <- sim_mvt_rf(n_pts = nrow(g), grid = g)

# n_pts = 100, n_knots = 30, n_draws = 20, gp_scale = 0.3,
# sigma_t = 0.2, mvt = TRUE, df = 3,
# grid = cbind(lon = runif(n_pts, 5, 15), lat = runif(n_pts, 5, 15)))

draws <- lapply(c(2, 1e9),
  function(x) {
    draws <- 5
    s <- sim_mvt_rf(df = x, grid = g, n_pts = nrow(g), seed = 29,
      n_draws = draws, gp_scale = 12, sigma_t = 0.2)
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
  theme(panel.margin = unit(0.001, "lines"))
ggsave("figs/nu-rf-illustration.pdf", width = 7, height = 3)

# for irregular spacing interpolation:
# ds <- akima::interp(x = lon, y = lat, z = v))
# with(fld, persp(x, y, z))

# df <- reshape2::melt(fld$z, na.rm = TRUE)
# df$lon <- fld$x[df$x]
# df$lat <- fld$y[df$y]
