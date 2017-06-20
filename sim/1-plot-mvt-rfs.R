library(dplyr)
library(ggplot2)
source("simulationTesting/sim_mvt_rf.R")
library(ggsidekick)
library(rrfields)

g <- expand.grid(lon = seq(1, 10, 0.5),
  lat = seq(1, 10, 0.5))

exp_cor <- function(delta_ij, phi) {
  exp(-phi * delta_ij)
}

draws <- lapply(c(2, 1e9),
  function(x) {
    g <- expand.grid(lon = seq(1, 10, length.out = 25),
      lat = seq(1, 10, length.out = 25))
    draws <- 3
    s <- rrfields::sim_rrfield(df = x, n_draws = draws,
      gp_scale = 1.6, gp_sigma = 0.3, n_knots = 30, seed = 9,
      g = g, n_data_points = nrow(g))
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
  nu_lab = c("MVT, v = 2", "MVN")
)

draws <- inner_join(draws, labels, by = "nu")

p <- draws %>%
  ggplot(aes(x = lon, y = lat, z = re, fill = re)) +
  facet_grid(nu_lab~i) +
  geom_raster() +
  viridis::scale_fill_viridis(option = "C") +
  theme_sleek() +
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
  theme(panel.spacing = unit(-0.15, "lines"))

ggsave("figs/nu-rf-illustration-small.pdf", width = 4.5, height = 3)

