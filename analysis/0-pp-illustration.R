library(rrfields)
library(ggsidekick)
library(tidyverse)
theme_set(theme_sleek())

set.seed(1)
x <- sim_rrfield(df = 99, n_data_points = 30, seed = NULL,
  n_draws = 1, n_knots = 8, gp_scale = 0.5, gp_sigma = 1,
  sd_obs = 0.01, obs_error = "normal")

# names(x)

initialize_plot <- function(title = "", cex = 0.7) {
  plot(x$dat$lon, x$dat$lat, type = "n", ann = FALSE, axes = FALSE)
  box(col = "grey50")
  mtext(title, cex = cex, line = 0.15, col = "grey30")
}
plot_data <- function(pch = 21, col = "grey30", cex = 1) {
  points(x$dat$lon, x$dat$lat, pch = pch, col = col, cex = cex)
}
plot_knots <- function() {
  points(x$knots[,1], x$knots[,2], col = "red", pch = 19, cex = 1.1)
}
plot_covariance <- function() {
  for (i in seq_len(nrow(x$knot))) {
    for (j in seq_len(nrow(x$dat))) {
      segments(x0 = x$knots[i,1], y0 = x$knots[i,2],
        x1 = x$dat$lon[j], y1 = x$dat$lat[j], col = "#00000020")
    }
  }
}


width <- 5.5
golden_ratio <- (1+sqrt(5))/2
pdf("figs/pp-illustration.pdf", width = width,
  height = width / golden_ratio)
par(mfrow = c(2, 3), mar = c(0, 0, 2.9, 0), oma = c(1, 3, 0, 1),
  cex = 0.6)

initialize_plot("Observe spatial data")
plot_data()
mtext("Before fitting model", side = 2, cex = 0.7, line = 0.3,
  col = "grey30")

initialize_plot("Select knot locations")
plot_data()
plot_knots()

initialize_plot("Calculate knot-data covariance")
plot_data()
plot_knots()
plot_covariance()

initialize_plot("Model knot values\nas random field")
plot_knots()
mtext("Model fitting process", side = 2, cex = 0.7, line = 0.3,
  col = "grey30")

initialize_plot("Project knot values to data\nlocations with knot-data covariance")
plot_data(pch = 4, cex = 1.3)
plot_knots()
plot_covariance()

initialize_plot("Evaluate likelihood")
plot_data(col = "#00000080")
plot_data(pch = 4, cex = 1.3, col = "#00000080")
dev.off()
