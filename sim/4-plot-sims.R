# This code creates a figure for both simulations
# The top fig is recapturing nu
# And the bottom fig is the model mismatch

# Recapture row:

library("beanplot")
library("dplyr")

axis_col <- "grey45"
cols <- RColorBrewer::brewer.pal(3, "YlOrRd")
margin_line <- 1.5
margin_color <- "grey45"
axis_lab_cex <- 0.8
mvn_bar_colour <- "#f1691390"
mvt_bar_colour <- "#2171b590"
axis_colour <- "grey45"
text_colour <- "grey45"
bar_colour <- "grey70"
bar_colour_dark <- "#2171b590"

add_label <- function(xfrac = 0, yfrac = 0.07, label = "", pos = 4, ...) {
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, ...)
}

out <- readRDS("sim/mvt-norm-sim-testing2.rds")
if ("gp_scale" %in% names(out)) out <- rename(out, gp_theta_est = gp_scale_est, gp_theta = gp_scale)

out_summary <- data.frame(out) %>%
  filter(rhat < 1.05, ess > 100) %>%
  select(-rhat, -ess)

out_summary <- mutate(out_summary, df_lab = paste0("nu==", df),
  df_lab = factor(df_lab, levels = c("nu==20", "nu==5", "nu==2.5")))

out_summary <- mutate(out_summary, n_draws_lab = paste0("Time~steps==", n_draws),
  n_draws_lab = factor(n_draws_lab, levels = c("Time~steps==25",
    "Time~steps==15", "Time~steps==5")))

plot_panel_base <- function(d, x, hlines = 2.5, col = rep(cols[[3]], 3), x_vals = c(1, 2, 3)) {
  col <- paste0(col, "")
  plot(1, 1, xlim = c(.6, 3.4), ylim = c(2, 30), type = "n",
    axes = FALSE, ann = FALSE, yaxs = "i")
  abline(h = hlines, lty = 2, col = "grey65")
  beanplot(as.formula(paste0("df_est ~ ", x)), data = d, what = c(0,1,0,0),
    log = "", col = list(col[1], col[2], col[3]), border = NA,
    add = TRUE, axes = FALSE, cutmin = 2)
  points(jitter(as.numeric(as.factor(d[,x])), amount = 0.09), d$df_est,
    col = "#00000020", cex = 0.8, pch = 20)
  axis(1, at = 1:3, labels = x_vals, col.axis = axis_col, col = axis_col, col.ticks = axis_col, las = 1)
  box(col = axis_col)
}

margin_line <- 1.5
margin_color <- "grey45"
pdf("figs/simulation-results.pdf", width = 6, height = 2.3)
par(mfrow = c(1, 3), mar = c(1.5, 0, .5, 0), oma = c(1.2, 3, 0, .5),
  cex = 0.8, tcl = -0.2, mgp = c(2, 0.4, 0))
filter(out_summary, sd_obs == "0.1", n_draws == 25) %>%
  plot_panel_base("df", hlines = c(2.5, 5, 20), col = rev(cols), x_vals = c(2.5, 5, 20))
mtext("Degrees of freedom parameter", 1, col = margin_color, line = margin_line, cex = axis_lab_cex)
mtext(expression(Estimated~nu), 2, col = margin_color, line = margin_line +.5, cex = axis_lab_cex)
axis(2, col = axis_col, col.ticks = axis_col, las = 1, col.axis = axis_col, at = c(2, 10, 20, 309))
x_text <- 0.5
cex_text <- 0.9
text(x_text, 28, "a)", pos = 4, cex = cex_text, col = margin_color)
text(x_text, 25, "Obs. CV = 0.1", pos = 4, cex = cex_text, col = margin_color)
text(x_text, 22, "25 time steps", pos = 4, cex = cex_text, col = margin_color)

filter(out_summary, sd_obs == "0.1", df == 2.5) %>%
  plot_panel_base("n_draws", x_vals = c(5, 15, 25))
mtext("Number of time steps", 1, col = margin_color, line = margin_line, cex = axis_lab_cex)
text(x_text, 28, "b)", pos = 4, cex = cex_text, col = margin_color)
text(x_text, 25, "Obs. CV = 0.1", pos = 4, cex = cex_text, col = margin_color)

filter(out_summary, df == "2.5", n_draws == 25) %>%
  plot_panel_base("sd_obs", x_vals = c(0.1, 0.6, 1.2))
text(x_text, 28, "c)", pos = 4, cex = cex_text, col = margin_color)
text(x_text, 25, "25 time steps", pos = 4, cex = cex_text, col = margin_color)
mtext("Observation error CV", 1, col = margin_color, line = margin_line, cex = axis_lab_cex)
dev.off()


# ----------------------------------------------------------------
gp_sigma <- 0.3
sigma <- 0.8
df <- 2
gp_theta <- 1.2
rmse <- readRDS("sim/match-mismatch-rmse.rds")
loo <- readRDS("sim/match-mismatch-loo.rds")
load("sim/match-mismatch-example.rda")
rmse <- mutate(rmse, perc_better=(rmse_wrong-rmse)/rmse*100)
# Plot example parameter estimates:

pdf("figs/simulation-perf.pdf", width = 6.3, height = 2.4)
par(mfrow = c(1, 3), mar = c(2.5, 2.8, .5, 0), oma = c(1.2, 0.7, 0, .5),
  cex = 0.8, tcl = -0.2, mgp = c(2, 0.4, 0))
op <- vector(mode = "list", length = 4)
for(i in seq_along(op)) {
  op[[i]]$post <- e1[[i]]
  op[[i]]$dens <- density(op[[i]]$post,
    from = quantile(op[[i]]$post, probs = 0.001)[[1]],
    to = quantile(op[[i]]$post, probs = 0.999)[[1]])
  op[[i]]$med_post <- median(op[[i]]$post)
  op[[i]]$med_dens_height <-
    op[[i]]$dens$y[max(which(op[[i]]$dens$x < op[[i]]$med_post))]
}
op_mvn <- vector(mode = "list", length = 4)
for(i in 2:length(op_mvn)) {
  op_mvn[[i]]$post <- e2[[i-1]]
  op_mvn[[i]]$dens <- density(op_mvn[[i]]$post,
    from = quantile(op_mvn[[i]]$post, probs = 0.001)[[1]],
    to = quantile(op_mvn[[i]]$post, probs = 0.999)[[1]])
  op_mvn[[i]]$med_post <- median(op_mvn[[i]]$post)
  op_mvn[[i]]$med_dens_height <-
    op_mvn[[i]]$dens$y[max(which(op_mvn[[i]]$dens$x < op_mvn[[i]]$med_post))]
}
par(tck = -0.02)
xlim <- c(0, 3.2)
params <- length(op)
plot(1, 1, xlim = xlim, ylim = c(params+.5, .8), type = "n",
  ylab = "", xlab = "", axes = FALSE, xaxs = "i")
abline(v = 0, lty = 2, col = "grey40", lwd = 0.6)
scaling_factor <- 85
gap <- 0.4
true <- c(df, gp_sigma, gp_theta, sigma)
names(e1)
for(i in 1:length(op)) {
  if (i == 4) scaling_factor <- scaling_factor * 1.5
  if (i > 1) {
    polygon(c(op_mvn[[i]]$dens$x, rev(op_mvn[[i]]$dens$x)),
      i + gap + c(op_mvn[[i]]$dens$y/scaling_factor, -rev(op_mvn[[i]]$dens$y/scaling_factor)),
      border = mvn_bar_colour, lwd = 0.5, col = mvn_bar_colour)
  }
  polygon(c(op[[i]]$dens$x, rev(op[[i]]$dens$x)),
    i + c(op[[i]]$dens$y/scaling_factor, -rev(op[[i]]$dens$y/scaling_factor)),
    border = mvt_bar_colour, lwd = 0.5, col = mvt_bar_colour)
  points(true[i], i+gap/2, col = "grey10", pch = 4)
  par(xpd = NA)
  par(xpd = FALSE)
}

text(1.7, 2.8+0.7, label = "MVN", cex = cex_text, col = mvn_bar_colour, pos = 4)
text(1.7, 2.8+0.2, label = "MVT", cex = cex_text, col = bar_colour_dark, pos = 4)

axis(2, at = seq_along(op)+gap/2,
  labels = c(expression(nu), expression(sigma[GP]), expression(theta[GP]), expression(sigma[obs])),
  las = 1, lwd = 0, line = 0.2, col.axis = margin_color, col = axis_colour)
axis(1, at = seq(0, 5, 1), mgp = c(2, 0.3, 0), col = axis_colour, col.axis = axis_colour)
mtext("Coefficient value", side = 1, line = margin_line, cex = axis_lab_cex, col = margin_color)
box(col = axis_colour)
add_label(label = "a)", cex = cex_text, col = margin_color)
mtext("Parameter", side = 2, line = 2.4, col = axis_colour, cex = 0.8)


h <- hist(rmse$perc_better, probability = TRUE,
  breaks = seq(range(rmse$perc_better)[1], range(rmse$perc_better)[2],
    length.out = 12), plot = FALSE,
  warn.unused = FALSE)
plot(0,0, xlim = range(rmse$perc_better), ylim = c(0, max(h$counts)), type = "n",
  xlab = "", axes = FALSE,
  ylab = "")
mtext("% worse RMSE\nwith MVN vs. MVT", side = 1, line = margin_line+1,
  cex = axis_lab_cex, col = margin_color)
for(j in seq_along(h$breaks)) {
  rect(h$breaks[j], 0, h$breaks[j+1], h$counts[j], border = "white",
    col = bar_colour, lwd = 1)
}
abline(v = 1, lty = 2, col = "grey30")
axis(1, col = axis_colour, col.axis = axis_colour)
box(col = axis_colour)
add_label(label = "b)", cex = cex_text, col = margin_color)
axis(2, col = axis_colour, col.axis = axis_colour, las = 1, cex.axis = 1)
mtext("Frequency", side = 2, line = 1.5, col = axis_colour, cex = 0.8)

h <- hist(loo$delta_loo, probability = TRUE,
  breaks = seq(range(loo$delta_loo)[1],
    range(loo$delta_loo)[2], length.out = 12), plot = FALSE,
  warn.unused = FALSE)
plot(0,0, xlim = range(loo$delta_loo), ylim = c(0, max(h$counts)), type = "n",
  xlab = "", axes = FALSE,
  ylab = "")
mtext(expression(Delta~LOOIC), side = 1, line = margin_line,
  cex = axis_lab_cex, col = margin_color)
mtext("(MVT - MVN)", side = 1, line = margin_line+1, cex = axis_lab_cex, col = margin_color)
for(j in seq_along(h$breaks)) {
  rect(h$breaks[j], 0, h$breaks[j+1], h$counts[j], border = "white",
    col = bar_colour, lwd = 1)
}
abline(v = 1, lty = 2, col = "grey30")
axis(1, col = axis_colour, col.axis = axis_colour)
box(col = axis_colour)
add_label(label = "c)", cex = cex_text, col = margin_color)
axis(2, col = axis_colour, col.axis = axis_colour, las = 1, cex.axis = 1)
mtext("Frequency", side = 2, line = 1.5, col = axis_colour, cex = 0.8)

dev.off()

