library(dplyr)
library(ggplot2)
library(rrfields)
library(rstan)
library(viridis)
library(assertthat)
library(ggsidekick) # devtools::install_github("seananderson/ggsidekick")

d <- readRDS("examples/beetles/mountain-pine-beetle-data.rds")
axis_colour <- "grey40"
text_colour <- "grey25"
bar_colour <- "#2171b570"
bar_colour_dark <- "#2171b590"

plot_jpeg = function(path, add=FALSE) {
  library('jpeg')
  jpg = readJPEG(path, native=TRUE) # read the file
  res = dim(jpg)[2:1] # get the resolution
  if (!add) # initialize an empty plot area if add==FALSE
    plot(1,1,xlim=c(1,res[1]),ylim=c(1,res[2]),asp=1,
      type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  rasterImage(jpg,1,1,res[1],res[2])
}

# --------------------
# This creates a map of Washington and Oregon State:

# Bring in the background map data

pdf("figs/beetle-performance.pdf", width = 6.0, height = 4.5)
layout(rbind(c(1, 1, 1, 2, 2, 2), c(3, 3, 4, 4, 5, 5)))
par(las = 1, mgp = c(1.8, 0.4, 0), xaxs = "i", yaxs = "i")
par(mar = c(3.0, 3.0, 0, 0), oma = c(1, 0.5, 0.5, 0.5))
par(cex = 0.75)
par(tcl = -0.2)
par(col.lab = text_colour)

library(mapdata)
library(maps)
library(rgdal)
mpc <- ggplot2::map_data("worldHires", "Canada")
mps <- ggplot2::map_data("state")
mpc$group <- mpc$group + max(mps$group)
mp <- rbind(mpc,mps)
mp2 <- mp
coordinates(mp2) <- c("long", "lat")
proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
proj4string(mp2) <- CRS("+proj=longlat +datum=NAD83")
mp2 <- spTransform(mp2, CRS(proj))
mp2 <- fortify(as.data.frame(mp2))

id <- "mountain-pine-beetle-pnw-raster"
nbin <- 500
d_highres <- readRDS(paste0("examples/beetles/",
  id, "-", "dataframe", "-", nbin, "x", nbin, "highres.rds"))

d_highres <- dplyr::filter(d_highres, cover > 0)

margin_gap_y <- c(-2, 2)
margin_gap_x <- c(-1, 2)

plot(0,0, xlim = range(d_highres$x)/1e5 + margin_gap_x,
  ylim = range(d_highres$y)/1e5 + margin_gap_y, type = "n", ann = TRUE,
  xlab = "", ylab = expression(10^5~UTM~North),
  asp = 1, axes = FALSE)
rect(-30, 20, -10, 35, col = "#e6e8ed")
for (i in unique(mp2$group)) {
  x <- filter(mp2, group == i)
  with(x, polygon(long/1e5, lat/1e5, border = "grey40", lwd = 0.7, col = "white"))
}

mtext(expression(10^5~UTM~West), side = 1, line = 1.5,
  cex = 0.75, col = text_colour)
with(filter(d_highres, year == 2010), points(x/1e5, y/1e5, cex = 0.1,
  col = "#c42e1708", pch = 19))
axis(1, col = axis_colour, col.axis = axis_colour)
axis(2, col = axis_colour, col.axis = axis_colour)
text(-17, 32, "Canada", col = text_colour, pos = 4)
text(-14, 28, "Montana", col = text_colour, pos = 4)
text(-22, 26, "Oregon", col = text_colour, pos = 4)
text(-21.3, 29, "Washington", col = text_colour, pos = 4)
# text(-23.7, 23, "California", col = text_colour, pos = 4)
text(-16, 24, "Idaho", col = text_colour, pos = 4)
text(-27, 29, "Pacific\nOcean", col = text_colour, pos = 4)
box(col = axis_colour)

plot_jpeg("examples/beetles/8621706223_2dbc7e1781_k-edited.jpg")
mtext("Pine beetle infested forest", side = 1, line = 1,
  cex = 0.75, col = text_colour)

# --------------------
# Create a panel showing the nu parameter from the mvt model
load("examples/beetles/mountain-pine-beetle-pnw-raster-mvn-lognormal-500x500.rda")
load("examples/beetles/mountain-pine-beetle-pnw-raster-mvt-lognormal-500x500.rda")
e <- extract(mvt$model)

prior <- data.frame(df = rgamma(1e6, shape = 2, rate = 0.1))  %>%
  filter(df > 2)

plot(0,0, xlim = c(2, 25), ylim = c(0, 0.58),type = "n",
  xlab = "", axes = FALSE,
  ylab = "Probability density")
mtext("Degrees of freedom parameter", side = 1, line = 1.5, cex = 0.75, col = text_colour)
mtext(expression(nu), side = 1, line = 2.5, cex = 0.75, col = text_colour)
h <- hist(e$df, probability = TRUE, breaks = seq(2, 40, 1.25), plot = FALSE,
  warn.unused = FALSE)
for(j in seq_along(h$breaks)) {
  rect(h$breaks[j], 0, h$breaks[j+1], h$density[j], border = "white",
    col = bar_colour, lwd = 1)
}
lines(density(prior$df, cut = 0, from = 2, to = 40), col = "grey40",
  lty = 2)
text(15, 0.07, "Prior", col = text_colour, pos = 4)
text(3.5, 0.3, "Posterior", col = bar_colour_dark, pos = 4)
axis(1, col = axis_colour, col.axis = axis_colour)
axis(2, col = axis_colour, col.axis = axis_colour)
box(col = axis_colour)

# ----------------
# What about the distribution of log predictive density on the held out data?
dat <-  dplyr::filter(d, hold_out)
p <- predict(mvt, type = "link", newdata =  dat, return_mcmc = TRUE)
pn <- predict(mvn, type = "link", newdata =  dat, return_mcmc = TRUE)
sigma <- extract(mvt$model)$sigma
sigman <- extract(mvn$model)$sigma

scores <- matrix(0, nrow = nrow(p), ncol = ncol(p))
for(draw in seq_len(ncol(p))){
  scores[, draw] <- dlnorm(dat$cover,
    mean = p[, draw], sd = sigma[draw], log = TRUE)
}

scoresn <- matrix(0, nrow = nrow(p), ncol = ncol(p))
for(draw in seq_len(ncol(pn))){
  scoresn[, draw] <- dlnorm(dat$cover,
    mean = pn[, draw], sd = sigman[draw], log = TRUE)
}

# Sum the log densities within each MCMC iteration:
sc <- data.frame(lp = apply(scores, 2, sum), model = "mvt")
sc2 <- data.frame(lp = apply(scoresn, 2, sum), model = "mvn")

# Make sure we summed across MCMC chains and not across rows of data:
assertthat::assert_that(identical(length(sigma), nrow(sc)))

all_samples <- c(sc$lp, sc2$lp)
l <- min(all_samples)
u <- max(all_samples)
h <- hist(sc$lp, probability = TRUE,
  breaks = seq(l, u, length.out = 25), plot = FALSE,
  warn.unused = FALSE)
h2 <- hist(sc2$lp, probability = TRUE,
  breaks = seq(l, u, length.out = 25), plot = FALSE,
  warn.unused = FALSE)
plot(0,0, xlim = c(l, u), ylim = c(0, max(h$density) * 1.02),type = "n",
  xlab = "", ylab = "", axes = FALSE)
mtext("Log predictive density", side = 1, line = 1.5, cex = 0.75, col = text_colour)
mtext("for held-out data", side = 1, line = 2.5, cex = 0.75, col = text_colour)
for(j in seq_along(h2$breaks)) {
  rect(h2$breaks[j], 0, h2$breaks[j+1], h2$density[j], border = "white",
    col = "#f1691370", lwd = 1)
}
for(j in seq_along(h$breaks)) {
  rect(h$breaks[j], 0, h$breaks[j+1], h$density[j], border = "white",
    col = bar_colour, lwd = 1)
}
text(-24, 0.08, "MVN", col = "#f1691390", pos = 4)
text(-12, 0.11, "MVT", col = bar_colour_dark, pos = 4)
axis(1, col = axis_colour, col.axis = axis_colour)
axis(2, col = axis_colour, col.axis = axis_colour)
box(col = axis_colour)

# --------------------------
# Look at the ratio of the credible intervals
if (!exists("pred")) {
  pred <- predict(mvt, interval = "confidence", conf_level = 0.95,
    newdata = d, type = "link")
  pred <- data.frame(pred, d)
  pred <- pred %>% mutate(residual = log(cover) - estimate)
}

if (!exists("pred_mvn")) {
  pred_mvn <- predict(mvn, interval = "confidence", conf_level = 0.95,
    newdata = d, type = "link")
  pred_mvn <- data.frame(pred_mvn, d)
  pred_mvn <- pred_mvn %>% mutate(residual = log(cover) - estimate)
}

combined <- bind_rows(
  mutate(pred, model = "mvt"),
  mutate(pred_mvn, model = "mvn"))

combined <- mutate(combined, id = paste(x, y))

cis <- mutate(combined, conf_width = exp(conf_high) - exp(conf_low)) %>%
  filter(!hold_out) %>%
  group_by(year, x, y) %>%
  summarize(conf_width_ratio =
      conf_width[model=="mvn"]/conf_width[model=="mvt"],
    conf_mvt = conf_width[model=="mvt"],
    conf_mvn = conf_width[model=="mvn"],
    hold_out = unique(hold_out)
  ) %>%
  ungroup()
h <- hist(cis$conf_width_ratio, probability = TRUE,
  breaks = seq(0, max(cis$conf_width_ratio)*1.001, length.out = 20), plot = FALSE,
  warn.unused = FALSE)
plot(0,0, xlim = c(0, max(cis$conf_width_ratio)*1.03),
  ylim = c(0, max(h$density) * 1.02), type = "n", axes = FALSE,
  xlab = "", ylab = "")
mtext("Ratio of credible intervals", side = 1, line = 1.5, cex = 0.75, col = text_colour)

mtext("(MVN/MVT)", side = 1, line = 2.5, cex = 0.75, col = text_colour)
for(j in seq_along(h$breaks)) {
  rect(h$breaks[j], 0, h$breaks[j+1], h$density[j], border = "white",
    col = "grey70", lwd = 1)
}
abline(v = 1, lty = 2, col = "grey30")
axis(1, col = axis_colour, col.axis = axis_colour)
axis(2, col = axis_colour, col.axis = axis_colour)
box(col = axis_colour)
text(0, 0.95, "MVN\nmore\nprecise", col = "#f1691390", pos = 4)
text(1.7, 0.95, "MVT\nmore\nprecise", col = bar_colour_dark, pos = 4)

dev.off()

# This creates a data frame to label the years:
yr <- unique(select(d, year))
yr <- mutate(yr, x = -17, y = 24)

g <- ggplot(dplyr::filter(pred),#, year %in% c(2002, 2006, 2010, 2014)),
  aes(x, y, colour = exp(estimate), fill = exp(estimate))) +
  geom_tile() +
  facet_wrap(~year) +
  theme_sleek() +
  scale_fill_viridis(option = "B", trans = "sqrt", breaks = seq(0.1, 1.3, 0.3)) +
  scale_color_viridis(option = "B", trans = "sqrt", breaks = seq(0.1, 1.3, 0.3)) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(-23, -17, 4)) +
  ylab(expression(10^5~UTM~North)) +
  xlab(expression(10^5~UTM~West)) +
  labs(color = "Beetle\n% cover", fill = "Beetle\n% cover") +
  theme(panel.spacing = unit(-0.1, "lines")) +
  theme(strip.background = element_blank(),
    strip.text.x = element_blank(), legend.position = "right") +
  geom_text(data = yr, aes(x = x, y = y, label = year), inherit.aes = FALSE, size = 3)

ggsave("figs/beetles-mvt-predictions.pdf", width = 5.4, height = 5)

# https://flic.kr/p/rj1GFA
# https://flic.kr/p/e8SuTi https://creativecommons.org/licenses/by/2.0/
