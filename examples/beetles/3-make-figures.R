library(dplyr)
library(ggplot2)
library(rrfields)
library(rstan)
library(viridis)
library(assertthat)
library(ggsidekick) # devtools::install_github("seananderson/ggsidekick")

d <- readRDS("examples/beetles/mountain-pine-beetle-data.rds")

# --------------------
# This creates a map of Washington and Oregon State:

# Bring in the background map data

par(las = 1, mgp = c(2, 0.6, 0), xaxs = "i", yaxs = "i")
par(mfrow = c(2, 2), mar = c(3, 3, 1, 1), oma = c(3, 3, 1, 1))

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
  xlab = expression(10^5~UTM~West), ylab = expression(10^5~UTM~North),
  asp = 1)
rect(-30, 20, -10, 35, col = "#e6e8ed")
for (i in unique(mp2$group)) {
  x <- filter(mp2, group == i)
  with(x, polygon(long/1e5, lat/1e5, border = "grey40", lwd = 0.7, col = "white"))
}
with(filter(d_highres, year == 2010), points(x/1e5, y/1e5, cex = 0.1,
  col = "#00000020", pch = 19))

# ggplot(filter(d_highres, year == 2010), aes(x/1e5, y/1e5)) +
#   geom_polygon(data = mp2, aes(long/1e5, lat/1e5, group = group),
#     inherit.aes = FALSE, fill = "white", col = "grey40", lwd = 0.3) +
#   geom_point(size = 0.1, alpha = 0.2) +
#   # facet_wrap(~year) +
#   theme_sleek()+
#   theme(panel.background = element_rect(fill="#e6e8ed")) +
#   coord_fixed(xlim = range(d_highres$x)/1e5 + margin_gap_x, ylim = range(d_highres$y)/1e5 + margin_gap_y) +
#   xlab(expression(10^5~UTM~West)) +
#   ylab(expression(10^5~UTM~North))
# ggsave("figs/map.pdf", width = 4, height = 5)

# --------------------
# Create a panel showing the nu parameter from the mvt model
load("examples/beetles/mountain-pine-beetle-pnw-raster-mvn-lognormal-500x500.rda")
load("examples/beetles/mountain-pine-beetle-pnw-raster-mvt-lognormal-500x500.rda")
e <- extract(mvt$model)

prior <- data.frame(df = rgamma(1e6, shape = 2, rate = 0.1))  %>%
  filter(df > 2)

# ggplot(data.frame(df = e$df), aes(df)) +
#   geom_density(fill = "grey60", colour = "grey60") +
#   theme_sleek() +
#   geom_density(data = prior, aes(df), inherit.aes = FALSE,
#     colour = "grey30", lty = 2) +
#   scale_x_continuous(limits = c(2, 30))

plot(0,0, xlim = c(2, 40), ylim = c(0, 0.7),type = "n",
  xlab = expression(nu~(degrees~of~freedom~parameter)),
  ylab = "Probability density")
h <- hist(e$df, probability = TRUE, breaks = seq(2, 40, 1.25), plot = FALSE,
  warn.unused = FALSE)
for(j in seq_along(h$breaks)) {
  rect(h$breaks[j], 0, h$breaks[j+1], h$density[j], border = "white",
    col = "grey70", lwd = 1)
}
lines(density(prior$df, cut = 0, from = 2, to = 40), col = "grey40",
  lty = 2)

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

# ggplot(sc, aes(lp)) +
  # geom_histogram(fill = "red", alpha = 0.5) +
  # geom_histogram(data = sc2, aes(lp), fill = "blue", alpha = 0.5)

all_samples <- c(sc$lp, sc2$lp)
l <- min(all_samples)
u <- max(all_samples)
h <- hist(sc$lp, probability = TRUE,
  breaks = seq(l, u, length.out = 25), plot = FALSE,
  warn.unused = FALSE)
h2 <- hist(sc2$lp, probability = TRUE,
  breaks = seq(l, u, length.out = 25), plot = FALSE,
  warn.unused = FALSE)
plot(0,0, xlim = c(l, u), ylim = c(0, max(h$counts) * 1.02),type = "n",
  xlab = "Log predictive density", ylab = "Count")
for(j in seq_along(h$breaks)) {
  rect(h$breaks[j], 0, h$breaks[j+1], h$counts[j], border = "white",
    col = "#FF000050", lwd = 1)
}
for(j in seq_along(h2$breaks)) {
  rect(h2$breaks[j], 0, h2$breaks[j+1], h2$counts[j], border = "white",
    col = "#00000050", lwd = 1)
}


# --------------------------
# Look at the ratio of the credible intervals
pred <- predict(mvt, interval = "confidence", conf_level = 0.95,
  newdata = d, type = "link")
pred <- data.frame(pred, d)
pred <- pred %>% mutate(residual = log(cover) - estimate)

pred_mvn <- predict(mvn, interval = "confidence", conf_level = 0.95,
  newdata = d, type = "link")
pred_mvn <- data.frame(pred_mvn, d)
pred_mvn <- pred_mvn %>% mutate(residual = log(cover) - estimate)

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
  ylim = c(0, max(h$counts) * 1.02), type = "n",
  xlab = "Ratio of credible intervals (MVN/MVT)", ylab = "Count")
for(j in seq_along(h$breaks)) {
  rect(h$breaks[j], 0, h$breaks[j+1], h$counts[j], border = "white",
    col = "grey70", lwd = 1)
}
abline(v = 1, lty = 2, col = "grey30")

# ggplot(cis, aes(log(conf_mvt))) +
  # geom_histogram(fill = "red", alpha = 0.5) +
  # geom_histogram(data = cis, aes(log(conf_mvn)), fill = "blue", alpha = 0.5)

# ----------------
# RMSE error plot, not going to show
# ############
# pred <- predict(mvt, interval = "prediction", conf_level = 0.95,
#   newdata = d, type = "response")
# pred <- data.frame(pred, d)
#
# pred_mvn <- predict(mvn, interval = "prediction", conf_level = 0.95,
#   newdata = d, type = "response")
# pred_mvn <- data.frame(pred_mvn, d)
#
# combined <- bind_rows(
#   mutate(pred, model = "mvt"),
#   mutate(pred_mvn, model = "mvn"))
# ##############
# combined <- mutate(combined, id = paste(x, y))
# combined %>%
#   filter(hold_out) %>%
#   group_by(model) %>%
#   mutate(included = cover > conf_low & cover < conf_high) %>%
#   summarise(coverage = sum(included) / n())
#
# rmse <- combined %>% filter(hold_out) %>%
#   group_by(model) %>%
#   mutate(error = log(estimate) - log(cover)) %>%
#   ungroup()
#
# ggplot(rmse, aes(abs(error), colour = model)) +
#   geom_freqpoly()
#
# r <- rmse %>%
#   group_by(model) %>%
#   summarise(rmse = sqrt(mean(error^2)))
#
# r <- as.data.frame(r)
# round(100 * (r[r$model=="mvn","rmse"] - r[r$model=="mvt","rmse"]) /
#     r[r$model=="mvt","rmse"], 1)

# ----------------------
# And the main spatiotemporal prediction plot:
# mean(holdout_data$scores - holdout_data$scoresn)

# sum(rowMeans(scores))/nrow(p)
# sum(rowMeans(scoresn))/nrow(pn)

# This creates a data frame to label the years:
yr <- unique(select(d, year))
yr <- mutate(yr, x = -17, y = 24)

g <- ggplot(dplyr::filter(pred),#, year %in% c(2002, 2006, 2010, 2014)),
  aes(x, y, colour = exp(estimate), fill = exp(estimate))) +
  # stat_summary_hex(binwidth =diff(sort(unique(d$y)))[[1]]) +
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

# print(g)
# ggsave("figs/beetles-mvt-predictions.pdf", width = 5.4, height = 5)


# https://flic.kr/p/rj1GFA
# https://flic.kr/p/e8SuTi https://creativecommons.org/licenses/by/2.0/
