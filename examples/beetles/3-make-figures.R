library(dplyr)
library(ggplot2)
library(rrfields)
library(rstan)
library(viridis)
library(assertthat)
library(ggsidekick) # devtools::install_github("seananderson/ggsidekick")

d <- readRDS("examples/beetles/mountain-pine-beetle-data.rds")

# Bring in the background map data

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

# --------------------
# This creates a map of Washington and Oregon State:
margin_gap_y <- c(-2, 2)
margin_gap_x <- c(-1, 2)

ggplot(filter(d_highres, year == 2010), aes(x/1e5, y/1e5)) +
  geom_polygon(data = mp2, aes(long/1e5, lat/1e5, group = group),
    inherit.aes = FALSE, fill = "white", col = "grey40", lwd = 0.3) +
  geom_point(size = 0.1, alpha = 0.2) +
  # facet_wrap(~year) +
  theme_sleek()+
  theme(panel.background = element_rect(fill="#e6e8ed")) +
  coord_fixed(xlim = range(d_highres$x)/1e5 + margin_gap_x, ylim = range(d_highres$y)/1e5 + margin_gap_y) +
  xlab(expression(10^5~UTM~West)) +
  ylab(expression(10^5~UTM~North))
ggsave("figs/map.pdf", width = 4, height = 5)

# --------------------
# Create a panel showing the nu parameter from the mvt model
load("examples/beetles/mountain-pine-beetle-pnw-raster-mvn-lognormal-500x500.rda")
load("examples/beetles/mountain-pine-beetle-pnw-raster-mvt-lognormal-500x500.rda")
e <- extract(mvt$model)

prior <- data.frame(df = rgamma(1e6, shape = 2, rate = 0.1))  %>%
  filter(df > 2)

ggplot(data.frame(df = e$df), aes(df)) +
  geom_density() +
  theme_sleek() +
  geom_density(data = prior, aes(df), inherit.aes = FALSE) +
  scale_x_continuous(limits = c(2, 20))


posterior <- density(e$df, cut = 0, from = 2, to = 20)
# hist(e$df, probability = TRUE)
par(las = 1, mgp = c(2, 0.6, 0), xaxs = "i")
plot(0,0, xlim = c(2, 20), ylim = c(0, 0.7),type = "n")
lines(posterior)
# hist(e$df, probability = TRUE)
lines(density(prior$df, cut = 0, from = 2, to = 20), col = "grey40",
  lty = 2)




# --------------------
# Make the main spatiotemporal prediction plot for the
# MVT random fields model
pred <- predict(mvt, interval = "confidence", conf_level = 0.95,
  newdata = d)
pred <- pred %>% mutate(x = d$x, y = d$y, observed = d$cover,
  year = d$year, residual = log(observed) - estimate)

pred_mvn <- predict(mvn, interval = "confidence", conf_level = 0.95,
  newdata = d)
pred_mvn <- pred_mvn %>% mutate(x = d$x, y = d$y, observed = d$cover,
  year = d$year, residual = log(observed) - estimate)

combined <- bind_rows(
  mutate(pred, model = "mvt"),
  mutate(pred_mvn, model = "mvn"))

combined <- mutate(combined, id = paste(x, y))
# combined$hold_out <- c(d$hold_out, d$hold_out)

cis <- mutate(combined, conf_width = conf_high - conf_low) %>%
  group_by(year, x, y) %>%
  summarize(conf_width_ratio =
      conf_width[model=="mvn"]/conf_width[model=="mvt"]
    # hold_out = unique(hold_out)
  ) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(median_ratio = median(conf_width_ratio))
hist(cis$conf_width_ratio)

########
# Let's try the log predictive posterior distribution
p <- predict(mvt, interval = "confidence", conf_level = 0.95,
  type = "link", newdata = filter(d, hold_out), return_mcmc = TRUE)
pn <- predict(mvn, interval = "confidence", conf_level = 0.95,
  type = "link", newdata = filter(d, hold_out), return_mcmc = TRUE)
sigma <- extract(mvt$model)$sigma
sigman <- extract(mvn$model)$sigma

holdout_data <- filter(d, hold_out)


scores <- matrix(0, nrow = nrow(p), ncol = ncol(p))
for(draw in seq_len(ncol(p))){
  scores[, draw] <- 1/dlnorm(holdout_data$cover,
    mean = p[, draw], sd = sigma[draw])
}

scoresn <- matrix(0, nrow = nrow(p), ncol = ncol(p))
for(draw in seq_len(ncol(p))){
  scoresn[, draw] <- 1/dlnorm(holdout_data$cover,
    mean = pn[, draw], sd = sigman[draw])
}

-mean(log(rowMeans(scores)))
-mean(log(rowMeans(scoresn)))

sum(rowMeans(scores))/nrow(p)
sum(rowMeans(scoresn))/nrow(pn)

lps <- -sum(log(rowMeans(scores)))/ length(TestIdx)

scores <- vector(mode = "numeric", length = nrow(p))
for (i in seq_len(nrow(p))) {
  scores[i] <- density(p[i,],
    from=log(holdout_data$cover[1]), to=log(holdout_data$cover[1]), n=1)$y
}
sum((scores))

scores <- vector(mode = "numeric", length = nrow(p))
for (i in seq_len(nrow(pt))) {
  scores[i] <- density(pt[i,],
    from=log(holdout_data$cover[1]), to=log(holdout_data$cover[1]), n=1)$y
}
sum((scores))

#################


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

print(g)
ggsave("figs/beetles-mvt-predictions.pdf", width = 5.4, height = 5)

