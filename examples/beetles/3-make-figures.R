library(dplyr)
library(ggplot2)
library(rrfields)
library(rstan)
library(viridis)
library(assertthat)
library(ggsidekick) # devtools::install_github("seananderson/ggsidekick")

id <- "mountain-pine-beetle-pnw-raster"
nbin <- 500
model_file <- paste0("examples/beetles/", id, "-",
  "mvt-normal", "-", nbin, "x", nbin, ".rda")
load(model_file)
model_file <- paste0("examples/beetles/", id, "-",
  "mvn-lognormal", "-", nbin, "x", nbin, ".rda")
load(model_file)

d <- readRDS(paste0("examples/beetles/",
  id, "-", "dataframe", "-", nbin, "x", nbin, ".rds"))
d_highres <- readRDS(paste0("examples/beetles/",
  id, "-", "dataframe", "-", nbin, "x", nbin, "highres.rds"))

d <- dplyr::filter(d, cover > 0)
d_highres <- dplyr::filter(d_highres, cover > 0)

d$x <- d$x/1e5
d$y <- d$y/1e5

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
plot(0,0, xlim = c(2, 10), ylim = c(0, 0.7),type = "n")
lines(posterior)
lines(density(prior$df, cut = 0, from = 2, to = 10), col = "grey40",
  lty = 2)




# --------------------
# Make the main spatiotemporal prediction plot for the 
# MVT random fields model 
pred <- predict(mvt, interval = "confidence", conf_level = 0.95,
  newdata = d)
pred <- pred %>% mutate(x = d$x, y = d$y, observed = d$cover,
  year = d$year, residual = log(observed) - estimate)

# This creates a data frame to label the years:
yr <- unique(select(d, year))
yr <- mutate(yr, x = -17, y = 24)

g <- ggplot(dplyr::filter(pred),#, year %in% c(2002, 2006, 2010, 2014)),
  aes(x, y, colour = exp(estimate), fill = exp(estimate))) +
  # stat_summary_hex(binwidth =diff(sort(unique(d$y)))[[1]]) +
  geom_tile() +
  facet_wrap(~year) +
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

