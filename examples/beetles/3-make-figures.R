library(dplyr)
library(ggplot2)
library(rrfields)
library(rstan)
library(viridis)
library(assertthat)

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

# ggplot(mp, aes(long,lat, group=group)) +
#   geom_polygon(fill = NA, colour = "grey50") +
#   coord_equal(xlim = c(-150, -80), ylim = c(20, 70))
mp2 <- mp
coordinates(mp2) <- c("long", "lat")
proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
proj4string(mp2) <- CRS("+proj=longlat +datum=NAD83")
proj4string(mp2)
mp2 <- spTransform(mp2, CRS(proj))
proj4string(mp2)

mp2 <- fortify(as.data.frame(mp2))

ggplot(filter(d_highres, year == 2010), aes(x, y)) +
  geom_polygon(data = mp2, aes(long, lat, group = group),
  inherit.aes = FALSE, fill = "white", col = "grey40", lwd = 0.3) +
  geom_point(size = 0.1, alpha = 0.1) +
  facet_wrap(~year) +
  theme(panel.background = element_rect(fill="#e6e8ed")) +
  coord_fixed(xlim = range(d_highres$x), ylim = range(d_highres$y))


pred <- predict(mvt, interval = "confidence", conf_level = 0.95,
  newdata = d)
pred <- pred %>% mutate(x = d$x, y = d$y, observed = d$cover,
  year = d$year, residual = log(observed) - estimate)

g <- ggplot(dplyr::filter(pred),
  aes(x, y, fill = exp(estimate))) +
  geom_tile() + facet_wrap(~year) +
  scale_fill_viridis(option = "D", trans = "sqrt", breaks = seq(0.1, 1.3, 0.3)) +
  coord_fixed() +
  ylab(expression(10^5~UTM~North)) +
  xlab(expression(10^5~UTM~West)) +
  labs(fill = "Beetle\n% cover")

print(g)

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

