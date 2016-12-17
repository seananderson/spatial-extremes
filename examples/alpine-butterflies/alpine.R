# Matter, S. F., N. Keyghobadi, and J. Roland. 2014. Ten years of abundance data within a spatial population network of the alpine butterfly, Parnassius smintheus. Ecology 95:2985–2985.

# see also:
# Caplins, S. A., K. J. Gilbert, C. Ciotir, J. Roland, S. F. Matter, and N. Keyghobadi. 2014. Landscape structure and the genetic effects of a population collapse. Proceedings of the Royal Society B: Biological Sciences 281:20141798–20141798.

library(tidyverse)
abund <- read_csv("examples/alpine-butterflies/data-raw/PopulationAbunTransect.csv") %>%
  rename(Meadow = Population) %>%
  mutate(Date = lubridate::dmy(Date))
sites <- read_csv("examples/alpine-butterflies/data-raw/Landscape.csv")

abund <- abund %>% inner_join(sites)
names(abund) <- tolower(names(abund))

abund <- abund %>% mutate(year = lubridate::year(date)) %>%
  rename(lat = `lat. (decimal degrees n)`, lon = `lon. (decimal degrees w)`)

abund %>% group_by(meadow, year) %>% summarise(n = n())

d <- abund %>% group_by(meadow, year, lon, lat) %>%
  summarise(n = max(observer1))

ggplot(d, aes(year, n)) + geom_point() +
  facet_wrap(~meadow, scales = "free_y")

ggplot(d, aes(lon, lat, colour = log(n+1))) + geom_point(size = 2) +
  facet_wrap(~year)

d$lon_scaled <- d$lon * 100
d$lat_scaled <- d$lat * 100

library(rrfields)
m <- rrfield(n ~ 1, data = d, time = "year", lon = "lon_scaled", lat = "lat_scaled",
  nknots = 13, obs_error = "nb2", station = "meadow",
  estimate_df = TRUE, estimate_ar = FALSE, algorithm = "sampling",
  chains = 4, iter = 2000,
  prior_intercept = student_t(3, 0, 10), prior_gp_sigma = half_t(1000, 0, 1),
  control = list(adapt_delta = 0.95), cores = 4)
m

# mark-recapture data

abund <- read_csv("examples/alpine-butterflies/data-raw/PopulationAbunMark-Recap.csv") %>%
  rename(Meadow = Popualtion) %>%
  mutate(Date = lubridate::dmy(Date))
sites <- read_csv("examples/alpine-butterflies/data-raw/Landscape.csv")

abund <- abund %>% inner_join(sites)
names(abund) <- tolower(names(abund))

abund <- abund %>% mutate(year = lubridate::year(date)) %>%
  rename(lat = `lat. (decimal degrees n)`, lon = `lon. (decimal degrees w)`)

abund %>% group_by(meadow, year) %>% summarise(n = n())

d <- abund %>% group_by(meadow, year, lon, lat) %>%
  summarise(n = mean(`est. n`))

ggplot(d, aes(year, n)) + geom_point() +
  facet_wrap(~meadow, scales = "free_y")

ggplot(d, aes(lon, lat, colour = log(n+1))) + geom_point(size = 2) +
  facet_wrap(~year)

d$lon_scaled <- d$lon * 100
d$lat_scaled <- d$lat * 100
d$n_int <- round(d$n)

library(rrfields)
m <- rrfield(n_int ~ 1, data = d, time = "year", lon = "lon_scaled", lat = "lat_scaled",
  nknots = 12, obs_error = "nb2", station = "meadow",
  estimate_df = TRUE, estimate_ar = FALSE, fixed_ar_value = 0, algorithm = "sampling",
  chains = 2, iter = 800,
  prior_intercept = student_t(3, 0, 6), prior_gp_sigma = half_t(1000, 0, 1),
  control = list(adapt_delta = 0.95), cores = 4)
m

pairs(m$model, pars = c("gp_sigma", "gp_scale", "df[1]", "nb2_phi[1]", "B[1]", "spatialEffectsKnots[1,1]"))
