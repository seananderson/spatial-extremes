library(dplyr)
library(ggplot2)
library(rrfields)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

gdb_folder <- "Pacific_Northwest"
id <- "mountain-pine-beetle-pnw-raster"
nbin <- 30

d <- readRDS(paste0(id, "-", "dataframe", "-", nbin, "x", nbin, ".rds"))

# rescale
d$x = d$x/1e5
d$y = d$y/1e5

d <- dplyr::filter(d, cover > 0)

ggplot(d, aes(x, y, fill = log10(cover))) + geom_tile() +
  facet_wrap(~year) +
  scale_fill_viridis(option = "D") +
  theme_light() +
  coord_fixed()

d <- mutate(d, station = paste(x, y))
# set.seed(1)
# holdout = sample(1:nrow(catch), size=round(nrow(catch)*0.1,0))

mean(log(d$cover))

mvt_lognormal <- rrfield(log(cover) ~ -1 + as.factor(year), data = d,
  time = "year", lon = "x", lat = "y",
  station = "station",
  nknots = 30,
  obs_error = "normal",
  prior_gp_sigma = half_t(7, 0, 2),
  prior_gp_scale = half_t(7, 0, 2),
  prior_intercept = student_t(100, 0, 10),
  prior_beta = student_t(100, 0, 10),
  prior_sigma = half_t(7, 0, 5),
  estimate_ar = FALSE,
  estimate_df = TRUE,
  chains = 2, iter = 700, control = list(adapt_delta = 0.95))

saveRDS(mvt_lognormal, paste0(id, "-", "lognormal", "-", nbin, "x", nbin, ".rds"))

plot(mvt_lognormal, type = "prediction")
plot(mvt_lognormal, type = "spatial-residual")
plot(mvt_lognormal, type = "residual-vs-fitted")

pred <- predict(mvt_lognormal, interval = "prediction", conf_level = 0.8)
pred$observed <- log(d$cover)
pred$year <- d$year
ggplot(pred, aes(observed, estimate)) +
  geom_point() +
  # geom_pointrange()
  theme_light() +
  coord_fixed() +
  facet_wrap(~year) +
  geom_abline(intercept = 0, slope = 1)
