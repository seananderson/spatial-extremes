library(tidyverse)
d <- readr::read_csv("examples/moths/moth-dat.csv") %>%
  left_join(readr::read_csv("examples/moths/moth-sites.csv"))
d <- d %>% mutate(easting = as.numeric(easting)/10, northing = as.numeric(northing)/10)
stopifnot({qq <- filter(d, is.na(easting));length(unique(qq$site)) == 0})

group_by(d, species) %>% summarize(n = n()) %>% arrange(n)

filter(d, species == "Xestia xanthographa") %>%
  ggplot(aes(easting, northing, colour = log(number+1))) + geom_point() +
  facet_wrap(~year)

library(rrfields)
m <- rrfield(number ~ 1, data = filter(d, species == "Xestia xanthographa"),
  time = "year", lon = "easting", lat = "northing",
  station = "site", estimate_ar = FALSE, estimate_df = TRUE, fixed_ar_value = 1,
  obs_error = "nb2",
  iter = 2000, chains = 4, cores = 4,
  prior_gp_sigma = half_t(1000, 0, 1), prior_intercept = student_t(1000, 0, 15),
  prior_sigma = half_t(1000, 0, 1), prior_gp_scale = half_t(1000, 0, 10), nknots = 10L,
  control = list(adapt_delta = 0.95))
m

plot(m)

library(bayesplot)
color_scheme_set("blue")
posterior <- rstan::extract(m$model, inc_warmup = FALSE, permuted = FALSE)
pars <- c("df[1]", "gp_sigma", "nb2_phi[1]", "gp_scale", "B[1]")
mcmc_trace(posterior,  pars = pars)

obs <- m$data$number
p <- predict(m)
plot(obs, p$estimate)
cor(obs, log(p$estimate))
segments(obs, p$conf_low, obs, p$conf_high, lwd = 0.5, col = "#00000020")
abline(a = 0, b = 1)
