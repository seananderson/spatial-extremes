# Plot the model coefficients from the 2 models

load("examples/beetles/mountain-pine-beetle-pnw-raster-mvn-lognormal-500x500.rda")
load("examples/beetles/mountain-pine-beetle-pnw-raster-mvt-lognormal-500x500.rda")

e <- extract(mvt$model)
names(mvt$model)
print(mvt)
library(dplyr)

library(broom)
posterior_t <- tidy(mvt$model, estimate.method = "median",
  conf.int = TRUE, rhat = TRUE, ess = TRUE) %>%
  filter(!(grepl("log_lik", term) | grepl("spatialEffects", term) | grepl("yearEffects", term)))

posterior_n <- tidy(mvn$model, estimate.method = "median",
  conf.int = TRUE, rhat = TRUE, ess = TRUE) %>%
  filter(!(grepl("log_lik", term) | grepl("spatialEffects", term) | grepl("yearEffects", term)))

posteriors <- bind_rows(data.frame(posterior_t, model = "MVT", stringsAsFactors = FALSE),
  data.frame(posterior_n, model = "MVN", stringsAsFactors = FALSE))
head(posteriors)

library(ggplot2)

print(posteriors)
ggplot(posteriors, aes(term, estimate, ymin = conf.low, ymax = conf.high, color = model)) +
  geom_pointrange()

