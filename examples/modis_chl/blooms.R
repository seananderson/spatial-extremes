# devtools::install_github("ropensci/rerddap")
# install.packages("rerddap")

library(dplyr)

rerddap::ed_search(query = "temp") %>% .$info

(out <- rerddap::info("osuBloomsModisChla"))

d <- list()
m <- c("02-01", "02-15", "03-01", "03-15",
  "04-01", "04-15", "05-01", "05-15",
  "06-01", "06-15", "07-01", "07-15",
  "08-01", "08-15", "09-01", "09-15")
t1 <- paste0("2014-", m, "T00:00:00Z")
t2 <- paste0("2014-", m, "T00:00:01Z")

for (i in seq_along(m)) {
  message(m[[i]])
  d[[i]] <- rerddap::griddap(out,
    time = c(t1[[i]], t2[[i]]),
    longitude = c(-126, -124),
    latitude = c(42, 48.5))$data
  d[[i]]$month <- m[[i]]
  d[[i]] <- as_data_frame(d[[i]])
}
d <- bind_rows(d)

x <- d$prc_chla+abs(min(d$prc_chla, na.rm = T))
x[x == 0] <- NA
hist(log(x))
d$prc_chla_trans <- log(x)

d <- d %>% group_by(month) %>%
  mutate(prc_chla_trans = prc_chla_trans - mean(prc_chla_trans, na.rm = TRUE)) %>%
  ungroup()

library(ggplot2)
filter(d, prc_chla_trans < 3, prc_chla_trans > -3) %>%
  ggplot(aes(lon, lat, fill = prc_chla_trans)) +
  geom_raster() +
  scale_fill_gradient2() +
  facet_wrap(~month)


# convert to UTMs
# sample grid - at ?
# fiddle with downsampling (500?), fiddle with knots (30-50)
