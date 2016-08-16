library(dplyr)
library(PBSmapping)
library(ggplot2)

dat = readRDS("data/triennial_annual_joined.rds")
dat$PID = 1
dat$POS = seq(1,dim(dat)[1])
dat$X = dat$haul_longitude_dd
dat$Y = dat$haul_latitude_dd
attr(dat, "projection") = "LL"
attr(dat, "zone") = 7
dat = convUL(dat)

data(nepacLL)
nepacLL = convUL(nepacLL)

group_by(dat, species) %>%
  summarize(n = length(which(!is.na(age_yrs)))) %>%
  arrange(-n)

group_by(dat, species, year) %>%
  summarize(n = length(which(!is.na(age_yrs) & age_yrs < 5))) %>%
  group_by(species) %>%
  summarize(n = length(n > 0)) %>%
    arrange(-n)


plot_species = function(x) {

p1 = ggplot(nepacLL, aes(X, Y, group = PID)) + geom_polygon() +
  geom_point(data = dat[dat$species==x,],
    aes(X, Y, colour = log(age_yrs), alpha = 0.02), size=0.5) +
  scale_colour_gradient(low = "blue") +
  coord_cartesian(xlim = c(1600,2400), ylim=c(4000,5600)) +
  xlab("Eastings") + ylab("Northings") + ggtitle(x)

p1

}

# JT's species:
# Bocaccio (5)
# canary (1)
# chilipepper (3)
# darkblotched (2)
# POP (4)
drkb = plot_species("darkblotched rockfish")
pop = plot_species("pacific ocean perch")
canary = plot_species("canary rockfish")
yellowtail = plot_species("yellowtail rockfish")
bocaccio = plot_species("bocaccio")


# There's a lot of weight but not age data from these spp -- can we boost sample
# size by doing some simple regressions?
# We have A - L - W data, so presumably can recover age from L or W

# L - W regressions work great, although may not be useful
# summary(lm(log(weight_kg) ~ log(length_cm), data = dat[dat$species=="darkblotched rockfish",]))

# ggplot(data = dat[dat$species=="darkblotched rockfish" & dat$age_yrs < 10,], aes(age_yrs, log(length_cm))) + geom_point()
# summary(lm(log(length_cm) ~ log(age_yrs), data = dat[dat$species=="darkblotched rockfish" & dat$age_yrs < 10,]))

# von B fir a number of different ways. implement it in optim, like this:
# http://rstudio-pubs-static.s3.amazonaws.com/1123_582a9931d78240c78181281fd831e0d0.html
SSQ_vonB <- function(theta) {
  Linf <- theta[1]
  K <- theta[2]
  t0 <- theta[3]
  epsilon <- rep(0, length(age))
  lpred <- rep(0, length(age))
  for (i in 1:length(age)) {
    lpred[i] <- Linf * (1 - exp(-K * (age[i] - t0)))
    epsilon[i] <- (lt[i] - lpred[i])^2
  }
  ssq <- sum(epsilon)
  return(ssq)
}

sub = dat[dat$species=="darkblotched rockfish" & !is.na(dat$age_yrs) & dat$age_yrs < 10,]
age = sub$age_yrs
lt = sub$length_cm
out <- optim(c(45, 0.2, -0.5), fn = SSQ_vonB, method = "BFGS", hessian = TRUE)

pred = data.frame(x = 1:10, y = out$par[1] * (1 - exp(-out$par[2] * (1:10 - out$par[3]))))

# This seems like predictions (predicting age from L) is highly variable
ggplot(data = sub,
  aes(age_yrs, length_cm)) + geom_point() +
  geom_line(data = pred, aes(x, y), colour = "blue")

# could prediction with rf be any better?
library(randomForest)

train = sample(seq(1, nrow(sub)), size=round(0.8*nrow(sub)))
sub$age_yrs = as.factor(sub$age_yrs)
sub$len2 = sub$length_cm ^ 2
rf = randomForest(age_yrs ~ length_cm, data=sub, type="classification")

