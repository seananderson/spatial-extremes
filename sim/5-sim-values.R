# This code generates 'magic number' values to insert into the results section

library(dplyr)

write_tex <- function(x, macro, ...) {
  out <- paste0("\\newcommand{\\", macro, "}{", x, "}")
  cat(out, file = zz)
  cat("\n", file = zz)
}
zz <- file("text/sim-values.tex", "w") # open .tex file to write to throughout

sim1 <- readRDS("sim/mvt-norm-sim-testing2.rds")
rmse <- readRDS("sim/match-mismatch-rmse.rds")
loo <- readRDS("sim/match-mismatch-loo.rds")
load("sim/match-mismatch-example.rda")

# Values for paragraph on simulation recovery testing
x0 <- filter(sim1, n_draws == 25, sd_obs == 0.1, df == 2.5)
mape0 <- median(abs(x0$df_est - x0$df)/x0$df)
write_tex(round(mape0, 2), "mapeTwentyFive")

x1 <- filter(sim1, n_draws == 15, sd_obs == 0.1, df == 2.5)
mape1 <- median(abs(x1$df_est - x1$df)/x1$df)
write_tex(round(mape1/mape0, 0), "mapeFifteenIncFold")

x2 <- filter(sim1, n_draws == 5, sd_obs == 0.1, df == 2.5)
mape2 <- median(abs(x2$df_est - x2$df)/x2$df)
write_tex(round(mape2/mape0, 0), "mapeFiveIncFold")

# Values for paragraph on match mismatch simulation testing
write_tex(round(median(rmse$rmse_wrong), 2), "rmseWrong")
write_tex(round(median(rmse$rmse), 2), "rmseRight")

# rmse
write_tex(round(median(rmse$perc_better_ci), 0), "medianMedianCIWiderSim")

write_tex(round(100 - sum(!loo$delta_loo < 0)/nrow(loo) * 100, 0), "looCorrectSim")

# Beetle

load("examples/beetles/mountain-pine-beetle-pnw-raster-mvn-lognormal-500x500.rda")
load("examples/beetles/mountain-pine-beetle-pnw-raster-mvt-lognormal-500x500.rda")
et <- rstan::extract(mvt$model)
en <- rstan::extract(mvn$model)
write_tex(sprintf("%0.1f", round(median(et$df), 1)), "nuBeetleMedian")
write_tex(sprintf("%0.1f", round(quantile(et$df, probs = 0.025)[[1]], 1)), "nuBeetleLower")
write_tex(sprintf("%0.1f", round(quantile(et$df, probs = 0.975)[[1]], 1)), "nuBeetleUpper")

cis <- readRDS("examples/beetles/beetle-cis.rds")
write_tex(round((mean(cis$conf_width_ratio) - 1) * 100, 0), "meanPercLargerCIsBeetles")
write_tex(round((1 - median(1/cis$conf_width_ratio)) * 100, 0), "medianPercSmallerCIsBeetles")

close(zz)
