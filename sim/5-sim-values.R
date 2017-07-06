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
rmse <- mutate(rmse, perc_better=(rmse_wrong-rmse)/rmse*100)
loo <- readRDS("sim/match-mismatch-loo.rds")
load("sim/match-mismatch-example.rda")

# Values for paragraph on simulation recovery testing
x0 <- filter(sim1, n_draws == 25, sd_obs == 0.1, df == 2.5)
mape0 <- median(abs(x0$df_est - x0$df)/x0$df)

x1 <- filter(sim1, n_draws == 15, sd_obs == 0.1, df == 2.5)
mape1 <- median(abs(x1$df_est - x1$df)/x1$df)

x2 <- filter(sim1, n_draws == 5, sd_obs == 0.1, df == 2.5)
mape2 <- median(abs(x2$df_est - x2$df)/x2$df)

sim1_mape0 <- round(mape0, 2)
mape1
mape2
sim1_mape_inc1 <- round((mape1 / mape0), 1)
sim1_mape_inc2 <- round((mape2 / mape0), 1)

# Values for paragraph on match mismatch simulation testing
median(rmse$rmse_wrong)
median(rmse$rmse)

write_tex(round(mean(x), 0) , "basePriorMean")

close(zz)
