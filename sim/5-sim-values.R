# This code generates 'magic number' values to insert into the results section

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



write_tex(round(mean(x), 0) , "basePriorMean")

close(zz)
