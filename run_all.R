# Runs the full pipeline in order. Each step reads the previous steps'
# outputs from data/, so steps can also be re-run individually with source().
# Expensive results (MC critical values, GARCH fits/forecasts, ML forecasts)
# are cached in data/; delete the .rds file to force re-estimation.

steps <- c("01_data.R", "02_bubble_detection.R", "03_regime_distributions.R",
           "04_garch.R", "05_ml_quantile.R", "06_backtesting.R", "07_early_warning.R")

for (s in steps) {
  message("\n==== ", s, " ====")
  source(s, local = new.env())
}
