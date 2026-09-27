rm(list = ls())

cat("🚀 Starting Project...\n")

source("scripts/01_load_data.R")
cat("✔ Data loaded\n")

source("scripts/02_clean_data.R")
cat("✔ Data cleaned\n")

source("scripts/03_exploratory_analysis.R")
cat("✔ Exploratory analysis done\n")

source("scripts/04_spatial_analysis.R")
cat("✔ Spatial analysis done\n")

source("scripts/06_spatial_autocorrelation.R")
cat("✔ Spatial autocorrelation (Moran's I) done\n")

source("scripts/07_sensitivity_analysis.R")
cat("✔ Sensitivity analysis (VAT shocks, weights, LISA) done\n")

# NOTE: scripts/08_spatial_regression.R (SAR/SEM) is intentionally NOT run
# as part of the main pipeline. With no regional covariates available,
# those models are intercept-only and reduce to a statistical exercise
# rather than an explanatory model (see paper.qmd, Sensitivity Analysis
# section). The script is kept on disk for future use once regional
# covariates are added. Run it manually if needed:
#   source("scripts/08_spatial_regression.R")

source("scripts/05_plots.R")
cat("✔ Plots generated\n")

cat("🎉 All done! Check results/ folder\n")
