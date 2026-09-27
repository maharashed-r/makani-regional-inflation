# ================================
# Formal spatial econometric model (reviewer-requested)
#
# Global and Local Moran's I establish that regional mean inflation is
# spatially structured, but they are diagnostic statistics, not models.
# This script fits actual spatial regression models -- Spatial Lag (SAR)
# and Spatial Error (SEM) -- against a non-spatial OLS baseline, to
# formalise the spatial dependence documented above.
#
# There are no additional regional covariates in this dataset (only CPI
# itself), so all three models are intercept-only (mean_inflation ~ 1).
# The quantity of interest is each spatial model's dependence coefficient
# (rho for SAR, lambda for SEM), whether it is significantly different
# from zero (likelihood-ratio test against the OLS baseline), and whether
# AIC favours the spatial specification over plain OLS. This is a minimal,
# single-model formalisation of the spatial dependence already documented
# via Moran's I -- not a full explanatory model of what drives regional
# inflation, which would require additional covariates (population, trade
# exposure, etc.) not available here.
#
# Depends on: inflation_sf, lw (04_spatial_analysis.R, 06_spatial_autocorrelation.R)
# ================================

library(spatialreg)
library(spdep)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

ols_model <- lm(mean_inflation ~ 1, data = inflation_sf)

sar_model <- lagsarlm(
  mean_inflation ~ 1, data = inflation_sf, listw = lw, zero.policy = TRUE
)
sem_model <- errorsarlm(
  mean_inflation ~ 1, data = inflation_sf, listw = lw, zero.policy = TRUE
)

sar_summary <- summary(sar_model)
sem_summary <- summary(sem_model)

# Likelihood-ratio test p-value for rho = 0 / lambda = 0. spatialreg's
# summary.sarlm normally exposes this as $LR1$p.value; fall back to
# LR.Sarlm() directly against the OLS model if that field is ever missing,
# and to NA (with a warning) rather than letting the pipeline fail.
get_lr_p <- function(model, summary_obj) {
  p <- tryCatch(summary_obj$LR1$p.value, error = function(e) NULL)
  if (is.null(p) || length(p) == 0 || is.na(p)) {
    p <- tryCatch(spatialreg::LR.Sarlm(model, ols_model)$p.value, error = function(e) NA_real_)
  }
  as.numeric(p)
}

sar_lr_p <- get_lr_p(sar_model, sar_summary)
sem_lr_p <- get_lr_p(sem_model, sem_summary)

spatial_model_comparison <- data.frame(
  Model     = c("OLS (no spatial term)", "Spatial Lag (SAR)", "Spatial Error (SEM)"),
  Parameter = c("-", "rho", "lambda"),
  Estimate  = c(NA_real_, round(unname(sar_model$rho), 3), round(unname(sem_model$lambda), 3)),
  LR_p_value = c(NA_real_, signif(sar_lr_p, 3), signif(sem_lr_p, 3)),
  AIC       = round(c(AIC(ols_model), AIC(sar_model), AIC(sem_model)), 2)
)

write.csv(spatial_model_comparison, "results/tables/spatial_model_comparison.csv", row.names = FALSE)
message("Spatial regression model comparison (OLS vs SAR vs SEM) written.")
print(spatial_model_comparison)

message("---- Full SAR summary ----")
print(sar_summary)
message("---- Full SEM summary ----")
print(sem_summary)
