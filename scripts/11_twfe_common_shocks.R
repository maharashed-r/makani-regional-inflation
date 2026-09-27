# ================================================================
# MAKANI v2 -- Stage 1 econometrics: two-way fixed effects (TWFE)
# benchmark models for regional inflation, primary outcome inflation_mom.
#
# NEW v2 FILE. Reads only the validated, frozen v2 panel:
#   data_v2/saudi_cpi_panel_v2.csv
# Does not read, source, or modify anything under baseline/, data_clean/,
# results/, or scripts/01-10. Writes only under results_v2/.
#
# Research question at this stage (descriptive, NOT causal): after
# removing (1) time-invariant regional heterogeneity (region fixed
# effects) and (2) common month-specific national shocks (month fixed
# effects), how much of the variation in regional inflation is
# attributable to each component, and what residual is left over for
# scripts/12_residual_spatial_tests.R to test spatially?
#
# No specialised panel/FE package (fixest, plm, lfe, ...) is used --
# none is in the locked renv environment and none was added, per
# instruction not to alter the package environment to obtain a
# preferred estimator. N=13, T=157 (MoM) / T=146 (YoY) is small enough
# that plain lm() with factor() dummies, and/or a direct group-mean
# decomposition (exact for a BALANCED two-way panel, which this is), is
# both correct and computationally trivial.
#
# Run manually: Rscript scripts/11_twfe_common_shocks.R
# ================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

dir.create("results_v2/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results_v2/figures", recursive = TRUE, showWarnings = FALSE)

panel <- read.csv("data_v2/saudi_cpi_panel_v2.csv", stringsAsFactors = FALSE)
panel$date <- as.Date(panel$date)

## ================================================================
## PRIMARY OUTCOME: inflation_mom
## ================================================================

## ---- 1. Estimation sample verification --------------------------------

mom_full   <- panel
mom_sample <- panel %>% filter(!is.na(inflation_mom))

n_regions_mom <- length(unique(mom_sample$region_id))
n_months_mom  <- length(unique(mom_sample$date))
n_obs_mom     <- nrow(mom_sample)
obs_per_region_mom <- mom_sample %>% count(region_id, name = "n_obs")
balanced_mom  <- length(unique(obs_per_region_mom$n_obs)) == 1

dropped_mom <- mom_full %>% filter(is.na(inflation_mom))
n_dropped_mom <- nrow(dropped_mom)
dropped_reason_mom <- if (n_dropped_mom > 0) {
  all_first_month <- all(dropped_mom$date == min(panel$date))
  if (all_first_month) {
    paste0("All ", n_dropped_mom, " dropped rows are the structural first-month-per-region ",
           "observation (", format(min(panel$date)), "), where inflation_mom is undefined ",
           "(no prior month exists) -- expected, not a data defect.")
  } else {
    "Dropped rows do NOT all correspond to the structural first month -- investigate."
  }
} else "No rows dropped."

mom_sample_summary <- data.frame(
  n_regions          = n_regions_mom,
  n_months           = n_months_mom,
  total_obs          = n_obs_mom,
  obs_per_region_min = min(obs_per_region_mom$n_obs),
  obs_per_region_max = max(obs_per_region_mom$n_obs),
  balanced           = balanced_mom,
  first_date         = format(min(mom_sample$date)),
  last_date          = format(max(mom_sample$date)),
  n_dropped          = n_dropped_mom,
  dropped_reason     = dropped_reason_mom
)
write.csv(mom_sample_summary, "results_v2/tables/twfe_estimation_sample_mom.csv", row.names = FALSE)
message("---- MoM estimation sample ----")
message("Regions: ", n_regions_mom, " | Months: ", n_months_mom, " | Total obs: ", n_obs_mom,
        " | Balanced: ", balanced_mom, " | Dropped: ", n_dropped_mom)
message(dropped_reason_mom)

if (!balanced_mom) {
  stop("MoM estimation sample is not balanced -- this script assumes a balanced panel ",
       "for the group-mean variance decomposition below; investigate before proceeding.")
}

## ---- 2. M0-M3 model fit -------------------------------------------------
## region_id and date are both used as unordered factors (dummy variables).
## Reference-category coding (default contr.treatment) is used only for
## the model-fit statistics (R^2, AIC, BIC, F-test) -- NOT for the
## variance decomposition or the saved fitted/residual values, which use
## the exact group-mean ("effects coding") decomposition in step 3, the
## numerically clean and unambiguous choice for a balanced panel.

mom_sample$region_id_f <- factor(mom_sample$region_id)
mom_sample$date_f      <- factor(mom_sample$date)

m0 <- lm(inflation_mom ~ 1, data = mom_sample)
m1 <- lm(inflation_mom ~ region_id_f, data = mom_sample)
m2 <- lm(inflation_mom ~ date_f, data = mom_sample)
m3 <- lm(inflation_mom ~ region_id_f + date_f, data = mom_sample)

fit_stats <- function(model, label) {
  s <- summary(model)
  data.frame(
    model        = label,
    n_obs        = length(residuals(model)),
    n_params     = length(coef(model)),
    r_squared    = s$r.squared,
    adj_r_squared = s$adj.r.squared,
    resid_se     = s$sigma,
    aic          = AIC(model),
    bic          = BIC(model),
    logLik       = as.numeric(logLik(model))
  )
}

model_fit_mom <- bind_rows(
  fit_stats(m0, "M0: pooled (mu only)"),
  fit_stats(m1, "M1: region FE only"),
  fit_stats(m2, "M2: month FE only"),
  fit_stats(m3, "M3: TWFE (region + month)")
)
write.csv(model_fit_mom, "results_v2/tables/twfe_model_fit_mom.csv", row.names = FALSE)
message("\n---- M0-M3 model fit (MoM) ----")
print(model_fit_mom)

## ---- 3. Variance decomposition (exact, balanced two-way ANOVA identity) --
## For a fully balanced N x T panel (verified above), the classical
## two-way decomposition is EXACT, not approximate:
##   pi_it = grand_mean + alpha_i + tau_t + resid_it
## where alpha_i = region_i mean - grand_mean, tau_t = month_t mean -
## grand_mean, and SS_total = SS_region + SS_month + SS_residual exactly
## (no interaction term, no omitted cross-term, verified numerically
## below). This is computed directly from group means rather than read
## off the lm() dummy-coefficient parameterisation, to avoid any
## reference-category ambiguity in what "the region effect" means.

grand_mean_mom <- mean(mom_sample$inflation_mom)

region_means_mom <- mom_sample %>%
  group_by(region_id, Region) %>%
  summarise(region_mean = mean(inflation_mom), .groups = "drop") %>%
  mutate(alpha_i = region_mean - grand_mean_mom)

month_means_mom <- mom_sample %>%
  group_by(date) %>%
  summarise(month_mean = mean(inflation_mom), .groups = "drop") %>%
  mutate(tau_t = month_mean - grand_mean_mom)

mom_decomp <- mom_sample %>%
  left_join(region_means_mom %>% select(region_id, alpha_i), by = "region_id") %>%
  left_join(month_means_mom  %>% select(date, tau_t),        by = "date") %>%
  mutate(
    fitted_twfe   = grand_mean_mom + alpha_i + tau_t,
    twfe_residual = inflation_mom - fitted_twfe
  )

ss_total_mom    <- sum((mom_decomp$inflation_mom - grand_mean_mom)^2)
ss_region_mom   <- n_months_mom  * sum(region_means_mom$alpha_i^2)
ss_month_mom    <- n_regions_mom * sum(month_means_mom$tau_t^2)
ss_residual_mom <- sum(mom_decomp$twfe_residual^2)

additivity_check_mom <- abs((ss_region_mom + ss_month_mom + ss_residual_mom) - ss_total_mom)
if (additivity_check_mom / ss_total_mom > 1e-8) {
  stop("Variance decomposition additivity check failed for MoM: SS_region + SS_month + ",
       "SS_residual does not equal SS_total within tolerance (diff = ", additivity_check_mom, "). ",
       "The balanced-panel identity this decomposition relies on may not hold -- investigate ",
       "before trusting the shares below.")
}

variance_decomp_mom <- data.frame(
  component     = c("Region (alpha_i)", "Month (tau_t)", "Residual", "Total"),
  sum_sq        = c(ss_region_mom, ss_month_mom, ss_residual_mom, ss_total_mom),
  share_of_total = c(ss_region_mom, ss_month_mom, ss_residual_mom, ss_total_mom) / ss_total_mom
)
write.csv(variance_decomp_mom, "results_v2/tables/twfe_variance_decomposition_mom.csv", row.names = FALSE)
message("\n---- Variance decomposition (MoM, exact for this balanced panel) ----")
message("Additivity check (should be ~0): ", signif(additivity_check_mom, 3))
print(variance_decomp_mom)
message(
  "\nNOTE (interpretation, not causal): this decomposition partitions variance, it does not ",
  "identify causal effects. A large month-effect share means regional inflation co-moves ",
  "strongly with common national-level shocks in a given month; it does not by itself imply ",
  "any region caused inflation in another."
)

## Region FE and month FE tables, saved separately for later plotting (fig B).
region_fe_mom <- region_means_mom %>% arrange(desc(alpha_i))
write.csv(region_fe_mom, "results_v2/tables/twfe_region_fe_mom.csv", row.names = FALSE)

month_fe_mom <- month_means_mom %>% arrange(date)
write.csv(month_fe_mom, "results_v2/tables/twfe_month_fe_mom.csv", row.names = FALSE)

## ---- 4. Save TWFE (M3) residuals, with validation -----------------------

residuals_mom <- mom_decomp %>%
  transmute(region_id, Region, date, inflation_mom, fitted = fitted_twfe, twfe_residual) %>%
  arrange(region_id, date)

write.csv(residuals_mom, "results_v2/tables/twfe_residuals_mom.csv", row.names = FALSE)

## Validation
resid_checks <- list()
resid_checks$n_residuals_matches_sample <- nrow(residuals_mom) == n_obs_mom
resid_checks$no_unexpected_missing      <- sum(is.na(residuals_mom$twfe_residual)) == 0

overall_mean_resid <- mean(residuals_mom$twfe_residual)
resid_checks$overall_mean_near_zero <- abs(overall_mean_resid) < 1e-8

by_month_mean <- residuals_mom %>% group_by(date) %>%
  summarise(mean_resid = mean(twfe_residual), .groups = "drop")
resid_checks$month_means_near_zero <- all(abs(by_month_mean$mean_resid) < 1e-8)

by_region_mean <- residuals_mom %>% group_by(region_id) %>%
  summarise(mean_resid = mean(twfe_residual), .groups = "drop")
resid_checks$region_means_near_zero <- all(abs(by_region_mean$mean_resid) < 1e-8)

resid_validation <- data.frame(
  check = names(resid_checks),
  pass  = unlist(resid_checks)
)
write.csv(resid_validation, "results_v2/tables/twfe_residual_validation_mom.csv", row.names = FALSE)

message("\n---- TWFE residual validation (MoM) ----")
message("n residuals = ", nrow(residuals_mom), " (expected ", n_obs_mom, ")")
message("overall mean residual = ", format(overall_mean_resid, scientific = TRUE))
message("max |monthly mean residual| = ", format(max(abs(by_month_mean$mean_resid)), scientific = TRUE))
message("max |regional mean residual| = ", format(max(abs(by_region_mean$mean_resid)), scientific = TRUE))
print(resid_validation)
if (!all(unlist(resid_checks))) {
  stop("TWFE residual validation failed -- see results_v2/tables/twfe_residual_validation_mom.csv")
}

## ---- 5. Figure A: regional monthly inflation over time -------------------

p_a <- ggplot(mom_sample, aes(x = date, y = inflation_mom, color = Region)) +
  geom_line(linewidth = 0.3, alpha = 0.8) +
  labs(title = "Regional Monthly Inflation (MoM) Over Time",
       subtitle = "13 Saudi regions, Feb 2013 - Feb 2026",
       x = NULL, y = "MoM inflation (%)") +
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_blank())
ggsave("results_v2/figures/fig_A_regional_inflation_mom_over_time.png", p_a,
       width = 10, height = 5.5, dpi = 300)

## ---- 6. Figure B: estimated common month effect tau_t over time ---------

p_b <- ggplot(month_fe_mom, aes(x = date, y = tau_t)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(color = "#2c3e50") +
  labs(title = "Estimated Common Month Effect (tau_t) Over Time",
       subtitle = "M3 two-way FE model, MoM inflation -- additive shock common to all 13 regions in a given month",
       x = NULL, y = expression(tau[t])) +
  theme_minimal()
ggsave("results_v2/figures/fig_B_month_fe_tau_over_time.png", p_b,
       width = 10, height = 5, dpi = 300)

message("\nWrote MoM TWFE outputs to results_v2/tables/ and figures A, B to results_v2/figures/.")

## ================================================================
## SECONDARY ROBUSTNESS OUTCOME: inflation_yoy
## Kept entirely separate from the MoM primary inference above -- no
## table mixes MoM and YoY results. Same TWFE framework (M0-M3 fit +
## exact balanced-panel variance decomposition); the monthly spatial
## Moran analysis in scripts/12 is run on the MoM residuals only (the
## primary specification), not repeated here.
## ================================================================

yoy_sample <- panel %>% filter(!is.na(inflation_yoy))

n_regions_yoy <- length(unique(yoy_sample$region_id))
n_months_yoy  <- length(unique(yoy_sample$date))
n_obs_yoy     <- nrow(yoy_sample)
obs_per_region_yoy <- yoy_sample %>% count(region_id, name = "n_obs")
balanced_yoy  <- length(unique(obs_per_region_yoy$n_obs)) == 1

dropped_yoy <- panel %>% filter(is.na(inflation_yoy))
n_dropped_yoy <- nrow(dropped_yoy)
first_12_months <- sort(unique(panel$date))[1:12]
dropped_reason_yoy <- if (n_dropped_yoy > 0) {
  all_first12 <- all(dropped_yoy$date %in% first_12_months) &&
    all(sapply(split(dropped_yoy$date, dropped_yoy$region_id), function(d) setequal(d, first_12_months)))
  if (all_first12) {
    paste0("All ", n_dropped_yoy, " dropped rows are each region's structural first-12-months ",
           "(", format(min(first_12_months)), " to ", format(max(first_12_months)), "), where ",
           "inflation_yoy is undefined (no observation 12 months earlier exists) -- expected.")
  } else {
    "Dropped rows do NOT all correspond to the structural first 12 months -- investigate."
  }
} else "No rows dropped."

yoy_sample_summary <- data.frame(
  n_regions = n_regions_yoy, n_months = n_months_yoy, total_obs = n_obs_yoy,
  obs_per_region_min = min(obs_per_region_yoy$n_obs), obs_per_region_max = max(obs_per_region_yoy$n_obs),
  balanced = balanced_yoy, first_date = format(min(yoy_sample$date)), last_date = format(max(yoy_sample$date)),
  n_dropped = n_dropped_yoy, dropped_reason = dropped_reason_yoy
)
write.csv(yoy_sample_summary, "results_v2/tables/twfe_estimation_sample_yoy.csv", row.names = FALSE)
message("\n---- YoY estimation sample (secondary robustness) ----")
message("Regions: ", n_regions_yoy, " | Months: ", n_months_yoy, " | Total obs: ", n_obs_yoy,
        " | Balanced: ", balanced_yoy, " | Dropped: ", n_dropped_yoy)

if (balanced_yoy) {
  yoy_sample$region_id_f <- factor(yoy_sample$region_id)
  yoy_sample$date_f      <- factor(yoy_sample$date)

  y0 <- lm(inflation_yoy ~ 1, data = yoy_sample)
  y1 <- lm(inflation_yoy ~ region_id_f, data = yoy_sample)
  y2 <- lm(inflation_yoy ~ date_f, data = yoy_sample)
  y3 <- lm(inflation_yoy ~ region_id_f + date_f, data = yoy_sample)

  model_fit_yoy <- bind_rows(
    fit_stats(y0, "M0: pooled (mu only)"),
    fit_stats(y1, "M1: region FE only"),
    fit_stats(y2, "M2: month FE only"),
    fit_stats(y3, "M3: TWFE (region + month)")
  )
  write.csv(model_fit_yoy, "results_v2/tables/twfe_model_fit_yoy.csv", row.names = FALSE)

  grand_mean_yoy <- mean(yoy_sample$inflation_yoy)
  region_means_yoy <- yoy_sample %>% group_by(region_id, Region) %>%
    summarise(region_mean = mean(inflation_yoy), .groups = "drop") %>%
    mutate(alpha_i = region_mean - grand_mean_yoy)
  month_means_yoy <- yoy_sample %>% group_by(date) %>%
    summarise(month_mean = mean(inflation_yoy), .groups = "drop") %>%
    mutate(tau_t = month_mean - grand_mean_yoy)
  yoy_decomp <- yoy_sample %>%
    left_join(region_means_yoy %>% select(region_id, alpha_i), by = "region_id") %>%
    left_join(month_means_yoy  %>% select(date, tau_t),        by = "date") %>%
    mutate(fitted_twfe = grand_mean_yoy + alpha_i + tau_t,
           twfe_residual = inflation_yoy - fitted_twfe)

  ss_total_yoy    <- sum((yoy_decomp$inflation_yoy - grand_mean_yoy)^2)
  ss_region_yoy   <- n_months_yoy  * sum(region_means_yoy$alpha_i^2)
  ss_month_yoy    <- n_regions_yoy * sum(month_means_yoy$tau_t^2)
  ss_residual_yoy <- sum(yoy_decomp$twfe_residual^2)
  additivity_check_yoy <- abs((ss_region_yoy + ss_month_yoy + ss_residual_yoy) - ss_total_yoy)

  variance_decomp_yoy <- data.frame(
    component = c("Region (alpha_i)", "Month (tau_t)", "Residual", "Total"),
    sum_sq    = c(ss_region_yoy, ss_month_yoy, ss_residual_yoy, ss_total_yoy),
    share_of_total = c(ss_region_yoy, ss_month_yoy, ss_residual_yoy, ss_total_yoy) / ss_total_yoy
  )
  write.csv(variance_decomp_yoy, "results_v2/tables/twfe_variance_decomposition_yoy.csv", row.names = FALSE)

  message("YoY additivity check (should be ~0): ", signif(additivity_check_yoy, 3))
  message("\n---- M0-M3 model fit (YoY, secondary robustness) ----")
  print(model_fit_yoy)
  message("\n---- Variance decomposition (YoY, secondary robustness) ----")
  print(variance_decomp_yoy)
} else {
  warning("YoY sample is not balanced -- skipping YoY TWFE (secondary robustness only, ",
          "does not affect the primary MoM analysis above).")
}

message("\n11_twfe_common_shocks.R complete.")
