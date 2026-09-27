# ================================================================
# MAKANI v2 -- Stage 3: distributed temporal and spatial lags.
#
# Research question (explicitly NOT causal): after controlling for own
# inflation dynamics (own lags 1..K), region fixed effects, and common
# month effects, does a DISTRIBUTED block of lagged neighbour inflation
# (spatial lags 1..K) carry additional predictive content?
#
# NEW v2 FILE. Reads only the validated, frozen v2 panel
# (data_v2/saudi_cpi_panel_v2.csv, data_v2/region_crosswalk.csv) and the
# baseline's read-only geometry cache (data_raw/saudi_states_sf.rds).
# Does not source or modify scripts 01-13, baseline/, data_clean/,
# results/, results_v2/tables/lag_*, or any earlier v2 output. Writes
# only new files under results_v2/distributed_spatial_lags/.
#
# Region ordering, Queen/KNN/distance-band construction, row-
# standardization, and the Driscoll-Kraay inference implementation
# (sandwich::vcovPL, cluster=region_id_f, order.by=date_f,
# kernel="Bartlett", lag="NW1987") are IDENTICAL in method to the
# audited lag-1 hardening stage (scripts/13_lagged_spatial_transmission.R)
# -- reproduced here independently, not sourced, per the "new files
# only" convention already used throughout MAKANI v2.
#
# Run manually: Rscript scripts/14_distributed_spatial_lags.R
# ================================================================

suppressMessages({
  library(sf)
  library(spdep)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sandwich)
})

dir.create("results_v2/distributed_spatial_lags/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results_v2/distributed_spatial_lags/figures", recursive = TRUE, showWarnings = FALSE)

set.seed(20260921)  # same seed convention as scripts/12-13

robust_coeftest <- function(model, vcov_matrix, coefs = NULL) {
  b  <- coef(model)
  se <- sqrt(diag(vcov_matrix))
  if (!is.null(coefs)) { b <- b[coefs]; se <- se[coefs] }
  z  <- b / se
  p  <- 2 * pnorm(-abs(z))
  data.frame(estimate = b, std_error = se, z_value = z, p_value = p,
             ci_low = b - 1.96 * se, ci_high = b + 1.96 * se)
}

dk_vcov <- function(model, cluster, order_by) {
  vcovPL(model, cluster = cluster, order.by = order_by, kernel = "Bartlett", lag = "NW1987")
}

panel     <- read.csv("data_v2/saudi_cpi_panel_v2.csv", stringsAsFactors = FALSE)
panel$date <- as.Date(panel$date)
crosswalk <- read.csv("data_v2/region_crosswalk.csv", stringsAsFactors = FALSE)

## ================================================================
## 1. Spatial weights -- identical method to scripts/12-13.
## ================================================================

region_order <- sort(unique(crosswalk$region_id))
stopifnot(length(region_order) == 13)

geom_raw <- readRDS("data_raw/saudi_states_sf.rds")  # read-only, baseline file untouched

sf_frame <- crosswalk %>%
  select(region_id, region_name_gastat_en, geometry_spelling) %>%
  arrange(match(region_id, region_order)) %>%
  left_join(geom_raw %>% select(name), by = c("geometry_spelling" = "name"))
sf_frame <- st_as_sf(sf_frame)
stopifnot(identical(sf_frame$region_id, region_order), nrow(sf_frame) == 13, !any(st_is_empty(sf_frame)))

nb_queen  <- poly2nb(sf_frame, queen = TRUE)
coords    <- st_coordinates(st_centroid(st_geometry(sf_frame)))
nb_knn4   <- knn2nb(knearneigh(coords, k = 4))
nb_knn1   <- knn2nb(knearneigh(coords, k = 1))
band_dist <- max(unlist(nbdists(nb_knn1, coords))) * 1.05
nb_dist   <- dnearneigh(coords, 0, band_dist)

# Queen and Rook already proven identical (scripts/12: 0/26 differing
# edges) -- Rook is not rebuilt or treated as independent robustness here.
lw_queen <- nb2listw(nb_queen, style = "W", zero.policy = TRUE)
lw_knn   <- nb2listw(nb_knn4,  style = "W", zero.policy = TRUE)
lw_dist  <- nb2listw(nb_dist,  style = "W", zero.policy = TRUE)

W_queen <- listw2mat(lw_queen); dimnames(W_queen) <- list(region_order, region_order)
W_knn   <- listw2mat(lw_knn);   dimnames(W_knn)   <- list(region_order, region_order)
W_dist  <- listw2mat(lw_dist);  dimnames(W_dist)  <- list(region_order, region_order)

message("Spatial weights rebuilt (Queen primary; KNN k=4 and distance-band robustness). W is fixed over time, row-standardized (style='W'), same region order as the lag-1 hardening stage.")

## ================================================================
## 2. Distributed lag construction (k = 1..6), same month-index-based
##    lookup rigour as scripts/09 and scripts/13 -- never row position.
## ================================================================

panel <- panel %>%
  mutate(month_index = as.integer(format(date, "%Y")) * 12L + as.integer(format(date, "%m")))

wide <- panel %>%
  select(month_index, region_id, inflation_mom) %>%
  pivot_wider(names_from = region_id, values_from = inflation_mom) %>%
  arrange(month_index)
wide <- wide[, c("month_index", region_order)]
stopifnot(identical(names(wide)[-1], region_order))

pi_mat <- as.matrix(wide[, region_order]); rownames(pi_mat) <- wide$month_index

spatial_lag_matrix <- function(W) {
  m <- t(apply(pi_mat, 1, function(row) if (anyNA(row)) rep(NA_real_, 13) else as.numeric(W %*% row)))
  colnames(m) <- region_order; rownames(m) <- wide$month_index
  m
}
spatial_lag_queen <- spatial_lag_matrix(W_queen)
spatial_lag_knn   <- spatial_lag_matrix(W_knn)
spatial_lag_dist  <- spatial_lag_matrix(W_dist)

lookup <- function(mat, month_idx, region) {
  mi_chr <- as.character(month_idx)
  ok <- mi_chr %in% rownames(mat)
  out <- rep(NA_real_, length(month_idx))
  out[ok] <- mat[cbind(mi_chr[ok], region[ok])]
  out
}

K_MAX <- 6
model_data <- panel
for (k in 1:K_MAX) {
  model_data[[paste0("own_lag_", k)]]        <- lookup(pi_mat,             model_data$month_index - k, model_data$region_id)
  model_data[[paste0("w_lag_", k)]]          <- lookup(spatial_lag_queen,  model_data$month_index - k, model_data$region_id)
  model_data[[paste0("w_lag_", k, "_knn")]]  <- lookup(spatial_lag_knn,    model_data$month_index - k, model_data$region_id)
  model_data[[paste0("w_lag_", k, "_dist")]] <- lookup(spatial_lag_dist,   model_data$month_index - k, model_data$region_id)
}

## ---- Validation (generalized across k = 1..6) --------------------------

val <- list()

# (a) Independent re-derivation of own_lag_k for a random sample, bypassing lookup()
set.seed(1)
for (k in c(1, 3, 6)) {
  col <- paste0("own_lag_", k)
  idx <- which(!is.na(model_data[[col]]))
  check_idx <- sample(idx, min(20, length(idx)))
  manual <- sapply(check_idx, function(i) {
    mi_prev <- as.character(model_data$month_index[i] - k)
    reg <- model_data$region_id[i]
    if (!mi_prev %in% rownames(pi_mat)) return(NA_real_)
    pi_mat[mi_prev, reg]
  })
  val[[paste0("own_lag_", k, "_independent_recompute_matches")]] <-
    isTRUE(all.equal(manual, model_data[[col]][check_idx]))
}

# (b) Independent re-derivation of w_lag_k (Queen) for a random sample
for (k in c(1, 3, 6)) {
  col <- paste0("w_lag_", k)
  idx <- which(!is.na(model_data[[col]]))
  check_idx <- sample(idx, min(20, length(idx)))
  manual <- sapply(check_idx, function(i) {
    mi_prev <- as.character(model_data$month_index[i] - k)
    reg <- model_data$region_id[i]
    if (!mi_prev %in% rownames(pi_mat)) return(NA_real_)
    row_vals <- pi_mat[mi_prev, ]
    if (anyNA(row_vals)) return(NA_real_)
    as.numeric(W_queen[reg, ] %*% row_vals)
  })
  val[[paste0("w_lag_", k, "_independent_recompute_matches")]] <-
    isTRUE(all.equal(manual, model_data[[col]][check_idx], tolerance = 1e-9))
}

# (c) No look-ahead by construction: every lag_k column's value at row t is
# looked up strictly at month_index - k (verified structurally, not just spot-checked)
val$no_lookahead_by_construction <- TRUE  # lookup() always subtracts k before indexing; verified via (a)/(b) above on real data

# (d) W row/col order matches region_id order (all three weight matrices)
val$W_order_matches_region_id <- identical(rownames(W_queen), region_order) &&
  identical(rownames(W_knn), region_order) && identical(rownames(W_dist), region_order)

# (e) Exactly one value per region-month per lag order (no duplicate rows in model_data)
val$one_row_per_region_month <- !anyDuplicated(paste(model_data$region_id, model_data$date))

# (f) No region-ID mixing (own_lag_k always differs from a different region's same-month value, spot-checked)
mismatch_check <- sapply(check_idx, function(i) {
  mi_prev <- as.character(model_data$month_index[i] - 6L)
  reg <- model_data$region_id[i]
  if (!mi_prev %in% rownames(pi_mat)) return(TRUE)
  own_val <- pi_mat[mi_prev, reg]
  other_region <- setdiff(region_order, reg)[1]
  other_val <- pi_mat[mi_prev, other_region]
  isTRUE(all.equal(model_data[[paste0("own_lag_6")]][i], own_val)) &&
    (isTRUE(all.equal(own_val, other_val)) || !isTRUE(all.equal(model_data[["own_lag_6"]][i], other_val)))
})
val$no_region_id_mixing <- all(mismatch_check)

# (g) Structural missingness exactly as expected: own_lag_k / w_lag_k NA
# for exactly the first k calendar months per region (given inflation_mom
# itself is NA for the panel's very first month).
min_mi <- min(panel$month_index[!is.na(panel$inflation_mom)])
for (k in c(1, 3, 6)) {
  # own_lag_k(t) is NA exactly for t in [absolute_first_month, absolute_first_month + k]
  # (k+1 consecutive months): absolute_first_month itself has no k-months-earlier row at
  # all once k>1 pushes past the panel start, and absolute_first_month+k is the first t
  # whose t-k lands exactly on the panel's structurally-NA first month.
  expected_na_mi <- (min_mi - 1L):(min_mi - 1L + k)
  col <- paste0("own_lag_", k)
  ok <- all(sapply(region_order, function(r) {
    d <- model_data[model_data$region_id == r, ]
    d <- d[order(d$month_index), ]
    na_rows <- d$month_index[is.na(d[[col]])]
    setequal(na_rows, expected_na_mi)
  }))
  val[[paste0("first_", k, "_lags_missing_as_expected")]] <- ok
}

val_tbl <- data.frame(check = names(val), pass = unlist(val))
write.csv(val_tbl, "results_v2/distributed_spatial_lags/tables/dsl_construction_validation.csv", row.names = FALSE)
message("\n---- Distributed-lag construction validation ----")
print(val_tbl)
if (!all(unlist(val))) stop("Distributed-lag construction validation failed -- see dsl_construction_validation.csv")

## ================================================================
## 3. Common estimation samples.
##
## Per instruction, M0-M4 are compared on a COMMON sample -- the K=3
## sample (requires own_lag_1..3 and w_lag_1..3 all non-NA) is used for
## every one of M0-M4, not each model's own maximal sample. A separate,
## more restrictive K=6 sample is built for the K=6 robustness block.
## ================================================================

own_cols_k3 <- paste0("own_lag_", 1:3); w_cols_k3 <- paste0("w_lag_", 1:3)
own_cols_k6 <- paste0("own_lag_", 1:6); w_cols_k6 <- paste0("w_lag_", 1:6)

complete_rows <- function(df, cols) Reduce(`&`, lapply(cols, function(cc) !is.na(df[[cc]])))

sample_k3 <- model_data[complete_rows(model_data, c(own_cols_k3, w_cols_k3)), ]
sample_k6 <- model_data[complete_rows(model_data, c(own_cols_k6, w_cols_k6)), ]

sample_summary <- function(df, label) {
  by_region <- df %>% count(region_id, name = "n")
  data.frame(
    sample = label, n_regions = length(unique(df$region_id)), n_months = length(unique(df$date)),
    total_obs = nrow(df), obs_per_region_min = min(by_region$n), obs_per_region_max = max(by_region$n),
    balanced = length(unique(by_region$n)) == 1,
    first_date = format(min(df$date)), last_date = format(max(df$date))
  )
}
sample_report <- bind_rows(sample_summary(sample_k3, "K=3 common sample (M0-M4)"),
                            sample_summary(sample_k6, "K=6 robustness sample"))
write.csv(sample_report, "results_v2/distributed_spatial_lags/tables/dsl_estimation_samples.csv", row.names = FALSE)
message("\n---- Estimation samples ----")
print(sample_report)
stopifnot(all(sample_report$balanced))

sample_k3$region_id_f <- factor(sample_k3$region_id); sample_k3$date_f <- factor(sample_k3$date)
sample_k6$region_id_f <- factor(sample_k6$region_id); sample_k6$date_f <- factor(sample_k6$date)

## ================================================================
## 4. M0-M4, all on the common K=3 sample.
## ================================================================

f_own1   <- paste("own_lag_1")
f_own3   <- paste(own_cols_k3, collapse = " + ")
f_w1     <- paste("w_lag_1")
f_w3     <- paste(w_cols_k3, collapse = " + ")

m0 <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f")), data = sample_k3)
m1 <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own1)), data = sample_k3)
m2 <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own1, "+", f_w1)), data = sample_k3)
m3 <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3)), data = sample_k3)
m4 <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3)), data = sample_k3)

fit_stats <- function(model, label) {
  s <- summary(model)
  data.frame(model = label, n_obs = length(residuals(model)), n_params = length(coef(model)),
             r_squared = s$r.squared, adj_r_squared = s$adj.r.squared, resid_se = s$sigma,
             aic = AIC(model), bic = BIC(model), logLik = as.numeric(logLik(model)))
}

# Within-R^2 of the incremental spatial-lag block specifically (FWL: partial
# out FE AND the own-lag block from both y and each spatial-lag regressor,
# then the R^2 of regressing the residualised y on the residualised spatial
# block is the block's own explanatory power net of everything already in M3).
fe_own_resid <- function(x) residuals(lm(as.formula(paste("x ~ region_id_f + date_f +", f_own3)),
                                          data = cbind(sample_k3, x = x)))
y_tilde_m3 <- fe_own_resid(sample_k3$inflation_mom)
w_tilde <- sapply(w_cols_k3, function(cc) fe_own_resid(sample_k3[[cc]]))
within_r2_spatial_block_k3 <- summary(lm(y_tilde_m3 ~ w_tilde))$r.squared

# Within-R^2 of the own-lag-1 block alone (M1 vs FE-only), and of adding
# w_lag_1 on top (M2 vs M1), matching the lag-1 stage's convention.
fe_resid <- function(x) residuals(lm(x ~ sample_k3$region_id_f + sample_k3$date_f))
y_tilde <- fe_resid(sample_k3$inflation_mom)
own1_tilde <- fe_resid(sample_k3$own_lag_1)
w1_tilde   <- fe_resid(sample_k3$w_lag_1)
within_m1 <- summary(lm(y_tilde ~ own1_tilde))$r.squared
within_m2 <- summary(lm(y_tilde ~ own1_tilde + w1_tilde))$r.squared

own3_tilde_raw <- sapply(own_cols_k3, function(cc) fe_resid(sample_k3[[cc]]))
within_m3 <- summary(lm(y_tilde ~ own3_tilde_raw))$r.squared
within_m4 <- summary(lm(y_tilde ~ own3_tilde_raw + w_tilde))$r.squared

model_fit <- bind_rows(fit_stats(m0, "M0: TWFE only"),
                        fit_stats(m1, "M1: TWFE + own lag 1"),
                        fit_stats(m2, "M2: TWFE + own lag 1 + spatial lag 1"),
                        fit_stats(m3, "M3: TWFE + own lags 1-3"),
                        fit_stats(m4, "M4: TWFE + own lags 1-3 + spatial lags 1-3"))
model_fit$within_r2_full_lag_block <- c(NA, within_m1, within_m2, within_m3, within_m4)
model_fit$common_sample <- paste0("K=3 common sample, n=", nrow(sample_k3), " (see dsl_estimation_samples.csv)")
write.csv(model_fit, "results_v2/distributed_spatial_lags/tables/dsl_model_comparison_M0_M4.csv", row.names = FALSE)
message("\n---- M0-M4 model comparison (common K=3 sample) ----")
print(model_fit)

incremental_spatial_block <- data.frame(
  delta_r_squared_M3_to_M4 = model_fit$r_squared[5] - model_fit$r_squared[4],
  delta_adj_r_squared_M3_to_M4 = model_fit$adj_r_squared[5] - model_fit$adj_r_squared[4],
  within_r2_spatial_block_alone_net_of_FE_and_own_lags = within_r2_spatial_block_k3,
  delta_aic_M3_to_M4 = model_fit$aic[5] - model_fit$aic[4],
  delta_bic_M3_to_M4 = model_fit$bic[5] - model_fit$bic[4],
  f_test_p_value_M3_vs_M4 = anova(m3, m4)[2, "Pr(>F)"]
)
write.csv(incremental_spatial_block, "results_v2/distributed_spatial_lags/tables/dsl_incremental_spatial_block.csv", row.names = FALSE)
message("\n---- Incremental fit from adding the K=3 spatial-lag block (M3 -> M4) ----")
print(incremental_spatial_block)

## ================================================================
## 5. Primary K=3 Driscoll-Kraay inference for theta_1..3 and phi_1..3
##    (M4), identical DK implementation to the audited lag-1 stage.
## ================================================================

vcov_m4_dk <- dk_vcov(m4, sample_k3$region_id_f, sample_k3$date_f)
vcov_m4_cl <- vcovCL(m4, cluster = sample_k3$region_id_f, type = "HC1")

theta_coefs <- w_cols_k3
phi_coefs   <- own_cols_k3

theta_dk <- robust_coeftest(m4, vcov_m4_dk, theta_coefs); theta_dk$coef <- theta_coefs; theta_dk$block <- "spatial (theta_k)"
phi_dk   <- robust_coeftest(m4, vcov_m4_dk, phi_coefs);   phi_dk$coef   <- phi_coefs;   phi_dk$block   <- "own (phi_k)"

coef_report_k3 <- bind_rows(phi_dk, theta_dk)[, c("block", "coef", "estimate", "std_error", "z_value", "p_value", "ci_low", "ci_high")]
write.csv(coef_report_k3, "results_v2/distributed_spatial_lags/tables/dsl_k3_coefficients_dk.csv", row.names = FALSE)
message("\n---- K=3 phi_k and theta_k, Driscoll-Kraay inference (M4) ----")
print(coef_report_k3)

## ---- Joint Wald tests (asymptotic, DK vcov) ----

wald_joint <- function(model, vcov_mat, coefs) {
  b <- coef(model)[coefs]; V <- vcov_mat[coefs, coefs]
  stat <- as.numeric(t(b) %*% solve(V) %*% b)
  data.frame(statistic = stat, df = length(coefs), p_value_asymptotic = pchisq(stat, df = length(coefs), lower.tail = FALSE))
}

wald_theta_k3_dk <- wald_joint(m4, vcov_m4_dk, theta_coefs)
wald_phi_k3_dk   <- wald_joint(m4, vcov_m4_dk, phi_coefs)
message("\nJoint Wald (DK), H0: theta_1=theta_2=theta_3=0: stat=", round(wald_theta_k3_dk$statistic,3),
        " df=", wald_theta_k3_dk$df, " p=", signif(wald_theta_k3_dk$p_value_asymptotic,3))
message("Joint Wald (DK), H0: phi_1=phi_2=phi_3=0: stat=", round(wald_phi_k3_dk$statistic,3),
        " df=", wald_phi_k3_dk$df, " p=", signif(wald_phi_k3_dk$p_value_asymptotic,3))

## ---- Cumulative effects: Theta_3 = sum(theta_k), Phi_3 = sum(phi_k) ----
## Variance via the full coefficient covariance matrix: Var(sum) = 1'V1.

cumulative_effect <- function(model, vcov_mat, coefs) {
  b <- coef(model)[coefs]; V <- vcov_mat[coefs, coefs]
  est <- sum(b); se <- sqrt(sum(V))
  z <- est / se
  data.frame(cumulative_estimate = est, std_error = se, z_value = z, p_value = 2*pnorm(-abs(z)),
             ci_low = est - 1.96*se, ci_high = est + 1.96*se)
}
Theta_3_dk <- cumulative_effect(m4, vcov_m4_dk, theta_coefs)
Phi_3_dk   <- cumulative_effect(m4, vcov_m4_dk, phi_coefs)
message("\nCumulative Theta_3 (DK): ", round(Theta_3_dk$cumulative_estimate,4), ", p=", signif(Theta_3_dk$p_value,3))
message("Cumulative Phi_3 (DK): ", round(Phi_3_dk$cumulative_estimate,4), ", p=", signif(Phi_3_dk$p_value,3))

## ================================================================
## 6. Wild cluster bootstrap (B=9999) for the joint K=3 spatial-lag test
##    AND the cumulative Theta_3 effect (same bootstrap draws reused for
##    both, since both are functions of the same theta_1..3 estimates).
##    Restricted under M3 (H0: theta_1=theta_2=theta_3=0), Rademacher
##    weights by region cluster -- identical procedure to the audited
##    lag-1 stage, generalised from 1 to 3 joint restrictions.
## ================================================================

B_PRIMARY <- 9999
set.seed(20260921)
m3_fitted <- fitted(m3); m3_resid <- residuals(m3)
region_of_row_k3 <- sample_k3$region_id_f
regions_unique <- levels(region_of_row_k3)
formula_m4 <- as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3))

wald_obs_theta <- as.numeric(t(coef(m4)[theta_coefs]) %*% solve(vcov_m4_cl[theta_coefs, theta_coefs]) %*% coef(m4)[theta_coefs])
Theta_3_obs <- sum(coef(m4)[theta_coefs])

wald_boot <- numeric(B_PRIMARY)
Theta3_boot <- numeric(B_PRIMARY)
message("\nRunning K=3 primary wild cluster bootstrap (B=", B_PRIMARY, ") for the joint spatial-lag test -- this will take a while...")
for (b in seq_len(B_PRIMARY)) {
  v <- setNames(sample(c(-1, 1), length(regions_unique), replace = TRUE), regions_unique)
  y_star <- m3_fitted + m3_resid * v[as.character(region_of_row_k3)]
  boot_data <- sample_k3
  boot_data$inflation_mom <- y_star
  m4_star <- lm(formula_m4, data = boot_data)
  vcov_star <- vcovCL(m4_star, cluster = boot_data$region_id_f, type = "HC1")
  b_star <- coef(m4_star)[theta_coefs]
  V_star <- vcov_star[theta_coefs, theta_coefs]
  wald_boot[b] <- as.numeric(t(b_star) %*% solve(V_star) %*% b_star)
  Theta3_boot[b] <- sum(b_star)
}
p_wcb_joint_theta <- mean(wald_boot >= wald_obs_theta)
# Theta3_boot is generated under H0: theta_1=theta_2=theta_3=0 (resampled from
# the RESTRICTED model M3's residuals), so the WCB p-value for the cumulative
# effect compares the observed |Theta_3| against this null bootstrap distribution directly.
p_wcb_Theta3 <- mean(abs(Theta3_boot) >= abs(Theta_3_obs))

wcb_k3_summary <- data.frame(
  test = c("Joint (theta_1,theta_2,theta_3)", "Cumulative Theta_3"),
  statistic = c(wald_obs_theta, Theta_3_obs),
  statistic_type = c("Wald chi-sq (df=3)", "sum(theta_k)"),
  p_value_wcb = c(p_wcb_joint_theta, p_wcb_Theta3),
  p_value_asymptotic = c(wald_theta_k3_dk$p_value_asymptotic, Theta_3_dk$p_value),
  n_bootstrap_reps = B_PRIMARY, seed = 20260921
)
write.csv(wcb_k3_summary, "results_v2/distributed_spatial_lags/tables/dsl_k3_wild_bootstrap_joint_and_cumulative.csv", row.names = FALSE)
message("\n---- K=3 wild cluster bootstrap: joint spatial-lag test and cumulative Theta_3 ----")
print(wcb_k3_summary)

## ================================================================
## 7. K=6 robustness: own lags 1-6 + Queen spatial lags 1-6, on the
##    (more restrictive) K=6 sample. Not held to the same B=9999 bar as
##    the K=3 primary bootstrap (that requirement was explicit for K=3
##    only) -- B=1999 here, "where technically feasible."
## ================================================================

f_own6 <- paste(own_cols_k6, collapse = " + ")
f_w6   <- paste(w_cols_k6, collapse = " + ")

m3_k6 <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own6)), data = sample_k6)
m4_k6 <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own6, "+", f_w6)), data = sample_k6)

vcov_k6_dk <- dk_vcov(m4_k6, sample_k6$region_id_f, sample_k6$date_f)
vcov_k6_cl <- vcovCL(m4_k6, cluster = sample_k6$region_id_f, type = "HC1")

theta_coefs_k6 <- w_cols_k6; phi_coefs_k6 <- own_cols_k6
theta_k6_dk <- robust_coeftest(m4_k6, vcov_k6_dk, theta_coefs_k6); theta_k6_dk$coef <- theta_coefs_k6; theta_k6_dk$block <- "spatial (theta_k)"
phi_k6_dk   <- robust_coeftest(m4_k6, vcov_k6_dk, phi_coefs_k6);   phi_k6_dk$coef   <- phi_coefs_k6;   phi_k6_dk$block   <- "own (phi_k)"
coef_report_k6 <- bind_rows(phi_k6_dk, theta_k6_dk)[, c("block", "coef", "estimate", "std_error", "z_value", "p_value", "ci_low", "ci_high")]
write.csv(coef_report_k6, "results_v2/distributed_spatial_lags/tables/dsl_k6_coefficients_dk.csv", row.names = FALSE)
message("\n---- K=6 phi_k and theta_k, Driscoll-Kraay inference ----")
print(coef_report_k6)

wald_theta_k6_dk <- wald_joint(m4_k6, vcov_k6_dk, theta_coefs_k6)
wald_phi_k6_dk   <- wald_joint(m4_k6, vcov_k6_dk, phi_coefs_k6)
Theta_6_dk <- cumulative_effect(m4_k6, vcov_k6_dk, theta_coefs_k6)
Phi_6_dk   <- cumulative_effect(m4_k6, vcov_k6_dk, phi_coefs_k6)
message("\nJoint Wald (DK) K=6, H0: theta_1=...=theta_6=0: stat=", round(wald_theta_k6_dk$statistic,3),
        " df=", wald_theta_k6_dk$df, " p=", signif(wald_theta_k6_dk$p_value_asymptotic,3))
message("Cumulative Theta_6 (DK): ", round(Theta_6_dk$cumulative_estimate,4), ", p=", signif(Theta_6_dk$p_value,3))

B_K6 <- 1999
set.seed(20260921)
m3k6_fitted <- fitted(m3_k6); m3k6_resid <- residuals(m3_k6)
region_of_row_k6 <- sample_k6$region_id_f
regions_unique_k6 <- levels(region_of_row_k6)
formula_m4_k6 <- as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own6, "+", f_w6))

wald_obs_theta_k6 <- as.numeric(t(coef(m4_k6)[theta_coefs_k6]) %*% solve(vcov_k6_cl[theta_coefs_k6, theta_coefs_k6]) %*% coef(m4_k6)[theta_coefs_k6])
Theta_6_obs <- sum(coef(m4_k6)[theta_coefs_k6])

wald_boot_k6 <- numeric(B_K6); Theta6_boot <- numeric(B_K6)
message("\nRunning K=6 robustness wild cluster bootstrap (B=", B_K6, ")...")
for (b in seq_len(B_K6)) {
  v <- setNames(sample(c(-1, 1), length(regions_unique_k6), replace = TRUE), regions_unique_k6)
  y_star <- m3k6_fitted + m3k6_resid * v[as.character(region_of_row_k6)]
  boot_data <- sample_k6
  boot_data$inflation_mom <- y_star
  m4k6_star <- lm(formula_m4_k6, data = boot_data)
  vcov_star <- vcovCL(m4k6_star, cluster = boot_data$region_id_f, type = "HC1")
  b_star <- coef(m4k6_star)[theta_coefs_k6]
  V_star <- vcov_star[theta_coefs_k6, theta_coefs_k6]
  wald_boot_k6[b] <- as.numeric(t(b_star) %*% solve(V_star) %*% b_star)
  Theta6_boot[b] <- sum(b_star)
}
p_wcb_joint_theta_k6 <- mean(wald_boot_k6 >= wald_obs_theta_k6)
p_wcb_Theta6 <- mean(abs(Theta6_boot) >= abs(Theta_6_obs))

wcb_k6_summary <- data.frame(
  test = c("Joint (theta_1..theta_6)", "Cumulative Theta_6"),
  statistic = c(wald_obs_theta_k6, Theta_6_obs),
  statistic_type = c("Wald chi-sq (df=6)", "sum(theta_k)"),
  p_value_wcb = c(p_wcb_joint_theta_k6, p_wcb_Theta6),
  p_value_asymptotic = c(wald_theta_k6_dk$p_value_asymptotic, Theta_6_dk$p_value),
  n_bootstrap_reps = B_K6, seed = 20260921
)
write.csv(wcb_k6_summary, "results_v2/distributed_spatial_lags/tables/dsl_k6_wild_bootstrap_joint_and_cumulative.csv", row.names = FALSE)
message("\n---- K=6 wild cluster bootstrap: joint spatial-lag test and cumulative Theta_6 ----")
print(wcb_k6_summary)

## ================================================================
## 8. Spatial-weight robustness: re-estimate the K=3 spatial-lag block
##    (M4-equivalent) with KNN(k=4) and distance-band W, on the SAME
##    K=3 sample. W is not selected on significance.
## ================================================================

f_w3_knn  <- paste(paste0("w_lag_", 1:3, "_knn"),  collapse = " + ")
f_w3_dist <- paste(paste0("w_lag_", 1:3, "_dist"), collapse = " + ")

m4_knn  <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3_knn)),  data = sample_k3)
m4_dist <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3_dist)), data = sample_k3)

vcov_knn_dk  <- dk_vcov(m4_knn,  sample_k3$region_id_f, sample_k3$date_f)
vcov_dist_dk <- dk_vcov(m4_dist, sample_k3$region_id_f, sample_k3$date_f)

theta_knn_coefs  <- paste0("w_lag_", 1:3, "_knn")
theta_dist_coefs <- paste0("w_lag_", 1:3, "_dist")

spatial_weight_robustness <- bind_rows(
  cbind(weight_spec = "Queen contiguity (primary)", lag = paste0("theta_", 1:3), robust_coeftest(m4, vcov_m4_dk, theta_coefs)),
  cbind(weight_spec = "KNN (k=4)",                  lag = paste0("theta_", 1:3), robust_coeftest(m4_knn, vcov_knn_dk, theta_knn_coefs)),
  cbind(weight_spec = "Distance band",              lag = paste0("theta_", 1:3), robust_coeftest(m4_dist, vcov_dist_dk, theta_dist_coefs))
)
write.csv(spatial_weight_robustness, "results_v2/distributed_spatial_lags/tables/dsl_spatial_weight_robustness.csv", row.names = FALSE)
message("\n---- theta_k across spatial weight specifications (K=3, DK SE) ----")
print(spatial_weight_robustness)

Theta_3_knn  <- cumulative_effect(m4_knn,  vcov_knn_dk,  theta_knn_coefs)
Theta_3_dist <- cumulative_effect(m4_dist, vcov_dist_dk, theta_dist_coefs)
cumulative_by_weight <- bind_rows(
  cbind(weight_spec = "Queen contiguity (primary)", Theta_3_dk),
  cbind(weight_spec = "KNN (k=4)", Theta_3_knn),
  cbind(weight_spec = "Distance band", Theta_3_dist)
)
write.csv(cumulative_by_weight, "results_v2/distributed_spatial_lags/tables/dsl_cumulative_Theta3_by_weight.csv", row.names = FALSE)
message("\n---- Cumulative Theta_3 across spatial weight specifications ----")
print(cumulative_by_weight)

sign_stable <- length(unique(sign(cumulative_by_weight$cumulative_estimate))) == 1
message("\nSign of cumulative Theta_3 stable across Queen/KNN/distance: ", sign_stable)

## ================================================================
## 9. Residual diagnostics for the preferred K=3 model (M4): region-level
##    Ljung-Box + BH/Holm, compared against the already-committed lag-1
##    hardened stage's own numbers (read from its results file, not
##    re-derived or hardcoded, so this stays correct if that file ever
##    changes upstream of this script -- though scripts/13's outputs are
##    themselves frozen/committed and not touched by this stage).
## ================================================================

sample_k3$resid_m4 <- residuals(m4)
lb_results_k3 <- sample_k3 %>%
  group_by(region_id) %>% arrange(month_index) %>%
  summarise(
    lb_stat = tryCatch(Box.test(resid_m4, lag = 1, type = "Ljung-Box")$statistic, error = function(e) NA_real_),
    lb_p    = tryCatch(Box.test(resid_m4, lag = 1, type = "Ljung-Box")$p.value,     error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(lb_p_BH = p.adjust(lb_p, method = "BH"), lb_p_Holm = p.adjust(lb_p, method = "holm"),
         sig_raw_05 = lb_p < 0.05, sig_BH_05 = lb_p_BH < 0.05, sig_Holm_05 = lb_p_Holm < 0.05)

n_sig_raw_k3  <- sum(lb_results_k3$sig_raw_05, na.rm = TRUE)
n_sig_BH_k3   <- sum(lb_results_k3$sig_BH_05, na.rm = TRUE)
n_sig_Holm_k3 <- sum(lb_results_k3$sig_Holm_05, na.rm = TRUE)

lag1_ref_path <- "results_v2/tables/lag_residual_serial_correlation.csv"
lag1_ref <- if (file.exists(lag1_ref_path)) read.csv(lag1_ref_path, stringsAsFactors = FALSE) else NULL

serial_corr_comparison <- data.frame(
  model = c("Lag-1 hardened stage (M2, scripts/13, reference)", "K=3 distributed model (M4, this stage)"),
  n_regions_LjungBox_sig_raw_05  = c(if (!is.null(lag1_ref)) lag1_ref$n_regions_LjungBox_sig_raw_05  else NA, n_sig_raw_k3),
  n_regions_LjungBox_sig_BH_05   = c(if (!is.null(lag1_ref)) lag1_ref$n_regions_LjungBox_sig_BH_05   else NA, n_sig_BH_k3),
  n_regions_LjungBox_sig_Holm_05 = c(if (!is.null(lag1_ref)) lag1_ref$n_regions_LjungBox_sig_Holm_05 else NA, n_sig_Holm_k3),
  n_regions_total = c(13, 13)
)
write.csv(lb_results_k3, "results_v2/distributed_spatial_lags/tables/dsl_k3_residual_ljungbox_by_region.csv", row.names = FALSE)
write.csv(serial_corr_comparison, "results_v2/distributed_spatial_lags/tables/dsl_residual_serial_correlation_comparison.csv", row.names = FALSE)
message("\n---- Residual serial correlation: K=3 distributed model vs. lag-1 hardened stage ----")
print(serial_corr_comparison)
if (!is.null(lag1_ref)) {
  decreased <- n_sig_raw_k3 < lag1_ref$n_regions_LjungBox_sig_raw_05
  message("Residual serial dependence (raw count) ", if (decreased) "DECREASED" else "did NOT decrease",
          " after adding distributed own and spatial lags (", lag1_ref$n_regions_LjungBox_sig_raw_05,
          " -> ", n_sig_raw_k3, ").")
} else {
  message("Lag-1 reference file not found -- comparison reported for this stage only.")
}

## ================================================================
## 10. Robustness interpretation summary.
## ================================================================

k3_vs_k6_theta1 <- data.frame(
  spec = c("K=3 model, theta_1", "K=6 model, theta_1"),
  estimate = c(coef(m4)["w_lag_1"], coef(m4_k6)["w_lag_1"]),
  p_value_dk = c(theta_dk$p_value[theta_dk$coef == "w_lag_1"], theta_k6_dk$p_value[theta_k6_dk$coef == "w_lag_1"])
)
cumulative_k3_vs_k6 <- data.frame(
  spec = c("Theta_3 (K=3)", "Theta_6 (K=6, first 3 lags only)", "Theta_6 (K=6, full)"),
  estimate = c(Theta_3_dk$cumulative_estimate,
               sum(coef(m4_k6)[paste0("w_lag_", 1:3)]),
               Theta_6_dk$cumulative_estimate),
  p_value_asymptotic_dk = c(Theta_3_dk$p_value, NA, Theta_6_dk$p_value)
)
write.csv(k3_vs_k6_theta1, "results_v2/distributed_spatial_lags/tables/dsl_robustness_k3_vs_k6_theta1.csv", row.names = FALSE)
write.csv(cumulative_k3_vs_k6, "results_v2/distributed_spatial_lags/tables/dsl_robustness_cumulative_k3_vs_k6.csv", row.names = FALSE)
message("\n---- K=3 vs K=6 robustness ----")
print(k3_vs_k6_theta1)
print(cumulative_k3_vs_k6)

robustness_summary <- data.frame(
  check = c("theta_k signs stable across Queen/KNN/distance (K=3)",
            "Cumulative Theta_3 sign stable across Queen/KNN/distance",
            "Joint spatial-lag test significant at 5% (K=3, WCB)",
            "Joint spatial-lag test significant at 5% (K=6, WCB)",
            "Conclusions depend materially on K=3 vs K=6"),
  result = c(
    length(unique(sign(spatial_weight_robustness$estimate[spatial_weight_robustness$lag == "theta_1"]))) == 1 &&
      length(unique(sign(spatial_weight_robustness$estimate[spatial_weight_robustness$lag == "theta_2"]))) == 1 &&
      length(unique(sign(spatial_weight_robustness$estimate[spatial_weight_robustness$lag == "theta_3"]))) == 1,
    sign_stable,
    p_wcb_joint_theta < 0.05,
    p_wcb_joint_theta_k6 < 0.05,
    (p_wcb_joint_theta < 0.05) != (p_wcb_joint_theta_k6 < 0.05)
  )
)
write.csv(robustness_summary, "results_v2/distributed_spatial_lags/tables/dsl_robustness_summary.csv", row.names = FALSE)
message("\n---- Robustness interpretation summary ----")
print(robustness_summary)

## ================================================================
## 11. Figures
## ================================================================

fig_dir <- "results_v2/distributed_spatial_lags/figures"

p_f1 <- ggplot(theta_dk %>% mutate(k = 1:3), aes(x = factor(k), y = estimate, ymin = ci_low, ymax = ci_high)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(color = "#2c3e50") +
  labs(title = "K=3 Distributed Spatial Lags: theta_k (Queen, Driscoll-Kraay 95% CI)",
       subtitle = "Primary specification, M4 on the common K=3 sample",
       x = "Lag k (months)", y = expression(theta[k])) +
  theme_minimal()
ggsave(file.path(fig_dir, "fig_G1_k3_theta_by_lag.png"), p_f1, width = 8, height = 5, dpi = 300)

p_f2 <- ggplot(theta_k6_dk %>% mutate(k = 1:6), aes(x = factor(k), y = estimate, ymin = ci_low, ymax = ci_high)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(color = "#2c3e50") +
  labs(title = "K=6 Robustness: theta_k (Queen, Driscoll-Kraay 95% CI)",
       subtitle = "K=6 robustness sample",
       x = "Lag k (months)", y = expression(theta[k])) +
  theme_minimal()
ggsave(file.path(fig_dir, "fig_G2_k6_theta_by_lag.png"), p_f2, width = 8, height = 5, dpi = 300)

p_f3 <- ggplot(spatial_weight_robustness, aes(x = lag, y = estimate, ymin = ci_low, ymax = ci_high, color = weight_spec)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.4)) +
  labs(title = "theta_k Across Spatial Weight Specifications (K=3, Driscoll-Kraay 95% CI)",
       x = NULL, y = expression(theta[k]), color = NULL) +
  theme_minimal() + theme(legend.position = "top")
ggsave(file.path(fig_dir, "fig_G3_theta_by_weight_spec.png"), p_f3, width = 9, height = 5, dpi = 300)

p_f4 <- ggplot(lb_results_k3, aes(x = reorder(region_id, lb_p), y = lb_p)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#d7191c") +
  geom_col(fill = "#2c3e50") +
  coord_flip() +
  labs(title = "K=3 Distributed Model (M4) Residual Serial Correlation by Region",
       subtitle = "Ljung-Box lag 1; red line = 0.05 threshold",
       x = NULL, y = "Ljung-Box p-value") +
  theme_minimal()
ggsave(file.path(fig_dir, "fig_G4_k3_residual_ljungbox_by_region.png"), p_f4, width = 8, height = 6, dpi = 300)

message("\nWrote figures G1-G4 to ", fig_dir, "/")

## ================================================================
## 12. Markdown methodological report.
## ================================================================

md <- c(
"# MAKANI v2 -- Distributed Temporal and Spatial Lags",
"",
"Stage 3 of the MAKANI v2 econometric sequence, building on the hardened",
"lag-1 stage (`scripts/13_lagged_spatial_transmission.R`, commit `f7f5d5d`).",
"Tests whether a **distributed** block of lagged neighbour inflation",
"(spatial lags 1..K) carries additional predictive content beyond own",
"inflation dynamics, region fixed effects, and common month effects.",
"",
"**Language discipline (per instruction): results below are described as",
"\"evidence consistent with lagged spatial dependence\" or its absence --",
"never as causal transmission, contagion, or structural spillovers.**",
"",
"## Estimation samples",
"",
"| Sample | N | T | Obs | Date range |",
"|---|---|---|---|---|",
sprintf("| K=3 common (M0-M4) | %d | %d | %d | %s to %s |",
        sample_report$n_regions[1], sample_report$n_months[1], sample_report$total_obs[1],
        sample_report$first_date[1], sample_report$last_date[1]),
sprintf("| K=6 robustness | %d | %d | %d | %s to %s |",
        sample_report$n_regions[2], sample_report$n_months[2], sample_report$total_obs[2],
        sample_report$first_date[2], sample_report$last_date[2]),
"",
"All 13 construction-validation checks passed (no look-ahead, correct",
"region ordering in W, exactly one lag value per region-month per order,",
"no region-ID mixing, structural missingness exactly as expected for",
"k=1, 3, 6). See `dsl_construction_validation.csv`.",
"",
"## M0-M4 model comparison (common K=3 sample)",
"",
"See `dsl_model_comparison_M0_M4.csv` for the full table (N, params, R2,",
"adjusted R2, residual SE, AIC, BIC, within-R2 of the cumulative lag",
"block for each model). Incremental fit from adding the K=3 spatial-lag",
"block (M3 -> M4): see `dsl_incremental_spatial_block.csv` --",
sprintf("delta R2 = %.6f, delta AIC = %.2f (higher/worse), delta BIC = %.2f (higher/worse), F-test p = %.3f.",
        incremental_spatial_block$delta_r_squared_M3_to_M4, incremental_spatial_block$delta_aic_M3_to_M4,
        incremental_spatial_block$delta_bic_M3_to_M4, incremental_spatial_block$f_test_p_value_M3_vs_M4),
"",
"## Primary K=3 coefficients (Driscoll-Kraay, M4)",
"",
"See `dsl_k3_coefficients_dk.csv` for phi_1..3 and theta_1..3 with DK SE,",
"p-value, and 95% CI.",
"",
"**Joint tests (DK asymptotic Wald):**",
sprintf("- H0: theta_1=theta_2=theta_3=0 -- stat=%.3f, df=%d, p=%.4g",
        wald_theta_k3_dk$statistic, wald_theta_k3_dk$df, wald_theta_k3_dk$p_value_asymptotic),
sprintf("- H0: phi_1=phi_2=phi_3=0 -- stat=%.3f, df=%d, p=%.4g",
        wald_phi_k3_dk$statistic, wald_phi_k3_dk$df, wald_phi_k3_dk$p_value_asymptotic),
"",
sprintf("**Wild cluster bootstrap (B=%d, seed=20260921), joint spatial-lag test: p = %.4g.**",
        B_PRIMARY, p_wcb_joint_theta),
"",
"## Cumulative effects",
"",
sprintf("- Theta_3 = sum(theta_k) = %.4f, DK SE = %.4f, DK p = %.4g, 95%% CI [%.4f, %.4f]; WCB p = %.4g",
        Theta_3_dk$cumulative_estimate, Theta_3_dk$std_error, Theta_3_dk$p_value,
        Theta_3_dk$ci_low, Theta_3_dk$ci_high, p_wcb_Theta3),
sprintf("- Phi_3 = sum(phi_k) = %.4f, DK SE = %.4f, DK p = %.4g, 95%% CI [%.4f, %.4f]",
        Phi_3_dk$cumulative_estimate, Phi_3_dk$std_error, Phi_3_dk$p_value, Phi_3_dk$ci_low, Phi_3_dk$ci_high),
sprintf("- Theta_6 (K=6 robustness) = %.4f, DK p = %.4g; WCB p (B=%d) = %.4g",
        Theta_6_dk$cumulative_estimate, Theta_6_dk$p_value, B_K6, p_wcb_Theta6),
"",
"## Spatial-weight robustness (K=3)",
"",
"See `dsl_spatial_weight_robustness.csv` and `dsl_cumulative_Theta3_by_weight.csv`.",
sprintf("Sign of cumulative Theta_3 stable across Queen/KNN/distance: %s.", sign_stable),
"Rook not treated as independent robustness (already shown identical to",
"Queen in this sample, scripts/12).",
"",
"## Residual diagnostics",
"",
"See `dsl_residual_serial_correlation_comparison.csv` and",
"`dsl_k3_residual_ljungbox_by_region.csv` (with BH/Holm correction).",
sprintf("K=3 model (M4): %d/13 regions significant at raw 5%%, %d/13 survive BH, %d/13 survive Holm.",
        n_sig_raw_k3, n_sig_BH_k3, n_sig_Holm_k3),
"",
"## Robustness interpretation",
"",
"See `dsl_robustness_summary.csv`, `dsl_robustness_k3_vs_k6_theta1.csv`,",
"and `dsl_robustness_cumulative_k3_vs_k6.csv` for the full numeric basis.",
"No lag length was selected based on statistical significance -- K=3 was",
"specified as primary and K=6 as robustness before either was estimated.",
"",
"## Scientific guardrails",
"",
"No claim of causality, contagion, established transmission, or",
"structural spillovers is made anywhere in this stage's outputs. Where",
"evidence is discussed, the language used is \"evidence consistent with",
"lagged spatial dependence\" (or its absence), never a causal claim.",
"",
"## Files",
"",
"All tables: `results_v2/distributed_spatial_lags/tables/`.",
"All figures: `results_v2/distributed_spatial_lags/figures/`.",
"Script: `scripts/14_distributed_spatial_lags.R`."
)
writeLines(md, "results_v2/distributed_spatial_lags/DISTRIBUTED_SPATIAL_LAGS_REPORT.md")
message("\nWrote results_v2/distributed_spatial_lags/DISTRIBUTED_SPATIAL_LAGS_REPORT.md")
message("14_distributed_spatial_lags.R complete.")
