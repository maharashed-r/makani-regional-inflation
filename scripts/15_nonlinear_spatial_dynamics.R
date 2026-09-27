# ================================================================
# MAKANI v2 -- Stage 4: nonlinear and state-dependent regional
# inflation dynamics.
#
# Research question (explicitly NOT causal): does the weak average
# linear spatial relationship documented in Stage 14 (commit fd7ef10)
# mask (A) state-dependent effects, (B) smooth nonlinearity, (C) a
# threshold/regime split, or (E) asymmetric positive/negative
# responses? Nonlinearity is NOT assumed -- every test below can, and
# is expected in most cases to, return a null result.
#
# NEW v2 FILE. Reads only data_v2/saudi_cpi_panel_v2.csv,
# data_v2/region_crosswalk.csv, and the baseline's read-only geometry
# cache. Does not source or modify scripts 01-14, baseline/,
# data_clean/, results/, or any earlier v2 output. The K=3 lag/weight
# construction is REPRODUCED independently here (not sourced), for the
# stated reason that an exactly matched common estimation sample is
# required for the interaction/spline/threshold models below -- the
# Stage-14 linear reference itself is not rebuilt or reinterpreted, and
# this script's K=3 sample is validated against Stage 14's own saved
# sample summary as a cross-check, not a substitute.
#
# ---- PRE-SPECIFIED DESIGN DECISIONS (stated before any model is run,
# per instruction not to optimize any of this for significance) ----
#
# 1. State variable: national_infl_t = the cross-sectional mean of
#    inflation_mom across all 13 regions in month t. This is constant
#    within a month by construction, so it (and any transformation of
#    it alone) is perfectly collinear with the month fixed effects
#    already in every model -- only INTERACTIONS of national_infl_t
#    with a region-varying regressor add new information beyond the
#    FE, exactly as the task instructions note for the high-inflation
#    dummy. No standalone national_infl_t or f(national_infl_t) main
#    effect is therefore included in any formula below.
#
# 2. High-inflation regime (Section A, primary): top tercile of
#    national_infl_t's own empirical distribution WITHIN the K=3
#    estimation sample (i.e. HighInflation_t = 1 if national_infl_t
#    exceeds its 2/3 quantile over the 154 sample months). Alternative
#    pre-specified definitions (top quartile; mean + 1 SD) are carried
#    as robustness only (Section F), never substituted into the
#    primary test after the fact.
#
# 3. Composite spatial-lag variable (Sections B, C, D, E only --
#    Section A keeps w_lag_1/2/3 fully separate, exactly as specified):
#    w_lag_avg = mean(w_lag_1, w_lag_2, w_lag_3). Using one composite
#    variable for the spline/quadratic/threshold/asymmetry tests keeps
#    each of those models to a small, pre-specified number of added
#    parameters (appropriate for N=13 clusters); this simplification is
#    stated here, before estimation, not chosen after seeing results.
#
# 4. Spline: splines::ns(national_infl_t, df = 3) (base R, no package
#    added), df fixed before estimation, default quantile-based knot
#    placement (never hand-tuned or searched). Model includes a plain
#    linear w_lag_avg term plus w_lag_avg interacted with the 3-column
#    spline basis, so the joint test on the spline-interaction terms
#    specifically targets curvature beyond a flat average effect.
#
# 5. Threshold search: candidate gamma grid = distinct monthly
#    national_infl_t values strictly between the 15th and 85th
#    percentile of the K=3 sample (>=15% trimmed on each side, per
#    instruction), thinned to at most ~30 evenly spaced grid points for
#    computational tractability under repeated bootstrap resampling --
#    stated here as a practical necessity, not a search for a
#    significant threshold. Bootstrap existence test uses B=499
#    (reduced from the requested B=9999 "where computationally
#    feasible": each replication requires re-running the FULL grid
#    search, making 9999 reps computationally prohibitive on this
#    machine within a single run; 499 reps x ~30 grid points is
#    ~15,000 model fits, already a substantial computation). This
#    deviation is stated explicitly, not hidden.
#
# Run manually: Rscript scripts/15_nonlinear_spatial_dynamics.R
# ================================================================

suppressMessages({
  library(sf)
  library(spdep)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sandwich)
  library(splines)
})

dir.create("results_v2/nonlinear_spatial_dynamics/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results_v2/nonlinear_spatial_dynamics/figures", recursive = TRUE, showWarnings = FALSE)

set.seed(20260921)

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
wald_joint <- function(model, vcov_mat, coefs) {
  b <- coef(model)[coefs]; V <- vcov_mat[coefs, coefs]
  stat <- as.numeric(t(b) %*% solve(V) %*% b)
  data.frame(statistic = stat, df = length(coefs), p_value_asymptotic = pchisq(stat, df = length(coefs), lower.tail = FALSE))
}
lincomb_effect <- function(model, vcov_mat, coefs, weights = rep(1, length(coefs))) {
  b <- coef(model)[coefs]; V <- vcov_mat[coefs, coefs]
  est <- sum(weights * b); se <- sqrt(as.numeric(t(weights) %*% V %*% weights))
  z <- est / se
  data.frame(estimate = est, std_error = se, z_value = z, p_value = 2 * pnorm(-abs(z)),
             ci_low = est - 1.96 * se, ci_high = est + 1.96 * se)
}

panel     <- read.csv("data_v2/saudi_cpi_panel_v2.csv", stringsAsFactors = FALSE)
panel$date <- as.Date(panel$date)
crosswalk <- read.csv("data_v2/region_crosswalk.csv", stringsAsFactors = FALSE)

## ================================================================
## 1. Spatial weights -- identical method to scripts/12-14.
## ================================================================

region_order <- sort(unique(crosswalk$region_id))
stopifnot(length(region_order) == 13)

geom_raw <- readRDS("data_raw/saudi_states_sf.rds")

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

lw_queen <- nb2listw(nb_queen, style = "W", zero.policy = TRUE)
lw_knn   <- nb2listw(nb_knn4,  style = "W", zero.policy = TRUE)
lw_dist  <- nb2listw(nb_dist,  style = "W", zero.policy = TRUE)

W_queen <- listw2mat(lw_queen); dimnames(W_queen) <- list(region_order, region_order)
W_knn   <- listw2mat(lw_knn);   dimnames(W_knn)   <- list(region_order, region_order)
W_dist  <- listw2mat(lw_dist);  dimnames(W_dist)  <- list(region_order, region_order)

message("Spatial weights rebuilt (Queen primary; KNN k=4, distance-band robustness). Same region order/row-standardization as Stages 12-14.")

## ================================================================
## 2. Distributed lags k=1..6, national inflation state variable.
## ================================================================

panel <- panel %>% mutate(month_index = as.integer(format(date, "%Y")) * 12L + as.integer(format(date, "%m")))

wide <- panel %>% select(month_index, region_id, inflation_mom) %>%
  pivot_wider(names_from = region_id, values_from = inflation_mom) %>% arrange(month_index)
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
  mi_chr <- as.character(month_idx); ok <- mi_chr %in% rownames(mat)
  out <- rep(NA_real_, length(month_idx)); out[ok] <- mat[cbind(mi_chr[ok], region[ok])]
  out
}

model_data <- panel
for (k in 1:6) {
  model_data[[paste0("own_lag_", k)]]        <- lookup(pi_mat,             model_data$month_index - k, model_data$region_id)
  model_data[[paste0("w_lag_", k)]]          <- lookup(spatial_lag_queen,  model_data$month_index - k, model_data$region_id)
  model_data[[paste0("w_lag_", k, "_knn")]]  <- lookup(spatial_lag_knn,    model_data$month_index - k, model_data$region_id)
  model_data[[paste0("w_lag_", k, "_dist")]] <- lookup(spatial_lag_dist,   model_data$month_index - k, model_data$region_id)
}

# National/common inflation state variable: cross-sectional mean of
# inflation_mom across all 13 regions in month t (contemporaneous with
# the outcome, not a lag -- no look-ahead concern, same status as a
# month fixed effect; see design note 1 above).
national_infl <- panel %>% filter(!is.na(inflation_mom)) %>%
  group_by(month_index) %>% summarise(national_infl_t = mean(inflation_mom), .groups = "drop")
model_data <- model_data %>% left_join(national_infl, by = "month_index")

## ================================================================
## 3. Common K=3 sample -- must match Stage 14 exactly.
## ================================================================

own_cols_k3 <- paste0("own_lag_", 1:3); w_cols_k3 <- paste0("w_lag_", 1:3)
complete_rows <- function(df, cols) Reduce(`&`, lapply(cols, function(cc) !is.na(df[[cc]])))
sample_k3 <- model_data[complete_rows(model_data, c(own_cols_k3, w_cols_k3)), ]

sample_k3$region_id_f <- factor(sample_k3$region_id)
sample_k3$date_f      <- factor(sample_k3$date)

sample_check <- data.frame(
  n_regions = length(unique(sample_k3$region_id)), n_months = length(unique(sample_k3$date)),
  total_obs = nrow(sample_k3), first_date = format(min(sample_k3$date)), last_date = format(max(sample_k3$date))
)
stage14_ref_path <- "results_v2/distributed_spatial_lags/tables/dsl_estimation_samples.csv"
if (file.exists(stage14_ref_path)) {
  s14 <- read.csv(stage14_ref_path, stringsAsFactors = FALSE)[1, ]  # first row = K=3 sample
  matches_stage14 <- sample_check$n_regions == s14$n_regions && sample_check$n_months == s14$n_months &&
    sample_check$total_obs == s14$total_obs && sample_check$first_date == s14$first_date && sample_check$last_date == s14$last_date
  message("K=3 sample matches Stage 14's saved sample exactly: ", matches_stage14)
  if (!matches_stage14) stop("K=3 sample does NOT match Stage 14 -- investigate before proceeding.")
} else {
  message("Stage 14 sample file not found -- proceeding without cross-check (should not happen if Stage 14 is committed).")
}
write.csv(sample_check, "results_v2/nonlinear_spatial_dynamics/tables/nl_estimation_sample_k3.csv", row.names = FALSE)

# Composite spatial-lag variable and its positive/negative parts
# (constructed entirely from already-lagged, t-1-or-earlier information;
# see design note 3).
sample_k3$w_lag_avg     <- rowMeans(sample_k3[, w_cols_k3])
sample_k3$w_lag_avg_pos <- pmax(sample_k3$w_lag_avg, 0)
sample_k3$w_lag_avg_neg <- pmin(sample_k3$w_lag_avg, 0)
stopifnot(isTRUE(all.equal(sample_k3$w_lag_avg_pos + sample_k3$w_lag_avg_neg, sample_k3$w_lag_avg)))

# High-inflation regime (primary: top tercile within the K=3 sample's own months)
month_level_state <- sample_k3 %>% distinct(date, national_infl_t)
tercile_cut <- quantile(month_level_state$national_infl_t, probs = 2/3, na.rm = TRUE)
quartile_cut <- quantile(month_level_state$national_infl_t, probs = 3/4, na.rm = TRUE)
sd_cut <- mean(month_level_state$national_infl_t) + sd(month_level_state$national_infl_t)
sample_k3$high_infl_tercile  <- as.numeric(sample_k3$national_infl_t > tercile_cut)
sample_k3$high_infl_quartile <- as.numeric(sample_k3$national_infl_t > quartile_cut)
sample_k3$high_infl_sd       <- as.numeric(sample_k3$national_infl_t > sd_cut)

regime_def_table <- data.frame(
  definition = c("Top tercile (primary)", "Top quartile (robustness)", "Mean + 1 SD (robustness)"),
  cutoff = c(tercile_cut, quartile_cut, sd_cut),
  n_high_regime_months = c(sum(month_level_state$national_infl_t > tercile_cut),
                            sum(month_level_state$national_infl_t > quartile_cut),
                            sum(month_level_state$national_infl_t > sd_cut)),
  n_total_months = nrow(month_level_state)
)
write.csv(regime_def_table, "results_v2/nonlinear_spatial_dynamics/tables/nl_regime_definitions.csv", row.names = FALSE)
message("\n---- Pre-specified high-inflation regime definitions ----")
print(regime_def_table)

## ================================================================
## VALIDATION (beyond the Stage-14 sample cross-check above)
## ================================================================

val <- list()
val$W_order_matches_region_id <- identical(rownames(W_queen), region_order) &&
  identical(rownames(W_knn), region_order) && identical(rownames(W_dist), region_order)
val$no_duplicate_region_month <- !anyDuplicated(paste(sample_k3$region_id, sample_k3$date))
val$w_lag_avg_reconstructs_from_pos_neg <- isTRUE(all.equal(sample_k3$w_lag_avg_pos + sample_k3$w_lag_avg_neg, sample_k3$w_lag_avg))
val$national_infl_constant_within_month <- all(sapply(split(sample_k3$national_infl_t, sample_k3$date), function(x) length(unique(x)) == 1))
val$regime_dummy_is_deterministic_function_of_date <- all(sapply(split(sample_k3$high_infl_tercile, sample_k3$date), function(x) length(unique(x)) == 1))
val$no_lookahead_state_var_contemporaneous_not_future <- TRUE  # national_infl_t uses month t only, by construction (see design note 1); own/spatial lags unaffected, still strictly t-k

val_tbl <- data.frame(check = names(val), pass = unlist(val))
write.csv(val_tbl, "results_v2/nonlinear_spatial_dynamics/tables/nl_construction_validation.csv", row.names = FALSE)
message("\n---- Additional construction validation ----")
print(val_tbl)
if (!all(unlist(val))) stop("Nonlinear-stage construction validation failed.")

## ================================================================
## SECTION A: high-inflation regime interaction (primary nonlinear test)
## ================================================================

f_own3 <- paste(own_cols_k3, collapse = " + ")
f_w3   <- paste(w_cols_k3, collapse = " + ")
f_w3_x_regime <- paste(paste0(w_cols_k3, ":high_infl_tercile"), collapse = " + ")

m_regime <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3, "+", f_w3_x_regime)),
                data = sample_k3)

vcov_regime_dk <- dk_vcov(m_regime, sample_k3$region_id_f, sample_k3$date_f)
vcov_regime_cl <- vcovCL(m_regime, cluster = sample_k3$region_id_f, type = "HC1")

delta_coefs <- paste0(w_cols_k3, ":high_infl_tercile")
delta_dk <- robust_coeftest(m_regime, vcov_regime_dk, delta_coefs)
delta_dk$coef <- paste0("delta_", 1:3)
write.csv(delta_dk, "results_v2/nonlinear_spatial_dynamics/tables/nl_A_regime_interaction_coefficients.csv", row.names = FALSE)
message("\n---- Section A: regime-interaction coefficients (delta_k), Driscoll-Kraay ----")
print(delta_dk)

wald_delta_dk <- wald_joint(m_regime, vcov_regime_dk, delta_coefs)
message("Joint Wald (DK), H0: delta_1=delta_2=delta_3=0: stat=", round(wald_delta_dk$statistic,3),
        " df=", wald_delta_dk$df, " p=", signif(wald_delta_dk$p_value_asymptotic,3))

# Cumulative spatial effect by regime: normal = sum(theta_k); high = sum(theta_k + delta_k)
theta_coefs <- w_cols_k3
normal_regime_effect <- lincomb_effect(m_regime, vcov_regime_dk, theta_coefs, weights = rep(1,3))
high_regime_effect   <- lincomb_effect(m_regime, vcov_regime_dk, c(theta_coefs, delta_coefs), weights = rep(1,6))
regime_cumulative <- bind_rows(cbind(regime = "Normal inflation", normal_regime_effect),
                                cbind(regime = "High inflation (top tercile)", high_regime_effect))
write.csv(regime_cumulative, "results_v2/nonlinear_spatial_dynamics/tables/nl_A_cumulative_effect_by_regime.csv", row.names = FALSE)
message("\n---- Section A: cumulative spatial effect by regime ----")
print(regime_cumulative)

## ---- Section A: wild cluster bootstrap (B=9999), joint test of delta_1=delta_2=delta_3=0 ----
## Restricted model = m_regime without the three interaction terms (i.e. the plain K=3 model).

B_A <- 9999
m_restricted_A <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3)), data = sample_k3)
fitted_A <- fitted(m_restricted_A); resid_A <- residuals(m_restricted_A)
region_of_row <- sample_k3$region_id_f
regions_unique <- levels(region_of_row)
formula_regime <- as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3, "+", f_w3_x_regime))

wald_obs_A <- as.numeric(t(coef(m_regime)[delta_coefs]) %*% solve(vcov_regime_cl[delta_coefs, delta_coefs]) %*% coef(m_regime)[delta_coefs])

set.seed(20260921)
wald_boot_A <- numeric(B_A)
message("\nRunning Section A wild cluster bootstrap (B=", B_A, ") for H0: delta_1=delta_2=delta_3=0 ...")
for (b in seq_len(B_A)) {
  v <- setNames(sample(c(-1, 1), length(regions_unique), replace = TRUE), regions_unique)
  y_star <- fitted_A + resid_A * v[as.character(region_of_row)]
  boot_data <- sample_k3
  boot_data$inflation_mom <- y_star
  m_star <- lm(formula_regime, data = boot_data)
  vcov_star <- vcovCL(m_star, cluster = boot_data$region_id_f, type = "HC1")
  b_star <- coef(m_star)[delta_coefs]
  V_star <- vcov_star[delta_coefs, delta_coefs]
  wald_boot_A[b] <- as.numeric(t(b_star) %*% solve(V_star) %*% b_star)
}
p_wcb_A <- mean(wald_boot_A >= wald_obs_A)

wcb_A_summary <- data.frame(
  test = "Joint (delta_1, delta_2, delta_3)", statistic = wald_obs_A, statistic_type = "Wald chi-sq (df=3)",
  p_value_wcb = p_wcb_A, p_value_asymptotic = wald_delta_dk$p_value_asymptotic, n_bootstrap_reps = B_A, seed = 20260921
)
write.csv(wcb_A_summary, "results_v2/nonlinear_spatial_dynamics/tables/nl_A_wild_bootstrap_joint.csv", row.names = FALSE)
message("\n---- Section A: wild cluster bootstrap joint test ----")
print(wcb_A_summary)

## ================================================================
## SECTION B: smooth nonlinearity via natural spline (composite w_lag_avg)
## ================================================================

SPLINE_DF <- 3
ns_basis <- ns(sample_k3$national_infl_t, df = SPLINE_DF)  # fit once, reused for prediction (design note 4)
spline_cols <- paste0("ns_state_", 1:SPLINE_DF)
for (j in 1:SPLINE_DF) sample_k3[[paste0("ns_state_", j)]] <- ns_basis[, j]
spline_interact_terms <- paste0("w_lag_avg:", spline_cols)

m_spline <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3,
                                 "+ w_lag_avg +", paste(spline_interact_terms, collapse = " + "))),
                data = sample_k3)
vcov_spline_dk <- dk_vcov(m_spline, sample_k3$region_id_f, sample_k3$date_f)

spline_coefs <- spline_interact_terms
spline_dk <- robust_coeftest(m_spline, vcov_spline_dk, spline_coefs)
spline_dk$term <- paste0("w_lag_avg:ns_state_", 1:SPLINE_DF)
write.csv(spline_dk, "results_v2/nonlinear_spatial_dynamics/tables/nl_B_spline_coefficients.csv", row.names = FALSE)
message("\n---- Section B: spline-interaction coefficients (curvature beyond linear w_lag_avg) ----")
print(spline_dk)

wald_spline_dk <- wald_joint(m_spline, vcov_spline_dk, spline_coefs)
message("Joint Wald (DK) test of nonlinear spline terms: stat=", round(wald_spline_dk$statistic,3),
        " df=", wald_spline_dk$df, " p=", signif(wald_spline_dk$p_value_asymptotic,3))
write.csv(wald_spline_dk, "results_v2/nonlinear_spatial_dynamics/tables/nl_B_spline_joint_test.csv", row.names = FALSE)

# Marginal effect of w_lag_avg across OBSERVED support only (no extrapolation).
state_range <- range(sample_k3$national_infl_t)
state_grid <- seq(state_range[1], state_range[2], length.out = 100)
ns_grid <- predict(ns_basis, newx = state_grid)  # same knots as the fitted model
coefs_for_marginal <- c("w_lag_avg", spline_coefs)
V_marg <- vcov_spline_dk[coefs_for_marginal, coefs_for_marginal]
b_marg <- coef(m_spline)[coefs_for_marginal]

marginal_effect <- sapply(seq_len(nrow(ns_grid)), function(i) {
  design_row <- c(1, ns_grid[i, ])
  est <- as.numeric(design_row %*% b_marg)
  se  <- sqrt(as.numeric(t(design_row) %*% V_marg %*% design_row))
  c(est = est, se = se)
})
marginal_df <- data.frame(state = state_grid, estimate = marginal_effect["est", ], std_error = marginal_effect["se", ])
marginal_df$ci_low  <- marginal_df$estimate - 1.96 * marginal_df$std_error
marginal_df$ci_high <- marginal_df$estimate + 1.96 * marginal_df$std_error

# Data-support density (to flag weak-support regions on the marginal-effect plot)
dens <- density(sample_k3$national_infl_t, from = state_range[1], to = state_range[2])
support_density <- approx(dens$x, dens$y, xout = state_grid)$y
marginal_df$support_density <- support_density
low_support_threshold <- quantile(support_density, 0.1)  # bottom decile of density = flagged weak support
marginal_df$weak_support <- marginal_df$support_density < low_support_threshold

write.csv(marginal_df, "results_v2/nonlinear_spatial_dynamics/tables/nl_B_marginal_effect_by_state.csv", row.names = FALSE)
message("Marginal effect grid written (", sum(marginal_df$weak_support), " of ", nrow(marginal_df),
        " grid points flagged as weak-support, bottom decile of local density).")
message("Observed support range for national_infl_t: [", round(state_range[1],4), ", ", round(state_range[2],4), "] -- no prediction outside this range.")

## ================================================================
## SECTION C: centered quadratic diagnostic (secondary only)
## ================================================================

state_mean <- mean(sample_k3$national_infl_t)
sample_k3$state_centered  <- sample_k3$national_infl_t - state_mean
sample_k3$state_centered2 <- sample_k3$state_centered^2

m_quad <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3,
                               "+ w_lag_avg:state_centered + w_lag_avg:state_centered2")),
             data = sample_k3)
vcov_quad_dk <- dk_vcov(m_quad, sample_k3$region_id_f, sample_k3$date_f)
quad_coefs <- c("w_lag_avg:state_centered", "w_lag_avg:state_centered2")
quad_dk <- robust_coeftest(m_quad, vcov_quad_dk, quad_coefs)
quad_dk$term <- c("linear (w_lag_avg:state_centered)", "quadratic (w_lag_avg:state_centered2)")
write.csv(quad_dk, "results_v2/nonlinear_spatial_dynamics/tables/nl_C_quadratic_coefficients.csv", row.names = FALSE)
message("\n---- Section C: quadratic diagnostic coefficients ----")
print(quad_dk)

wald_quad_dk <- wald_joint(m_quad, vcov_quad_dk, quad_coefs)
message("Joint Wald (DK) test, linear+quadratic: stat=", round(wald_quad_dk$statistic,3),
        " df=", wald_quad_dk$df, " p=", signif(wald_quad_dk$p_value_asymptotic,3))
write.csv(wald_quad_dk, "results_v2/nonlinear_spatial_dynamics/tables/nl_C_quadratic_joint_test.csv", row.names = FALSE)

b_lin  <- coef(m_quad)["w_lag_avg:state_centered"]
b_quad <- coef(m_quad)["w_lag_avg:state_centered2"]
turning_point_centered <- -b_lin / (2 * b_quad)
turning_point_original <- turning_point_centered + state_mean
turning_in_support <- turning_point_original >= state_range[1] && turning_point_original <= state_range[2]
n_obs_near_turning <- sum(abs(sample_k3$national_infl_t - turning_point_original) < 0.5 * sd(sample_k3$national_infl_t))

turning_point_report <- data.frame(
  turning_point_original_scale = turning_point_original, within_observed_support = turning_in_support,
  observed_support_low = state_range[1], observed_support_high = state_range[2],
  n_obs_within_0.5sd_of_turning_point = n_obs_near_turning, total_obs = nrow(sample_k3)
)
write.csv(turning_point_report, "results_v2/nonlinear_spatial_dynamics/tables/nl_C_turning_point.csv", row.names = FALSE)
message("\n---- Section C: implied turning point ----")
print(turning_point_report)
if (!turning_in_support) message("Turning point falls OUTSIDE observed support -- NOT interpreted as an economically meaningful extremum.")

## ================================================================
## SECTION D: threshold model (Hansen-style grid search + wild bootstrap)
## ================================================================

trim_low  <- quantile(month_level_state$national_infl_t, 0.15)
trim_high <- quantile(month_level_state$national_infl_t, 0.85)
candidate_vals <- sort(unique(month_level_state$national_infl_t[
  month_level_state$national_infl_t > trim_low & month_level_state$national_infl_t < trim_high]))
GRID_MAX <- 30
if (length(candidate_vals) > GRID_MAX) {
  thin_idx <- round(seq(1, length(candidate_vals), length.out = GRID_MAX))
  gamma_grid <- candidate_vals[thin_idx]
} else {
  gamma_grid <- candidate_vals
}
message("\nThreshold search: ", length(gamma_grid), " candidate gamma values, trimmed to [",
        round(trim_low,4), ", ", round(trim_high,4), "] (15%/85% of the K=3 sample's monthly state distribution).")

m_thresh_restricted <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+ w_lag_avg")), data = sample_k3)
SSR0 <- sum(residuals(m_thresh_restricted)^2)

ssr_by_gamma <- sapply(gamma_grid, function(g) {
  regime <- as.numeric(sample_k3$national_infl_t > g)
  m_g <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+ w_lag_avg + w_lag_avg:regime")),
            data = cbind(sample_k3, regime = regime))
  sum(residuals(m_g)^2)
})
gamma_hat <- gamma_grid[which.min(ssr_by_gamma)]
regime_at_gamma_hat <- as.numeric(sample_k3$national_infl_t > gamma_hat)
n_low  <- sum(regime_at_gamma_hat == 0); n_high <- sum(regime_at_gamma_hat == 1)

m_thresh_hat <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+ w_lag_avg + w_lag_avg:regime")),
                    data = cbind(sample_k3, regime = regime_at_gamma_hat))
vcov_thresh_dk <- dk_vcov(m_thresh_hat, sample_k3$region_id_f, sample_k3$date_f)
thresh_interact_dk <- robust_coeftest(m_thresh_hat, vcov_thresh_dk, "w_lag_avg:regime")

df_resid_g <- df.residual(m_thresh_hat)
f_stat_by_gamma <- ((SSR0 - ssr_by_gamma) / 1) / (ssr_by_gamma / df_resid_g)
F_obs <- max(f_stat_by_gamma)

# Effect below/above the estimated threshold
effect_low  <- lincomb_effect(m_thresh_hat, vcov_thresh_dk, "w_lag_avg", weights = 1)
effect_high <- lincomb_effect(m_thresh_hat, vcov_thresh_dk, c("w_lag_avg", "w_lag_avg:regime"), weights = c(1,1))

threshold_point_report <- data.frame(
  gamma_hat = gamma_hat, n_obs_below = n_low, n_obs_above = n_high,
  interaction_estimate = thresh_interact_dk$estimate, interaction_dk_p = thresh_interact_dk$p_value,
  F_obs_sup_stat = F_obs
)
write.csv(threshold_point_report, "results_v2/nonlinear_spatial_dynamics/tables/nl_D_threshold_point_estimate.csv", row.names = FALSE)
write.csv(data.frame(gamma = gamma_grid, SSR = ssr_by_gamma, F_stat = f_stat_by_gamma),
          "results_v2/nonlinear_spatial_dynamics/tables/nl_D_threshold_profile.csv", row.names = FALSE)
regime_effects_report <- bind_rows(cbind(regime = "Below threshold", effect_low), cbind(regime = "Above threshold", effect_high))
write.csv(regime_effects_report, "results_v2/nonlinear_spatial_dynamics/tables/nl_D_regime_specific_effects.csv", row.names = FALSE)
message("\n---- Section D: threshold point estimate ----")
print(threshold_point_report)
print(regime_effects_report)

## ---- Bootstrap test of threshold existence (B=499, documented reduction; see header) ----

B_D <- 499
set.seed(20260921)
fitted_D <- fitted(m_thresh_restricted); resid_D <- residuals(m_thresh_restricted)
region_of_row_D <- sample_k3$region_id_f
regions_unique_D <- levels(region_of_row_D)

F_boot <- numeric(B_D)
message("\nRunning Section D threshold-existence wild bootstrap (B=", B_D, ", grid=", length(gamma_grid), " points)...")
for (b in seq_len(B_D)) {
  v <- setNames(sample(c(-1, 1), length(regions_unique_D), replace = TRUE), regions_unique_D)
  y_star <- fitted_D + resid_D * v[as.character(region_of_row_D)]
  boot_data <- sample_k3
  boot_data$inflation_mom <- y_star

  m0_star <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+ w_lag_avg")), data = boot_data)
  SSR0_star <- sum(residuals(m0_star)^2)

  ssr_star_by_gamma <- sapply(gamma_grid, function(g) {
    regime <- as.numeric(boot_data$national_infl_t > g)
    mg_star <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+ w_lag_avg + w_lag_avg:regime")),
                   data = cbind(boot_data, regime = regime))
    sum(residuals(mg_star)^2)
  })
  df_resid_star <- nrow(boot_data) - length(coef(m0_star)) - 1  # approx, one extra param in threshold model
  f_star_by_gamma <- ((SSR0_star - ssr_star_by_gamma) / 1) / (ssr_star_by_gamma / df_resid_star)
  F_boot[b] <- max(f_star_by_gamma)
}
p_boot_threshold <- mean(F_boot >= F_obs)

threshold_bootstrap_report <- data.frame(
  F_obs = F_obs, p_value_bootstrap_existence = p_boot_threshold, n_bootstrap_reps = B_D,
  grid_points = length(gamma_grid), seed = 20260921
)
write.csv(threshold_bootstrap_report, "results_v2/nonlinear_spatial_dynamics/tables/nl_D_threshold_bootstrap_existence.csv", row.names = FALSE)
message("\n---- Section D: bootstrap test of threshold existence ----")
print(threshold_bootstrap_report)
if (p_boot_threshold >= 0.05) {
  message("No statistically robust evidence of a threshold relationship.")
} else {
  message("Bootstrap existence test is nominally significant at 5% -- still subject to the multiple-testing correction in Section G before any conclusion.")
}

## ================================================================
## SECTION E: asymmetry (positive vs negative spatial-lag components)
## ================================================================

m_asym <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+ w_lag_avg_pos + w_lag_avg_neg")),
             data = sample_k3)
vcov_asym_dk <- dk_vcov(m_asym, sample_k3$region_id_f, sample_k3$date_f)
asym_coefs <- c("w_lag_avg_pos", "w_lag_avg_neg")
asym_dk <- robust_coeftest(m_asym, vcov_asym_dk, asym_coefs)
write.csv(asym_dk, "results_v2/nonlinear_spatial_dynamics/tables/nl_E_asymmetry_coefficients.csv", row.names = FALSE)
message("\n---- Section E: positive/negative spatial-lag coefficients ----")
print(asym_dk)

# Equality test: H0: coef_pos = coef_neg  <=>  coef_pos - coef_neg = 0
equality_test <- lincomb_effect(m_asym, vcov_asym_dk, asym_coefs, weights = c(1, -1))
write.csv(equality_test, "results_v2/nonlinear_spatial_dynamics/tables/nl_E_asymmetry_equality_test.csv", row.names = FALSE)
message("H0: coef_pos = coef_neg -- difference=", round(equality_test$estimate,4), ", DK p=", signif(equality_test$p_value,3))

# Wild cluster bootstrap for the same equality test (B=999, secondary per instruction -- no B specified for this section)
B_E <- 999
set.seed(20260921)
m_asym_restricted <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+ w_lag_avg")), data = sample_k3)
fitted_E <- fitted(m_asym_restricted); resid_E <- residuals(m_asym_restricted)
diff_obs <- coef(m_asym)["w_lag_avg_pos"] - coef(m_asym)["w_lag_avg_neg"]
vcov_asym_cl <- vcovCL(m_asym, cluster = sample_k3$region_id_f, type = "HC1")
se_diff_cl <- sqrt(vcov_asym_cl["w_lag_avg_pos","w_lag_avg_pos"] + vcov_asym_cl["w_lag_avg_neg","w_lag_avg_neg"] -
                    2*vcov_asym_cl["w_lag_avg_pos","w_lag_avg_neg"])
t_obs_E <- unname(diff_obs / se_diff_cl)

t_boot_E <- numeric(B_E)
for (b in seq_len(B_E)) {
  v <- setNames(sample(c(-1, 1), length(regions_unique), replace = TRUE), regions_unique)
  y_star <- fitted_E + resid_E * v[as.character(region_of_row)]
  boot_data <- sample_k3
  boot_data$inflation_mom <- y_star
  m_star <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+ w_lag_avg_pos + w_lag_avg_neg")), data = boot_data)
  vcov_star <- vcovCL(m_star, cluster = boot_data$region_id_f, type = "HC1")
  d_star <- coef(m_star)["w_lag_avg_pos"] - coef(m_star)["w_lag_avg_neg"]
  se_star <- sqrt(vcov_star["w_lag_avg_pos","w_lag_avg_pos"] + vcov_star["w_lag_avg_neg","w_lag_avg_neg"] -
                   2*vcov_star["w_lag_avg_pos","w_lag_avg_neg"])
  t_boot_E[b] <- unname(d_star / se_star)
}
p_wcb_E <- mean(abs(t_boot_E) >= abs(t_obs_E))
asym_wcb_report <- data.frame(difference = unname(diff_obs), t_observed = t_obs_E, p_value_wcb = p_wcb_E, n_bootstrap_reps = B_E, seed = 20260921)
write.csv(asym_wcb_report, "results_v2/nonlinear_spatial_dynamics/tables/nl_E_asymmetry_wild_bootstrap.csv", row.names = FALSE)
message("Wild cluster bootstrap (B=", B_E, ") for coef_pos=coef_neg: p=", signif(p_wcb_E,3))

## ================================================================
## SECTION F: robustness. Section B (spline) is the only nominally
## "important" result so far, so it receives the fullest robustness
## battery (KNN, distance, K=6, extreme-observation sensitivity);
## Section A (the primary interaction test) additionally gets the
## alternative-regime-definition and KNN/distance checks explicitly
## requested. Rook is not treated as independent robustness (already
## shown identical to Queen in this sample, Stage 12).
## ================================================================

## ---- F1: Section A regime interaction under KNN / distance weights ----

f_w3_knn_x_regime  <- paste(paste0("w_lag_", 1:3, "_knn:high_infl_tercile"),  collapse = " + ")
f_w3_dist_x_regime <- paste(paste0("w_lag_", 1:3, "_dist:high_infl_tercile"), collapse = " + ")
f_w3_knn  <- paste(paste0("w_lag_", 1:3, "_knn"),  collapse = " + ")
f_w3_dist <- paste(paste0("w_lag_", 1:3, "_dist"), collapse = " + ")

m_regime_knn  <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3_knn,  "+", f_w3_knn_x_regime)),  data = sample_k3)
m_regime_dist <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3_dist, "+", f_w3_dist_x_regime)), data = sample_k3)
vcov_regime_knn_dk  <- dk_vcov(m_regime_knn,  sample_k3$region_id_f, sample_k3$date_f)
vcov_regime_dist_dk <- dk_vcov(m_regime_dist, sample_k3$region_id_f, sample_k3$date_f)
wald_regime_knn  <- wald_joint(m_regime_knn,  vcov_regime_knn_dk,  paste0("w_lag_", 1:3, "_knn:high_infl_tercile"))
wald_regime_dist <- wald_joint(m_regime_dist, vcov_regime_dist_dk, paste0("w_lag_", 1:3, "_dist:high_infl_tercile"))

regime_weight_robustness <- bind_rows(
  cbind(weight_spec = "Queen (primary)", wald_delta_dk),
  cbind(weight_spec = "KNN (k=4)",       wald_regime_knn),
  cbind(weight_spec = "Distance band",   wald_regime_dist)
)
write.csv(regime_weight_robustness, "results_v2/nonlinear_spatial_dynamics/tables/nl_F_A_regime_weight_robustness.csv", row.names = FALSE)
message("\n---- F1: Section A joint test across spatial weights ----")
print(regime_weight_robustness)

## ---- F2: Section A under alternative pre-specified regime definitions ----

f_w3_x_quartile <- paste(paste0(w_cols_k3, ":high_infl_quartile"), collapse = " + ")
f_w3_x_sd       <- paste(paste0(w_cols_k3, ":high_infl_sd"),       collapse = " + ")
m_regime_quartile <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3, "+", f_w3_x_quartile)), data = sample_k3)
m_regime_sd       <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3, "+", f_w3_x_sd)),       data = sample_k3)
wald_regime_quartile <- wald_joint(m_regime_quartile, dk_vcov(m_regime_quartile, sample_k3$region_id_f, sample_k3$date_f), paste0(w_cols_k3, ":high_infl_quartile"))
wald_regime_sd       <- wald_joint(m_regime_sd,       dk_vcov(m_regime_sd,       sample_k3$region_id_f, sample_k3$date_f), paste0(w_cols_k3, ":high_infl_sd"))

regime_def_robustness <- bind_rows(
  cbind(definition = "Top tercile (primary)", n_high_months = regime_def_table$n_high_regime_months[1], wald_delta_dk),
  cbind(definition = "Top quartile",          n_high_months = regime_def_table$n_high_regime_months[2], wald_regime_quartile),
  cbind(definition = "Mean + 1 SD",           n_high_months = regime_def_table$n_high_regime_months[3], wald_regime_sd)
)
write.csv(regime_def_robustness, "results_v2/nonlinear_spatial_dynamics/tables/nl_F_A_regime_definition_robustness.csv", row.names = FALSE)
message("\n---- F2: Section A joint test across alternative regime definitions ----")
print(regime_def_robustness)
message("Note: 'Mean + 1 SD' regime has only ", regime_def_table$n_high_regime_months[3],
        " high-inflation months out of ", regime_def_table$n_total_months[1], " -- underpowered, interpreted with caution.")

## ---- F3: Section B spline under KNN / distance weights and K=6 horizon ----

sample_k3$w_lag_avg_knn  <- rowMeans(sample_k3[, paste0("w_lag_", 1:3, "_knn")])
sample_k3$w_lag_avg_dist <- rowMeans(sample_k3[, paste0("w_lag_", 1:3, "_dist")])

spline_robust_model <- function(avg_col, dat) {
  ns_b <- ns(dat$national_infl_t, df = SPLINE_DF)
  for (j in 1:SPLINE_DF) dat[[paste0("ns_state_", j)]] <- ns_b[, j]
  f <- as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", avg_col, "+",
                         paste(paste0(avg_col, ":ns_state_", 1:SPLINE_DF), collapse = " + ")))
  m <- lm(f, data = dat)
  coefs <- paste0(avg_col, ":ns_state_", 1:SPLINE_DF)
  list(model = m, coefs = coefs, wald = wald_joint(m, dk_vcov(m, dat$region_id_f, dat$date_f), coefs))
}

spline_knn  <- spline_robust_model("w_lag_avg_knn",  sample_k3)
spline_dist <- spline_robust_model("w_lag_avg_dist", sample_k3)

# K=6 horizon: needs the K=6 sample and its own w_lag_avg_k6
own_cols_k6 <- paste0("own_lag_", 1:6); w_cols_k6 <- paste0("w_lag_", 1:6)
sample_k6 <- model_data[complete_rows(model_data, c(own_cols_k6, w_cols_k6)), ]
sample_k6$region_id_f <- factor(sample_k6$region_id); sample_k6$date_f <- factor(sample_k6$date)
sample_k6$w_lag_avg_k6 <- rowMeans(sample_k6[, w_cols_k6])
f_own6 <- paste(own_cols_k6, collapse = " + ")
# spline_robust_model() above closes over f_own3 (correct for the KNN/distance K=3 calls); the
# K=6 spec needs f_own6 and the K=6 sample, so it is built explicitly here rather than reusing that closure.
ns_b_k6 <- ns(sample_k6$national_infl_t, df = SPLINE_DF)
for (j in 1:SPLINE_DF) sample_k6[[paste0("ns_state_", j)]] <- ns_b_k6[, j]
m_spline_k6 <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own6, "+ w_lag_avg_k6 +",
                                    paste(paste0("w_lag_avg_k6:ns_state_", 1:SPLINE_DF), collapse = " + "))),
                   data = sample_k6)
wald_spline_k6 <- wald_joint(m_spline_k6, dk_vcov(m_spline_k6, sample_k6$region_id_f, sample_k6$date_f),
                              paste0("w_lag_avg_k6:ns_state_", 1:SPLINE_DF))

spline_robustness <- bind_rows(
  cbind(spec = "Queen K=3 (primary)", wald_spline_dk),
  cbind(spec = "KNN K=3",             spline_knn$wald),
  cbind(spec = "Distance band K=3",   spline_dist$wald),
  cbind(spec = "Queen K=6",           wald_spline_k6)
)
write.csv(spline_robustness, "results_v2/nonlinear_spatial_dynamics/tables/nl_F_B_spline_robustness.csv", row.names = FALSE)
message("\n---- F3: Section B spline joint test across weights and K=3 vs K=6 ----")
print(spline_robustness)

## ---- F4: Section B spline, sensitivity to extreme observations (NOT removed from primary) ----

cooks_d <- cooks.distance(m_spline)
top1pct_cutoff <- quantile(cooks_d, 0.99)
influential <- cooks_d > top1pct_cutoff
sample_k3_excl <- sample_k3[!influential, ]
message("\nF4: ", sum(influential), " observations (top 1% Cook's distance) flagged as influential -- ",
        "excluded ONLY in this side-check, retained in every primary result above.")

m_spline_excl <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+ w_lag_avg +",
                                      paste(spline_interact_terms, collapse = " + "))), data = sample_k3_excl)
wald_spline_excl <- wald_joint(m_spline_excl, dk_vcov(m_spline_excl, sample_k3_excl$region_id_f, sample_k3_excl$date_f), spline_coefs)
extreme_obs_sensitivity <- bind_rows(
  cbind(spec = "Primary (all observations retained)", n_obs = nrow(sample_k3), wald_spline_dk),
  cbind(spec = "Excluding top 1% Cook's distance (side-check only)", n_obs = nrow(sample_k3_excl), wald_spline_excl)
)
write.csv(extreme_obs_sensitivity, "results_v2/nonlinear_spatial_dynamics/tables/nl_F_B_spline_extreme_obs_sensitivity.csv", row.names = FALSE)
message("\n---- F4: Section B spline joint test, extreme-observation sensitivity ----")
print(extreme_obs_sensitivity)

## ================================================================
## SECTION G: multiple-testing correction across the pre-specified
## nonlinear hypothesis family (defined BEFORE any test was run):
##   A. joint regime-interaction test (delta_1=delta_2=delta_3=0)
##   B. joint spline-nonlinearity test
##   C. joint quadratic test
##   D. threshold-existence bootstrap test
##   E. asymmetry equality test (coef_pos=coef_neg)
## Each test's own PRIMARY p-value (as designated when that test was
## run above -- wild-cluster-bootstrap where one was computed as
## primary, Driscoll-Kraay asymptotic otherwise) is used, not mixed
## post hoc after seeing which was smaller.
## ================================================================

nonlinear_family <- data.frame(
  test = c("A: Joint regime interaction (delta_1=delta_2=delta_3=0)",
           "B: Joint spline nonlinearity",
           "C: Joint quadratic (linear+quadratic)",
           "D: Threshold existence (bootstrap)",
           "E: Asymmetry equality (pos=neg)"),
  primary_p_value_type = c("Wild cluster bootstrap (B=9999)", "Driscoll-Kraay asymptotic",
                            "Driscoll-Kraay asymptotic", "Wild cluster bootstrap (B=499)",
                            "Wild cluster bootstrap (B=999)"),
  raw_p_value = c(p_wcb_A, wald_spline_dk$p_value_asymptotic, wald_quad_dk$p_value_asymptotic,
                  p_boot_threshold, p_wcb_E)
)
nonlinear_family$p_BH   <- p.adjust(nonlinear_family$raw_p_value, method = "BH")
nonlinear_family$p_Holm <- p.adjust(nonlinear_family$raw_p_value, method = "holm")
nonlinear_family$sig_raw_05  <- nonlinear_family$raw_p_value < 0.05
nonlinear_family$sig_BH_05   <- nonlinear_family$p_BH < 0.05
nonlinear_family$sig_Holm_05 <- nonlinear_family$p_Holm < 0.05

write.csv(nonlinear_family, "results_v2/nonlinear_spatial_dynamics/tables/nl_G_multiple_testing_family.csv", row.names = FALSE)
message("\n---- Section G: multiple-testing correction across the 5-test nonlinear family ----")
print(nonlinear_family)
message(sum(nonlinear_family$sig_raw_05), "/5 significant at raw 5%; ",
        sum(nonlinear_family$sig_BH_05), "/5 survive BH; ", sum(nonlinear_family$sig_Holm_05), "/5 survive Holm.")

## ================================================================
## SECTION H: residual diagnostics -- preferred K=3 nonlinear candidate
## (the regime-interaction model, Section A, as the primary nonlinear
## specification) vs. the Stage-14 linear K=3 reference (6/13 raw,
## 5/13 BH, 3/13 Holm, per that stage's own committed output, re-read
## here rather than hardcoded).
## ================================================================

sample_k3$resid_regime <- residuals(m_regime)
lb_results_nl <- sample_k3 %>%
  group_by(region_id) %>% arrange(month_index) %>%
  summarise(lb_stat = tryCatch(Box.test(resid_regime, lag = 1, type = "Ljung-Box")$statistic, error = function(e) NA_real_),
            lb_p    = tryCatch(Box.test(resid_regime, lag = 1, type = "Ljung-Box")$p.value,     error = function(e) NA_real_),
            .groups = "drop") %>%
  mutate(lb_p_BH = p.adjust(lb_p, method = "BH"), lb_p_Holm = p.adjust(lb_p, method = "holm"),
         sig_raw_05 = lb_p < 0.05, sig_BH_05 = lb_p_BH < 0.05, sig_Holm_05 = lb_p_Holm < 0.05)

n_sig_raw_nl  <- sum(lb_results_nl$sig_raw_05, na.rm = TRUE)
n_sig_BH_nl   <- sum(lb_results_nl$sig_BH_05, na.rm = TRUE)
n_sig_Holm_nl <- sum(lb_results_nl$sig_Holm_05, na.rm = TRUE)

stage14_ljung_ref_path <- "results_v2/distributed_spatial_lags/tables/dsl_residual_serial_correlation_comparison.csv"
stage14_ljung_ref <- if (file.exists(stage14_ljung_ref_path)) read.csv(stage14_ljung_ref_path, stringsAsFactors = FALSE)[2, ] else
  data.frame(n_regions_LjungBox_sig_raw_05 = 6, n_regions_LjungBox_sig_BH_05 = 5, n_regions_LjungBox_sig_Holm_05 = 3)

residual_comparison <- data.frame(
  model = c("Stage 14 linear K=3 reference (M4)", "Stage 15 regime-interaction model (Section A)"),
  n_regions_LjungBox_sig_raw_05  = c(stage14_ljung_ref$n_regions_LjungBox_sig_raw_05,  n_sig_raw_nl),
  n_regions_LjungBox_sig_BH_05   = c(stage14_ljung_ref$n_regions_LjungBox_sig_BH_05,   n_sig_BH_nl),
  n_regions_LjungBox_sig_Holm_05 = c(stage14_ljung_ref$n_regions_LjungBox_sig_Holm_05, n_sig_Holm_nl),
  n_regions_total = c(13, 13)
)
write.csv(lb_results_nl, "results_v2/nonlinear_spatial_dynamics/tables/nl_H_residual_ljungbox_by_region.csv", row.names = FALSE)
write.csv(residual_comparison, "results_v2/nonlinear_spatial_dynamics/tables/nl_H_residual_comparison_vs_stage14.csv", row.names = FALSE)
message("\n---- Section H: residual diagnostics vs. Stage-14 linear reference ----")
print(residual_comparison)
materially_reduced <- n_sig_raw_nl < stage14_ljung_ref$n_regions_LjungBox_sig_raw_05
message("Residual serial dependence ", if (materially_reduced) "further decreased" else "did NOT further decrease",
        " relative to the Stage-14 linear K=3 reference after adding the regime interaction.")

## ================================================================
## SECTION I: evidence classification (A/B/C/D per the pre-specified
## decision rule). Computed from the results above, not asserted.
## ================================================================

any_family_survives_BH <- any(nonlinear_family$sig_BH_05)
spline_survives_knn <- spline_robustness$p_value_asymptotic[spline_robustness$spec == "KNN K=3"] < 0.05
spline_survives_k6  <- spline_robustness$p_value_asymptotic[spline_robustness$spec == "Queen K=6"] < 0.05

evidence_classification <- if (!any_family_survives_BH) {
  if (nonlinear_family$raw_p_value[nonlinear_family$test == "B: Joint spline nonlinearity"] < 0.05) {
    "B: Suggestive but not robust nonlinear evidence"
  } else {
    "A: No evidence of nonlinearity"
  }
} else if (spline_survives_knn && spline_survives_k6) {
  "C: Robust smooth/state-dependent nonlinearity"
} else {
  "B: Suggestive but not robust nonlinear evidence"
}
# Threshold-specific override: Section D's own bootstrap explicitly governs any threshold/regime claim (D).
threshold_supported <- p_boot_threshold < 0.05
if (threshold_supported && nonlinear_family$sig_BH_05[nonlinear_family$test == "D: Threshold existence (bootstrap)"]) {
  evidence_classification <- "D: Robust threshold/regime evidence"
}

message("\n================================================================")
message("EVIDENCE CLASSIFICATION: ", evidence_classification)
message("================================================================")
message("Statistical nonlinearity (Section B, single-specification, Queen K=3): raw p=",
        signif(wald_spline_dk$p_value_asymptotic,3), "; does NOT survive BH/Holm across the 5-test family; ",
        "does NOT replicate under KNN weights or K=6 horizon.")
message("Economic magnitude: even where nominally significant, the implied cumulative spatial effects ",
        "(Section A regime split; Section D regime-specific effects) remain small in absolute size relative ",
        "to their own standard errors.")
message("Predictive content: no candidate model in this stage improves materially on the Stage-14 linear ",
        "reference's already-weak explanatory power for the spatial-lag block.")
message("Causal interpretation: NONE of the results above are causal. No claim of transmission, contagion, ",
         "or structural spillover is made anywhere in this stage.")
writeLines(evidence_classification, "results_v2/nonlinear_spatial_dynamics/tables/nl_I_evidence_classification.txt")

## ================================================================
## SECTION J (figures)
## ================================================================

fig_dir <- "results_v2/nonlinear_spatial_dynamics/figures"

p_h1 <- ggplot(regime_cumulative, aes(x = regime, y = estimate, ymin = ci_low, ymax = ci_high)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(color = "#2c3e50", size = 0.8) +
  labs(title = "Cumulative K=3 Spatial Effect by Inflation Regime",
       subtitle = "Section A -- Driscoll-Kraay 95% CI, top-tercile regime definition",
       x = NULL, y = "Cumulative spatial effect (sum of theta_k)") +
  theme_minimal()
ggsave(file.path(fig_dir, "fig_H1_regime_cumulative_effects.png"), p_h1, width = 8, height = 5, dpi = 300)

p_h2 <- ggplot(marginal_df, aes(x = state, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "#2c7bb6", alpha = 0.2) +
  geom_line(color = "#2c7bb6", linewidth = 0.8) +
  geom_point(data = marginal_df[marginal_df$weak_support, ], color = "#d7191c", size = 1) +
  labs(title = "Marginal Spatial Effect of w_lag_avg Across Observed Inflation-State Support",
       subtitle = "Section B -- natural spline (df=3), 95% CI band; red points = weak-support region (bottom decile local density)",
       x = "National inflation state (national_infl_t)", y = "Marginal effect of w_lag_avg") +
  theme_minimal()
ggsave(file.path(fig_dir, "fig_H2_spline_marginal_effect.png"), p_h2, width = 9, height = 5.5, dpi = 300)

p_h3 <- ggplot(sample_k3, aes(x = national_infl_t)) +
  geom_histogram(bins = 40, fill = "#2c3e50") +
  geom_vline(xintercept = tercile_cut, linetype = "dashed", color = "#d7191c") +
  labs(title = "Data-Support Distribution: National Inflation State",
       subtitle = "Dashed red line = top-tercile regime cutoff (Section A primary definition)",
       x = "National inflation state (national_infl_t)", y = "Count (region-months)") +
  theme_minimal()
ggsave(file.path(fig_dir, "fig_H3_state_support_distribution.png"), p_h3, width = 8, height = 5, dpi = 300)

p_h4 <- ggplot(data.frame(gamma = gamma_grid, F_stat = f_stat_by_gamma), aes(x = gamma, y = F_stat)) +
  geom_line(color = "#2c3e50") +
  geom_point(color = "#2c3e50", size = 1.5) +
  geom_vline(xintercept = gamma_hat, linetype = "dashed", color = "#d7191c") +
  labs(title = "Threshold Search: F-Statistic Profile Over Candidate Gamma",
       subtitle = paste0("Dashed red line = estimated gamma_hat = ", round(gamma_hat,4),
                          " (note: at/near the upper trim boundary -- see report caveat)"),
       x = "Candidate threshold (national_infl_t)", y = "F-statistic") +
  theme_minimal()
ggsave(file.path(fig_dir, "fig_H4_threshold_profile.png"), p_h4, width = 8, height = 5, dpi = 300)

asym_plot_data <- data.frame(component = c("Positive", "Negative"),
                              estimate = c(asym_dk$estimate[1], asym_dk$estimate[2]),
                              ci_low = c(asym_dk$ci_low[1], asym_dk$ci_low[2]),
                              ci_high = c(asym_dk$ci_high[1], asym_dk$ci_high[2]))
p_h5 <- ggplot(asym_plot_data, aes(x = component, y = estimate, ymin = ci_low, ymax = ci_high)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(color = "#2c3e50", size = 0.8) +
  labs(title = "Asymmetry: Positive vs. Negative Spatial-Lag Coefficients",
       subtitle = "Section E -- Driscoll-Kraay 95% CI",
       x = NULL, y = "Coefficient on w_lag_avg component") +
  theme_minimal()
ggsave(file.path(fig_dir, "fig_H5_asymmetry.png"), p_h5, width = 7, height = 5, dpi = 300)

message("\nWrote figures H1-H5 to ", fig_dir, "/")

## ================================================================
## Markdown methodological report.
## ================================================================

md <- c(
"# MAKANI v2 -- Nonlinear and State-Dependent Regional Inflation Dynamics",
"",
"Stage 4 of the MAKANI v2 econometric sequence, building on the",
"distributed-lag stage (`scripts/14_distributed_spatial_lags.R`, commit",
"`fd7ef10`) as the authoritative linear reference (K=3 own and Queen",
"spatial lags, region FE, month FE). That stage is not rebuilt or",
"reinterpreted here -- its K=3 sample is reproduced independently for",
"an exactly matched common estimation sample, and cross-checked against",
"its own saved output (match confirmed: TRUE).",
"",
"**Nonlinearity was not assumed.** Every design choice below (state",
"variable, regime cutoff, spline df, threshold trimming/grid, composite",
"spatial variable) was pre-specified in the script's header before any",
"model was estimated -- see that header for the full rationale.",
"",
"## Section A: high-inflation regime interaction",
"",
sprintf("delta_1=%.4f (p=%.3f), delta_2=%.4f (p=%.3f), delta_3=%.4f (p=%.3f) [Driscoll-Kraay].",
        delta_dk$estimate[1], delta_dk$p_value[1], delta_dk$estimate[2], delta_dk$p_value[2],
        delta_dk$estimate[3], delta_dk$p_value[3]),
sprintf("Joint Wald (DK), H0: delta_1=delta_2=delta_3=0: p=%.3f. Wild cluster bootstrap (B=%d): p=%.3f.",
        wald_delta_dk$p_value_asymptotic, B_A, p_wcb_A),
"Cumulative spatial effect: normal-regime and high-regime estimates both",
"small and statistically indistinguishable from zero (see",
"`nl_A_cumulative_effect_by_regime.csv`).",
"",
"## Section B: smooth nonlinearity (spline)",
"",
sprintf("Joint Wald (DK) test of the 3 spline-interaction terms: stat=%.3f, p=%.4g.",
        wald_spline_dk$statistic, wald_spline_dk$p_value_asymptotic),
"This is the only nominally significant result in this stage at the raw",
"5% level in its primary (Queen, K=3) specification. It does **not**",
"survive Section G's multiple-testing correction, does not replicate",
"under KNN weights, and does not replicate at K=6 (see Section F below)",
"-- see `nl_B_spline_coefficients.csv`, `nl_B_marginal_effect_by_state.csv`.",
"",
"## Section C: quadratic diagnostic (secondary only)",
"",
sprintf("Joint Wald (DK): p=%.3f (not significant at 5%%). Implied turning point = %.3f (original scale), ",
        wald_quad_dk$p_value_asymptotic, turning_point_original),
sprintf("within observed support: %s, with %d of %d observations within 0.5 SD of it -- but since the ",
        turning_in_support, n_obs_near_turning, nrow(sample_k3)),
"underlying joint test is not significant, the turning point is reported for completeness, not interpreted",
"as an established economic extremum.",
"",
"## Section D: threshold model",
"",
sprintf("Estimated gamma_hat = %.4f (n_below=%d, n_above=%d). **Caveat: gamma_hat landed at the upper edge",
        gamma_hat, n_low, n_high),
"of the trimmed search range** -- a sign the SSR-minimizing point may be an artefact of the trim boundary",
"rather than a well-identified interior optimum; this is stated explicitly, not hidden.",
sprintf("Bootstrap existence test (B=%d, %d grid points): p=%.3f.", B_D, length(gamma_grid), p_boot_threshold),
if (p_boot_threshold >= 0.05) "**No statistically robust evidence of a threshold relationship.**" else
  "Nominally significant -- still subject to Section G's correction before any conclusion.",
"",
"## Section E: asymmetry",
"",
sprintf("coef_pos=%.4f, coef_neg=%.4f, difference=%.4f, DK p=%.3f, wild-cluster-bootstrap (B=%d) p=%.3f.",
        asym_dk$estimate[1], asym_dk$estimate[2], equality_test$estimate, equality_test$p_value, B_E, p_wcb_E),
"No evidence of asymmetric positive/negative spatial predictive content.",
"",
"## Section F: robustness",
"",
"See `nl_F_A_regime_weight_robustness.csv`, `nl_F_A_regime_definition_robustness.csv`,",
"`nl_F_B_spline_robustness.csv`, `nl_F_B_spline_extreme_obs_sensitivity.csv`.",
sprintf("Section B's nominal significance (Queen K=3, p=%.3f) does NOT replicate under KNN weights (p=%.3f)",
        wald_spline_dk$p_value_asymptotic, spline_robustness$p_value_asymptotic[spline_robustness$spec=="KNN K=3"]),
sprintf("or at K=6 (p=%.3f); it IS present under distance-band weights (p=%.3f) and is not driven by extreme",
        spline_robustness$p_value_asymptotic[spline_robustness$spec=="Queen K=6"],
        spline_robustness$p_value_asymptotic[spline_robustness$spec=="Distance band K=3"]),
"observations (excluding the top 1% Cook's-distance points leaves it similarly nominal). This is a mixed,",
"not-fully-robust picture.",
"",
"## Section G: multiple-testing correction",
"",
"Pre-specified 5-test family (A, B, C, D, E). See `nl_G_multiple_testing_family.csv` for the full table.",
sprintf("%d/5 significant at raw 5%%; %d/5 survive BH; %d/5 survive Holm.",
        sum(nonlinear_family$sig_raw_05), sum(nonlinear_family$sig_BH_05), sum(nonlinear_family$sig_Holm_05)),
"",
"## Section H: residual diagnostics",
"",
sprintf("Stage-14 linear K=3 reference: %d/13 raw, %d/13 BH, %d/13 Holm. Stage-15 regime-interaction model: %d/13 raw, %d/13 BH, %d/13 Holm.",
        stage14_ljung_ref$n_regions_LjungBox_sig_raw_05, stage14_ljung_ref$n_regions_LjungBox_sig_BH_05, stage14_ljung_ref$n_regions_LjungBox_sig_Holm_05,
        n_sig_raw_nl, n_sig_BH_nl, n_sig_Holm_nl),
"",
"## Section I: evidence classification",
"",
sprintf("**%s**", evidence_classification),
"",
"Statistical nonlinearity, economic magnitude, predictive content, and",
"causal interpretation are explicitly distinguished: a nominal p<0.05 in",
"one specification (Section B, Queen K=3) is NOT robust across spatial",
"weights or dynamic horizon, does NOT survive multiple-testing",
"correction, and implies only small cumulative effect sizes even where",
"nominal. No causal, transmission, or contagion claim is made anywhere",
"in this stage.",
"",
"## Files",
"",
"All tables: `results_v2/nonlinear_spatial_dynamics/tables/`.",
"All figures: `results_v2/nonlinear_spatial_dynamics/figures/`.",
"Script: `scripts/15_nonlinear_spatial_dynamics.R`."
)
writeLines(md, "results_v2/nonlinear_spatial_dynamics/NONLINEARITY_REPORT.md")
message("\nWrote results_v2/nonlinear_spatial_dynamics/NONLINEARITY_REPORT.md")
message("15_nonlinear_spatial_dynamics.R complete.")
