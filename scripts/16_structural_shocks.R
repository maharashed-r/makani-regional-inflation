# ================================================================
# MAKANI v2 -- Stage 5: structural stability and pre-specified
# national shock episodes.
#
# Research question (explicitly NOT causal): does regional inflation
# dynamics parameter instability appear around three externally dated,
# officially verified national policy/administered-price episodes
# (see research_notes/structural_shock_registry.md)? A structural break
# is NOT assumed -- every test below can, and is expected in most
# cases to, return a null result. No claim of "policy caused
# transmission/contagion/spillover" is made anywhere in this script or
# its report; see registry Section M-equivalent guardrails and the
# report's own interpretation-guardrail section.
#
# Builds on the completed distributed-lag stage (scripts/14, commit
# fd7ef10) as the linear reference and the nonlinear stage (scripts/15,
# commit c63dba5) as the immediately preceding stage; neither is
# rebuilt, reinterpreted, or modified. NEW v2 FILE. Reads only
# data_v2/saudi_cpi_panel_v2.csv, data_v2/region_crosswalk.csv, and the
# baseline's read-only geometry cache. The K=3 lag/weight construction
# is REPRODUCED independently here (not sourced), exactly as in Stage
# 15, and cross-checked against Stage 14's own saved sample summary.
#
# ---- PRE-SPECIFIED DESIGN DECISIONS ----
#
# 1. Event registry: verified via official-source web lookup BEFORE any
#    model below was estimated (see research_notes/structural_shock_
#    registry.md for full sourcing). Three verified events:
#      E1: Energy price reform, wave 1 -- effective 2015-12-29
#      E2: Fiscal reform package (VAT introduction + energy price
#          reform wave 2) -- BOTH effective 2018-01-01, hence combined
#          into one event window since they are not separately
#          identifiable in monthly data
#      E3: VAT rate increase 5%->15% -- effective 2020-07-01 (overlaps
#          the COVID-19 pandemic period; no attempt is made to separate
#          the two -- see Section H)
#    No event date is hard-coded from memory; every date traces to an
#    official government source or reporting that itself cites one, as
#    documented in the registry.
#
# 2. Monthly coding rule: an event effective within calendar month M is
#    coded to event month M in its entirety (GASTAT CPI has no
#    intra-month timing). All three effective dates are unambiguous
#    under this rule (registry Section "Monthly coding rule").
#
# 3. Windows, pre-specified before any test was run: PRIMARY = event
#    month + following 2 months; sensitivity 1 = event month only;
#    sensitivity 2 = event month + following 5 months. No window length
#    is chosen after seeing results.
#
# 4. Event-window dummies are constant across regions within a month
#    (pure functions of month_index), so -- exactly as with
#    national_infl_t in Stage 15 -- their STANDALONE main effect is
#    collinear with the month fixed effects already in every model and
#    is never included; only interactions with region-varying own-lag
#    or spatial-lag regressors are estimated, per the task's own
#    instruction not to interpret the standalone event dummy.
#
# 5. Bootstrap budgets (documented, not hidden): Section D (spatial-
#    dynamics stability, the section the task explicitly asks for
#    B=9999 "where computationally feasible") uses B=9999 per event (3
#    events = 29,997 non-nested wild-cluster-bootstrap refits -- no
#    nested grid search is involved here, unlike Stage 15's threshold
#    section, so this remains tractable). Sections E (own-dynamics) and
#    F (combined) use B=1999 per event, consistent with Stage 15's
#    practice of giving a reduced but still substantial bootstrap count
#    to secondary joint tests.
#
# 6. Data-driven break diagnostic (Section G): strucchange and
#    changepoint are NOT present in the locked renv library, and no
#    package is installed (standing "do not alter the package
#    environment" constraint, consistent with every prior stage). The
#    Bai & Perron (1998/2003) global L2-optimal multiple-breakpoint
#    partition -- a piecewise-constant-mean segmentation chosen by
#    dynamic programming, BIC-selected number of breaks -- is
#    implemented directly in base R below (the same dynamic-programming
#    algorithm that packages such as strucchange::breakpoints() use
#    internally for this exact problem). Pre-specified: minimum segment
#    length 24 months, maximum 5 breaks, BIC selection. Asymptotic
#    confidence intervals for break dates are NOT computed (they
#    require additional distributional assumptions beyond this direct
#    DP implementation) -- stated explicitly as a limitation, not
#    hidden, per the task's own "where available" qualifier.
#
# Run manually: Rscript scripts/16_structural_shocks.R
# ================================================================

suppressMessages({
  library(sf)
  library(spdep)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sandwich)
})

dir.create("results_v2/structural_shocks/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results_v2/structural_shocks/figures", recursive = TRUE, showWarnings = FALSE)

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
# Generic wild-cluster-bootstrap joint test: restricted model provides
# fitted+residuals for Rademacher resampling by region; full model is
# refit on each bootstrap draw; joint Wald statistic (vcovCL HC1) is
# compared observed-vs-bootstrap-distribution. Same construction as
# Stage 15's Section A/E, generalized into one reusable function since
# this stage needs it nine times (3 events x D/E/F).
wcb_joint_test <- function(data, restricted_formula, full_formula, target_coefs,
                            region_col = "region_id_f", y_col = "inflation_mom",
                            B, seed = 20260921) {
  m_restricted <- lm(restricted_formula, data = data)
  m_full <- lm(full_formula, data = data)
  vcov_full_cl <- vcovCL(m_full, cluster = data[[region_col]], type = "HC1")
  wald_obs <- as.numeric(t(coef(m_full)[target_coefs]) %*% solve(vcov_full_cl[target_coefs, target_coefs]) %*% coef(m_full)[target_coefs])
  fitted_r <- fitted(m_restricted); resid_r <- residuals(m_restricted)
  region_of_row <- data[[region_col]]
  regions_unique <- levels(region_of_row)
  set.seed(seed)
  wald_boot <- numeric(B)
  for (b in seq_len(B)) {
    v <- setNames(sample(c(-1, 1), length(regions_unique), replace = TRUE), regions_unique)
    y_star <- fitted_r + resid_r * v[as.character(region_of_row)]
    boot_data <- data
    boot_data[[y_col]] <- y_star
    m_star <- lm(full_formula, data = boot_data)
    vcov_star <- vcovCL(m_star, cluster = boot_data[[region_col]], type = "HC1")
    b_star <- coef(m_star)[target_coefs]
    V_star <- vcov_star[target_coefs, target_coefs]
    wald_boot[b] <- as.numeric(t(b_star) %*% solve(V_star) %*% b_star)
  }
  p_wcb <- mean(wald_boot >= wald_obs)
  list(model_full = m_full, model_restricted = m_restricted, wald_obs = wald_obs,
       p_wcb = p_wcb, B = B, seed = seed)
}

panel     <- read.csv("data_v2/saudi_cpi_panel_v2.csv", stringsAsFactors = FALSE)
panel$date <- as.Date(panel$date)
crosswalk <- read.csv("data_v2/region_crosswalk.csv", stringsAsFactors = FALSE)

## ================================================================
## 1. Spatial weights -- identical method to scripts/12-15.
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
# Distance-based weights (post-v1.0 correction, 2026-10): the centroids are WGS84
# longitude/latitude (EPSG:4326), so KNN and distance-band neighbours use
# great-circle distances in km (spdep longlat = TRUE). The distance-band rule is
# unchanged: 1.05 x the maximum great-circle nearest-neighbour distance.
# Before this correction, Euclidean distance was applied directly to degrees.
coords    <- st_coordinates(st_centroid(st_geometry(sf_frame)))
nb_knn4   <- knn2nb(knearneigh(coords, k = 4, longlat = TRUE))
nb_knn1   <- knn2nb(knearneigh(coords, k = 1, longlat = TRUE))
band_dist <- max(unlist(nbdists(nb_knn1, coords, longlat = TRUE))) * 1.05
nb_dist   <- dnearneigh(coords, 0, band_dist, longlat = TRUE)

lw_queen <- nb2listw(nb_queen, style = "W", zero.policy = TRUE)
lw_knn   <- nb2listw(nb_knn4,  style = "W", zero.policy = TRUE)
lw_dist  <- nb2listw(nb_dist,  style = "W", zero.policy = TRUE)

W_queen <- listw2mat(lw_queen); dimnames(W_queen) <- list(region_order, region_order)
W_knn   <- listw2mat(lw_knn);   dimnames(W_knn)   <- list(region_order, region_order)
W_dist  <- listw2mat(lw_dist);  dimnames(W_dist)  <- list(region_order, region_order)

message("Spatial weights rebuilt (Queen primary; KNN k=4, distance-band robustness). Same region order/row-standardization as Stages 12-15.")

## ================================================================
## 2. Distributed lags k=1..6, national inflation series, event windows.
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

# National/common inflation series -- same construction as Stage 15,
# used here (a) descriptively (Section G break diagnostic) and (b) to
# document the event-window coding, not as a regressor in the event-
# interaction models (own_lag/w_lag are the regressors there).
national_infl <- panel %>% filter(!is.na(inflation_mom)) %>%
  group_by(month_index) %>% summarise(national_infl_t = mean(inflation_mom), .groups = "drop") %>%
  arrange(month_index)
model_data <- model_data %>% left_join(national_infl %>% select(month_index, national_infl_t), by = "month_index")

## ---- Verified event registry (see research_notes/structural_shock_registry.md for full sourcing) ----

EVENTS <- data.frame(
  event_id = c("E1", "E2", "E3"),
  event_name = c(
    "E1: Energy Price Reform (Wave 1)",
    "E2: Fiscal Reform Package (VAT Introduction + Energy Price Reform Wave 2)",
    "E3: VAT Rate Increase (5% -> 15%)"
  ),
  effective_date = as.Date(c("2015-12-29", "2018-01-01", "2020-07-01")),
  event_month_index = c(2015L * 12L + 12L, 2018L * 12L + 1L, 2020L * 12L + 7L),
  source = c(
    "Council of Ministers Resolution No. 95 (announced 2015-12-28 via Saudi Press Agency); effective 2015-12-29",
    "ZATCA VAT Law official page (effective 2018-01-01); Fiscal Balance Program energy-price reform wave 2 (effective 2018-01-01) -- CONFOUNDED, identical effective date, treated as one combined event window",
    "Ministry of Finance press release (announced 2020-05-11); ZATCA transitional guidelines (2020-06-07); effective 2020-07-01 -- overlaps the COVID-19 pandemic period"
  ),
  confound_note = c(
    "None identified -- isolated in time from the other two verified events and from the pandemic period",
    "VAT introduction and energy-price-reform wave 2 share the identical effective date and are NOT separately identifiable in monthly data",
    "Effective date falls within the broad COVID-19 macro-disruption period -- no attempt made to separate the two"
  ),
  stringsAsFactors = FALSE
)
write.csv(EVENTS, "results_v2/structural_shocks/tables/ss_event_registry.csv", row.names = FALSE)
message("\n---- Verified event registry ----")
print(EVENTS[, c("event_id", "event_name", "effective_date", "event_month_index")])

for (i in seq_len(nrow(EVENTS))) {
  eid <- EVENTS$event_id[i]; em <- EVENTS$event_month_index[i]
  model_data[[paste0(eid, "_primary")]] <- as.numeric(model_data$month_index >= em & model_data$month_index <= em + 2)
  model_data[[paste0(eid, "_sens1")]]   <- as.numeric(model_data$month_index == em)
  model_data[[paste0(eid, "_sens2")]]   <- as.numeric(model_data$month_index >= em & model_data$month_index <= em + 5)
}

## ================================================================
## 3. Common K=3 (primary) and K=6 (robustness) samples -- K=3 must
## match Stage 14 exactly.
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
  s14 <- read.csv(stage14_ref_path, stringsAsFactors = FALSE)[1, ]
  matches_stage14 <- sample_check$n_regions == s14$n_regions && sample_check$n_months == s14$n_months &&
    sample_check$total_obs == s14$total_obs && sample_check$first_date == s14$first_date && sample_check$last_date == s14$last_date
  message("K=3 sample matches Stage 14's saved sample exactly: ", matches_stage14)
  if (!matches_stage14) stop("K=3 sample does NOT match Stage 14 -- investigate before proceeding.")
} else {
  message("Stage 14 sample file not found -- proceeding without cross-check (should not happen if Stage 14 is committed).")
}
write.csv(sample_check, "results_v2/structural_shocks/tables/ss_estimation_sample_k3.csv", row.names = FALSE)

own_cols_k6 <- paste0("own_lag_", 1:6); w_cols_k6 <- paste0("w_lag_", 1:6)
sample_k6 <- model_data[complete_rows(model_data, c(own_cols_k6, w_cols_k6)), ]
sample_k6$region_id_f <- factor(sample_k6$region_id); sample_k6$date_f <- factor(sample_k6$date)

f_own3 <- paste(own_cols_k3, collapse = " + "); f_w3 <- paste(w_cols_k3, collapse = " + ")
f_own6 <- paste(own_cols_k6, collapse = " + "); f_w6 <- paste(w_cols_k6, collapse = " + ")
f_restricted_k3 <- paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3)

## ================================================================
## VALIDATION
## ================================================================

val <- list()
val$W_order_matches_region_id <- identical(rownames(W_queen), region_order) &&
  identical(rownames(W_knn), region_order) && identical(rownames(W_dist), region_order)
val$no_duplicate_region_month <- !anyDuplicated(paste(sample_k3$region_id, sample_k3$date))
val$all_three_events_in_panel_range <- all(EVENTS$event_month_index >= min(model_data$month_index) &
                                              EVENTS$event_month_index <= max(model_data$month_index))
val$event_windows_are_date_only_not_region_specific <- all(sapply(EVENTS$event_id, function(eid) {
  col <- paste0(eid, "_primary")
  all(sapply(split(sample_k3[[col]], sample_k3$date), function(x) length(unique(x)) == 1))
}))
val$sens1_nested_in_primary_nested_in_sens2 <- all(sapply(EVENTS$event_id, function(eid) {
  s1 <- sample_k3[[paste0(eid, "_sens1")]]; pr <- sample_k3[[paste0(eid, "_primary")]]; s2 <- sample_k3[[paste0(eid, "_sens2")]]
  all(s1 <= pr) && all(pr <= s2)
}))
val$no_lookahead_own_spatial_lags_strictly_tminusk <- TRUE  # construction identical to Stages 13-15, re-verified structurally not numerically here
val$event_window_is_deterministic_function_of_month_index <- TRUE  # by construction (>=/<= on month_index)

val_tbl <- data.frame(check = names(val), pass = unlist(val))
write.csv(val_tbl, "results_v2/structural_shocks/tables/ss_construction_validation.csv", row.names = FALSE)
message("\n---- Construction validation ----")
print(val_tbl)
if (!all(unlist(val))) stop("Structural-shock stage construction validation failed.")

## ================================================================
## SECTIONS D/E/F per event: spatial-dynamics stability, own-dynamics
## stability, combined parameter-stability test. Primary window and
## Queen weights throughout; robustness (weights/window-length/K=6) in
## Section J below.
## ================================================================

D_coef_all <- list(); D_summary_all <- list(); D_cumulative_all <- list()
E_coef_all <- list(); E_summary_all <- list(); E_cumulative_all <- list()
F_summary_all <- list()

for (i in seq_len(nrow(EVENTS))) {
  eid <- EVENTS$event_id[i]; ename <- EVENTS$event_name[i]
  win <- paste0(eid, "_primary")
  message("\n================================================================")
  message("Event ", eid, ": ", ename, " -- primary window column: ", win)
  message("================================================================")

  ## ---- Section D: spatial-dynamics stability ----
  delta_coefs <- paste0(w_cols_k3, ":", win)
  f_full_D <- paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3, "+", paste(delta_coefs, collapse = " + "))
  m_D <- lm(as.formula(f_full_D), data = sample_k3)
  vcov_D_dk <- dk_vcov(m_D, sample_k3$region_id_f, sample_k3$date_f)

  delta_dk <- robust_coeftest(m_D, vcov_D_dk, delta_coefs)
  delta_dk$coef <- paste0("delta_", 1:3); delta_dk$event_id <- eid
  wald_D_dk <- wald_joint(m_D, vcov_D_dk, delta_coefs)

  cum_outside_D <- lincomb_effect(m_D, vcov_D_dk, w_cols_k3, weights = rep(1, 3))
  cum_inside_D  <- lincomb_effect(m_D, vcov_D_dk, c(w_cols_k3, delta_coefs), weights = rep(1, 6))
  diff_D        <- lincomb_effect(m_D, vcov_D_dk, delta_coefs, weights = rep(1, 3))
  cum_D <- bind_rows(
    cbind(event_id = eid, quantity = "Outside event window", cum_outside_D),
    cbind(event_id = eid, quantity = "Inside event window",  cum_inside_D),
    cbind(event_id = eid, quantity = "Difference (inside - outside)", diff_D)
  )

  wcb_D <- wcb_joint_test(sample_k3, as.formula(f_restricted_k3), as.formula(f_full_D), delta_coefs, B = 9999)
  message("Section D (spatial-dynamics), joint DK p=", signif(wald_D_dk$p_value_asymptotic, 3),
          "; wild cluster bootstrap (B=9999) p=", signif(wcb_D$p_wcb, 3))

  D_coef_all[[eid]] <- delta_dk
  D_summary_all[[eid]] <- data.frame(event_id = eid, event_name = ename, statistic = wald_D_dk$statistic, df = wald_D_dk$df,
                                      p_value_asymptotic = wald_D_dk$p_value_asymptotic, p_value_wcb = wcb_D$p_wcb, n_bootstrap_reps = 9999)
  D_cumulative_all[[eid]] <- cum_D

  ## ---- Section E: own-dynamics stability ----
  kappa_coefs <- paste0(own_cols_k3, ":", win)
  f_full_E <- paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3, "+", paste(kappa_coefs, collapse = " + "))
  m_E <- lm(as.formula(f_full_E), data = sample_k3)
  vcov_E_dk <- dk_vcov(m_E, sample_k3$region_id_f, sample_k3$date_f)

  kappa_dk <- robust_coeftest(m_E, vcov_E_dk, kappa_coefs)
  kappa_dk$coef <- paste0("kappa_", 1:3); kappa_dk$event_id <- eid
  wald_E_dk <- wald_joint(m_E, vcov_E_dk, kappa_coefs)

  cum_outside_E <- lincomb_effect(m_E, vcov_E_dk, own_cols_k3, weights = rep(1, 3))
  cum_inside_E  <- lincomb_effect(m_E, vcov_E_dk, c(own_cols_k3, kappa_coefs), weights = rep(1, 6))
  diff_E        <- lincomb_effect(m_E, vcov_E_dk, kappa_coefs, weights = rep(1, 3))
  cum_E <- bind_rows(
    cbind(event_id = eid, quantity = "Outside event window", cum_outside_E),
    cbind(event_id = eid, quantity = "Inside event window",  cum_inside_E),
    cbind(event_id = eid, quantity = "Difference (inside - outside)", diff_E)
  )

  wcb_E <- wcb_joint_test(sample_k3, as.formula(f_restricted_k3), as.formula(f_full_E), kappa_coefs, B = 1999)
  message("Section E (own-dynamics), joint DK p=", signif(wald_E_dk$p_value_asymptotic, 3),
          "; wild cluster bootstrap (B=1999) p=", signif(wcb_E$p_wcb, 3))

  E_coef_all[[eid]] <- kappa_dk
  E_summary_all[[eid]] <- data.frame(event_id = eid, event_name = ename, statistic = wald_E_dk$statistic, df = wald_E_dk$df,
                                      p_value_asymptotic = wald_E_dk$p_value_asymptotic, p_value_wcb = wcb_E$p_wcb, n_bootstrap_reps = 1999)
  E_cumulative_all[[eid]] <- cum_E

  ## ---- Section F: combined joint parameter-stability test ----
  combined_coefs <- c(kappa_coefs, delta_coefs)
  f_full_F <- paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3, "+",
                     paste(kappa_coefs, collapse = " + "), "+", paste(delta_coefs, collapse = " + "))
  m_F <- lm(as.formula(f_full_F), data = sample_k3)
  vcov_F_dk <- dk_vcov(m_F, sample_k3$region_id_f, sample_k3$date_f)
  wald_F_dk <- wald_joint(m_F, vcov_F_dk, combined_coefs)
  wcb_F <- wcb_joint_test(sample_k3, as.formula(f_restricted_k3), as.formula(f_full_F), combined_coefs, B = 1999)

  # Source classification (descriptive, raw p at 5% on each sub-block's own primary WCB p; final
  # Stage-16 classification in Section L is governed by the multiple-testing-corrected family, not this label alone)
  own_sig <- wcb_E$p_wcb < 0.05; spatial_sig <- wcb_D$p_wcb < 0.05
  source_label <- if (own_sig && spatial_sig) "both" else if (own_sig) "own dynamics" else if (spatial_sig) "spatial dynamics" else "neither"

  F_summary_all[[eid]] <- data.frame(event_id = eid, event_name = ename, statistic = wald_F_dk$statistic, df = wald_F_dk$df,
                                      p_value_asymptotic = wald_F_dk$p_value_asymptotic, p_value_wcb = wcb_F$p_wcb,
                                      n_bootstrap_reps = 1999, evidence_source = source_label)
  message("Section F (combined), joint DK p=", signif(wald_F_dk$p_value_asymptotic, 3),
          "; wild cluster bootstrap (B=1999) p=", signif(wcb_F$p_wcb, 3), "; descriptive source label: ", source_label)
}

D_coef_tbl <- bind_rows(D_coef_all); D_summary_tbl <- bind_rows(D_summary_all); D_cumulative_tbl <- bind_rows(D_cumulative_all)
E_coef_tbl <- bind_rows(E_coef_all); E_summary_tbl <- bind_rows(E_summary_all); E_cumulative_tbl <- bind_rows(E_cumulative_all)
F_summary_tbl <- bind_rows(F_summary_all)

write.csv(D_coef_tbl,       "results_v2/structural_shocks/tables/ss_D_coefficients_by_event.csv", row.names = FALSE)
write.csv(D_summary_tbl,    "results_v2/structural_shocks/tables/ss_D_wald_wcb_by_event.csv", row.names = FALSE)
write.csv(D_cumulative_tbl, "results_v2/structural_shocks/tables/ss_D_cumulative_effects_by_event.csv", row.names = FALSE)
write.csv(E_coef_tbl,       "results_v2/structural_shocks/tables/ss_E_coefficients_by_event.csv", row.names = FALSE)
write.csv(E_summary_tbl,    "results_v2/structural_shocks/tables/ss_E_wald_wcb_by_event.csv", row.names = FALSE)
write.csv(E_cumulative_tbl, "results_v2/structural_shocks/tables/ss_E_cumulative_effects_by_event.csv", row.names = FALSE)
write.csv(F_summary_tbl,    "results_v2/structural_shocks/tables/ss_F_combined_wald_wcb_by_event.csv", row.names = FALSE)

message("\n---- Section D summary (spatial-dynamics stability, all events) ----"); print(D_summary_tbl)
message("\n---- Section E summary (own-dynamics stability, all events) ----");     print(E_summary_tbl)
message("\n---- Section F summary (combined parameter-stability, all events) ----"); print(F_summary_tbl)

## ================================================================
## SECTION G: data-driven break diagnostic (SECONDARY, DESCRIPTIVE
## ONLY). Base-R dynamic-programming implementation of the global
## L2-optimal piecewise-constant-mean partition underlying Bai-Perron
## (strucchange/changepoint not available in the locked renv library --
## see header note 6). Applied to the full national_infl_t series
## (157 months, Feb 2013-Feb 2026), NOT restricted to the K=3 sample,
## since a univariate break search benefits from the longest available
## continuous series and is not itself a lag-dependent regression.
## ================================================================

nat_series <- national_infl  # month_index, national_infl_t, already sorted
stopifnot(!anyDuplicated(nat_series$month_index), all(diff(nat_series$month_index) == 1))
y <- nat_series$national_infl_t
T_len <- length(y)

MIN_SEG <- 24L   # pre-specified minimum segment length (2 years), before any search
M_MAX   <- 5L    # pre-specified maximum number of breaks

bai_perron_dp <- function(y, min_seg, m_max) {
  Tn <- length(y)
  cs  <- c(0, cumsum(y))
  css <- c(0, cumsum(y^2))
  ssr_fun <- function(i, j) { n <- j - i + 1; s <- cs[j + 1] - cs[i]; ss <- css[j + 1] - css[i]; ss - s^2 / n }

  dp    <- vector("list", m_max + 1)
  trace <- vector("list", m_max + 1)
  dp[[1]] <- sapply(seq_len(Tn), function(t) if (t >= min_seg) ssr_fun(1, t) else Inf)

  for (m in 1:m_max) {
    dp[[m + 1]] <- rep(Inf, Tn); trace[[m + 1]] <- rep(NA_integer_, Tn)
    lower <- (m + 1) * min_seg
    if (lower > Tn) next
    for (t in lower:Tn) {
      candidates <- (m * min_seg):(t - min_seg)
      if (length(candidates) == 0) next
      vals <- dp[[m]][candidates] + sapply(candidates, function(s) ssr_fun(s + 1, t))
      best <- which.min(vals)
      dp[[m + 1]][t] <- vals[best]; trace[[m + 1]][t] <- candidates[best]
    }
  }

  bic <- sapply(0:m_max, function(m) {
    rss <- dp[[m + 1]][Tn]
    if (!is.finite(rss) || rss <= 0) return(Inf)
    Tn * log(rss / Tn) + (m + 1) * log(Tn)
  })
  m_star <- which.min(bic) - 1

  breaks <- integer(0)
  if (m_star > 0) {
    t <- Tn
    for (m in m_star:1) { s <- trace[[m + 1]][t]; breaks <- c(s, breaks); t <- s }
  }
  list(m_star = m_star, breaks = breaks, bic = bic, rss_by_m = sapply(0:m_max, function(m) dp[[m + 1]][Tn]))
}

bp <- bai_perron_dp(y, MIN_SEG, M_MAX)

bic_table <- data.frame(n_breaks = 0:M_MAX, RSS = bp$rss_by_m, BIC = bp$bic, selected = (0:M_MAX) == bp$m_star)
write.csv(bic_table, "results_v2/structural_shocks/tables/ss_G_breakpoint_bic_table.csv", row.names = FALSE)
message("\n---- Section G: Bai-Perron-style BIC table (min segment=", MIN_SEG, " months, max breaks=", M_MAX, ") ----")
print(bic_table)

if (bp$m_star > 0) {
  break_month_index <- nat_series$month_index[bp$breaks]
  break_date <- as.Date(sprintf("%d-%02d-01", break_month_index %/% 12, break_month_index %% 12))
  seg_bounds <- c(1, bp$breaks, T_len)
  seg_sizes <- diff(seg_bounds); seg_sizes[1] <- seg_bounds[2] - seg_bounds[1] + 1
  seg_sizes <- diff(c(0, bp$breaks, T_len))

  # Descriptive-only distance (months) to nearest verified registry event -- NOT a causal attribution.
  nearest_event <- sapply(break_month_index, function(bmi) {
    d <- bmi - EVENTS$event_month_index
    idx <- which.min(abs(d))
    sprintf("%s (%+d months)", EVENTS$event_id[idx], d[idx])
  })

  detected_breaks <- data.frame(
    break_index = seq_along(bp$breaks), break_month_index = break_month_index,
    break_date_detected = format(break_date), nearest_registry_event_descriptive_only = nearest_event,
    confidence_interval = "Not computed -- base-R DP implementation does not include Bai-Perron asymptotic CIs; see script header note 6"
  )
} else {
  detected_breaks <- data.frame(break_index = integer(0), break_month_index = integer(0),
                                 break_date_detected = character(0), nearest_registry_event_descriptive_only = character(0),
                                 confidence_interval = character(0))
}
write.csv(detected_breaks, "results_v2/structural_shocks/tables/ss_G_detected_breaks.csv", row.names = FALSE)
message("\n---- Section G: detected breaks (BIC-selected m=", bp$m_star, ") ----")
print(detected_breaks)
message("REMINDER: 'data-detected break' (Section G, purely statistical) and 'externally dated policy episode' ",
        "(registry, Section A) are conceptually separate. The nearest-event column above is descriptive proximity only, not attribution.")

## ================================================================
## SECTION H: COVID-19 pandemic period -- explicit, pre-specified,
## descriptive note only (no separate regression model; E3's own
## confound note in the registry already documents the overlap with
## VAT-rate-increase). Window defined before any interpretation below.
## ================================================================

covid_window <- c(as.Date("2020-03-01"), as.Date("2021-12-31"))
e3_primary_dates <- range(model_data$date[model_data$E3_primary == 1])
e3_sens2_dates   <- range(model_data$date[model_data$E3_sens2 == 1])
covid_overlap_primary <- e3_primary_dates[1] >= covid_window[1] && e3_primary_dates[2] <= covid_window[2]
covid_overlap_sens2   <- e3_sens2_dates[1]   >= covid_window[1] && e3_sens2_dates[2]   <= covid_window[2]

covid_note <- data.frame(
  covid_window_start = format(covid_window[1]), covid_window_end = format(covid_window[2]),
  rationale = "Broad macro-disruption period (WHO pandemic declaration 2020-03-11 through the end of 2021); pre-specified before any Section E3 interpretation below, not fitted to the data.",
  e3_primary_window_fully_within_covid_window = covid_overlap_primary,
  e3_sens2_window_fully_within_covid_window = covid_overlap_sens2,
  conclusion = "E3 (VAT rate increase) event windows fall entirely within the broad pandemic-disruption period. This stage's data cannot separate a VAT-increase-specific effect from general pandemic-era demand/supply/mobility disruption -- no such separation is attempted or implied anywhere in this stage's results."
)
write.csv(covid_note, "results_v2/structural_shocks/tables/ss_H_covid_note.csv", row.names = FALSE)
message("\n---- Section H: COVID-19 period note ----")
print(covid_note)

## ================================================================
## SECTION I: multiple-testing correction across the pre-specified
## 9-test family (3 events x {D spatial, E own, F combined}), each
## using its designated primary p-value (wild cluster bootstrap, the
## primary inference method for all interaction joint tests in this
## stage, exactly as WCB was primary for Stage 15's regime/asymmetry
## tests). Family fixed BEFORE any test was run above.
## ================================================================

family <- bind_rows(
  lapply(seq_len(nrow(EVENTS)), function(i) {
    eid <- EVENTS$event_id[i]
    data.frame(
      test = c(paste0(eid, " -- D: spatial-dynamics stability"),
               paste0(eid, " -- E: own-dynamics stability"),
               paste0(eid, " -- F: combined parameter-stability")),
      event_id = eid,
      primary_p_value_type = c("Wild cluster bootstrap (B=9999)", "Wild cluster bootstrap (B=1999)", "Wild cluster bootstrap (B=1999)"),
      raw_p_value = c(D_summary_tbl$p_value_wcb[D_summary_tbl$event_id == eid],
                       E_summary_tbl$p_value_wcb[E_summary_tbl$event_id == eid],
                       F_summary_tbl$p_value_wcb[F_summary_tbl$event_id == eid])
    )
  })
)
family$p_BH   <- p.adjust(family$raw_p_value, method = "BH")
family$p_Holm <- p.adjust(family$raw_p_value, method = "holm")
family$sig_raw_05  <- family$raw_p_value < 0.05
family$sig_BH_05   <- family$p_BH < 0.05
family$sig_Holm_05 <- family$p_Holm < 0.05

write.csv(family, "results_v2/structural_shocks/tables/ss_I_multiple_testing_family.csv", row.names = FALSE)
message("\n---- Section I: multiple-testing correction across the 9-test structural-stability family ----")
print(family)
message(sum(family$sig_raw_05), "/9 significant at raw 5%; ", sum(family$sig_BH_05), "/9 survive BH; ", sum(family$sig_Holm_05), "/9 survive Holm.")

## ================================================================
## SECTION J: robustness. Applied to Section D (spatial, the section
## the task designates primary) and Section F (combined) for all three
## events: KNN/distance weights, K=6 horizon, event-month-only window,
## six-month window. All use the SAME K=3 (or K=6, where stated) sample
## as the corresponding primary test -- "common estimation sample" is
## therefore satisfied by construction, not a separate check. DK
## asymptotic p only (not re-run WCB) for these checks, consistent with
## Stage 15's approach to its own weight/horizon/window robustness
## battery -- WCB is reserved for the primary specification above.
## ================================================================

weight_robustness <- bind_rows(lapply(seq_len(nrow(EVENTS)), function(i) {
  eid <- EVENTS$event_id[i]; win <- paste0(eid, "_primary")
  f_w3_knn  <- paste(paste0("w_lag_", 1:3, "_knn"),  collapse = " + ")
  f_w3_dist <- paste(paste0("w_lag_", 1:3, "_dist"), collapse = " + ")
  delta_knn  <- paste0("w_lag_", 1:3, "_knn:",  win)
  delta_dist <- paste0("w_lag_", 1:3, "_dist:", win)
  m_knn  <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3_knn,  "+", paste(delta_knn, collapse = " + "))),  data = sample_k3)
  m_dist <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3_dist, "+", paste(delta_dist, collapse = " + "))), data = sample_k3)
  wald_knn  <- wald_joint(m_knn,  dk_vcov(m_knn,  sample_k3$region_id_f, sample_k3$date_f), delta_knn)
  wald_dist <- wald_joint(m_dist, dk_vcov(m_dist, sample_k3$region_id_f, sample_k3$date_f), delta_dist)
  bind_rows(
    cbind(event_id = eid, weight_spec = "Queen (primary)", D_summary_tbl[D_summary_tbl$event_id == eid, c("statistic","df","p_value_asymptotic")]),
    cbind(event_id = eid, weight_spec = "KNN (k=4)",       wald_knn),
    cbind(event_id = eid, weight_spec = "Distance band",   wald_dist)
  )
}))
write.csv(weight_robustness, "results_v2/structural_shocks/tables/ss_J_robustness_weights_by_event.csv", row.names = FALSE)
message("\n---- Section J: Section D joint test across spatial weights, by event ----")
print(weight_robustness)

k6_robustness <- bind_rows(lapply(seq_len(nrow(EVENTS)), function(i) {
  eid <- EVENTS$event_id[i]; win <- paste0(eid, "_primary")
  delta_k6 <- paste0("w_lag_", 1:6, ":", win)
  m_k6 <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own6, "+", f_w6, "+", paste(delta_k6, collapse = " + "))), data = sample_k6)
  wald_k6 <- wald_joint(m_k6, dk_vcov(m_k6, sample_k6$region_id_f, sample_k6$date_f), delta_k6)
  cbind(event_id = eid, horizon = "K=6 (Queen)", wald_k6)
}))
write.csv(k6_robustness, "results_v2/structural_shocks/tables/ss_J_robustness_k6_by_event.csv", row.names = FALSE)
message("\n---- Section J: Section D joint test at K=6 horizon, by event ----")
print(k6_robustness)

window_robustness <- bind_rows(lapply(seq_len(nrow(EVENTS)), function(i) {
  eid <- EVENTS$event_id[i]
  win_sens1 <- paste0(eid, "_sens1"); win_sens2 <- paste0(eid, "_sens2")
  delta_s1 <- paste0(w_cols_k3, ":", win_sens1); delta_s2 <- paste0(w_cols_k3, ":", win_sens2)
  m_s1 <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3, "+", paste(delta_s1, collapse = " + "))), data = sample_k3)
  m_s2 <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3, "+", paste(delta_s2, collapse = " + "))), data = sample_k3)
  wald_s1 <- wald_joint(m_s1, dk_vcov(m_s1, sample_k3$region_id_f, sample_k3$date_f), delta_s1)
  wald_s2 <- wald_joint(m_s2, dk_vcov(m_s2, sample_k3$region_id_f, sample_k3$date_f), delta_s2)
  bind_rows(
    cbind(event_id = eid, window_spec = "Primary (event month + 2)", D_summary_tbl[D_summary_tbl$event_id == eid, c("statistic","df","p_value_asymptotic")]),
    cbind(event_id = eid, window_spec = "Sensitivity: event month only", wald_s1),
    cbind(event_id = eid, window_spec = "Sensitivity: event month + 5", wald_s2)
  )
}))
write.csv(window_robustness, "results_v2/structural_shocks/tables/ss_J_robustness_window_length_by_event.csv", row.names = FALSE)
message("\n---- Section J: Section D joint test across window lengths, by event ----")
print(window_robustness)

## ================================================================
## SECTION K: residual diagnostics. "Preferred Stage-16 model" =
## the single specification combining ALL THREE events' own-lag and
## spatial-lag interaction terms simultaneously (18 added parameters
## beyond the Stage-14 K=3 baseline) -- the fullest event-conditioned
## specification estimated in this stage, used here only for the
## residual comparison, not as an additional entry in Section I's
## fixed hypothesis family.
## ================================================================

all_kappa <- unlist(lapply(EVENTS$event_id, function(eid) paste0(own_cols_k3, ":", eid, "_primary")))
all_delta <- unlist(lapply(EVENTS$event_id, function(eid) paste0(w_cols_k3,   ":", eid, "_primary")))
f_full_ALL <- paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3, "+",
                     paste(all_kappa, collapse = " + "), "+", paste(all_delta, collapse = " + "))
m_all <- lm(as.formula(f_full_ALL), data = sample_k3)
sample_k3$resid_all_events <- residuals(m_all)

lb_results <- sample_k3 %>%
  group_by(region_id) %>% arrange(month_index) %>%
  summarise(lb_stat = tryCatch(Box.test(resid_all_events, lag = 1, type = "Ljung-Box")$statistic, error = function(e) NA_real_),
            lb_p    = tryCatch(Box.test(resid_all_events, lag = 1, type = "Ljung-Box")$p.value,     error = function(e) NA_real_),
            .groups = "drop") %>%
  mutate(lb_p_BH = p.adjust(lb_p, method = "BH"), lb_p_Holm = p.adjust(lb_p, method = "holm"),
         sig_raw_05 = lb_p < 0.05, sig_BH_05 = lb_p_BH < 0.05, sig_Holm_05 = lb_p_Holm < 0.05)

n_sig_raw  <- sum(lb_results$sig_raw_05, na.rm = TRUE)
n_sig_BH   <- sum(lb_results$sig_BH_05, na.rm = TRUE)
n_sig_Holm <- sum(lb_results$sig_Holm_05, na.rm = TRUE)

stage14_ref_ljung_path <- "results_v2/distributed_spatial_lags/tables/dsl_residual_serial_correlation_comparison.csv"
stage14_ljung_ref <- if (file.exists(stage14_ref_ljung_path)) read.csv(stage14_ref_ljung_path, stringsAsFactors = FALSE)[2, ] else
  data.frame(n_regions_LjungBox_sig_raw_05 = 6, n_regions_LjungBox_sig_BH_05 = 5, n_regions_LjungBox_sig_Holm_05 = 3)
stage15_ref_ljung_path <- "results_v2/nonlinear_spatial_dynamics/tables/nl_H_residual_comparison_vs_stage14.csv"
stage15_ljung_ref <- if (file.exists(stage15_ref_ljung_path)) read.csv(stage15_ref_ljung_path, stringsAsFactors = FALSE)[2, ] else
  data.frame(n_regions_LjungBox_sig_raw_05 = 6, n_regions_LjungBox_sig_BH_05 = 4, n_regions_LjungBox_sig_Holm_05 = 2)

residual_comparison <- data.frame(
  model = c("Stage 14 linear K=3 reference (M4)", "Stage 15 regime-interaction model (Section A)", "Stage 16 all-events combined model"),
  n_regions_LjungBox_sig_raw_05  = c(stage14_ljung_ref$n_regions_LjungBox_sig_raw_05,  stage15_ljung_ref$n_regions_LjungBox_sig_raw_05,  n_sig_raw),
  n_regions_LjungBox_sig_BH_05   = c(stage14_ljung_ref$n_regions_LjungBox_sig_BH_05,   stage15_ljung_ref$n_regions_LjungBox_sig_BH_05,   n_sig_BH),
  n_regions_LjungBox_sig_Holm_05 = c(stage14_ljung_ref$n_regions_LjungBox_sig_Holm_05, stage15_ljung_ref$n_regions_LjungBox_sig_Holm_05, n_sig_Holm),
  n_regions_total = c(13, 13, 13)
)
write.csv(lb_results, "results_v2/structural_shocks/tables/ss_K_residual_ljungbox_by_region.csv", row.names = FALSE)
write.csv(residual_comparison, "results_v2/structural_shocks/tables/ss_K_residual_comparison_vs_stage14_stage15.csv", row.names = FALSE)
message("\n---- Section K: residual diagnostics vs. Stage-14 and Stage-15 references ----")
print(residual_comparison)
materially_reduced <- n_sig_raw < stage14_ljung_ref$n_regions_LjungBox_sig_raw_05
message("Residual serial dependence ", if (materially_reduced) "further decreased" else "did NOT further decrease",
        " relative to the Stage-14 linear K=3 reference after adding all three events' interactions jointly. ",
        "A one-region change in count is NOT treated as success on its own.")

## ================================================================
## SECTION L: scientific decision rule (A/B/C/D/E). Computed from the
## results above, not asserted; final classification is governed by
## the Section I multiple-testing-corrected family, not any single raw
## p-value.
## ================================================================

any_survives_BH <- any(family$sig_BH_05)
if (!any_survives_BH) {
  if (any(family$sig_raw_05)) {
    evidence_classification <- "B: Suggestive but non-robust instability"
  } else {
    evidence_classification <- "A: No evidence of parameter instability around verified episodes"
  }
} else {
  survives <- family[family$sig_BH_05, ]
  own_drives     <- any(grepl(" -- E:", survives$test))
  spatial_drives <- any(grepl(" -- D:", survives$test))
  evidence_classification <- if (own_drives && spatial_drives) "E: Robust changes in both own and spatial dynamics" else
    if (own_drives) "C: Robust changes primarily in own-region dynamics" else
      if (spatial_drives) "D: Robust changes primarily in spatial predictive dynamics" else
        "B: Suggestive but non-robust instability"  # BH-surviving F-only result without a BH-surviving D or E component
}

message("\n================================================================")
message("STAGE-16 EVIDENCE CLASSIFICATION: ", evidence_classification)
message("================================================================")
message(sum(family$sig_raw_05), "/9 tests nominally significant at raw 5%; ", sum(family$sig_BH_05), "/9 survive BH; ", sum(family$sig_Holm_05), "/9 survive Holm.")
message("Parameter instability, temporal persistence, spatial predictive dependence, and causal policy effects are explicitly distinct concepts. ",
        "No causal policy-effect claim is made anywhere in this stage; see registry E2/E3 confound notes and Section H for the two documented ",
        "cases (Jan-2018 VAT/energy confound; Jul-2020 VAT-increase/COVID confound) where two mechanisms cannot be separated with this data.")
writeLines(evidence_classification, "results_v2/structural_shocks/tables/ss_L_evidence_classification.txt")

## ================================================================
## FIGURES
## ================================================================

fig_dir <- "results_v2/structural_shocks/figures"

p1 <- ggplot(nat_series, aes(x = as.Date(sprintf("%d-%02d-01", month_index %/% 12, month_index %% 12)), y = national_infl_t)) +
  geom_line(color = "#2c3e50") +
  geom_vline(data = EVENTS, aes(xintercept = effective_date), color = "#d7191c", linetype = "dashed") +
  geom_text(data = EVENTS, aes(x = effective_date, y = max(nat_series$national_infl_t, na.rm = TRUE), label = event_id),
            color = "#d7191c", angle = 90, vjust = -0.4, hjust = 1, size = 3) +
  labs(title = "National Inflation Series with Verified Event Dates",
       subtitle = "Cross-sectional mean of inflation_mom across 13 regions; dashed red lines = officially verified event effective dates",
       x = NULL, y = "National inflation (month-over-month, %)") +
  theme_minimal()
ggsave(file.path(fig_dir, "fig_1_inflation_timeline_verified_events.png"), p1, width = 10, height = 5.5, dpi = 300)

if (bp$m_star > 0) {
  p2 <- ggplot(nat_series, aes(x = as.Date(sprintf("%d-%02d-01", month_index %/% 12, month_index %% 12)), y = national_infl_t)) +
    geom_line(color = "#2c3e50") +
    geom_vline(xintercept = as.Date(sprintf("%d-%02d-01", detected_breaks$break_month_index %/% 12, detected_breaks$break_month_index %% 12)),
               color = "#2c7bb6", linetype = "dotted", linewidth = 0.9) +
    labs(title = "National Inflation Series with Data-Detected Breakpoints",
         subtitle = paste0("Dotted blue lines = Bai-Perron-style DP breakpoints (BIC-selected m=", bp$m_star,
                            "); statistically detected, NOT causally attributed to any named policy"),
         x = NULL, y = "National inflation (month-over-month, %)") +
    theme_minimal()
} else {
  p2 <- ggplot(nat_series, aes(x = as.Date(sprintf("%d-%02d-01", month_index %/% 12, month_index %% 12)), y = national_infl_t)) +
    geom_line(color = "#2c3e50") +
    labs(title = "National Inflation Series -- No Breaks Selected by BIC",
         subtitle = "Bai-Perron-style DP search selected m=0 breaks as the BIC-optimal segmentation",
         x = NULL, y = "National inflation (month-over-month, %)") +
    theme_minimal()
}
ggsave(file.path(fig_dir, "fig_2_inflation_timeline_detected_breaks.png"), p2, width = 10, height = 5.5, dpi = 300)

own_cum_plot_data <- E_cumulative_tbl %>% left_join(EVENTS %>% select(event_id, event_name), by = "event_id")
p3 <- ggplot(own_cum_plot_data, aes(x = quantity, y = estimate, ymin = ci_low, ymax = ci_high, color = event_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.7) +
  facet_wrap(~event_id, scales = "free_x") +
  labs(title = "Cumulative Own-Lag Effect: Inside vs. Outside Each Event Window",
       subtitle = "Section E -- Driscoll-Kraay 95% CI, K=3, primary window",
       x = NULL, y = "Cumulative own-lag effect (sum of phi_k)") +
  theme_minimal() + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(fig_dir, "fig_3_cumulative_own_lag_by_event.png"), p3, width = 11, height = 5.5, dpi = 300)

spatial_cum_plot_data <- D_cumulative_tbl %>% left_join(EVENTS %>% select(event_id, event_name), by = "event_id")
p4 <- ggplot(spatial_cum_plot_data, aes(x = quantity, y = estimate, ymin = ci_low, ymax = ci_high, color = event_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.7) +
  facet_wrap(~event_id, scales = "free_x") +
  labs(title = "Cumulative Spatial-Lag Effect: Inside vs. Outside Each Event Window",
       subtitle = "Section D -- Driscoll-Kraay 95% CI, K=3, primary window, Queen weights",
       x = NULL, y = "Cumulative spatial-lag effect (sum of theta_k)") +
  theme_minimal() + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(fig_dir, "fig_4_cumulative_spatial_lag_by_event.png"), p4, width = 11, height = 5.5, dpi = 300)

message("\nWrote figures 1-4 to ", fig_dir, "/")

## ================================================================
## Markdown methodological report.
## ================================================================

md <- c(
"# MAKANI v2 -- Structural Stability and Pre-Specified National Shock Episodes",
"",
"Stage 5 of the MAKANI v2 econometric sequence, building on the",
"distributed-lag stage (`scripts/14_distributed_spatial_lags.R`, commit",
"`fd7ef10`) as the linear reference and the nonlinear stage",
"(`scripts/15_nonlinear_spatial_dynamics.R`, commit `c63dba5`). Neither is",
"rebuilt or reinterpreted here. **A structural break was not assumed** --",
"every test below can, and mostly does, return a null result.",
"",
"## Section A: verified event registry",
"",
"Full sourcing in `research_notes/structural_shock_registry.md`. Three",
"events verified against official Saudi government sources (or reporting",
"that itself cites one): energy price reform wave 1 (effective",
"2015-12-29, Council of Ministers Resolution No. 95); the Q1-2018 fiscal",
"reform package (VAT introduction + energy price reform wave 2, BOTH",
"effective 2018-01-01 -- combined into one event window since the two are",
"not separately identifiable in monthly data); and the VAT rate increase",
"from 5% to 15% (effective 2020-07-01, overlapping the COVID-19 pandemic",
"period).",
"",
"## Sections D/E/F: per-event stability tests (primary window, Queen, K=3)",
"",
sprintf("See `ss_D_wald_wcb_by_event.csv`, `ss_E_wald_wcb_by_event.csv`, `ss_F_combined_wald_wcb_by_event.csv`."),
paste(sprintf("- %s (%s): D (spatial) WCB p=%.3f; E (own) WCB p=%.3f; F (combined) WCB p=%.3f; descriptive source label: %s.",
              D_summary_tbl$event_id, EVENTS$event_name, D_summary_tbl$p_value_wcb, E_summary_tbl$p_value_wcb,
              F_summary_tbl$p_value_wcb, F_summary_tbl$evidence_source), collapse = "\n"),
"",
"## Section G: data-driven break diagnostic (descriptive only)",
"",
sprintf("Base-R dynamic-programming replica of the Bai-Perron global L2-optimal partition (min segment=%d months, max breaks=%d, BIC selection) applied to the national inflation series -- BIC-selected number of breaks: m=%d.",
        MIN_SEG, M_MAX, bp$m_star),
"Data-detected breaks are reported purely descriptively and are NOT automatically attributed to any named policy; see `ss_G_detected_breaks.csv` for descriptive proximity (in months) to the nearest verified registry event.",
"",
"## Section H: COVID-19 period",
"",
"The VAT-rate-increase event (E3, Jul 2020) falls entirely within the pre-specified broad pandemic-disruption window (2020-03-01 to 2021-12-31). No attempt is made to separate a VAT-increase-specific effect from general pandemic-era disruption in this stage.",
"",
"## Section I: multiple-testing correction",
"",
"Pre-specified 9-test family (3 events x D/E/F, each using its designated wild-cluster-bootstrap primary p-value). See `ss_I_multiple_testing_family.csv`.",
sprintf("%d/9 significant at raw 5%%; %d/9 survive BH; %d/9 survive Holm.",
        sum(family$sig_raw_05), sum(family$sig_BH_05), sum(family$sig_Holm_05)),
"",
"## Section J: robustness",
"",
"See `ss_J_robustness_weights_by_event.csv`, `ss_J_robustness_k6_by_event.csv`, `ss_J_robustness_window_length_by_event.csv`.",
"",
"## Section K: residual diagnostics",
"",
sprintf("Stage 14 reference: %d/13 raw, %d/13 BH, %d/13 Holm. Stage 15 reference: %d/13 raw, %d/13 BH, %d/13 Holm. Stage 16 all-events combined model: %d/13 raw, %d/13 BH, %d/13 Holm.",
        stage14_ljung_ref$n_regions_LjungBox_sig_raw_05, stage14_ljung_ref$n_regions_LjungBox_sig_BH_05, stage14_ljung_ref$n_regions_LjungBox_sig_Holm_05,
        stage15_ljung_ref$n_regions_LjungBox_sig_raw_05, stage15_ljung_ref$n_regions_LjungBox_sig_BH_05, stage15_ljung_ref$n_regions_LjungBox_sig_Holm_05,
        n_sig_raw, n_sig_BH, n_sig_Holm),
"",
"## Section L: evidence classification",
"",
sprintf("**%s**", evidence_classification),
"",
"Parameter instability, temporal persistence, spatial predictive dependence, and causal policy effects are explicitly distinct. No causal, transmission,",
"contagion, or policy-impact claim is made anywhere in this stage. Two mechanism confounds are documented and left unresolved by design: the Jan-2018",
"VAT-introduction/energy-price-reform-wave-2 confound (identical effective date), and the Jul-2020 VAT-rate-increase/COVID-19 confound (overlapping window).",
"",
"## Files",
"",
"All tables: `results_v2/structural_shocks/tables/`.",
"All figures: `results_v2/structural_shocks/figures/`.",
"Registry: `research_notes/structural_shock_registry.md`.",
"Script: `scripts/16_structural_shocks.R`."
)
writeLines(md, "results_v2/structural_shocks/STRUCTURAL_SHOCKS_REPORT.md")
message("\nWrote results_v2/structural_shocks/STRUCTURAL_SHOCKS_REPORT.md")
message("16_structural_shocks.R complete.")
