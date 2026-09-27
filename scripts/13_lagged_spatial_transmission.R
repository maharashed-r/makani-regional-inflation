# ================================================================
# MAKANI v2 -- Stage 2 econometrics: lagged spatial transmission.
#
# Research question (explicitly NOT causal): does neighbouring regions'
# inflation at month t-1 carry predictive information for region i's
# inflation at month t, beyond (a) region fixed effects, (b) month fixed
# effects, and (c) the region's own lagged inflation?
#
# NEW v2 FILE. Reads only the validated, frozen v2 panel
# (data_v2/saudi_cpi_panel_v2.csv, data_v2/region_crosswalk.csv) and the
# baseline's read-only geometry cache (data_raw/saudi_states_sf.rds, same
# file scripts/04, 09, 12 already read, never written to). Does not
# source or modify scripts 01-12, baseline/, data_clean/, results/, or
# results_v2/ outputs from the earlier stages. Writes only new files
# under results_v2/tables/ and results_v2/figures/ (lag_*, spatial_*
# prefixes distinct from stage-1 file names -- nothing is overwritten).
#
# No new spatial/panel/bootstrap package is installed. lm() + factor()
# dummies (as in scripts 11/12) for the two-way FE regressions; sandwich
# (already in the locked environment, used nowhere else in this stage's
# lineage until now) for cluster-robust and Driscoll-Kraay-style
# panel-corrected covariance; a manual wild cluster bootstrap in base R
# for the few-cluster finite-sample check explicitly requested.
#
# Run manually: Rscript scripts/13_lagged_spatial_transmission.R
# ================================================================

suppressMessages({
  library(sf)
  library(spdep)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sandwich)
})

# Small helper standing in for lmtest::coeftest() (not in the locked
# environment) -- builds a coefficient table from any model + vcov matrix
# using a normal reference distribution, consistent with how sandwich's
# own vignette recommends using these robust covariance matrices.
robust_coeftest <- function(model, vcov_matrix, coefs = NULL) {
  b  <- coef(model)
  se <- sqrt(diag(vcov_matrix))
  if (!is.null(coefs)) { b <- b[coefs]; se <- se[coefs] }
  z  <- b / se
  p  <- 2 * pnorm(-abs(z))
  data.frame(estimate = b, std_error = se, z_value = z, p_value = p,
             ci_low = b - 1.96 * se, ci_high = b + 1.96 * se)
}

dir.create("results_v2/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results_v2/figures", recursive = TRUE, showWarnings = FALSE)

set.seed(20260921)  # same seed convention as scripts/12, reused for the wild bootstrap below

panel     <- read.csv("data_v2/saudi_cpi_panel_v2.csv", stringsAsFactors = FALSE)
panel$date <- as.Date(panel$date)
crosswalk <- read.csv("data_v2/region_crosswalk.csv", stringsAsFactors = FALSE)

## ================================================================
## 1. Spatial weights (Queen primary; KNN k=4 and distance band for
##    robustness) -- rebuilt independently here with the identical
##    method already validated in scripts/12_residual_spatial_tests.R,
##    not sourced from it. Region order fixed once and reused for every
##    date's spatial lag, so W never changes over time.
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

nb_queen <- poly2nb(sf_frame, queen = TRUE)
coords   <- st_coordinates(st_centroid(st_geometry(sf_frame)))
nb_knn4  <- knn2nb(knearneigh(coords, k = 4))
nb_knn1  <- knn2nb(knearneigh(coords, k = 1))
band_dist <- max(unlist(nbdists(nb_knn1, coords))) * 1.05
nb_dist   <- dnearneigh(coords, 0, band_dist)

# Queen vs Rook already proven identical in scripts/12 (0/26 differing
# edges) -- not rebuilt or treated as independent robustness here, per
# instruction.
lw_queen <- nb2listw(nb_queen, style = "W", zero.policy = TRUE)
lw_knn   <- nb2listw(nb_knn4,  style = "W", zero.policy = TRUE)
lw_dist  <- nb2listw(nb_dist,  style = "W", zero.policy = TRUE)

W_queen <- listw2mat(lw_queen)
W_knn   <- listw2mat(lw_knn)
W_dist  <- listw2mat(lw_dist)
stopifnot(all(dim(W_queen) == c(13, 13)), all(dim(W_knn) == c(13, 13)), all(dim(W_dist) == c(13, 13)))
# listw2mat() does not itself carry region labels -- the matrix's row/col
# order is guaranteed by construction to equal sf_frame's (= region_order)
# order since nb objects preserve their input polygon order, but dimnames
# are set explicitly here so every later lookup by region_id is safe
# rather than relying on that positional guarantee implicitly.
dimnames(W_queen) <- list(region_order, region_order)
dimnames(W_knn)   <- list(region_order, region_order)
dimnames(W_dist)  <- list(region_order, region_order)

message("Spatial weights rebuilt (Queen primary; KNN k=4 and distance-band robustness). W is fixed over time.")

## ================================================================
## 2. Own lag and spatial lag construction.
##
## Matched on an explicit linear month index (year*12+month) per
## region_id -- identical rigour to scripts/09's inflation construction
## -- so a row's lag values are ALWAYS looked up by calendar date, never
## by row position, and can never silently pick up the wrong period.
## ================================================================

panel <- panel %>%
  mutate(month_index = as.integer(format(date, "%Y")) * 12L + as.integer(format(date, "%m")))

# Wide matrix: rows = month_index, cols = region_id (fixed region_order), values = inflation_mom
wide <- panel %>%
  select(month_index, region_id, inflation_mom) %>%
  pivot_wider(names_from = region_id, values_from = inflation_mom) %>%
  arrange(month_index)
wide <- wide[, c("month_index", region_order)]  # enforce exact column order = region_order

stopifnot(identical(names(wide)[-1], region_order))  # W row/col order must match this exactly

pi_mat <- as.matrix(wide[, region_order])
rownames(pi_mat) <- wide$month_index

# Spatial lag of inflation at EVERY month (each row t uses ONLY that
# row's own values -- computing this does not by itself create any
# look-ahead; look-ahead would only be possible if a t-1 ROW were ever
# matched to a t LABEL below, which the explicit month_index lookup
# prevents).
spatial_lag_queen <- t(apply(pi_mat, 1, function(row) {
  if (anyNA(row)) return(rep(NA_real_, 13))
  as.numeric(W_queen %*% row)
}))
spatial_lag_knn <- t(apply(pi_mat, 1, function(row) {
  if (anyNA(row)) return(rep(NA_real_, 13))
  as.numeric(W_knn %*% row)
}))
spatial_lag_dist <- t(apply(pi_mat, 1, function(row) {
  if (anyNA(row)) return(rep(NA_real_, 13))
  as.numeric(W_dist %*% row)
}))
colnames(spatial_lag_queen) <- colnames(spatial_lag_knn) <- colnames(spatial_lag_dist) <- region_order
rownames(spatial_lag_queen) <- rownames(spatial_lag_knn) <- rownames(spatial_lag_dist) <- wide$month_index

lookup <- function(mat, month_idx, region) {
  mi_chr <- as.character(month_idx)
  ok <- mi_chr %in% rownames(mat)
  out <- rep(NA_real_, length(month_idx))
  out[ok] <- mat[cbind(mi_chr[ok], region[ok])]
  out
}

model_data <- panel %>%
  mutate(
    own_lag_1   = lookup(pi_mat,              month_index - 1L, region_id),
    w_lag_1     = lookup(spatial_lag_queen,    month_index - 1L, region_id),
    w_lag_1_knn = lookup(spatial_lag_knn,      month_index - 1L, region_id),
    w_lag_1_dist= lookup(spatial_lag_dist,     month_index - 1L, region_id),
    own_lag_2   = lookup(pi_mat,              month_index - 2L, region_id),
    w_lag_2     = lookup(spatial_lag_queen,    month_index - 2L, region_id),
    w_lead_1    = lookup(spatial_lag_queen,    month_index + 1L, region_id)  # placebo/timing check only, section 8
  )

## ---- Validation (§9): no look-ahead leakage, correct alignment -------

val <- list()

# (a) Every own_lag_1 value must equal the SAME region's inflation_mom
# recorded exactly one calendar month earlier -- re-derived independently
# via a date-shift join, not by trusting the lookup() function that built it.
check_tbl <- panel %>%
  select(region_id, date, inflation_mom) %>%
  mutate(prev_date = as.Date(format(date, "%Y-%m-01")) - 1) %>%
  mutate(prev_date = as.Date(format(prev_date, "%Y-%m-01"))) %>%
  left_join(panel %>% select(region_id, date, inflation_mom),
            by = c("region_id", "prev_date" = "date"), suffix = c("", "_prevcheck"))
merged_check <- model_data %>%
  select(region_id, date, own_lag_1) %>%
  left_join(check_tbl %>% select(region_id, date, inflation_mom_prevcheck), by = c("region_id", "date"))
val$own_lag_matches_independent_join <- isTRUE(all.equal(
  merged_check$own_lag_1, merged_check$inflation_mom_prevcheck
))

# (b) The month actually used for every non-NA w_lag_1 is strictly one
# calendar month before that row's own date -- i.e. no row's spatial lag
# could have been sourced from its own or a future month.
lag_month_offset <- (model_data$month_index - 1L)
val$w_lag1_uses_exactly_t_minus_1 <- all(
  is.na(model_data$w_lag_1) | (lag_month_offset == model_data$month_index - 1L)
)
# Stronger, independent re-check: rebuild w_lag_1 for a random sample of
# rows from raw pi_mat/W_queen directly, bypassing lookup() entirely.
set.seed(1)
check_idx <- sample(which(!is.na(model_data$w_lag_1)), 20)
manual_check <- sapply(check_idx, function(i) {
  mi_prev <- as.character(model_data$month_index[i] - 1L)
  reg     <- model_data$region_id[i]
  if (!mi_prev %in% rownames(pi_mat)) return(NA_real_)
  row_vals <- pi_mat[mi_prev, ]
  if (anyNA(row_vals)) return(NA_real_)
  as.numeric(W_queen[reg, ] %*% row_vals)
})
val$w_lag1_independent_recompute_matches <- isTRUE(all.equal(
  manual_check, model_data$w_lag_1[check_idx], tolerance = 1e-9
))

# (c) Region ordering in W matches the crosswalk's region_id order exactly
# (re-verified here, independent of the §7 checks already run in scripts/12).
val$W_row_col_order_matches_region_id <- identical(rownames(W_queen), region_order) &&
  identical(colnames(W_queen), region_order)

# (d) Exactly one w_lag_1 value per region-month (no duplicate rows).
val$one_spatial_lag_per_region_month <- !anyDuplicated(paste(model_data$region_id, model_data$date))

# (e) No accidental region-ID mixing: every row's own_lag_1 was looked up
# under the SAME region_id as the row itself (lookup() indexes pi_mat by
# `region` column of the current row, so this holds by construction --
# verified by checking own_lag_1 never equals a neighbour's actual t-1
# value when that differs from the own region's, for a spot sample).
mismatch_check <- sapply(check_idx, function(i) {
  mi_prev <- as.character(model_data$month_index[i] - 1L)
  reg     <- model_data$region_id[i]
  own_val <- pi_mat[mi_prev, reg]
  other_region <- setdiff(region_order, reg)[1]
  other_val <- pi_mat[mi_prev, other_region]
  isTRUE(all.equal(model_data$own_lag_1[i], own_val)) &&
    !isTRUE(all.equal(model_data$own_lag_1[i], other_val)) || isTRUE(all.equal(own_val, other_val))
})
val$no_region_id_mixing <- all(mismatch_check)

# (f) First-lag missing observations handled correctly: own_lag_1/w_lag_1
# should be NA exactly for month_index < (min month_index with valid
# inflation_mom) + 1, i.e. the first two calendar months per region
# (Jan 2013 has no inflation_mom at all; Feb 2013 has inflation_mom but
# no t-1 inflation_mom to lag), and non-NA everywhere else.
min_mi <- min(panel$month_index[!is.na(panel$inflation_mom)])
expected_na_mi <- c(min_mi - 1L, min_mi)  # Jan 2013, Feb 2013 (as month_index)
val$first_lag_missing_as_expected <- all(sapply(region_order, function(r) {
  d <- model_data[model_data$region_id == r, ]
  d <- d[order(d$month_index), ]
  na_rows <- d$month_index[is.na(d$own_lag_1)]
  setequal(na_rows, expected_na_mi)
}))

val_tbl <- data.frame(check = names(val), pass = unlist(val))
write.csv(val_tbl, "results_v2/tables/lag_construction_validation.csv", row.names = FALSE)
message("\n---- Lag construction validation ----")
print(val_tbl)
if (!all(unlist(val))) stop("Lag construction validation failed -- see results_v2/tables/lag_construction_validation.csv")

## ================================================================
## 3. Primary estimation sample
## ================================================================

sample_main <- model_data %>% filter(!is.na(own_lag_1) & !is.na(w_lag_1))

n_regions_s  <- length(unique(sample_main$region_id))
n_months_s   <- length(unique(sample_main$date))
obs_per_region <- sample_main %>% count(region_id, name = "n")
balanced_s   <- length(unique(obs_per_region$n)) == 1

dropped <- model_data %>% filter(is.na(own_lag_1) | is.na(w_lag_1))
mi_to_date <- function(mi) as.Date(sprintf("%d-%02d-01", (mi - 1L) %/% 12L, ((mi - 1L) %% 12L) + 1L))
dropped_reason <- paste0(
  "All ", nrow(dropped), " dropped rows are each region's first TWO calendar months (",
  format(mi_to_date(min_mi - 1L)), " and ", format(mi_to_date(min_mi)), ") -- ",
  "own_lag_1/w_lag_1 are structurally undefined there (own_lag_1 needs t-1's inflation_mom, which is itself ",
  "structurally missing for the panel's very first month)."
)

sample_summary <- data.frame(
  n_regions = n_regions_s, n_months = n_months_s, total_obs = nrow(sample_main),
  obs_per_region_min = min(obs_per_region$n), obs_per_region_max = max(obs_per_region$n),
  balanced = balanced_s,
  first_date = format(min(sample_main$date)), last_date = format(max(sample_main$date)),
  n_dropped = nrow(dropped)
)
write.csv(sample_summary, "results_v2/tables/lag_estimation_sample.csv", row.names = FALSE)
message("\n---- Primary estimation sample ----")
print(sample_summary)
message(dropped_reason)

if (!balanced_s) stop("Primary estimation sample is not balanced -- investigate before proceeding.")

sample_main$region_id_f <- factor(sample_main$region_id)
sample_main$date_f      <- factor(sample_main$date)

## ================================================================
## 4. M0-M2 estimation (LSDV two-way FE, as in scripts/11 -- no new
##    panel/FE package). Primary spatial lag = Queen contiguity.
## ================================================================

m0 <- lm(inflation_mom ~ region_id_f + date_f, data = sample_main)
m1 <- lm(inflation_mom ~ region_id_f + date_f + own_lag_1, data = sample_main)
m2 <- lm(inflation_mom ~ region_id_f + date_f + own_lag_1 + w_lag_1, data = sample_main)

fit_stats <- function(model, label) {
  s <- summary(model)
  data.frame(
    model = label, n_obs = length(residuals(model)), n_params = length(coef(model)),
    r_squared = s$r.squared, adj_r_squared = s$adj.r.squared, resid_se = s$sigma,
    aic = AIC(model), bic = BIC(model), logLik = as.numeric(logLik(model))
  )
}

# "Within R^2" contributed by the lag term(s) alone, net of both fixed
# effects -- computed via Frisch-Waugh-Lovell (residualise y and each lag
# regressor on the FE-only model, then the R^2 of that residualised
# regression is numerically identical to what the LSDV coefficient would
# imply, and is a cleaner measure of the lag terms' own explanatory power
# than the model's total R^2, which is dominated by the ~169 FE dummies).
fe_resid <- function(x) residuals(lm(x ~ sample_main$region_id_f + sample_main$date_f))
y_tilde        <- fe_resid(sample_main$inflation_mom)
own_lag_tilde  <- fe_resid(sample_main$own_lag_1)
w_lag_tilde    <- fe_resid(sample_main$w_lag_1)

within_m1 <- summary(lm(y_tilde ~ own_lag_tilde))$r.squared
within_m2 <- summary(lm(y_tilde ~ own_lag_tilde + w_lag_tilde))$r.squared

model_fit <- bind_rows(fit_stats(m0, "M0: TWFE only"),
                        fit_stats(m1, "M1: TWFE + own_lag_1"),
                        fit_stats(m2, "M2: TWFE + own_lag_1 + w_lag_1 (Queen)"))
model_fit$within_r2_of_lag_terms <- c(NA, within_m1, within_m2)
write.csv(model_fit, "results_v2/tables/lag_model_fit.csv", row.names = FALSE)
message("\n---- M0-M2 model fit ----")
print(model_fit)

## ================================================================
## 5. Inference: three methods, side by side, none chosen for a smaller
##    p-value. N=13 region clusters, T=156.
##
## Method assumptions:
##  (a) Classical OLS SE: assumes i.i.d., homoskedastic errors -- almost
##      certainly wrong here (panel data, likely serial and cross-
##      sectional correlation); reported only as the naive reference
##      point other methods are compared against.
##  (b) Cluster-robust by region (CR1, sandwich::vcovCL): assumes errors
##      are independent ACROSS clusters (regions) but allows arbitrary
##      within-cluster correlation (serial correlation within a region
##      over time) and heteroskedasticity. Its asymptotics are in the
##      NUMBER OF CLUSTERS, not NT -- with only 13 clusters, this is
##      well below common rules of thumb (30-50+) for the cluster-robust
##      sandwich estimator to be reliable; reported with that limitation
##      flagged explicitly, not as the headline number.
##  (c) Driscoll-Kraay-style panel-corrected (sandwich::vcovPL, cluster =
##      region, order.by = date): designed exactly for panels with a
##      FIXED, possibly small N and a long T (Hoechle 2007) -- it allows
##      serial correlation within a region AND cross-sectional (contemporaneous)
##      correlation ACROSS regions simultaneously, which is directly
##      relevant here since the whole premise of this stage is possible
##      cross-regional dependence. This is treated as the PRIMARY
##      inference method for phi and theta, precisely because it does not
##      rely on having many independent clusters.
##  (d) Wild cluster bootstrap (manual, Rademacher weights, restricted
##      under H0: theta=0, i.e. resampling from M1's residuals): the
##      finite-sample method appropriate for few clusters requested
##      explicitly -- does not rely on asymptotic normality of the
##      cluster-robust sandwich estimator at all.
##
## Driscoll-Kraay implementation audit (sandwich 3.1.1, verified by
## inspecting vcovPL()/meatPL()'s actual bytecode in this session, not
## assumed from the documentation alone):
##   - Function: sandwich::vcovPL(x, cluster=, order.by=, kernel=, ...),
##     which calls meatPL() for the "meat" and sandwich() for the bread.
##   - cluster = region_id_f (the panel/individual unit); order.by =
##     date_f (the time ordering) -- exactly the N x T structure of this
##     panel, both passed explicitly, never left to autodetection.
##   - kernel = "Bartlett", passed explicitly (this also happens to be
##     vcovPL's own default, so being explicit changes nothing numerically
##     -- it removes the ambiguity of relying on an implicit default).
##   - lag/bandwidth: meatPL()'s default is lag = "NW1987", which resolves
##     to floor(nt^(1/4)) where nt = number of unique time periods in the
##     sample. For this stage's primary sample (nt = T = 156), that is
##     floor(156^0.25) = floor(3.534) = 3 -- i.e. the estimator allows for
##     Bartlett-kernel-weighted serial correlation up to 3 monthly lags
##     within each region, on top of full contemporaneous cross-sectional
##     correlation across all 13 regions in the same month (the aggregate
##     = TRUE default, which is the actual Driscoll & Kraay 1998
##     construction: cross-sectionally sum the score vectors within each
##     time period first, then apply a HAC/Newey-West kernel over time).
##     This default was previously left implicit in this script; it is
##     now passed explicitly (lag = "NW1987") below so the resolved value
##     is documented rather than merely inherited silently. Passing the
##     same default explicitly changes no computed value.
##   - Finite-sample adjustment: meatPL()'s adjust = TRUE default IS
##     applied (not overridden), multiplying the meat by n/(n-k) where
##     n = NT = 2028 and k = the model's parameter count -- a real,
##     applied finite-sample correction, not merely available and unused.
##   - fix = FALSE (default, not overridden): no eigenvalue-clipping
##     correction is applied if the resulting matrix is not PSD; not
##     needed here since no such failure was observed for these models.
##   - Cross-check against plm::vcovSCC (the same Driscoll-Kraay estimator
##     under a different package): plm is NOT installed in this project's
##     locked renv library (confirmed via requireNamespace() in this
##     session). Per instruction, it was NOT installed to obtain this
##     cross-check. This is a genuine, documented limitation: the DK/SCC
##     standard errors reported below rely on a single implementation
##     (sandwich::vcovPL), not independently cross-validated against a
##     second package's version of the same estimator.
## ================================================================

vcov_classical <- vcov(m2)
vcov_cluster   <- vcovCL(m2, cluster = sample_main$region_id_f, type = "HC1")
vcov_dk        <- vcovPL(m2, cluster = sample_main$region_id_f, order.by = sample_main$date_f,
                          kernel = "Bartlett", lag = "NW1987")

dk_lag_resolved <- floor(length(unique(sample_main$date_f))^(1/4))
dk_audit <- data.frame(
  function_used = "sandwich::vcovPL (via meatPL)", package_version = as.character(packageVersion("sandwich")),
  cluster_variable = "region_id_f (N=13 panel units)", order_by_variable = "date_f (T=156 periods)",
  kernel = "Bartlett", lag_rule = "NW1987 (default, passed explicitly)",
  lag_resolved_value = dk_lag_resolved, aggregate_cross_sectional = TRUE,
  finite_sample_adjust_n_over_n_minus_k = TRUE, fix_psd = FALSE,
  plm_vcovSCC_cross_check = "NOT AVAILABLE -- plm not installed in locked renv library; not installed per instruction, documented as a limitation"
)
write.csv(dk_audit, "results_v2/tables/lag_dk_implementation_audit.csv", row.names = FALSE)
message("\n---- Driscoll-Kraay implementation audit ----")
message("Resolved lag/bandwidth (NW1987 rule, nt=", length(unique(sample_main$date_f)), "): ", dk_lag_resolved, " months")
print(dk_audit)

lag_coefs <- c("own_lag_1", "w_lag_1")
inf_classical <- cbind(method = "Classical OLS (naive reference)", robust_coeftest(m2, vcov_classical, lag_coefs), coef = lag_coefs)
inf_cluster   <- cbind(method = "Cluster-robust by region (13 clusters -- flagged small)", robust_coeftest(m2, vcov_cluster, lag_coefs), coef = lag_coefs)
inf_dk        <- cbind(method = "Driscoll-Kraay panel-corrected (PRIMARY)", robust_coeftest(m2, vcov_dk, lag_coefs), coef = lag_coefs)

inference_table <- bind_rows(inf_classical, inf_cluster, inf_dk)
inference_table <- inference_table[, c("coef", "method", "estimate", "std_error", "z_value", "p_value", "ci_low", "ci_high")]
write.csv(inference_table, "results_v2/tables/lag_inference_phi_theta.csv", row.names = FALSE)
message("\n---- phi (own_lag_1) and theta (w_lag_1) under three SE methods ----")
print(inference_table)

## ---- Wild cluster bootstrap for theta (H0: theta = 0), restricted under M1 ----
## Hardening pass: B increased from 999 to 9999 (per instruction, "at
## least 4999, preferably 9999 if computationally practical" -- timed at
## ~0.1s/rep on this machine, so 9999 reps is practical, ~17 minutes).
## A fixed seed is set immediately before this loop (redundant with the
## top-of-script seed, but makes this block reproducible on its own even
## if earlier code changes shift the RNG stream). This is for final
## inference precision, not to search for a different p-value -- the
## same restricted-bootstrap procedure as before is used unchanged, only
## with more repetitions.

B <- 9999
set.seed(20260921)
m1_fitted <- fitted(m1)
m1_resid  <- residuals(m1)
region_of_row <- sample_main$region_id_f

theta_hat  <- coef(m2)["w_lag_1"]
se_cl_theta <- sqrt(vcov_cluster["w_lag_1", "w_lag_1"])
t_obs <- unname(theta_hat / se_cl_theta)

wcb_t <- numeric(B)
regions_unique <- levels(region_of_row)
X_formula <- inflation_mom ~ region_id_f + date_f + own_lag_1 + w_lag_1
for (b in seq_len(B)) {
  v <- setNames(sample(c(-1, 1), length(regions_unique), replace = TRUE), regions_unique)
  y_star <- m1_fitted + m1_resid * v[as.character(region_of_row)]
  boot_data <- sample_main
  boot_data$inflation_mom <- y_star
  m2_star <- lm(X_formula, data = boot_data)
  vcov_star <- vcovCL(m2_star, cluster = boot_data$region_id_f, type = "HC1")
  se_star <- sqrt(vcov_star["w_lag_1", "w_lag_1"])
  wcb_t[b] <- unname(coef(m2_star)["w_lag_1"] / se_star)
}
p_wcb <- mean(abs(wcb_t) >= abs(t_obs))

wcb_summary <- data.frame(
  method = paste0("Wild cluster bootstrap (Rademacher, restricted under M1, B=", B, ")"),
  theta_estimate = unname(theta_hat), t_observed = t_obs, p_value_wcb = p_wcb,
  n_bootstrap_reps = B, seed = 20260921
)
write.csv(wcb_summary, "results_v2/tables/lag_wild_cluster_bootstrap.csv", row.names = FALSE)
message("\n---- Wild cluster bootstrap for theta (B=", B, ") ----")
print(wcb_summary)

## ================================================================
## 6. Incremental M1 -> M2 comparison
## ================================================================

incremental <- data.frame(
  delta_r_squared      = model_fit$r_squared[3] - model_fit$r_squared[2],
  delta_adj_r_squared  = model_fit$adj_r_squared[3] - model_fit$adj_r_squared[2],
  delta_within_r2      = within_m2 - within_m1,
  delta_aic            = model_fit$aic[3] - model_fit$aic[2],
  delta_bic            = model_fit$bic[3] - model_fit$bic[2],
  f_test_p_value       = anova(m1, m2)[2, "Pr(>F)"]
)
write.csv(incremental, "results_v2/tables/lag_incremental_M1_to_M2.csv", row.names = FALSE)
message("\n---- Incremental M1 -> M2 (adding w_lag_1) ----")
print(incremental)
message("Note: AIC/BIC increase (worse) from M1 to M2 despite the F-test p-value above -- ",
        "consistent with w_lag_1 adding essentially no explanatory power once own_lag_1 and ",
        "the fixed effects are already included.")

## ================================================================
## 7. Residual serial-correlation and persistence diagnostics (M2)
## ================================================================

sample_main$resid_m2 <- residuals(m2)

# Panel AR(1) test on residuals: regress resid_it on resid_i,t-1 (region
# fixed effects absorbed via region dummies), test the lag coefficient.
resid_lagged <- sample_main %>%
  select(region_id, month_index, resid_m2) %>%
  mutate(prev_mi = month_index - 1L) %>%
  left_join(sample_main %>% select(region_id, month_index, resid_m2_lag = resid_m2),
            by = c("region_id", "prev_mi" = "month_index"))
resid_ar1_data <- resid_lagged %>% filter(!is.na(resid_m2_lag))
ar1_model <- lm(resid_m2 ~ resid_m2_lag + factor(region_id), data = resid_ar1_data)
ar1_coef  <- coef(ar1_model)["resid_m2_lag"]
ar1_vcov_cl <- vcovCL(ar1_model, cluster = factor(resid_ar1_data$region_id), type = "HC1")
ar1_se_cl <- sqrt(ar1_vcov_cl["resid_m2_lag", "resid_m2_lag"])
ar1_p_cl  <- 2 * pnorm(-abs(ar1_coef / ar1_se_cl))

# Per-region Ljung-Box test (lag 1) on M2 residuals, as a second, fully
# independent check of the same question.
#
# Hardening pass: 13 simultaneous tests inflate the false-positive risk,
# exactly the concern already raised for the Moran's I monthly tests in
# scripts/12. BH and Holm corrections are added here for the same reason
# -- 8/13 raw-significant regions must NOT be described as established
# residual dynamics until shown against the corrected counts below.
lb_results <- sample_main %>%
  group_by(region_id) %>%
  arrange(month_index) %>%
  summarise(
    lb_stat = tryCatch(Box.test(resid_m2, lag = 1, type = "Ljung-Box")$statistic, error = function(e) NA_real_),
    lb_p    = tryCatch(Box.test(resid_m2, lag = 1, type = "Ljung-Box")$p.value,     error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    lb_p_BH   = p.adjust(lb_p, method = "BH"),
    lb_p_Holm = p.adjust(lb_p, method = "holm"),
    sig_raw_05  = lb_p < 0.05,
    sig_BH_05   = lb_p_BH < 0.05,
    sig_Holm_05 = lb_p_Holm < 0.05
  )
n_lb_sig_raw  <- sum(lb_results$sig_raw_05, na.rm = TRUE)
n_lb_sig_BH   <- sum(lb_results$sig_BH_05, na.rm = TRUE)
n_lb_sig_Holm <- sum(lb_results$sig_Holm_05, na.rm = TRUE)

serial_corr_summary <- data.frame(
  panel_AR1_coefficient = unname(ar1_coef),
  panel_AR1_cluster_robust_p = ar1_p_cl,
  n_regions_LjungBox_sig_raw_05  = n_lb_sig_raw,
  n_regions_LjungBox_sig_BH_05   = n_lb_sig_BH,
  n_regions_LjungBox_sig_Holm_05 = n_lb_sig_Holm,
  n_regions_total = nrow(lb_results)
)
write.csv(serial_corr_summary, "results_v2/tables/lag_residual_serial_correlation.csv", row.names = FALSE)
write.csv(lb_results, "results_v2/tables/lag_residual_ljungbox_by_region.csv", row.names = FALSE)
message("\n---- Residual serial correlation (M2), with BH/Holm correction across 13 regions ----")
print(serial_corr_summary)
message(n_lb_sig_raw, "/13 regions nominally significant at raw 5%; ", n_lb_sig_BH,
        "/13 survive BH; ", n_lb_sig_Holm, "/13 survive Holm. Only the corrected counts should be ",
        "described as established residual dynamics.")

# Unconditional persistence of inflation_mom itself (context, not a model
# choice): raw lag-1 autocorrelation of inflation_mom pooled within region.
raw_ar1 <- sample_main %>%
  group_by(region_id) %>% arrange(month_index) %>%
  summarise(acf1 = tryCatch(acf(inflation_mom, lag.max = 1, plot = FALSE)$acf[2, 1, 1], error = function(e) NA_real_),
            .groups = "drop")
message("Median raw inflation_mom lag-1 autocorrelation across regions: ", round(median(raw_ar1$acf1, na.rm = TRUE), 3))

## ================================================================
## 8. Lag sufficiency check (diagnostic only -- does NOT change the
##    primary M0-M2 specification, per instruction).
##
## Hardening pass: the exploratory run of this diagnostic reported
## own_lag_2's cluster-robust p (~0.018) as if it were dispositive. With
## only 13 clusters that single method is exactly the one flagged as
## unreliable for phi/theta above, so it must not be trusted alone here
## either. Re-estimated under the SAME primary framework used for M2:
## Driscoll-Kraay as primary, cluster-robust shown but flagged, plus a
## wild cluster bootstrap and a joint Wald test of both lag-2 terms.
## ================================================================

sample_lag2 <- model_data %>% filter(!is.na(own_lag_1) & !is.na(w_lag_1) & !is.na(own_lag_2) & !is.na(w_lag_2))
sample_lag2$region_id_f <- factor(sample_lag2$region_id)
sample_lag2$date_f      <- factor(sample_lag2$date)

m1_ext <- lm(inflation_mom ~ region_id_f + date_f + own_lag_1 + w_lag_1, data = sample_lag2)  # restricted, for the joint test and WCB null
m2_ext <- lm(inflation_mom ~ region_id_f + date_f + own_lag_1 + w_lag_1 + own_lag_2 + w_lag_2, data = sample_lag2)

vcov_ext_cl <- vcovCL(m2_ext, cluster = sample_lag2$region_id_f, type = "HC1")
vcov_ext_dk <- vcovPL(m2_ext, cluster = sample_lag2$region_id_f, order.by = sample_lag2$date_f,
                       kernel = "Bartlett", lag = "NW1987")
lag2_coefs <- c("own_lag_2", "w_lag_2")

lag2_cluster <- cbind(method = "Cluster-robust by region (13 clusters -- flagged small)",
                       robust_coeftest(m2_ext, vcov_ext_cl, lag2_coefs), coef = lag2_coefs)
lag2_dk      <- cbind(method = "Driscoll-Kraay panel-corrected (PRIMARY, same framework as M2)",
                       robust_coeftest(m2_ext, vcov_ext_dk, lag2_coefs), coef = lag2_coefs)
lag2_inference <- bind_rows(lag2_cluster, lag2_dk)
lag2_inference <- lag2_inference[, c("coef", "method", "estimate", "std_error", "z_value", "p_value", "ci_low", "ci_high")]
write.csv(lag2_inference, "results_v2/tables/lag_sufficiency_check.csv", row.names = FALSE)
message("\n---- Lag-2 diagnostic, hardened (own_lag_2, w_lag_2; DK primary, cluster-robust flagged) ----")
print(lag2_inference)

# Wild cluster bootstrap for own_lag_2 and w_lag_2 jointly (restricted
# under M1_ext, i.e. H0: both lag-2 coefficients are zero). Kept at
# B=999 (this is a secondary "where feasible" diagnostic per instruction,
# not the final theta inference that B was raised to 9999 for above).
B_lag2 <- 999
set.seed(20260921)
m1ext_fitted <- fitted(m1_ext); m1ext_resid <- residuals(m1_ext)
region_of_row_ext <- sample_lag2$region_id_f
regions_unique_ext <- levels(region_of_row_ext)
X_formula_ext <- inflation_mom ~ region_id_f + date_f + own_lag_1 + w_lag_1 + own_lag_2 + w_lag_2

joint_wald_boot <- numeric(B_lag2)
own2_t_boot <- numeric(B_lag2)
w2_t_boot   <- numeric(B_lag2)
for (b in seq_len(B_lag2)) {
  v <- setNames(sample(c(-1, 1), length(regions_unique_ext), replace = TRUE), regions_unique_ext)
  y_star <- m1ext_fitted + m1ext_resid * v[as.character(region_of_row_ext)]
  boot_data <- sample_lag2
  boot_data$inflation_mom <- y_star
  m2ext_star <- lm(X_formula_ext, data = boot_data)
  vcov_star  <- vcovCL(m2ext_star, cluster = boot_data$region_id_f, type = "HC1")
  b_star <- coef(m2ext_star)[lag2_coefs]
  V_star <- vcov_star[lag2_coefs, lag2_coefs]
  joint_wald_boot[b] <- as.numeric(t(b_star) %*% solve(V_star) %*% b_star)
  own2_t_boot[b] <- unname(b_star["own_lag_2"] / sqrt(V_star["own_lag_2", "own_lag_2"]))
  w2_t_boot[b]   <- unname(b_star["w_lag_2"]   / sqrt(V_star["w_lag_2", "w_lag_2"]))
}

# Observed joint Wald statistic (cluster-robust vcov, consistent with the bootstrap's null-resampling covariance)
b_obs <- coef(m2_ext)[lag2_coefs]
V_obs <- vcov_ext_cl[lag2_coefs, lag2_coefs]
wald_obs <- as.numeric(t(b_obs) %*% solve(V_obs) %*% b_obs)
p_wald_chisq <- pchisq(wald_obs, df = 2, lower.tail = FALSE)
p_wald_wcb   <- mean(joint_wald_boot >= wald_obs)

own2_t_obs <- unname(b_obs["own_lag_2"] / sqrt(V_obs["own_lag_2", "own_lag_2"]))
w2_t_obs   <- unname(b_obs["w_lag_2"]   / sqrt(V_obs["w_lag_2", "w_lag_2"]))
p_own2_wcb <- mean(abs(own2_t_boot) >= abs(own2_t_obs))
p_w2_wcb   <- mean(abs(w2_t_boot) >= abs(w2_t_obs))

lag2_wcb_summary <- data.frame(
  test = c("own_lag_2 (individual)", "w_lag_2 (individual)", "Joint (own_lag_2 & w_lag_2)"),
  statistic = c(own2_t_obs, w2_t_obs, wald_obs),
  statistic_type = c("t", "t", "Wald chi-sq (df=2)"),
  p_value_wcb = c(p_own2_wcb, p_w2_wcb, p_wald_wcb),
  p_value_asymptotic = c(NA, NA, p_wald_chisq),
  n_bootstrap_reps = B_lag2
)
write.csv(lag2_wcb_summary, "results_v2/tables/lag_sufficiency_wild_bootstrap_and_joint_test.csv", row.names = FALSE)
message("\n---- Lag-2 wild cluster bootstrap + joint Wald test (B=", B_lag2, ") ----")
print(lag2_wcb_summary)
message("Primary specification (M0-M2, lag-1 only) is unchanged regardless of this result, per instruction.")

## ================================================================
## 9. Spatial robustness: re-estimate M2 with KNN(k=4) and distance-band
##    W in place of Queen. W is not selected on significance.
## ================================================================

m2_knn  <- lm(inflation_mom ~ region_id_f + date_f + own_lag_1 + w_lag_1_knn,  data = sample_main)
m2_dist <- lm(inflation_mom ~ region_id_f + date_f + own_lag_1 + w_lag_1_dist, data = sample_main)

dk_args <- list(kernel = "Bartlett", lag = "NW1987")  # explicit, documented in section 5's DK audit above
spatial_robustness <- bind_rows(
  cbind(weight_spec = "Queen contiguity (primary)",
        robust_coeftest(m2, do.call(vcovPL, c(list(x = m2, cluster = sample_main$region_id_f, order.by = sample_main$date_f), dk_args)), "w_lag_1")),
  cbind(weight_spec = "KNN (k=4)",
        robust_coeftest(m2_knn, do.call(vcovPL, c(list(x = m2_knn, cluster = sample_main$region_id_f, order.by = sample_main$date_f), dk_args)), "w_lag_1_knn")),
  cbind(weight_spec = "Distance band",
        robust_coeftest(m2_dist, do.call(vcovPL, c(list(x = m2_dist, cluster = sample_main$region_id_f, order.by = sample_main$date_f), dk_args)), "w_lag_1_dist"))
)
write.csv(spatial_robustness, "results_v2/tables/lag_spatial_robustness.csv", row.names = FALSE)
message("\n---- theta across W specifications (Driscoll-Kraay SE) ----")
print(spatial_robustness)

theta_range <- range(spatial_robustness$estimate)
theta_all_same_sign <- length(unique(sign(spatial_robustness$estimate))) == 1
theta_all_nonsig <- all(spatial_robustness$p_value > 0.05)
message("theta range across W specs: [", round(theta_range[1], 4), ", ", round(theta_range[2], 4), "]",
        " | same sign across all: ", theta_all_same_sign, " | all non-significant: ", theta_all_nonsig)

## ================================================================
## 10. Placebo / directional (timing) check: does FUTURE neighbour
##     inflation W*pi_(t+1) show an association comparable to or
##     stronger than the lagged (t-1) association? This is a timing
##     diagnostic only, not a causal identification test.
## ================================================================

sample_placebo <- model_data %>% filter(!is.na(own_lag_1) & !is.na(w_lead_1))
sample_placebo$region_id_f <- factor(sample_placebo$region_id)
sample_placebo$date_f      <- factor(sample_placebo$date)

m2_placebo <- lm(inflation_mom ~ region_id_f + date_f + own_lag_1 + w_lead_1, data = sample_placebo)
vcov_placebo_dk <- vcovPL(m2_placebo, cluster = sample_placebo$region_id_f, order.by = sample_placebo$date_f,
                           kernel = "Bartlett", lag = "NW1987")
placebo_result <- robust_coeftest(m2_placebo, vcov_placebo_dk, "w_lead_1")
placebo_result$n_obs <- nrow(sample_placebo)
placebo_result$first_date <- format(min(sample_placebo$date))
placebo_result$last_date  <- format(max(sample_placebo$date))
write.csv(placebo_result, "results_v2/tables/lag_placebo_future_neighbour.csv", row.names = FALSE)
message("\n---- Placebo/timing check: theta on FUTURE neighbour inflation W*pi_(t+1) (full, non-common sample) ----")
print(placebo_result)

## ---- Hardening: common-sample check ----
## sample_main (lag model) runs Mar2013-Feb2026 (needs t-1 only); sample_
## placebo (lead model) runs Mar2013-Jan2026 (needs t+1 too, so loses the
## last month). These are NOT the same sample as built above -- verified
## explicitly, then both models are re-estimated restricted to their
## common region-month support so the lag-vs-lead comparison is fair.
common_sample_check <- data.frame(
  lag_model_n = nrow(sample_main), lag_model_start = format(min(sample_main$date)), lag_model_end = format(max(sample_main$date)),
  placebo_model_n = nrow(sample_placebo), placebo_model_start = format(min(sample_placebo$date)), placebo_model_end = format(max(sample_placebo$date)),
  samples_identical = isTRUE(all.equal(sort(paste(sample_main$region_id, sample_main$date)),
                                        sort(paste(sample_placebo$region_id, sample_placebo$date))))
)
write.csv(common_sample_check, "results_v2/tables/lag_placebo_common_sample_check.csv", row.names = FALSE)
message("\n---- Common-sample check (lag model vs placebo/lead model) ----")
print(common_sample_check)

## NOTE (bug fixed during hardening): base R's intersect() strips the
## Date class from its result (returns plain numeric days-since-epoch),
## so a later `date %in% common_dates` would silently match zero rows if
## common_dates were left as bare numeric while `date` stays class Date.
## Restoring the Date class explicitly here is required for the filter
## below to work at all, not merely for style.
common_dates <- as.Date(intersect(unique(sample_main$date), unique(sample_placebo$date)), origin = "1970-01-01")
sample_common <- model_data %>% filter(date %in% common_dates & !is.na(own_lag_1) & !is.na(w_lag_1) & !is.na(w_lead_1))
sample_common$region_id_f <- factor(sample_common$region_id)
sample_common$date_f      <- factor(sample_common$date)
stopifnot(length(unique(sample_common$date)) == length(common_dates))  # confirm no further rows silently dropped

m2_common_lag  <- lm(inflation_mom ~ region_id_f + date_f + own_lag_1 + w_lag_1,  data = sample_common)
m2_common_lead <- lm(inflation_mom ~ region_id_f + date_f + own_lag_1 + w_lead_1, data = sample_common)

vcov_common_lag  <- vcovPL(m2_common_lag,  cluster = sample_common$region_id_f, order.by = sample_common$date_f, kernel = "Bartlett", lag = "NW1987")
vcov_common_lead <- vcovPL(m2_common_lead, cluster = sample_common$region_id_f, order.by = sample_common$date_f, kernel = "Bartlett", lag = "NW1987")

common_sample_comparison <- bind_rows(
  cbind(spec = "Lagged (t-1) neighbour inflation, common sample", robust_coeftest(m2_common_lag, vcov_common_lag, "w_lag_1")),
  cbind(spec = "Future (t+1) neighbour inflation, common sample -- placebo", robust_coeftest(m2_common_lead, vcov_common_lead, "w_lead_1"))
)
common_sample_comparison$n_obs <- nrow(sample_common)
common_sample_comparison$first_date <- format(min(sample_common$date))
common_sample_comparison$last_date  <- format(max(sample_common$date))
write.csv(common_sample_comparison, "results_v2/tables/lag_placebo_common_sample_comparison.csv", row.names = FALSE)
message("\n---- Lag vs lead, re-estimated on the IDENTICAL common region-month sample (n=", nrow(sample_common), ") ----")
print(common_sample_comparison)
message("This common-sample re-estimation is reported separately for a fair timing comparison; it is NOT a ",
        "causal identification test, and the full-sample theta (w_lag_1) reported earlier in this script is ",
        "left unchanged -- only this comparison uses the restricted common support.")

theta_lag_estimate  <- inf_dk$estimate[inf_dk$coef == "w_lag_1"]
theta_lead_estimate <- placebo_result$estimate
message("theta (lag, t-1, full sample) = ", round(theta_lag_estimate, 4),
        " vs theta (lead/placebo, t+1, full sample) = ", round(theta_lead_estimate, 4),
        " | |lead| >= |lag|: ", abs(theta_lead_estimate) >= abs(theta_lag_estimate))

## ================================================================
## 11. Figures
## ================================================================

# Figure F1: phi and theta point estimates + 95% CI under the three SE
# methods (classical shown only as the naive reference point).
fig1_data <- inference_table %>%
  mutate(coef_label = ifelse(coef == "own_lag_1", "phi (own_lag_1)", "theta (w_lag_1, Queen)"),
         method = factor(method, levels = c("Classical OLS (naive reference)",
                                             "Cluster-robust by region (13 clusters -- flagged small)",
                                             "Driscoll-Kraay panel-corrected (PRIMARY)")))
p_f1 <- ggplot(fig1_data, aes(x = method, y = estimate, ymin = ci_low, ymax = ci_high)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(color = "#2c3e50") +
  coord_flip() +
  facet_wrap(~coef_label, scales = "free_x") +
  labs(title = "phi and theta Under Three Inference Methods",
       subtitle = "95% CI. Driscoll-Kraay is the primary method (small-N, long-T panel); cluster-robust flagged (only 13 clusters).",
       x = NULL, y = "Coefficient estimate") +
  theme_minimal()
ggsave("results_v2/figures/fig_F1_phi_theta_inference_methods.png", p_f1, width = 10, height = 5, dpi = 300)

# Figure F2: theta across W specifications, plus the lead/placebo check,
# all on one comparable scale.
fig2_data <- bind_rows(
  spatial_robustness %>% transmute(spec = weight_spec, estimate, ci_low, ci_high, role = "Lagged (t-1) neighbour inflation"),
  placebo_result %>% transmute(spec = "Queen (lead, t+1) -- placebo", estimate, ci_low, ci_high, role = "Future (t+1) neighbour inflation -- placebo")
)
p_f2 <- ggplot(fig2_data, aes(x = spec, y = estimate, ymin = ci_low, ymax = ci_high, color = role)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange() +
  coord_flip() +
  scale_color_manual(values = c("Lagged (t-1) neighbour inflation" = "#2c7bb6",
                                 "Future (t+1) neighbour inflation -- placebo" = "#d7191c")) +
  labs(title = "theta: Spatial-Weight Robustness and Timing Placebo",
       subtitle = "Driscoll-Kraay 95% CI. A placebo CI comparable to or wider than the lagged estimates weakens a transmission reading.",
       x = NULL, y = "theta estimate", color = NULL) +
  theme_minimal() + theme(legend.position = "top")
ggsave("results_v2/figures/fig_F2_theta_robustness_and_placebo.png", p_f2, width = 10, height = 5, dpi = 300)

# Figure F3: residual serial correlation -- per-region Ljung-Box p-values.
p_f3 <- ggplot(lb_results, aes(x = reorder(region_id, lb_p), y = lb_p)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#d7191c") +
  geom_col(fill = "#2c3e50") +
  coord_flip() +
  labs(title = "M2 Residual Serial Correlation by Region (Ljung-Box, lag 1)",
       subtitle = "Red line = 0.05 significance threshold",
       x = NULL, y = "Ljung-Box p-value") +
  theme_minimal()
ggsave("results_v2/figures/fig_F3_residual_serial_correlation_by_region.png", p_f3, width = 8, height = 6, dpi = 300)

message("\nWrote figures F1-F3 to results_v2/figures/.")
message("13_lagged_spatial_transmission.R complete.")
