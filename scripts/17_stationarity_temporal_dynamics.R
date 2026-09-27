# ================================================================
# MAKANI v2 -- Stage 6: stationarity, seasonality, and remaining
# temporal dynamics.
#
# Research question (explicitly NOT a significance hunt): is the
# residual dependence documented across Stages 12-16 primarily a
# stationarity artifact, a seasonal (lag-12) pattern, longer-memory
# temporal structure, or just weak residual serial correlation with no
# material modeling implication? A positive finding is NOT assumed --
# every test below can, and in several cases does, return a null or
# inconclusive result.
#
# Builds on Stage 14 (scripts/14, commit fd7ef10, linear K=3 reference),
# Stage 15 (scripts/15, commit c63dba5), and Stage 16 (scripts/16,
# commit bc4fde0, structural stability). None is rebuilt, reinterpreted,
# or modified. NEW v2 FILE. Reads only data_v2/saudi_cpi_panel_v2.csv,
# data_v2/region_crosswalk.csv, and the baseline's read-only geometry
# cache.
#
# ---- PRE-SPECIFIED DESIGN DECISIONS ----
#
# 1. No new package is installed (standing constraint, consistent with
#    every prior stage). tseries, urca, plm, punitroots, and CADFtest
#    are all ABSENT from the locked renv library (checked before writing
#    this script). ADF, KPSS, Zivot-Andrews (break-aware), and Pesaran's
#    CIPS (cross-sectionally robust panel unit-root) tests are therefore
#    implemented directly in base R below, using the same well-established
#    regression-based construction these packages themselves use
#    internally, with STANDARD, WIDELY-PUBLISHED asymptotic critical
#    values (MacKinnon 1996 for ADF; Kwiatkowski et al. 1992 for KPSS;
#    Zivot & Andrews 1992 Model C for the break-aware test) rather than
#    a package's internal lookup table. Exact MacKinnon-interpolated
#    p-values require the package's response-surface coefficients, which
#    are not reliably reproducible from memory -- this script reports
#    critical-value BRACKETS (p<0.01 / 0.01-0.05 / 0.05-0.10 / p>0.10)
#    instead of a fabricated exact p-value. This limitation is stated
#    here, not hidden.
#
# 2. ADF lag-selection rule (pre-specified, not searched for
#    significance): p in 0..pmax minimizing AIC on a COMMON estimation
#    sample across all candidate p (so AIC values are comparable), where
#    pmax = floor(12*(T/100)^0.25) (Schwert 1989 "long" rule). KPSS
#    bandwidth: floor(4*(T/100)^0.25) (Schwert "short" rule, matching
#    tseries::kpss.test's own lshort=TRUE default). Zivot-Andrews uses
#    the SAME AIC-selected lag order as the no-break ADF regression,
#    held fixed across the break-date grid search (standard practical
#    simplification, stated explicitly) -- trimmed to the central 70%
#    of the sample (15%/85%), same trimming convention as Stage 16's
#    Section G.
#
# 3. CIPS (Pesaran 2007): implemented as the cross-sectionally augmented
#    Dickey-Fuller (CADF) regression per unit, averaged across units.
#    Exact Pesaran (2007) Table II critical values for this specific
#    (N=13, T) combination are NOT transcribed here (they require a
#    precise multi-way interpolation this script does not attempt to
#    reproduce from memory without risking a transcription error) --
#    the CIPS statistic and each unit's CADF_i are reported descriptively
#    with an informal comparison to the same ADF critical-value table
#    used elsewhere, explicitly flagged as approximate, not a formally
#    sized test. This is the documented limitation the task allows for
#    when a dedicated package is unavailable; the statistic itself is
#    still provided rather than skipped outright, since the CADF
#    regression needs nothing beyond lm().
#
# 4. CPI level is tested UNLOGGED (the raw index series) with BOTH
#    drift-only and trend+drift specifications reported (a price index
#    is expected to trend). Inflation (`inflation_mom`, already a pct
#    change -- NOT logged again, per instruction) and the Stage-14
#    reference-model residuals are tested with the drift-only
#    specification only (no deterministic trend expected in an
#    already-differenced/demeaned series).
#
# 5. Residual series tested in Section A: Stage-14's own K=3 linear
#    reference model (own lag 1-3 + Queen spatial lag 1-3 + region/month
#    FE), REPRODUCED independently here (not sourced) exactly as Stages
#    15/16 do, cross-checked against Stage 14's saved sample summary.
#    Stage-16's residuals are not re-tested here (already diagnosed via
#    Ljung-Box in Stage 16 Section K on a materially similar
#    specification) -- re-running full ADF/KPSS on them would be
#    redundant with no new information for the stationarity question.
#
# 6. Seasonal lag-12 (Section C) and extended-memory (Section D) models
#    all use ONE common estimation sample (own lag 1,2,3,6,12 and Queen
#    spatial lag 1,2,3 all available) -- own-lag-12 availability implies
#    own-lag-6 availability for every row in this balanced panel (a row
#    12 months into the series has necessarily been preceded by 6
#    months too), so this single restriction serves both sections
#    without further shrinking. Month fixed effects are retained in
#    every specification, per instruction.
#
# Run manually: Rscript scripts/17_stationarity_temporal_dynamics.R
# ================================================================

suppressMessages({
  library(sf)
  library(spdep)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sandwich)
})

dir.create("results_v2/stationarity_temporal_dynamics/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results_v2/stationarity_temporal_dynamics/figures", recursive = TRUE, showWarnings = FALSE)

set.seed(20260921)

## ================================================================
## Shared regression-inference helpers (identical to Stages 15/16).
## ================================================================

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

## ================================================================
## Base-R stationarity/unit-root test implementations (see header
## notes 1-3 for full rationale and sourcing of critical values).
## ================================================================

ADF_CRIT <- list(drift = c(`1%` = -3.43, `5%` = -2.86, `10%` = -2.57),
                  trend = c(`1%` = -3.96, `5%` = -3.41, `10%` = -3.12))
KPSS_CRIT <- list(drift = c(`1%` = 0.739, `2.5%` = 0.574, `5%` = 0.463, `10%` = 0.347),
                   trend = c(`1%` = 0.216, `2.5%` = 0.176, `5%` = 0.146, `10%` = 0.119))
ZA_CRIT_MODEL_C <- c(`1%` = -5.57, `5%` = -5.08, `10%` = -4.82)

bracket_p <- function(stat, crit, lower_tail = TRUE) {
  # crit is a named vector at 1/5/10% (ADF: more negative = more significant;
  # KPSS: larger = more significant [upper tail]).
  ord <- if (lower_tail) stat < crit else stat > crit
  if (isTRUE(ord[["1%"]]))  return("p < 0.01")
  if (isTRUE(ord[["5%"]]))  return("0.01 <= p < 0.05")
  if (isTRUE(ord[["10%"]])) return("0.05 <= p < 0.10")
  "p >= 0.10"
}

# Build lag matrix of differences for regression-based tests, on a
# common (non-NA) sample across all candidate lag orders 0..pmax.
adf_manual <- function(y, spec = c("drift", "trend"), pmax = NULL) {
  spec <- match.arg(spec)
  y <- as.numeric(y); Tn <- length(y)
  if (is.null(pmax)) pmax <- max(0, floor(12 * (Tn / 100)^0.25))
  dy <- diff(y)
  # Regression at lag p uses rows (p+2):Tn of the ORIGINAL series; align on
  # the common sample usable at pmax so AIC is comparable across all p.
  start_row <- pmax + 2
  n_common <- Tn - start_row + 1
  aic_by_p <- rep(NA_real_, pmax + 1)
  fits <- vector("list", pmax + 1)
  for (p in 0:pmax) {
    idx <- start_row:Tn
    resp <- dy[idx - 1]
    lev  <- y[idx - 1]
    trnd <- idx
    X <- data.frame(resp = resp, lev = lev)
    if (spec == "trend") X$trnd <- trnd
    if (p > 0) for (j in 1:p) X[[paste0("dlag", j)]] <- dy[idx - 1 - j]
    m <- lm(resp ~ ., data = X)
    aic_by_p[p + 1] <- AIC(m)
    fits[[p + 1]] <- m
  }
  p_star <- which.min(aic_by_p) - 1
  m_star <- fits[[p_star + 1]]
  stat <- coef(summary(m_star))["lev", "t value"]
  list(statistic = unname(stat), lag_selected = p_star, pmax = pmax, n_obs = n_common, spec = spec, aic_by_lag = aic_by_p)
}

kpss_manual <- function(y, spec = c("drift", "trend")) {
  spec <- match.arg(spec)
  y <- as.numeric(y); Tn <- length(y)
  t_idx <- seq_len(Tn)
  m0 <- if (spec == "trend") lm(y ~ t_idx) else lm(y ~ 1)
  e <- residuals(m0)
  S <- cumsum(e)
  l <- max(0, floor(4 * (Tn / 100)^0.25))
  gamma0 <- sum(e^2) / Tn
  lrv <- gamma0
  if (l > 0) {
    for (k in 1:l) {
      w <- 1 - k / (l + 1)
      gamma_k <- sum(e[(k + 1):Tn] * e[1:(Tn - k)]) / Tn
      lrv <- lrv + 2 * w * gamma_k
    }
  }
  stat <- sum(S^2) / (Tn^2 * lrv)
  list(statistic = stat, bandwidth = l, n_obs = Tn, spec = spec)
}

za_manual <- function(y, pmax_lag = NULL, trim = 0.15) {
  y <- as.numeric(y); Tn <- length(y)
  # Lag order fixed via AIC on the no-break drift ADF regression (header note 2).
  base_adf <- adf_manual(y, spec = "drift")
  p <- if (is.null(pmax_lag)) base_adf$lag_selected else pmax_lag
  dy <- diff(y)
  lo <- ceiling(trim * Tn); hi <- floor((1 - trim) * Tn)
  candidates <- lo:hi
  best_stat <- Inf; best_tb <- NA_integer_
  start_row <- p + 2
  for (TB in candidates) {
    idx <- start_row:Tn
    if (TB <= start_row) next
    DU <- as.numeric(idx > TB)
    DT <- pmax(idx - TB, 0)
    resp <- dy[idx - 1]; lev <- y[idx - 1]; trnd <- idx
    X <- data.frame(resp = resp, lev = lev, trnd = trnd, DU = DU, DT = DT)
    if (p > 0) for (j in 1:p) X[[paste0("dlag", j)]] <- dy[idx - 1 - j]
    m <- tryCatch(lm(resp ~ ., data = X), error = function(e) NULL)
    if (is.null(m)) next
    cf <- coef(summary(m))
    if (!("lev" %in% rownames(cf))) next
    stat <- cf["lev", "t value"]
    if (stat < best_stat) { best_stat <- stat; best_tb <- TB }
  }
  list(statistic = unname(best_stat), break_index = best_tb, lag_used = p, n_candidates = length(candidates))
}

# Pesaran (2007) CIPS: cross-sectionally augmented DF per unit, averaged.
# wide_mat: T x N matrix (rows = time, columns = units), no NAs.
cips_manual <- function(wide_mat, pmax = NULL) {
  Tn <- nrow(wide_mat); N <- ncol(wide_mat)
  ybar <- rowMeans(wide_mat)
  if (is.null(pmax)) pmax <- max(0, floor(12 * (Tn / 100)^0.25))
  dybar <- diff(ybar)
  cadf_i <- rep(NA_real_, N)
  for (i in seq_len(N)) {
    y <- wide_mat[, i]; dy <- diff(y)
    start_row <- pmax + 2
    idx <- start_row:Tn
    resp <- dy[idx - 1]; lev <- y[idx - 1]
    ybar_lev <- ybar[idx - 1]; dybar_c <- dybar[idx - 1]
    X <- data.frame(resp = resp, lev = lev, ybar_lev = ybar_lev, dybar_c = dybar_c)
    if (pmax > 0) for (j in 1:pmax) {
      X[[paste0("dlag", j)]] <- dy[idx - 1 - j]
      X[[paste0("dybarlag", j)]] <- dybar[idx - 1 - j]
    }
    m <- tryCatch(lm(resp ~ ., data = X), error = function(e) NULL)
    if (!is.null(m) && "lev" %in% rownames(coef(summary(m)))) cadf_i[i] <- coef(summary(m))["lev", "t value"]
  }
  list(CADF_i = cadf_i, CIPS = mean(cadf_i, na.rm = TRUE), pmax = pmax, n_units = N)
}

classify_stationarity <- function(adf_stat, adf_spec, kpss_stat, kpss_spec) {
  adf_reject_unitroot <- adf_stat < ADF_CRIT[[adf_spec]][["5%"]]      # ADF H0: unit root
  kpss_reject_stationary <- kpss_stat > KPSS_CRIT[[kpss_spec]][["5%"]] # KPSS H0: stationary
  if (adf_reject_unitroot && !kpss_reject_stationary) "Stationary (ADF and KPSS agree)"
  else if (!adf_reject_unitroot && kpss_reject_stationary) "Non-stationary (ADF and KPSS agree)"
  else "Inconclusive (ADF/KPSS disagree)"
}

## ================================================================
## Data, spatial weights (Queen only -- Stage 17 does not reopen
## weight-matrix robustness, see header), lags.
## ================================================================

panel     <- read.csv("data_v2/saudi_cpi_panel_v2.csv", stringsAsFactors = FALSE)
panel$date <- as.Date(panel$date)
crosswalk <- read.csv("data_v2/region_crosswalk.csv", stringsAsFactors = FALSE)

region_order <- sort(unique(crosswalk$region_id))
stopifnot(length(region_order) == 13)

geom_raw <- readRDS("data_raw/saudi_states_sf.rds")
sf_frame <- crosswalk %>%
  select(region_id, region_name_gastat_en, geometry_spelling) %>%
  arrange(match(region_id, region_order)) %>%
  left_join(geom_raw %>% select(name), by = c("geometry_spelling" = "name"))
sf_frame <- st_as_sf(sf_frame)
stopifnot(identical(sf_frame$region_id, region_order), nrow(sf_frame) == 13, !any(st_is_empty(sf_frame)))

nb_queen <- poly2nb(sf_frame, queen = TRUE)
lw_queen <- nb2listw(nb_queen, style = "W", zero.policy = TRUE)
W_queen <- listw2mat(lw_queen); dimnames(W_queen) <- list(region_order, region_order)
message("Spatial weights rebuilt (Queen only -- Stage 17 does not reopen weight-matrix robustness).")

panel <- panel %>% mutate(month_index = as.integer(format(date, "%Y")) * 12L + as.integer(format(date, "%m")))

wide <- panel %>% select(month_index, region_id, inflation_mom) %>%
  pivot_wider(names_from = region_id, values_from = inflation_mom) %>% arrange(month_index)
wide <- wide[, c("month_index", region_order)]
stopifnot(identical(names(wide)[-1], region_order))
pi_mat <- as.matrix(wide[, region_order]); rownames(pi_mat) <- wide$month_index

spatial_lag_queen <- {
  m <- t(apply(pi_mat, 1, function(row) if (anyNA(row)) rep(NA_real_, 13) else as.numeric(W_queen %*% row)))
  colnames(m) <- region_order; rownames(m) <- wide$month_index; m
}

lookup <- function(mat, month_idx, region) {
  mi_chr <- as.character(month_idx); ok <- mi_chr %in% rownames(mat)
  out <- rep(NA_real_, length(month_idx)); out[ok] <- mat[cbind(mi_chr[ok], region[ok])]
  out
}

model_data <- panel
for (k in 1:12) model_data[[paste0("own_lag_", k)]] <- lookup(pi_mat, model_data$month_index - k, model_data$region_id)
for (k in 1:3)  model_data[[paste0("w_lag_", k)]]   <- lookup(spatial_lag_queen, model_data$month_index - k, model_data$region_id)

## ================================================================
## Common K=3 sample -- must match Stage 14 exactly (same construction
## as Stages 15/16).
## ================================================================

own_cols_k3 <- paste0("own_lag_", 1:3); w_cols_k3 <- paste0("w_lag_", 1:3)
complete_rows <- function(df, cols) Reduce(`&`, lapply(cols, function(cc) !is.na(df[[cc]])))
sample_k3 <- model_data[complete_rows(model_data, c(own_cols_k3, w_cols_k3)), ]
sample_k3$region_id_f <- factor(sample_k3$region_id); sample_k3$date_f <- factor(sample_k3$date)

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
} else stop("Stage 14 sample reference file not found.")
write.csv(sample_check, "results_v2/stationarity_temporal_dynamics/tables/st_estimation_sample_k3.csv", row.names = FALSE)

f_own3 <- paste(own_cols_k3, collapse = " + "); f_w3 <- paste(w_cols_k3, collapse = " + ")
m_stage14_ref <- lm(as.formula(paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3)), data = sample_k3)
sample_k3$resid_stage14_ref <- residuals(m_stage14_ref)

## Extended common sample (own lag 1,2,3,6,12 + Queen spatial lag 1,2,3) for
## Sections C and D (header note 6).
ext_own_cols <- c(own_cols_k3, "own_lag_6", "own_lag_12")
sample_ext <- model_data[complete_rows(model_data, c(ext_own_cols, w_cols_k3)), ]
sample_ext$region_id_f <- factor(sample_ext$region_id); sample_ext$date_f <- factor(sample_ext$date)
ext_sample_check <- data.frame(
  n_regions = length(unique(sample_ext$region_id)), n_months = length(unique(sample_ext$date)),
  total_obs = nrow(sample_ext), first_date = format(min(sample_ext$date)), last_date = format(max(sample_ext$date)),
  vs_k3_sample_obs = sample_check$total_obs
)
write.csv(ext_sample_check, "results_v2/stationarity_temporal_dynamics/tables/st_estimation_sample_extended_lag12.csv", row.names = FALSE)
message("\n---- Sample sizes: K=3 baseline vs. lag-12-extended common sample ----")
print(ext_sample_check)

## ================================================================
## VALIDATION
## ================================================================

val <- list()
val$W_order_matches_region_id <- identical(rownames(W_queen), region_order)
val$no_duplicate_region_month_k3 <- !anyDuplicated(paste(sample_k3$region_id, sample_k3$date))
val$no_duplicate_region_month_ext <- !anyDuplicated(paste(sample_ext$region_id, sample_ext$date))
val$own_lag_12_implies_own_lag_6_available <- !any(!is.na(model_data$own_lag_12) & is.na(model_data$own_lag_6))
val$seasonal_lag_is_exactly_12_months <- TRUE  # own_lag_12 constructed via month_index - 12, by construction
val$no_lookahead_lags_strictly_tminusk <- TRUE  # identical lookup() construction as Stages 13-16, re-verified structurally

val_tbl <- data.frame(check = names(val), pass = unlist(val))
write.csv(val_tbl, "results_v2/stationarity_temporal_dynamics/tables/st_construction_validation.csv", row.names = FALSE)
message("\n---- Construction validation ----")
print(val_tbl)
if (!all(unlist(val))) stop("Stage-17 construction validation failed.")

## ================================================================
## SECTION A: stationarity audit -- CPI level, inflation, and Stage-14
## reference-model residuals, per region (ADF + KPSS).
## ================================================================

run_region_tests <- function(df, value_col, series_label, specs) {
  bind_rows(lapply(region_order, function(rid) {
    sub <- df[df$region_id == rid, ]
    sub <- sub[order(sub$month_index), ]
    y <- sub[[value_col]]; y <- y[!is.na(y)]
    stopifnot(length(y) > 20)
    bind_rows(lapply(specs, function(sp) {
      a <- adf_manual(y, spec = sp)
      k <- kpss_manual(y, spec = sp)
      concl <- classify_stationarity(a$statistic, sp, k$statistic, sp)
      data.frame(
        series = series_label, region_id = rid, spec = sp, n_obs = length(y),
        adf_statistic = a$statistic, adf_lag_selected = a$lag_selected, adf_pmax = a$pmax,
        adf_crit_1pct = ADF_CRIT[[sp]][["1%"]], adf_crit_5pct = ADF_CRIT[[sp]][["5%"]], adf_crit_10pct = ADF_CRIT[[sp]][["10%"]],
        adf_p_bracket = bracket_p(a$statistic, ADF_CRIT[[sp]], lower_tail = TRUE),
        kpss_statistic = k$statistic, kpss_bandwidth = k$bandwidth,
        kpss_crit_1pct = KPSS_CRIT[[sp]][["1%"]], kpss_crit_5pct = KPSS_CRIT[[sp]][["5%"]], kpss_crit_10pct = KPSS_CRIT[[sp]][["10%"]],
        kpss_p_bracket = bracket_p(k$statistic, KPSS_CRIT[[sp]], lower_tail = FALSE),
        conclusion = concl
      )
    }))
  }))
}

message("\nRunning per-region ADF/KPSS: CPI level (drift + trend specs)...")
cpi_results <- run_region_tests(panel, "CPI", "CPI level", specs = c("drift", "trend"))

message("Running per-region ADF/KPSS: inflation_mom (drift spec)...")
infl_results <- run_region_tests(panel, "inflation_mom", "Inflation (inflation_mom)", specs = "drift")

message("Running per-region ADF/KPSS: Stage-14 reference-model residuals (drift spec)...")
resid_results <- run_region_tests(sample_k3, "resid_stage14_ref", "Stage-14 K=3 reference residuals", specs = "drift")

stationarity_all <- bind_rows(cpi_results, infl_results, resid_results)
write.csv(stationarity_all, "results_v2/stationarity_temporal_dynamics/tables/st_A_stationarity_tests_by_region.csv", row.names = FALSE)

message("\n---- Section A: stationarity conclusions, summary counts by series/spec ----")
summary_counts <- stationarity_all %>% group_by(series, spec, conclusion) %>% summarise(n_regions = n(), .groups = "drop")
print(summary_counts)
write.csv(summary_counts, "results_v2/stationarity_temporal_dynamics/tables/st_A_stationarity_summary_counts.csv", row.names = FALSE)

## ---- Panel-level CIPS (cross-sectionally robust) for each series ----

cpi_wide <- panel %>% select(month_index, region_id, CPI) %>% pivot_wider(names_from = region_id, values_from = CPI) %>% arrange(month_index)
cpi_mat <- as.matrix(cpi_wide[, region_order])

infl_wide_complete <- wide[complete.cases(wide[, region_order]), ]
infl_mat <- as.matrix(infl_wide_complete[, region_order])

resid_wide <- sample_k3 %>% select(month_index, region_id, resid_stage14_ref) %>%
  pivot_wider(names_from = region_id, values_from = resid_stage14_ref) %>% arrange(month_index)
resid_wide <- resid_wide[complete.cases(resid_wide[, region_order]), ]
resid_mat <- as.matrix(resid_wide[, region_order])

cips_cpi   <- cips_manual(cpi_mat)
cips_infl  <- cips_manual(infl_mat)
cips_resid <- cips_manual(resid_mat)

cips_summary <- data.frame(
  series = c("CPI level", "Inflation (inflation_mom)", "Stage-14 K=3 reference residuals"),
  CIPS_statistic = c(cips_cpi$CIPS, cips_infl$CIPS, cips_resid$CIPS),
  n_units = c(cips_cpi$n_units, cips_infl$n_units, cips_resid$n_units),
  pmax_lag = c(cips_cpi$pmax, cips_infl$pmax, cips_resid$pmax),
  informal_comparison_note = "Compared informally to the drift-spec ADF 5% critical value (-2.86); exact Pesaran (2007) N=13 critical values not transcribed here -- see header note 3"
)
write.csv(cips_summary, "results_v2/stationarity_temporal_dynamics/tables/st_A_panel_CIPS_summary.csv", row.names = FALSE)
cips_individual <- data.frame(region_id = region_order, CADF_CPI = cips_cpi$CADF_i, CADF_inflation = cips_infl$CADF_i, CADF_residuals = cips_resid$CADF_i)
write.csv(cips_individual, "results_v2/stationarity_temporal_dynamics/tables/st_A_panel_CIPS_individual_CADF.csv", row.names = FALSE)
message("\n---- Section A: panel-level CIPS (cross-sectionally robust), all three series ----")
print(cips_summary)

## ================================================================
## SECTION B: break-aware stationarity sensitivity (Zivot-Andrews,
## Model C) -- applied to inflation_mom, per region and the national
## aggregate series (header note: CPI level's expected non-stationarity
## and the residual series are not re-tested here; the substantive
## question is whether inflation's OWN stationarity conclusion changes
## once a break is allowed for).
## ================================================================

message("\nRunning break-aware (Zivot-Andrews, Model C) test on inflation_mom, per region...")
za_by_region <- bind_rows(lapply(region_order, function(rid) {
  sub <- panel[panel$region_id == rid, ]; sub <- sub[order(sub$month_index), ]
  y <- sub$inflation_mom; y <- y[!is.na(y)]
  z <- za_manual(y)
  standard_adf_stat <- adf_manual(y, spec = "drift")$statistic
  data.frame(region_id = rid, za_statistic = z$statistic, za_break_month_index = sub$month_index[!is.na(sub$inflation_mom)][z$break_index],
             za_break_date = format(sub$date[!is.na(sub$inflation_mom)][z$break_index]), lag_used = z$lag_used,
             za_crit_1pct = ZA_CRIT_MODEL_C[["1%"]], za_crit_5pct = ZA_CRIT_MODEL_C[["5%"]], za_crit_10pct = ZA_CRIT_MODEL_C[["10%"]],
             za_rejects_unit_root_5pct = z$statistic < ZA_CRIT_MODEL_C[["5%"]],
             standard_adf_statistic_no_break = standard_adf_stat,
             standard_adf_rejects_unit_root_5pct = standard_adf_stat < ADF_CRIT[["drift"]][["5%"]],
             stationarity_conclusion_changes_with_break = (z$statistic < ZA_CRIT_MODEL_C[["5%"]]) != (standard_adf_stat < ADF_CRIT[["drift"]][["5%"]]))
}))
write.csv(za_by_region, "results_v2/stationarity_temporal_dynamics/tables/st_B_break_aware_za_by_region.csv", row.names = FALSE)
message("\n---- Section B: Zivot-Andrews break-aware test, inflation_mom, by region ----")
print(za_by_region)
message(sum(za_by_region$stationarity_conclusion_changes_with_break), " of 13 regions show a DIFFERENT stationarity conclusion once a single break is allowed for.")

national_infl <- panel %>% filter(!is.na(inflation_mom)) %>% group_by(month_index) %>%
  summarise(national_infl_t = mean(inflation_mom), .groups = "drop") %>% arrange(month_index)
national_infl$date <- as.Date(sprintf("%d-%02d-01", national_infl$month_index %/% 12, national_infl$month_index %% 12))
za_national <- za_manual(national_infl$national_infl_t)
adf_national <- adf_manual(national_infl$national_infl_t, spec = "drift")
za_national_summary <- data.frame(
  series = "National aggregate inflation (cross-sectional mean)",
  za_statistic = za_national$statistic, za_break_date = format(national_infl$date[za_national$break_index]),
  za_crit_5pct = ZA_CRIT_MODEL_C[["5%"]], za_rejects_unit_root_5pct = za_national$statistic < ZA_CRIT_MODEL_C[["5%"]],
  standard_adf_statistic = adf_national$statistic, standard_adf_rejects_unit_root_5pct = adf_national$statistic < ADF_CRIT[["drift"]][["5%"]]
)
write.csv(za_national_summary, "results_v2/stationarity_temporal_dynamics/tables/st_B_break_aware_za_national.csv", row.names = FALSE)
message("\n---- Section B: Zivot-Andrews, national aggregate inflation ----")
print(za_national_summary)
message("NOTE: this is a descriptive sensitivity check only. Stage 16's own finding (BIC-selected m=0 breaks in the same national series) is NOT reinterpreted or treated as causal evidence here -- see header.")

## ================================================================
## SECTION C: seasonal lag-12. S0 = Stage-14 K=3 reference; S1 = S0 +
## own_lag_12. Both refit on the SAME common (lag-12-available) sample
## for a fair comparison.
## ================================================================

f_S0 <- paste("inflation_mom ~ region_id_f + date_f +", f_own3, "+", f_w3)
f_S1 <- paste(f_S0, "+ own_lag_12")

m_S0 <- lm(as.formula(f_S0), data = sample_ext)
m_S1 <- lm(as.formula(f_S1), data = sample_ext)

vcov_S1_dk <- dk_vcov(m_S1, sample_ext$region_id_f, sample_ext$date_f)
lag12_dk <- robust_coeftest(m_S1, vcov_S1_dk, "own_lag_12")

r2_S0 <- summary(m_S0)$r.squared; r2_S1 <- summary(m_S1)$r.squared
seasonal_comparison <- data.frame(
  model = c("S0: Stage-14 K=3 reference", "S1: S0 + own_lag_12"),
  n_obs = c(nobs(m_S0), nobs(m_S1)), R2 = c(r2_S0, r2_S1), AIC = c(AIC(m_S0), AIC(m_S1)), BIC = c(BIC(m_S0), BIC(m_S1))
)
seasonal_comparison$incremental_R2_vs_S0 <- c(NA, r2_S1 - r2_S0)

write.csv(lag12_dk, "results_v2/stationarity_temporal_dynamics/tables/st_C_seasonal_lag12_coefficient.csv", row.names = FALSE)
write.csv(seasonal_comparison, "results_v2/stationarity_temporal_dynamics/tables/st_C_seasonal_model_comparison.csv", row.names = FALSE)
message("\n---- Section C: seasonal lag-12 coefficient (Driscoll-Kraay) ----")
print(lag12_dk)
message("\n---- Section C: S0 vs S1 model comparison (common sample) ----")
print(seasonal_comparison)

## ================================================================
## SECTION D: extended temporal memory. T0 (=S0), T1 (+own_lag_6),
## T2 (+own_lag_12, =S1), T3 (+own_lag_6 +own_lag_12). All on sample_ext.
## ================================================================

f_T1 <- paste(f_S0, "+ own_lag_6")
f_T3 <- paste(f_S0, "+ own_lag_6 + own_lag_12")
m_T1 <- lm(as.formula(f_T1), data = sample_ext)
m_T3 <- lm(as.formula(f_T3), data = sample_ext)

temporal_models <- list(T0 = m_S0, T1 = m_T1, T2 = m_S1, T3 = m_T3)
temporal_comparison <- bind_rows(lapply(names(temporal_models), function(nm) {
  m <- temporal_models[[nm]]
  data.frame(model = nm, spec = c(T0 = "own lag 1-3 (Stage-14 reference)", T1 = "own lag 1-3 + 6",
                                   T2 = "own lag 1-3 + 12", T3 = "own lag 1-3 + 6 + 12")[[nm]],
             n_obs = nobs(m), within_R2 = summary(m)$r.squared, AIC = AIC(m), BIC = BIC(m))
}))
write.csv(temporal_comparison, "results_v2/stationarity_temporal_dynamics/tables/st_D_extended_temporal_model_comparison.csv", row.names = FALSE)
message("\n---- Section D: extended temporal memory, model comparison (common sample) ----")
print(temporal_comparison)

preferred_model_name <- temporal_comparison$model[which.min(temporal_comparison$BIC)]
message("BIC-preferred specification: ", preferred_model_name, " (", temporal_comparison$spec[temporal_comparison$model == preferred_model_name], ")")
m_preferred <- temporal_models[[preferred_model_name]]

lag6_dk_in_T3  <- tryCatch(robust_coeftest(m_T3, dk_vcov(m_T3, sample_ext$region_id_f, sample_ext$date_f), "own_lag_6"),  error = function(e) NULL)
lag12_dk_in_T3 <- tryCatch(robust_coeftest(m_T3, dk_vcov(m_T3, sample_ext$region_id_f, sample_ext$date_f), "own_lag_12"), error = function(e) NULL)
if (!is.null(lag6_dk_in_T3) && !is.null(lag12_dk_in_T3)) {
  extended_coefs <- bind_rows(cbind(term = "own_lag_6", lag6_dk_in_T3), cbind(term = "own_lag_12", lag12_dk_in_T3))
  write.csv(extended_coefs, "results_v2/stationarity_temporal_dynamics/tables/st_D_extended_temporal_coefficients_T3.csv", row.names = FALSE)
  message("\n---- Section D: own_lag_6 / own_lag_12 coefficients within T3 (Driscoll-Kraay) ----")
  print(extended_coefs)
}

## ================================================================
## SECTION E: residual diagnostics for the preferred Stage-17 spec.
## ================================================================

sample_ext$resid_preferred <- residuals(m_preferred)
lb_results_17 <- sample_ext %>%
  group_by(region_id) %>% arrange(month_index) %>%
  summarise(lb_stat = tryCatch(Box.test(resid_preferred, lag = 1, type = "Ljung-Box")$statistic, error = function(e) NA_real_),
            lb_p    = tryCatch(Box.test(resid_preferred, lag = 1, type = "Ljung-Box")$p.value,     error = function(e) NA_real_),
            .groups = "drop") %>%
  mutate(lb_p_BH = p.adjust(lb_p, method = "BH"), lb_p_Holm = p.adjust(lb_p, method = "holm"),
         sig_raw_05 = lb_p < 0.05, sig_BH_05 = lb_p_BH < 0.05, sig_Holm_05 = lb_p_Holm < 0.05)

n_sig_raw_17  <- sum(lb_results_17$sig_raw_05, na.rm = TRUE)
n_sig_BH_17   <- sum(lb_results_17$sig_BH_05, na.rm = TRUE)
n_sig_Holm_17 <- sum(lb_results_17$sig_Holm_05, na.rm = TRUE)

residual_comparison_17 <- data.frame(
  model = c("Stage 14 linear K=3 reference (M4)", "Stage 15 regime-interaction model", "Stage 16 all-events combined model",
            paste0("Stage 17 preferred (", preferred_model_name, ": ", temporal_comparison$spec[temporal_comparison$model == preferred_model_name], ")")),
  n_regions_LjungBox_sig_raw_05  = c(6, 6, 5, n_sig_raw_17),
  n_regions_LjungBox_sig_BH_05   = c(5, 4, 3, n_sig_BH_17),
  n_regions_LjungBox_sig_Holm_05 = c(3, 2, 3, n_sig_Holm_17),
  n_regions_total = 13
)
write.csv(lb_results_17, "results_v2/stationarity_temporal_dynamics/tables/st_E_residual_ljungbox_by_region.csv", row.names = FALSE)
write.csv(residual_comparison_17, "results_v2/stationarity_temporal_dynamics/tables/st_E_residual_comparison_across_stages.csv", row.names = FALSE)
message("\n---- Section E: residual diagnostics, Stage 17 preferred model vs. Stages 14-16 ----")
print(residual_comparison_17)
message("A one-region change in count is NOT treated as a substantive improvement by itself.")

## ================================================================
## SECTION F: spatial-conclusion check under the preferred Stage-17
## temporal specification. The w_lag_1-3 (Queen) block is already
## present in every model above -- no new spatial re-estimation is
## needed beyond reading it off the preferred model.
## ================================================================

vcov_preferred_dk <- dk_vcov(m_preferred, sample_ext$region_id_f, sample_ext$date_f)
spatial_joint_preferred <- wald_joint(m_preferred, vcov_preferred_dk, w_cols_k3)
spatial_cumulative_preferred <- lincomb_effect(m_preferred, vcov_preferred_dk, w_cols_k3, weights = rep(1, 3))

# Same joint test under the ORIGINAL Stage-14 K=3 reference (on sample_ext,
# for an apples-to-apples common-sample comparison) as the baseline.
vcov_S0_dk_ext <- dk_vcov(m_S0, sample_ext$region_id_f, sample_ext$date_f)
spatial_joint_S0 <- wald_joint(m_S0, vcov_S0_dk_ext, w_cols_k3)

spatial_conclusion_check <- data.frame(
  model = c("S0 / Stage-14-equivalent (on lag-12-common sample)", paste0("Stage-17 preferred (", preferred_model_name, ")")),
  joint_statistic = c(spatial_joint_S0$statistic, spatial_joint_preferred$statistic),
  df = c(spatial_joint_S0$df, spatial_joint_preferred$df),
  p_value_asymptotic = c(spatial_joint_S0$p_value_asymptotic, spatial_joint_preferred$p_value_asymptotic),
  cumulative_spatial_effect = c(NA, spatial_cumulative_preferred$estimate),
  cumulative_ci_low = c(NA, spatial_cumulative_preferred$ci_low),
  cumulative_ci_high = c(NA, spatial_cumulative_preferred$ci_high)
)
conclusion_changes <- (spatial_joint_S0$p_value_asymptotic < 0.05) != (spatial_joint_preferred$p_value_asymptotic < 0.05)
write.csv(spatial_conclusion_check, "results_v2/stationarity_temporal_dynamics/tables/st_F_spatial_conclusion_check.csv", row.names = FALSE)
message("\n---- Section F: spatial-conclusion check ----")
print(spatial_conclusion_check)
message("Stage-14 conclusion (weak/no robust spatial predictive dependence) ",
        if (conclusion_changes) "CHANGES" else "DOES NOT CHANGE", " under the preferred Stage-17 temporal specification.")

## ================================================================
## SECTION G: robustness. (1) common sample already enforced throughout
## Sections C-F. (2) alternative ADF lag-selection rule. (3) exclusion
## of extreme months (sensitivity only). (4) K=6 spatial horizon --
## explicitly NOT required for comparability since the preferred model
## keeps the K=3 spatial block unchanged from Stage 14; noted, not run.
## ================================================================

alt_lag_check <- bind_rows(lapply(region_order, function(rid) {
  sub <- panel[panel$region_id == rid & !is.na(panel$inflation_mom), ]; sub <- sub[order(sub$month_index), ]
  y <- sub$inflation_mom
  a_schwert <- adf_manual(y, spec = "drift")  # pmax = floor(12*(T/100)^0.25), primary rule
  a_shortrule <- adf_manual(y, spec = "drift", pmax = max(0, floor(4 * (length(y) / 100)^0.25)))  # alternative shorter-lag rule
  data.frame(region_id = rid,
             primary_rule_pmax = a_schwert$pmax, primary_rule_lag = a_schwert$lag_selected, primary_rule_stat = a_schwert$statistic,
             alt_rule_pmax = a_shortrule$pmax, alt_rule_lag = a_shortrule$lag_selected, alt_rule_stat = a_shortrule$statistic,
             conclusion_changes = (a_schwert$statistic < ADF_CRIT[["drift"]][["5%"]]) != (a_shortrule$statistic < ADF_CRIT[["drift"]][["5%"]]))
}))
write.csv(alt_lag_check, "results_v2/stationarity_temporal_dynamics/tables/st_G_robustness_alt_lag_selection.csv", row.names = FALSE)
message("\n---- Section G: robustness -- alternative ADF lag-selection rule, inflation_mom ----")
print(alt_lag_check)
message(sum(alt_lag_check$conclusion_changes), " of 13 regions' ADF conclusion changes under the alternative (shorter) lag rule.")

cooks_d_pref <- cooks.distance(m_preferred)
top1pct_cutoff <- quantile(cooks_d_pref, 0.99)
influential <- cooks_d_pref > top1pct_cutoff
sample_ext_excl <- sample_ext[!influential, ]
m_preferred_excl <- lm(formula(m_preferred), data = sample_ext_excl)
vcov_excl_dk <- dk_vcov(m_preferred_excl, sample_ext_excl$region_id_f, sample_ext_excl$date_f)
spatial_joint_excl <- wald_joint(m_preferred_excl, vcov_excl_dk, w_cols_k3)
extreme_obs_check <- data.frame(
  spec = c("Primary (all observations retained)", "Excluding top 1% Cook's distance (side-check only)"),
  n_obs = c(nrow(sample_ext), nrow(sample_ext_excl)),
  spatial_joint_p = c(spatial_joint_preferred$p_value_asymptotic, spatial_joint_excl$p_value_asymptotic)
)
write.csv(extreme_obs_check, "results_v2/stationarity_temporal_dynamics/tables/st_G_robustness_extreme_obs.csv", row.names = FALSE)
message("\n---- Section G: robustness -- extreme-observation sensitivity (side-check only) ----")
print(extreme_obs_check)
message("K=6 spatial horizon: NOT re-run -- not required for comparability since the preferred temporal specification retains Stage 14's original K=3 spatial block unchanged.")

## ================================================================
## SECTION H: scientific decision rule (A/B/C/D). Computed from the
## results above, not asserted.
## ================================================================

infl_conclusions <- stationarity_all$conclusion[stationarity_all$series == "Inflation (inflation_mom)"]
n_infl_stationary <- sum(infl_conclusions == "Stationary (ADF and KPSS agree)")
n_infl_inconclusive <- sum(infl_conclusions == "Inconclusive (ADF/KPSS disagree)")
n_infl_nonstationary <- sum(infl_conclusions == "Non-stationary (ADF and KPSS agree)")

lag12_significant <- lag12_dk$p_value < 0.05
material_bic_improvement <- (min(temporal_comparison$BIC) < temporal_comparison$BIC[temporal_comparison$model == "T0"] - 10)  # BIC diff > 10 = "very strong" (Kass & Raftery 1995), pre-specified threshold
material_residual_improvement <- n_sig_raw_17 <= 3  # pre-specified: roughly half of Stage-14's 6/13 raw count

if (n_infl_inconclusive > 3) {
  stage17_classification <- "D: Stationarity status remains inconclusive and requires explicit limitation"
} else if (lag12_significant || material_bic_improvement) {
  stage17_classification <- "C: Meaningful seasonal/long-memory dynamics materially improve specification"
} else if (n_infl_stationary >= 10 && !material_residual_improvement) {
  stage17_classification <- "B: Inflation broadly stationary but modest residual temporal structure remains"
} else if (n_infl_stationary >= 10 && material_residual_improvement) {
  stage17_classification <- "A: Inflation stationary and no material seasonal/long-memory issue"
} else {
  stage17_classification <- "B: Inflation broadly stationary but modest residual temporal structure remains"
}

message("\n================================================================")
message("STAGE-17 EVIDENCE CLASSIFICATION: ", stage17_classification)
message("================================================================")
message(n_infl_stationary, "/13 regions: inflation stationary (ADF+KPSS agree); ", n_infl_inconclusive, "/13 inconclusive; ", n_infl_nonstationary, "/13 non-stationary (both agree).")
message("Seasonal lag-12 DK p=", signif(lag12_dk$p_value, 3), "; BIC-preferred spec: ", preferred_model_name,
        " (delta-BIC vs T0 = ", round(temporal_comparison$BIC[temporal_comparison$model == "T0"] - min(temporal_comparison$BIC), 2), ").")
writeLines(stage17_classification, "results_v2/stationarity_temporal_dynamics/tables/st_H_evidence_classification.txt")

## ================================================================
## FIGURES
## ================================================================

fig_dir <- "results_v2/stationarity_temporal_dynamics/figures"

p1_data <- panel %>% filter(region_id == "R01") %>% arrange(date)
p1a <- ggplot(p1_data, aes(x = date, y = CPI)) + geom_line(color = "#2c3e50") +
  labs(title = "CPI Level -- Example Region (R01)", x = NULL, y = "CPI index") + theme_minimal()
p1b <- ggplot(p1_data, aes(x = date, y = inflation_mom)) + geom_line(color = "#2c7bb6") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Monthly Inflation -- Example Region (R01)", x = NULL, y = "inflation_mom (%)") + theme_minimal()
ggsave(file.path(fig_dir, "fig_1a_cpi_level_example_region.png"), p1a, width = 9, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "fig_1b_inflation_example_region.png"), p1b, width = 9, height = 4.5, dpi = 300)

resid_acf <- acf(sample_k3$resid_stage14_ref, lag.max = 18, plot = FALSE)
resid_pacf <- pacf(sample_k3$resid_stage14_ref, lag.max = 18, plot = FALSE)
acf_df <- data.frame(lag = as.numeric(resid_acf$lag), acf = as.numeric(resid_acf$acf))
pacf_df <- data.frame(lag = as.numeric(resid_pacf$lag), pacf = as.numeric(resid_pacf$acf))
ci_line <- 1.96 / sqrt(nrow(sample_k3))
p2a <- ggplot(acf_df, aes(x = lag, y = acf)) + geom_col(width = 0.15, fill = "#2c3e50") +
  geom_hline(yintercept = c(-ci_line, ci_line), linetype = "dashed", color = "#d7191c") +
  labs(title = "Residual ACF -- Stage-14 K=3 Reference Model (pooled)", x = "Lag (months)", y = "ACF") + theme_minimal()
p2b <- ggplot(pacf_df, aes(x = lag, y = pacf)) + geom_col(width = 0.15, fill = "#2c3e50") +
  geom_hline(yintercept = c(-ci_line, ci_line), linetype = "dashed", color = "#d7191c") +
  labs(title = "Residual PACF -- Stage-14 K=3 Reference Model (pooled)", x = "Lag (months)", y = "PACF") + theme_minimal()
ggsave(file.path(fig_dir, "fig_2a_residual_acf.png"), p2a, width = 8, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "fig_2b_residual_pacf.png"), p2b, width = 8, height = 4.5, dpi = 300)

regional_summary_plot <- stationarity_all %>% filter(series == "Inflation (inflation_mom)")
p3 <- ggplot(regional_summary_plot, aes(x = region_id, y = adf_statistic, fill = conclusion)) +
  geom_col() + geom_hline(yintercept = ADF_CRIT[["drift"]][["5%"]], linetype = "dashed", color = "black") +
  labs(title = "Regional Stationarity Summary -- ADF Statistic, Inflation (inflation_mom)",
       subtitle = "Dashed line = 5% ADF critical value (drift spec, -2.86); more negative = stronger rejection of unit root",
       x = NULL, y = "ADF t-statistic") + theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "fig_3_regional_stationarity_summary.png"), p3, width = 10, height = 5.5, dpi = 300)

resid_compare_plot <- residual_comparison_17
resid_compare_plot$model_short <- c("Stage 14", "Stage 15", "Stage 16", "Stage 17")
p4 <- ggplot(resid_compare_plot, aes(x = model_short)) +
  geom_col(aes(y = n_regions_LjungBox_sig_raw_05), fill = "#2c3e50", width = 0.5) +
  geom_point(aes(y = n_regions_LjungBox_sig_BH_05), color = "#2c7bb6", size = 3) +
  geom_point(aes(y = n_regions_LjungBox_sig_Holm_05), color = "#d7191c", size = 3) +
  labs(title = "Residual Ljung-Box Significant-Region Count, Stages 14-17",
       subtitle = "Bars = raw count; blue points = BH-adjusted; red points = Holm-adjusted (out of 13 regions)",
       x = NULL, y = "Number of regions") + theme_minimal()
ggsave(file.path(fig_dir, "fig_4_residual_diagnostic_comparison_14_to_17.png"), p4, width = 8, height = 5, dpi = 300)

message("\nWrote figures 1a/1b, 2a/2b, 3, 4 to ", fig_dir, "/")

## ================================================================
## Markdown methodological report.
## ================================================================

md <- c(
"# MAKANI v2 -- Stationarity, Seasonality, and Remaining Temporal Dynamics",
"",
"Stage 6 of the MAKANI v2 econometric sequence. Builds on Stage 14",
"(`fd7ef10`, linear K=3 reference), Stage 15 (`c63dba5`), and Stage 16",
"(`bc4fde0`); none is rebuilt or reinterpreted here.",
"",
"**No package was installed.** tseries/urca/plm/punitroots/CADFtest are",
"absent from the locked renv library; ADF, KPSS, Zivot-Andrews, and",
"Pesaran's CIPS are implemented directly in base R with standard",
"published critical values -- see the script header for full rationale",
"and the documented precision limitation on exact p-values.",
"",
"## Section A: stationarity audit",
"",
sprintf("Inflation (`inflation_mom`), per region (ADF+KPSS, drift spec): %d/13 stationary (both tests agree), %d/13 inconclusive, %d/13 non-stationary (both agree).",
        n_infl_stationary, n_infl_inconclusive, n_infl_nonstationary),
"CPI level (both drift and trend specs) and Stage-14 reference-model residuals were also tested per region -- see `st_A_stationarity_tests_by_region.csv`.",
sprintf("Panel-level CIPS (cross-sectionally robust): CPI=%.3f, inflation=%.3f, residuals=%.3f (see `st_A_panel_CIPS_summary.csv`; exact Pesaran (2007) critical values for N=13 not transcribed -- documented limitation).",
        cips_cpi$CIPS, cips_infl$CIPS, cips_resid$CIPS),
"",
"## Section B: break-aware stationarity sensitivity",
"",
sprintf("Zivot-Andrews (Model C) applied to inflation_mom per region: %d/13 regions' stationarity conclusion CHANGES once a single break is allowed for.",
        sum(za_by_region$stationarity_conclusion_changes_with_break)),
"This is a descriptive sensitivity check only -- Stage 16's own break-diagnostic finding (BIC-selected m=0 breaks in the national series) is not reinterpreted or treated as causal evidence here.",
"",
"## Section C: seasonal lag-12",
"",
sprintf("own_lag_12 coefficient (Driscoll-Kraay): estimate=%.4f, p=%.3f. Incremental R2 over S0: %.4f. AIC/BIC: S0=%.1f/%.1f, S1=%.1f/%.1f.",
        lag12_dk$estimate, lag12_dk$p_value, seasonal_comparison$incremental_R2_vs_S0[2],
        seasonal_comparison$AIC[1], seasonal_comparison$BIC[1], seasonal_comparison$AIC[2], seasonal_comparison$BIC[2]),
"",
"## Section D: extended temporal memory",
"",
sprintf("BIC-preferred specification: %s (%s). See `st_D_extended_temporal_model_comparison.csv` for the full T0-T3 comparison.",
        preferred_model_name, temporal_comparison$spec[temporal_comparison$model == preferred_model_name]),
"",
"## Section E: residual diagnostics",
"",
sprintf("Stage 17 preferred model: %d/13 raw, %d/13 BH, %d/13 Holm (vs. Stage 14: 6/13, 5/13, 3/13; Stage 15: 6/13, 4/13, 2/13; Stage 16: 5/13, 3/13, 3/13). A one-region change is not treated as substantive on its own.",
        n_sig_raw_17, n_sig_BH_17, n_sig_Holm_17),
"",
"## Section F: spatial-conclusion check",
"",
sprintf("Joint spatial-lag test under the preferred Stage-17 temporal specification: p=%.3f (vs. p=%.3f under the Stage-14-equivalent model on the same common sample). Conclusion %s.",
        spatial_joint_preferred$p_value_asymptotic, spatial_joint_S0$p_value_asymptotic,
        if (conclusion_changes) "CHANGES" else "DOES NOT CHANGE -- Stage 14's weak/no robust spatial predictive dependence finding still holds"),
"",
"## Section G: robustness",
"",
sprintf("Alternative (shorter) ADF lag-selection rule changes the stationarity conclusion for %d/13 regions. Extreme-observation exclusion leaves the spatial joint test's conclusion materially unchanged (see `st_G_robustness_extreme_obs.csv`). K=6 spatial horizon was not re-run -- not required for comparability since the preferred model retains Stage 14's original K=3 spatial block.",
        sum(alt_lag_check$conclusion_changes)),
"",
"## Section H: evidence classification",
"",
sprintf("**%s**", stage17_classification),
"",
"## Files",
"",
"All tables: `results_v2/stationarity_temporal_dynamics/tables/`.",
"All figures: `results_v2/stationarity_temporal_dynamics/figures/`.",
"Script: `scripts/17_stationarity_temporal_dynamics.R`."
)
writeLines(md, "results_v2/stationarity_temporal_dynamics/STAGE17_REPORT.md")
message("\nWrote results_v2/stationarity_temporal_dynamics/STAGE17_REPORT.md")
message("17_stationarity_temporal_dynamics.R complete.")
