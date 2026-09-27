# ================================================================
# MAKANI v2 -- Stage 8 (pre-submission refinement): simulation-based
# power / minimum-detectable-effect (MDE) analysis for the LOCKED
# Stage-14 K=3 spatial-lag design.
#
# This is a DESIGN-SENSITIVITY analysis, not a new empirical finding.
# It answers: "given N=13, T of the actual K=3 estimation sample, the
# observed Queen weights, the locked dynamic specification, and the
# observed residual dependence, what magnitude of cumulative 3-month
# spatial predictive effect (Theta_3 = theta_1+theta_2+theta_3) could
# this design plausibly detect?" It does NOT re-estimate, reinterpret,
# or alter any Stage 14-18 empirical value -- those are read read-only
# from their committed CSVs for calibration and as a validation
# reference point, never modified.
#
# ---- PRE-SPECIFIED DESIGN (stated before any simulation is run) ----
#
# 1. Estimand: cumulative Theta_3, at a pre-specified symmetric grid
#    {0, +/-0.025, +/-0.05, +/-0.075, +/-0.10, +/-0.15, +/-0.20} (13
#    points), chosen before any simulation and NOT adjusted afterward
#    regardless of where power turns out to fall on this grid -- if
#    MDE_80 exceeds the largest grid point, that is reported as such,
#    not used as a reason to re-run a wider grid.
#
# 2. Lag-shape: PRIMARY = equal allocation (theta_k = Theta_3/3 for
#    k=1,2,3). SECONDARY (sensitivity only) = a pre-specified
#    front-loaded profile, weights (0.6, 0.3, 0.1) applied to Theta_3,
#    i.e. most of the cumulative effect realized at lag 1 and declining
#    over lags 2-3 -- a plausible alternative to "equal", not tuned to
#    maximize detected power.
#
# 3. Data-generating process: a FIXED-REGRESSOR Monte Carlo design.
#    own_lag_1-3, w_lag_1-3 (Queen), region_id_f, and date_f are held at
#    their OBSERVED values from the actual K=3 estimation sample
#    (reproduced independently here, identical construction to Stages
#    14-18, cross-checked against Stage 14's saved sample). The
#    simulated outcome is
#        y*_it = fitted0_it + theta_1*w_lag_1_it + theta_2*w_lag_2_it
#                + theta_3*w_lag_3_it + e*_it
#    where fitted0 is the fitted value from the RESTRICTED model (own
#    lags + region/month FE, no spatial lag -- i.e. the "no true spatial
#    effect" baseline), and e* is a residual draw from a MOVING BLOCK
#    BOOTSTRAP (block length L=6 months, primary) of that restricted
#    model's OBSERVED residuals, resampled as full 13-region row-blocks
#    (preserving observed cross-sectional dependence within a block)
#    concatenated across time (preserving observed short-range temporal
#    dependence up to the block length). L=6 is close to the standard
#    T^(1/3) ~= 5.4 rule-of-thumb block-bootstrap length for T=154.
#    L=12 (one calendar year) is run as a documented sensitivity check
#    on a reduced grid. This is a simplification relative to a fully
#    recursive dynamic re-simulation of the panel (which would also
#    require regenerating w_lag_k from a newly simulated y each period);
#    it is stated here explicitly, not hidden, and is a standard,
#    well-established Monte Carlo power-analysis design (fixed
#    regressors, resampled residuals) for exactly this purpose.
#
# 4. Inference rule: PRIMARY = the same Driscoll-Kraay joint Wald test
#    (vcovPL, Bartlett kernel, NW1987 bandwidth) used as primary
#    inference throughout Stages 14-18, applied to every simulated
#    replication (computationally tractable, unlike a full B=9999 wild-
#    cluster bootstrap nested inside thousands of replications, which
#    header note 5 below explains is not attempted for the full grid).
#    VALIDATED against a smaller nested wild-cluster-bootstrap (B=199)
#    benchmark at 3 selected effect sizes (Theta_3 = 0, 0.10, 0.20),
#    100 replications each -- reported as an approximation check, not
#    presented as "exact wild-bootstrap power" for the full grid.
#
# 5. Replications: 2,000 per point for the primary (equal-shape, L=6)
#    grid (13 points = 26,000 model fits); 1,000 per point for the
#    front-loaded-shape sensitivity (13 points = 13,000 fits); 500 per
#    point for the L=12 block-length sensitivity (5-point reduced grid
#    = 2,500 fits); 100 outer x B=199 inner for the WCB validation (3
#    points). A single fixed seed (20260921, the seed used throughout
#    MAKANI) is set ONCE at the start of the script; no per-replication
#    reseeding, so the full sequence is reproducible end to end.
#
# Run manually: Rscript scripts/20_power_mde_presubmission.R
# ================================================================

suppressMessages({
  library(sf)
  library(spdep)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sandwich)
})

dir.create("results_v2/presubmission_refinement/power_mde", recursive = TRUE, showWarnings = FALSE)

GLOBAL_SEED <- 20260921
set.seed(GLOBAL_SEED)

## ================================================================
## Shared helpers (identical to Stages 15-18).
## ================================================================

dk_vcov <- function(model, cluster, order_by) {
  vcovPL(model, cluster = cluster, order.by = order_by, kernel = "Bartlett", lag = "NW1987")
}
wald_joint <- function(model, vcov_mat, coefs) {
  b <- coef(model)[coefs]; V <- vcov_mat[coefs, coefs]
  stat <- as.numeric(t(b) %*% solve(V) %*% b)
  data.frame(statistic = stat, df = length(coefs), p_value_asymptotic = pchisq(stat, df = length(coefs), lower.tail = FALSE))
}

## ================================================================
## Reproduce the locked Stage-14 K=3 sample (read-only calibration;
## no Stage 14-18 file is modified).
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
W_queen  <- listw2mat(lw_queen); dimnames(W_queen) <- list(region_order, region_order)

panel <- panel %>% mutate(month_index = as.integer(format(date, "%Y")) * 12L + as.integer(format(date, "%m")))
wide <- panel %>% select(month_index, region_id, inflation_mom) %>%
  pivot_wider(names_from = region_id, values_from = inflation_mom) %>% arrange(month_index)
wide <- wide[, c("month_index", region_order)]
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
for (k in 1:3) {
  model_data[[paste0("own_lag_", k)]] <- lookup(pi_mat, model_data$month_index - k, model_data$region_id)
  model_data[[paste0("w_lag_", k)]]   <- lookup(spatial_lag_queen, model_data$month_index - k, model_data$region_id)
}

own_cols_k3 <- paste0("own_lag_", 1:3); w_cols_k3 <- paste0("w_lag_", 1:3)
complete_rows <- function(df, cols) Reduce(`&`, lapply(cols, function(cc) !is.na(df[[cc]])))
sample_k3 <- model_data[complete_rows(model_data, c(own_cols_k3, w_cols_k3)), ]
sample_k3 <- sample_k3[order(sample_k3$month_index, sample_k3$region_id), ]
sample_k3$region_id_f <- factor(sample_k3$region_id, levels = region_order)
sample_k3$date_f      <- factor(sample_k3$date)

sample_check <- data.frame(n_regions = length(unique(sample_k3$region_id)), n_months = length(unique(sample_k3$date)),
                            total_obs = nrow(sample_k3), first_date = format(min(sample_k3$date)), last_date = format(max(sample_k3$date)))
s14 <- read.csv("results_v2/distributed_spatial_lags/tables/dsl_estimation_samples.csv", stringsAsFactors = FALSE)[1, ]
matches_stage14 <- sample_check$n_regions == s14$n_regions && sample_check$n_months == s14$n_months &&
  sample_check$total_obs == s14$total_obs && sample_check$first_date == s14$first_date && sample_check$last_date == s14$last_date
message("K=3 sample matches Stage 14's saved sample exactly: ", matches_stage14)
if (!matches_stage14) stop("K=3 sample does NOT match Stage 14 -- investigate before proceeding.")

f_own3 <- paste(own_cols_k3, collapse = " + "); f_w3 <- paste(w_cols_k3, collapse = " + ")
f_restricted <- paste("inflation_mom ~ region_id_f + date_f +", f_own3)
f_full       <- paste(f_restricted, "+", f_w3)

m0 <- lm(as.formula(f_restricted), data = sample_k3)
fitted0 <- fitted(m0)
resid0  <- residuals(m0)

## ================================================================
## Fixed-regressor DGP scaffolding: residual matrix (T x 13, ordered by
## date then the frozen region_order), and (date_idx, region_idx) map
## for each row of sample_k3 for fast vectorised lookup during simulation.
## ================================================================

dates_sorted <- sort(unique(sample_k3$date))
Tn <- length(dates_sorted)
resid_mat <- matrix(NA_real_, nrow = Tn, ncol = 13, dimnames = list(as.character(dates_sorted), region_order))
resid_mat[cbind(match(sample_k3$date, dates_sorted), match(sample_k3$region_id, region_order))] <- resid0
stopifnot(!anyNA(resid_mat))  # sample_k3 is a balanced 13 x Tn panel, so this must be fully populated

date_idx   <- match(sample_k3$date, dates_sorted)
region_idx <- match(sample_k3$region_id, region_order)

BLOCK_LEN_PRIMARY <- 6L   # ~ T^(1/3) = 5.36 for T=154; documented rule-of-thumb, not tuned
BLOCK_LEN_ALT     <- 12L  # sensitivity: one calendar year

moving_block_bootstrap <- function(resid_mat, block_len) {
  Tn <- nrow(resid_mat)
  n_blocks <- ceiling(Tn / block_len)
  starts <- sample(1:(Tn - block_len + 1), n_blocks, replace = TRUE)
  out <- do.call(rbind, lapply(starts, function(s) resid_mat[s:(s + block_len - 1), , drop = FALSE]))
  out[1:Tn, , drop = FALSE]
}

lag_shape_equal       <- function(Theta3) rep(Theta3 / 3, 3)
lag_shape_frontloaded <- function(Theta3) Theta3 * c(0.6, 0.3, 0.1)

## ================================================================
## Core simulation: one replication -> DK-based joint-test rejection
## indicator (and, optionally, the WCB rejection indicator too).
## ================================================================

simulate_one <- function(theta_true, block_len, run_wcb = FALSE, wcb_B = 199) {
  e_star_mat <- moving_block_bootstrap(resid_mat, block_len)
  e_star <- e_star_mat[cbind(date_idx, region_idx)]
  y_star <- fitted0 + theta_true[1] * sample_k3$w_lag_1 + theta_true[2] * sample_k3$w_lag_2 +
    theta_true[3] * sample_k3$w_lag_3 + e_star
  dat <- sample_k3
  dat$inflation_mom <- y_star
  m1 <- lm(as.formula(f_full), data = dat)
  vcov1 <- dk_vcov(m1, dat$region_id_f, dat$date_f)
  wald1 <- wald_joint(m1, vcov1, w_cols_k3)
  reject_dk <- as.integer(wald1$p_value_asymptotic < 0.05)

  reject_wcb <- NA_integer_
  if (run_wcb) {
    m1_restricted <- lm(as.formula(f_restricted), data = dat)
    fitted_r <- fitted(m1_restricted); resid_r <- residuals(m1_restricted)
    vcov1_cl <- vcovCL(m1, cluster = dat$region_id_f, type = "HC1")
    wald_obs <- as.numeric(t(coef(m1)[w_cols_k3]) %*% solve(vcov1_cl[w_cols_k3, w_cols_k3]) %*% coef(m1)[w_cols_k3])
    regions_unique <- levels(dat$region_id_f)
    wald_boot <- numeric(wcb_B)
    for (b in seq_len(wcb_B)) {
      v <- setNames(sample(c(-1, 1), length(regions_unique), replace = TRUE), regions_unique)
      yb <- fitted_r + resid_r * v[as.character(dat$region_id_f)]
      bd <- dat; bd$inflation_mom <- yb
      mb <- lm(as.formula(f_full), data = bd)
      vcb <- vcovCL(mb, cluster = bd$region_id_f, type = "HC1")
      bb <- coef(mb)[w_cols_k3]
      wald_boot[b] <- as.numeric(t(bb) %*% solve(vcb[w_cols_k3, w_cols_k3]) %*% bb)
    }
    p_wcb <- mean(wald_boot >= wald_obs)
    reject_wcb <- as.integer(p_wcb < 0.05)
  }
  c(reject_dk = reject_dk, reject_wcb = reject_wcb)
}

## ================================================================
## 3B/3C/3D/3E: primary power curve (equal lag-shape, block=6).
## ================================================================

GRID <- c(0, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20)
GRID_SIGNED <- sort(unique(c(-GRID, GRID)))  # 13 points, symmetric, includes 0 once
R_PRIMARY <- 2000L

message("\nRunning PRIMARY power curve: ", length(GRID_SIGNED), " grid points x ", R_PRIMARY,
        " reps (equal lag-shape, block length ", BLOCK_LEN_PRIMARY, ")...")
primary_results <- bind_rows(lapply(GRID_SIGNED, function(g) {
  theta_true <- lag_shape_equal(g)
  rej <- vapply(seq_len(R_PRIMARY), function(r) simulate_one(theta_true, BLOCK_LEN_PRIMARY)[["reject_dk"]], numeric(1))
  power_hat <- mean(rej)
  mc_se <- sqrt(power_hat * (1 - power_hat) / R_PRIMARY)
  message("  Theta_3=", g, ": power=", round(power_hat, 3), " (MC SE=", round(mc_se, 4), ")")
  data.frame(Theta_3 = g, lag_shape = "Equal (primary)", block_len = BLOCK_LEN_PRIMARY, n_reps = R_PRIMARY,
             power = power_hat, mc_se = mc_se)
}))

## ================================================================
## Secondary sensitivity: front-loaded lag-shape, same grid, fewer reps.
## ================================================================

R_SECONDARY <- 1000L
message("\nRunning SECONDARY power curve: front-loaded lag-shape, ", length(GRID_SIGNED), " points x ", R_SECONDARY, " reps...")
frontloaded_results <- bind_rows(lapply(GRID_SIGNED, function(g) {
  theta_true <- lag_shape_frontloaded(g)
  rej <- vapply(seq_len(R_SECONDARY), function(r) simulate_one(theta_true, BLOCK_LEN_PRIMARY)[["reject_dk"]], numeric(1))
  power_hat <- mean(rej)
  mc_se <- sqrt(power_hat * (1 - power_hat) / R_SECONDARY)
  message("  Theta_3=", g, ": power=", round(power_hat, 3))
  data.frame(Theta_3 = g, lag_shape = "Front-loaded (sensitivity)", block_len = BLOCK_LEN_PRIMARY, n_reps = R_SECONDARY,
             power = power_hat, mc_se = mc_se)
}))

## ================================================================
## Block-length sensitivity (L=12), reduced grid, equal shape.
## ================================================================

GRID_REDUCED <- c(-0.20, -0.10, 0, 0.10, 0.20)
R_BLOCKALT <- 500L
message("\nRunning block-length sensitivity (L=", BLOCK_LEN_ALT, "), ", length(GRID_REDUCED), " points x ", R_BLOCKALT, " reps...")
blockalt_results <- bind_rows(lapply(GRID_REDUCED, function(g) {
  theta_true <- lag_shape_equal(g)
  rej <- vapply(seq_len(R_BLOCKALT), function(r) simulate_one(theta_true, BLOCK_LEN_ALT)[["reject_dk"]], numeric(1))
  power_hat <- mean(rej)
  mc_se <- sqrt(power_hat * (1 - power_hat) / R_BLOCKALT)
  message("  Theta_3=", g, ": power=", round(power_hat, 3))
  data.frame(Theta_3 = g, lag_shape = "Equal (block-length sensitivity)", block_len = BLOCK_LEN_ALT, n_reps = R_BLOCKALT,
             power = power_hat, mc_se = mc_se)
}))

power_curve <- bind_rows(primary_results, frontloaded_results, blockalt_results)
write.csv(power_curve, "results_v2/presubmission_refinement/power_mde/power_curve.csv", row.names = FALSE)

## ================================================================
## 3F validation: DK-based vs. nested-WCB-based rejection rate at 3
## selected effect sizes (equal shape, block=6), smaller rep count.
## ================================================================

VALIDATION_POINTS <- c(0, 0.10, 0.20)
R_VALIDATION <- 100L
WCB_B_VALIDATION <- 199L
message("\nRunning DK-vs-WCB validation: ", length(VALIDATION_POINTS), " points x ", R_VALIDATION,
        " reps x WCB(B=", WCB_B_VALIDATION, ")...")
validation_results <- bind_rows(lapply(VALIDATION_POINTS, function(g) {
  theta_true <- lag_shape_equal(g)
  out <- vapply(seq_len(R_VALIDATION), function(r) simulate_one(theta_true, BLOCK_LEN_PRIMARY, run_wcb = TRUE, wcb_B = WCB_B_VALIDATION),
                numeric(2))
  power_dk  <- mean(out["reject_dk", ])
  power_wcb <- mean(out["reject_wcb", ])
  message("  Theta_3=", g, ": power_DK=", round(power_dk, 3), " power_WCB(B=", WCB_B_VALIDATION, ")=", round(power_wcb, 3))
  data.frame(Theta_3 = g, n_reps = R_VALIDATION, wcb_B = WCB_B_VALIDATION, power_DK = power_dk, power_WCB = power_wcb,
             agreement_within_2SE = abs(power_dk - power_wcb) <= 2 * sqrt(0.25 / R_VALIDATION))
}))
write.csv(validation_results, "results_v2/presubmission_refinement/power_mde/simulation_validation.csv", row.names = FALSE)
message("\n---- DK-vs-WCB validation summary ----")
print(validation_results)

## ================================================================
## 3G: MDE_80 / MDE_90, simulated (primary curve) + analytic sanity check.
## ================================================================

# Fold the signed primary grid into |Theta_3| by averaging +/- power at each magnitude (symmetric-by-construction joint test).
primary_by_magnitude <- primary_results %>%
  mutate(abs_theta = abs(Theta_3)) %>%
  group_by(abs_theta) %>%
  summarise(power = mean(power), n_reps = sum(n_reps), .groups = "drop") %>%
  arrange(abs_theta)

interpolate_mde <- function(df, target_power) {
  if (all(df$power < target_power)) return(NA_real_)  # never reached within the grid
  if (df$power[1] >= target_power) return(df$abs_theta[1])
  idx <- which(df$power >= target_power)[1]
  x0 <- df$abs_theta[idx - 1]; x1 <- df$abs_theta[idx]
  y0 <- df$power[idx - 1];     y1 <- df$power[idx]
  x0 + (target_power - y0) * (x1 - x0) / (y1 - y0)
}
MDE_80 <- interpolate_mde(primary_by_magnitude, 0.80)
MDE_90 <- interpolate_mde(primary_by_magnitude, 0.90)

# Analytic single-parameter sanity check, from the LOCKED Stage-18 observed cumulative-effect DK SE (read-only reference).
fr_final <- read.csv("results_v2/final_robustness/tables/fr_K_final_result_table.csv", stringsAsFactors = FALSE)
observed_theta3_se <- fr_final$dk_se[fr_final$parameter == "Cumulative spatial effect (Theta_3)"]
observed_theta3_est <- fr_final$estimate[fr_final$parameter == "Cumulative spatial effect (Theta_3)"]
analytic_MDE_80 <- (qnorm(0.975) + qnorm(0.80)) * observed_theta3_se
analytic_MDE_90 <- (qnorm(0.975) + qnorm(0.90)) * observed_theta3_se

max_grid_power <- max(primary_by_magnitude$power)
mde_summary <- data.frame(
  quantity = c("MDE_80 (simulated, primary curve)", "MDE_90 (simulated, primary curve)",
               "MDE_80 (analytic sanity check)", "MDE_90 (analytic sanity check)",
               "Max power achieved within the pre-specified grid (|Theta_3|<=0.20)",
               "Observed Theta_3 estimate (Stage 18, locked)", "Observed Theta_3 DK SE (Stage 18, locked)"),
  value = c(MDE_80, MDE_90, analytic_MDE_80, analytic_MDE_90, max_grid_power, observed_theta3_est, observed_theta3_se),
  note = c(
    if (is.na(MDE_80)) "80% power NOT reached within the pre-specified grid (|Theta_3] <= 0.20) -- reported as exceeding the grid, not extrapolated" else "linear interpolation between adjacent grid points",
    if (is.na(MDE_90)) "90% power NOT reached within the pre-specified grid" else "linear interpolation between adjacent grid points",
    "single-parameter z-approximation using the LOCKED observed DK SE; sanity check only, not a substitute for the simulation",
    "single-parameter z-approximation using the LOCKED observed DK SE; sanity check only",
    "", "read-only from fr_K_final_result_table.csv (Stage 18)", "read-only from fr_K_final_result_table.csv (Stage 18)"
  )
)
write.csv(mde_summary, "results_v2/presubmission_refinement/power_mde/mde_summary.csv", row.names = FALSE)
message("\n---- MDE summary ----")
print(mde_summary[, c("quantity", "value")])

## ================================================================
## 3I: power classification (P1/P2/P3), from simulation results only.
## ================================================================

power_classification <- if (!is.na(MDE_80) && MDE_80 <= 0.05) {
  "P1: Adequate power to detect small effects"
} else if (max_grid_power >= 0.80) {
  "P2: Adequate power for moderate/large effects but limited power for small effects"
} else {
  "P3: Low power even for moderate effects"
}
message("\n================================================================")
message("POWER CLASSIFICATION: ", power_classification)
message("================================================================")
writeLines(power_classification, "results_v2/presubmission_refinement/power_mde/power_classification.txt")

## ================================================================
## Simulation design documentation table.
## ================================================================

simulation_design <- data.frame(
  parameter = c("Estimand", "Effect-size grid (pre-specified)", "Primary lag-shape", "Secondary lag-shape",
                "DGP", "Primary block length", "Block-length sensitivity", "Primary replications per point",
                "Secondary-shape replications per point", "Block-sensitivity replications per point",
                "Primary inference rule", "Validation inference rule", "Validation replications", "Seed"),
  value = c(
    "Theta_3 = theta_1+theta_2+theta_3 (cumulative 3-month spatial predictive effect)",
    paste(GRID_SIGNED, collapse = ", "),
    "Equal: theta_k = Theta_3/3 for k=1,2,3",
    "Front-loaded: theta = Theta_3 x (0.6, 0.3, 0.1)",
    "Fixed-regressor (observed own/spatial lags, region/month FE); y* = fitted(restricted model) + true spatial effect + moving-block-bootstrapped residual",
    paste0(BLOCK_LEN_PRIMARY, " months (~T^(1/3)=", round(Tn^(1/3), 2), ")"),
    paste0(BLOCK_LEN_ALT, " months, reduced grid (", paste(GRID_REDUCED, collapse=", "), "), R=", R_BLOCKALT),
    as.character(R_PRIMARY), as.character(R_SECONDARY), as.character(R_BLOCKALT),
    "Driscoll-Kraay joint Wald test (vcovPL, Bartlett, NW1987), same as Stage 14/18 primary inference",
    paste0("Wild cluster bootstrap, B=", WCB_B_VALIDATION, ", at Theta_3 in {", paste(VALIDATION_POINTS, collapse=", "), "}"),
    as.character(R_VALIDATION), as.character(GLOBAL_SEED)
  )
)
write.csv(simulation_design, "results_v2/presubmission_refinement/power_mde/simulation_design.csv", row.names = FALSE)

## ================================================================
## Figures.
## ================================================================

fig_dir <- "results_v2/presubmission_refinement/power_mde"

p1 <- ggplot(primary_results, aes(x = Theta_3, y = power)) +
  geom_hline(yintercept = c(0.05, 0.80, 0.90), linetype = "dotted", color = "grey60") +
  geom_line(color = "#2c3e50") + geom_point(color = "#2c3e50") +
  geom_errorbar(aes(ymin = pmax(0, power - 1.96 * mc_se), ymax = pmin(1, power + 1.96 * mc_se)), width = 0.005, color = "#2c3e50") +
  labs(title = "Simulated Power Curve vs. True Cumulative Spatial Effect (Theta_3)",
       subtitle = paste0("Equal lag-shape, block length ", BLOCK_LEN_PRIMARY, " months, R=", R_PRIMARY, " reps/point; dotted lines = 5%/80%/90%"),
       x = "True Theta_3", y = "Estimated power (DK joint Wald test, alpha=0.05)") +
  theme_minimal()
ggsave(file.path(fig_dir, "fig_1_power_curve.png"), p1, width = 8, height = 5, dpi = 300)

p2 <- ggplot(bind_rows(primary_results, frontloaded_results), aes(x = Theta_3, y = power, color = lag_shape)) +
  geom_hline(yintercept = 0.80, linetype = "dotted", color = "grey60") +
  geom_line() + geom_point() +
  labs(title = "Power Curve by Assumed Lag-Shape",
       subtitle = "Power is conditional on the assumed distribution of Theta_3 across lags 1-3",
       x = "True Theta_3", y = "Estimated power", color = "Lag shape") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "fig_2_power_by_lag_shape.png"), p2, width = 8, height = 5, dpi = 300)

detect_plot_data <- data.frame(
  label = c("Observed estimate", "Simulated grid (|Theta_3|)"),
  x = c(observed_theta3_est, NA)
)
p3 <- ggplot() +
  geom_rect(aes(xmin = -max(GRID_SIGNED), xmax = max(GRID_SIGNED), ymin = 0, ymax = 1), fill = "#eef3f8") +
  geom_point(data = data.frame(x = GRID_SIGNED, y = 0.5), aes(x = x, y = y), shape = 4, color = "#2c7bb6", size = 2) +
  geom_vline(xintercept = observed_theta3_est, color = "#2c3e50", linewidth = 1) +
  geom_vline(xintercept = fr_final$ci_low[fr_final$parameter == "Cumulative spatial effect (Theta_3)"], linetype = "dashed", color = "#2c3e50") +
  geom_vline(xintercept = fr_final$ci_high[fr_final$parameter == "Cumulative spatial effect (Theta_3)"], linetype = "dashed", color = "#2c3e50") +
  geom_vline(xintercept = c(-1, 1) * ifelse(is.na(MDE_80), max(GRID_SIGNED), MDE_80), color = "#d7191c", linetype = "dotdash") +
  annotate("text", x = observed_theta3_est, y = 0.85, label = "Observed estimate", angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = ifelse(is.na(MDE_80), max(GRID_SIGNED), MDE_80), y = 0.15, label = ifelse(is.na(MDE_80), "MDE_80 not reached in grid", "MDE_80"), angle = 90, vjust = -0.5, size = 3, color = "#d7191c") +
  labs(title = "Observed Estimate and CI Against the Simulated Detectable-Effect Scale",
       subtitle = "Blue x's = simulation grid points; solid black = observed Theta_3; dashed black = 95% CI; red dot-dash = MDE_80 (or grid boundary if not reached)",
       x = "Theta_3", y = NULL) +
  theme_minimal() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave(file.path(fig_dir, "fig_3_observed_vs_detectable_scale.png"), p3, width = 9, height = 4.5, dpi = 300)

message("\nWrote figures 1-3 to ", fig_dir, "/")

## ================================================================
## POWER_MDE_REPORT.md
## ================================================================

md <- c(
"# MAKANI v2 -- Simulation-Based Power / Minimum-Detectable-Effect (MDE) Analysis",
"",
"**This is a design-sensitivity analysis, not a new empirical finding.** It",
"does not re-estimate, reinterpret, or alter any Stage 14-18 result. See",
"`scripts/20_power_mde_presubmission.R` for the full pre-specified design.",
"",
"## Question answered",
"",
"Given N=13, T=154 (the locked Stage-14 K=3 estimation sample), the",
"observed Queen spatial-weight structure, the locked dynamic",
"specification, and the observed residual dependence, what magnitude of",
"cumulative 3-month spatial predictive effect (Theta_3) could this design",
"plausibly detect at conventional power (80%/90%)?",
"",
"## Design summary",
"",
sprintf("- Effect-size grid (pre-specified, NOT adjusted after seeing results): %s", paste(GRID_SIGNED, collapse = ", ")),
"- Primary lag-shape: equal allocation across lags 1-3. Secondary (sensitivity): front-loaded (0.6/0.3/0.1).",
sprintf("- Fixed-regressor DGP with moving-block-bootstrapped residuals (block length %d months primary, %d months sensitivity).", BLOCK_LEN_PRIMARY, BLOCK_LEN_ALT),
sprintf("- Primary inference: Driscoll-Kraay joint Wald test (same as Stage 14/18). Validated against a B=%d wild-cluster bootstrap at 3 effect sizes.", WCB_B_VALIDATION),
sprintf("- Replications: %d (primary), %d (front-loaded sensitivity), %d (block-length sensitivity), %d x WCB(B=%d) (validation). Seed=%d, set once.",
        R_PRIMARY, R_SECONDARY, R_BLOCKALT, R_VALIDATION, WCB_B_VALIDATION, GLOBAL_SEED),
"",
"## Validation: DK asymptotic vs. nested wild-cluster bootstrap",
"",
"See `simulation_validation.csv`. The DK-based rejection rate is used as",
"the primary simulation inference rule throughout the full grid (running",
"a full B=9999 WCB inside thousands of replications is not computationally",
"tractable); this validation checks that DK-based power tracks WCB-based",
"power reasonably at 3 selected effect sizes, on a reduced replication",
"count -- it is an approximation check, not a claim that the full-grid",
"power curve IS wild-bootstrap power.",
"",
"## Power curve",
"",
"See `power_curve.csv` and Figures 1-2.",
sprintf("Maximum power achieved anywhere within the pre-specified grid (|Theta_3| <= 0.20): %.3f.", max_grid_power),
"",
"## MDE",
"",
"See `mde_summary.csv`.",
sprintf("- MDE_80 (simulated): %s", if (is.na(MDE_80)) "NOT reached within the pre-specified grid (|Theta_3| <= 0.20)" else sprintf("%.4f", MDE_80)),
sprintf("- MDE_90 (simulated): %s", if (is.na(MDE_90)) "NOT reached within the pre-specified grid" else sprintf("%.4f", MDE_90)),
sprintf("- MDE_80 (analytic sanity check, single-parameter z-approximation from the locked observed DK SE = %.4f): %.4f", observed_theta3_se, analytic_MDE_80),
sprintf("- MDE_90 (analytic sanity check): %.4f", analytic_MDE_90),
"The analytic approximation is a sanity check only and is not substituted for the simulation result.",
"",
"## Power classification",
"",
sprintf("**%s**", power_classification),
"",
if (power_classification != "P1: Adequate power to detect small effects") {
  "This design cannot reliably distinguish a true zero cumulative spatial effect from a small-to-moderate one within the pre-specified grid. Small spatial effects (well within the range of what would be economically unsurprising) CANNOT be ruled out by this study's null/weak finding -- absence of evidence is not evidence of absence here. Only effects at or above the simulated MDE (or, if MDE_80 was not reached within the grid, effects larger than the largest grid point tested, |Theta_3|=0.20) would have been reliably detected."
} else {
  "This design has adequate power to detect even small cumulative spatial effects, strengthening the interpretation of the null/weak empirical result as informative about the absence of even modest spatial predictive dependence."
},
"",
"## Limitations of this power analysis",
"",
"- Fixed-regressor design: own-lag and spatial-lag regressor VALUES are held at their observed values rather than being regenerated recursively from a fully simulated dynamic panel (which would itself require assuming a data-generating process for how a true spatial effect would alter neighbours' inflation period by period). This is a standard simplification for power analysis, stated explicitly.",
"- Moving-block-bootstrap residuals approximate, but do not perfectly reproduce, the true (unknown) joint temporal/cross-sectional dependence structure of the actual error process.",
"- The DK-based rejection rule is validated against WCB at only 3 effect sizes and 100 replications each -- a modest validation sample, not an exhaustive one; Monte Carlo standard errors for that validation are reported in `simulation_validation.csv`.",
"- Power is estimated only at the pre-specified grid points; values between them are linearly interpolated for MDE_80/MDE_90, not separately simulated.",
"- This analysis assumes the SAME lag horizon (K=3) and spatial-weight structure (Queen) as the locked design; it does not assess power for K=6, alternative weights, or nonlinear specifications, which were not reopened per this stage's guardrails.",
"",
"## Files",
"",
"`power_curve.csv`, `mde_summary.csv`, `simulation_design.csv`, `simulation_validation.csv`, `power_classification.txt`, `fig_1_power_curve.png`, `fig_2_power_by_lag_shape.png`, `fig_3_observed_vs_detectable_scale.png`.",
"Script: `scripts/20_power_mde_presubmission.R`."
)
writeLines(md, "results_v2/presubmission_refinement/power_mde/POWER_MDE_REPORT.md")
message("\nWrote results_v2/presubmission_refinement/power_mde/POWER_MDE_REPORT.md")
message("20_power_mde_presubmission.R complete.")
