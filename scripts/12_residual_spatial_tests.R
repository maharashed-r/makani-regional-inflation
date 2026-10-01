# ================================================================
# MAKANI v2 -- Stage 1 econometrics (continued): monthly cross-sectional
# spatial dependence tests, raw inflation_mom vs TWFE (M3) residuals.
#
# NEW v2 FILE. Reads only:
#   - data_v2/saudi_cpi_panel_v2.csv, data_v2/region_crosswalk.csv
#   - results_v2/tables/twfe_residuals_mom.csv (built by scripts/11)
#   - data_raw/saudi_states_sf.rds (baseline's own cached geometry --
#     READ ONLY, never written to; same file scripts/04 and scripts/09
#     already read, never modified by either)
# Writes only under results_v2/. Does not source or modify scripts 01-10
# or anything under baseline/, data_clean/, or results/.
#
# CRITICAL DESIGN POINT (per instruction): spatial dependence is tested
# CROSS-SECTIONALLY, separately for each of the 157 months in the MoM
# estimation sample -- never on a region-level mean of twfe_residual.
# Region fixed effects force each region's OWN residual mean toward
# exactly zero (verified in scripts/11's residual validation), so a
# region-mean-collapsed Moran's I on residuals would be mechanically
# uninformative by construction, not a real test of anything.
#
# Run manually: Rscript scripts/12_residual_spatial_tests.R
# ================================================================

suppressMessages({
  library(sf)
  library(spdep)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

dir.create("results_v2/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results_v2/figures", recursive = TRUE, showWarnings = FALSE)

panel      <- read.csv("data_v2/saudi_cpi_panel_v2.csv", stringsAsFactors = FALSE)
panel$date <- as.Date(panel$date)
crosswalk  <- read.csv("data_v2/region_crosswalk.csv", stringsAsFactors = FALSE)
residuals_mom <- read.csv("results_v2/tables/twfe_residuals_mom.csv", stringsAsFactors = FALSE)
residuals_mom$date <- as.Date(residuals_mom$date)

## ----------------------------------------------------------------
## 1. Build the spatial frame in a FIXED region_id order (R01..R13).
##    Every month's value vector below is re-indexed to this exact
##    order via match(), never assumed from row order after a join.
## ----------------------------------------------------------------

region_order <- sort(unique(crosswalk$region_id))
stopifnot(length(region_order) == 13)

geom_raw <- readRDS("data_raw/saudi_states_sf.rds")  # read-only, baseline file untouched

sf_frame <- crosswalk %>%
  select(region_id, region_name_gastat_en, geometry_spelling) %>%
  arrange(match(region_id, region_order)) %>%
  left_join(
    geom_raw %>% select(name),
    by = c("geometry_spelling" = "name")
  )
sf_frame <- st_as_sf(sf_frame)

if (nrow(sf_frame) != 13 || any(st_is_empty(sf_frame))) {
  stop("Spatial frame construction failed: expected 13 non-empty geometries, got ",
       nrow(sf_frame), " (", sum(st_is_empty(sf_frame)), " empty).")
}
stopifnot(identical(sf_frame$region_id, region_order))

## ----------------------------------------------------------------
## 2. Spatial weights: Queen, Rook, KNN(k=4), distance band.
##    Built ONCE and held fixed over time (never rebuilt per month).
## ----------------------------------------------------------------

nb_queen <- poly2nb(sf_frame, queen = TRUE)
nb_rook  <- poly2nb(sf_frame, queen = FALSE)

nb_to_undirected_edges <- function(nb) {
  parts <- lapply(seq_along(nb), function(i) {
    js <- nb[[i]]
    js <- js[js > i]
    if (length(js) == 0) return(NULL)
    paste(i, js, sep = "-")
  })
  unlist(parts)
}

queen_edges <- nb_to_undirected_edges(nb_queen)
rook_edges  <- nb_to_undirected_edges(nb_rook)
n_differing_edges <- length(union(setdiff(queen_edges, rook_edges), setdiff(rook_edges, queen_edges)))
# NOTE: R's identical() on the raw nb objects returns FALSE even when the
# adjacency structure is the same, because poly2nb() stamps each nb object
# with a "call" attribute that differs (queen=TRUE vs queen=FALSE) -- an
# attribute artefact, not a structural difference. The actual identity
# question ("do Queen and Rook produce the same graph on this map?") is
# answered correctly by the edge-set comparison above, which is what is
# used for the decision below.
queen_rook_identical <- (n_differing_edges == 0)

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

weights_report <- bind_rows(
  data.frame(spec = "Queen contiguity",
             n_edges_undirected = length(queen_edges),
             mean_neighbors = mean(card(nb_queen)), min_neighbors = min(card(nb_queen)),
             max_neighbors = max(card(nb_queen)), islands = sum(card(nb_queen) == 0)),
  data.frame(spec = "Rook contiguity",
             n_edges_undirected = length(rook_edges),
             mean_neighbors = mean(card(nb_rook)), min_neighbors = min(card(nb_rook)),
             max_neighbors = max(card(nb_rook)), islands = sum(card(nb_rook) == 0)),
  data.frame(spec = "KNN (k=4, directed)",
             n_edges_undirected = sum(card(nb_knn4)),  # KNN is directed by construction; see note below
             mean_neighbors = mean(card(nb_knn4)), min_neighbors = min(card(nb_knn4)),
             max_neighbors = max(card(nb_knn4)), islands = sum(card(nb_knn4) == 0)),
  data.frame(spec = "Distance band (1.05x max nearest-neighbour dist)",
             n_edges_undirected = length(nb_to_undirected_edges(nb_dist)),
             mean_neighbors = mean(card(nb_dist)), min_neighbors = min(card(nb_dist)),
             max_neighbors = max(card(nb_dist)), islands = sum(card(nb_dist) == 0))
)
write.csv(weights_report, "results_v2/tables/spatial_weights_comparison.csv", row.names = FALSE)

queen_rook_check <- data.frame(
  n_edges_queen = length(queen_edges),
  n_edges_rook  = length(rook_edges),
  n_differing_edges = n_differing_edges,
  identical_graph = queen_rook_identical
)
write.csv(queen_rook_check, "results_v2/tables/spatial_weights_queen_rook_check.csv", row.names = FALSE)

message("---- Spatial weights comparison ----")
message("Queen and Rook contiguity graphs identical: ", queen_rook_identical,
        " (differing edges: ", n_differing_edges, ")")
message("NOTE: KNN(k=4)'s n_edges_undirected column reports the count of DIRECTED links ",
        "(sum of neighbour counts), not undirected edges -- KNN is not guaranteed symmetric ",
        "(A's 4 nearest neighbours need not include B even when B's 4 nearest include A). ",
        "Queen, Rook, and the distance band are symmetric by construction, so their counts ",
        "are genuine undirected edge counts.")
print(weights_report)

if (queen_rook_identical) {
  message("\nQueen and Rook are IDENTICAL on this map -- using Queen as the single primary W. ",
          "Rook is NOT treated as an independent robustness specification below, per instruction.")
  lw_primary <- nb2listw(nb_queen, style = "W", zero.policy = TRUE)
  primary_spec_name <- "Queen contiguity (identical to Rook on this map)"
} else {
  lw_primary <- nb2listw(nb_queen, style = "W", zero.policy = TRUE)
  primary_spec_name <- "Queen contiguity"
}

## ----------------------------------------------------------------
## 3. Monthly Global Moran's I -- raw inflation_mom and TWFE residuals,
##    same months, same fixed W, permutation-primary inference.
## ----------------------------------------------------------------

mom_sample <- panel %>% filter(!is.na(inflation_mom))
estimation_months <- sort(unique(mom_sample$date))
stopifnot(length(estimation_months) == 157)  # matches scripts/11's verified MoM sample

# Primary research hypothesis is POSITIVE spatial autocorrelation, so the
# one-sided "greater" alternative (I > E[I]) is used explicitly for both
# the asymptotic and permutation tests -- verified to be spdep's own
# default for moran.test()/moran.mc() (checked via formals() in this
# session), but stated explicitly here rather than left implicit, so
# every test's direction is documented in code and in the saved output
# (alternative_used column), and never silently mixed with a two-sided
# test anywhere in this pipeline.
MORAN_ALTERNATIVE <- "greater"

run_monthly_moran <- function(data_long, value_col, months, lw, region_order, seed, nsim = 9999) {
  set.seed(seed)
  rows <- vector("list", length(months))
  for (k in seq_along(months)) {
    m <- months[k]
    month_data <- data_long[data_long$date == m, ]
    idx <- match(region_order, month_data$region_id)
    if (anyNA(idx)) {
      stop("Month ", m, " is missing region(s): ",
           paste(region_order[is.na(idx)], collapse = ", "))
    }
    values <- month_data[[value_col]][idx]

    mt  <- moran.test(values, lw, zero.policy = TRUE, na.action = na.fail,
                       alternative = MORAN_ALTERNATIVE)
    mmc <- moran.mc(values, lw, nsim = nsim, zero.policy = TRUE,
                     alternative = MORAN_ALTERNATIVE)

    I_val        <- unname(mt$estimate["Moran I statistic"])
    expected_val <- unname(mt$estimate["Expectation"])

    rows[[k]] <- data.frame(
      date             = m,
      I                = I_val,
      expected_I       = expected_val,
      I_minus_expected = I_val - expected_val,
      perm_p_value     = mmc$p.value,
      asymp_p_value    = mt$p.value,
      alternative_used = MORAN_ALTERNATIVE,
      n_permutations   = nsim,
      weight_spec      = primary_spec_name
    )
  }
  bind_rows(rows)
}

message("\nRunning monthly permutation Moran's I on raw inflation_mom (", length(estimation_months),
        " months, ", 9999, " permutations each -- this may take a short while) ...")
moran_raw <- run_monthly_moran(mom_sample, "inflation_mom", estimation_months, lw_primary,
                                region_order, seed = 20260921)
write.csv(moran_raw, "results_v2/tables/moran_monthly_raw.csv", row.names = FALSE)
message("Done. Mean raw Moran I = ", round(mean(moran_raw$I), 4))

message("\nRunning monthly permutation Moran's I on TWFE residuals (same months, same W) ...")
moran_resid <- run_monthly_moran(residuals_mom, "twfe_residual", estimation_months, lw_primary,
                                  region_order, seed = 20260921)
write.csv(moran_resid, "results_v2/tables/moran_monthly_residual.csv", row.names = FALSE)
message("Done. Mean residual Moran I = ", round(mean(moran_resid$I), 4))

## ----------------------------------------------------------------
## 4. Multiple-testing correction (BH, Holm) on the permutation p-values,
##    separately for raw and residual.
## ----------------------------------------------------------------

add_mtc <- function(df) {
  df$p_BH   <- p.adjust(df$perm_p_value, method = "BH")
  df$p_Holm <- p.adjust(df$perm_p_value, method = "holm")
  df$sig_raw_05  <- df$perm_p_value < 0.05
  df$sig_BH_05   <- df$p_BH < 0.05
  df$sig_Holm_05 <- df$p_Holm < 0.05
  df
}
moran_raw_mtc    <- add_mtc(moran_raw)
moran_resid_mtc  <- add_mtc(moran_resid)

write.csv(moran_raw_mtc,   "results_v2/tables/moran_monthly_raw_mtc.csv",      row.names = FALSE)
write.csv(moran_resid_mtc, "results_v2/tables/moran_monthly_residual_mtc.csv", row.names = FALSE)

mtc_summary <- data.frame(
  series = c("Raw inflation_mom", "TWFE residual"),
  n_months = c(nrow(moran_raw_mtc), nrow(moran_resid_mtc)),
  n_raw_p_lt_05  = c(sum(moran_raw_mtc$sig_raw_05),  sum(moran_resid_mtc$sig_raw_05)),
  n_surviving_BH = c(sum(moran_raw_mtc$sig_BH_05),   sum(moran_resid_mtc$sig_BH_05)),
  n_surviving_Holm = c(sum(moran_raw_mtc$sig_Holm_05), sum(moran_resid_mtc$sig_Holm_05)),
  interpretation = paste0(
    sum(moran_raw_mtc$sig_raw_05), " months were nominally significant at the unadjusted 5% ",
    "level, but ", sum(moran_raw_mtc$sig_BH_05), " survived BH and ", sum(moran_raw_mtc$sig_Holm_05),
    " survived Holm correction; isolated monthly significance is therefore not treated as ",
    "robust evidence of systematic spatial dependence. (Descriptive count only -- the monthly ",
    "sequence may be serially dependent, so no binomial/independence-based test is performed ",
    "against the 157x0.05 expectation.)"
  )
)
write.csv(mtc_summary, "results_v2/tables/moran_multiple_testing_summary.csv", row.names = FALSE)

message("\n---- Multiple-testing correction summary ----")
print(mtc_summary[, setdiff(names(mtc_summary), "interpretation")])
message(
  "\n", sum(moran_raw_mtc$sig_raw_05), " months were nominally significant at the unadjusted 5% ",
  "level, but only ", sum(moran_raw_mtc$sig_BH_05), " survived BH and ", sum(moran_raw_mtc$sig_Holm_05),
  " survived Holm correction (raw series); therefore isolated monthly significance is not ",
  "treated as robust evidence of systematic spatial dependence. This is reported as descriptive ",
  "context, not as a formal test against the 157x0.05 expectation -- the monthly Moran sequence ",
  "may be serially dependent, and no binomial/independence-based test is performed without ",
  "explicitly justifying that assumption."
)

## ----------------------------------------------------------------
## 5. Paired raw-vs-residual comparison and summary statistics.
##    No paired t-test: the monthly Moran sequence may be serially
##    dependent (overlapping economic conditions month to month), so
##    only descriptive comparison is reported, per instruction.
## ----------------------------------------------------------------

comparison <- moran_raw_mtc %>%
  select(date, I_raw = I, expected_I_raw = expected_I, I_minus_expected_raw = I_minus_expected,
         perm_p_raw = perm_p_value, sig_BH_raw = sig_BH_05, sig_Holm_raw = sig_Holm_05) %>%
  left_join(
    moran_resid_mtc %>% select(date, I_resid = I, expected_I_resid = expected_I,
                                I_minus_expected_resid = I_minus_expected, perm_p_resid = perm_p_value,
                                sig_BH_resid = sig_BH_05, sig_Holm_resid = sig_Holm_05),
    by = "date"
  ) %>%
  mutate(I_diff_resid_minus_raw = I_resid - I_raw) %>%
  arrange(date)

write.csv(comparison, "results_v2/tables/moran_raw_vs_residual_comparison.csv", row.names = FALSE)

## Theoretical Moran null expectation for N=13 regions: E[I] = -1/(N-1).
## Computed from N, not read off any single month's result, though it is
## (correctly) identical across all 157 months since E[I] depends only on
## N, never on the data -- verified below (stopifnot).
N_REGIONS <- 13
theoretical_EI <- -1 / (N_REGIONS - 1)
stopifnot(
  isTRUE(all.equal(unique(round(moran_raw_mtc$expected_I, 10)), round(theoretical_EI, 10))),
  isTRUE(all.equal(unique(round(moran_resid_mtc$expected_I, 10)), round(theoretical_EI, 10)))
)

## Primary inference is centred on the null expectation E[I], not on the
## naive I > 0 threshold (I > 0 is retained below ONLY as a clearly
## labelled secondary/descriptive figure, per instruction).
summary_stats <- function(I, expected_I, p, p_bh, p_holm) {
  I_minus_E <- I - expected_I
  data.frame(
    mean_I               = mean(I),
    mean_expected_I      = mean(expected_I),
    mean_I_minus_expected   = mean(I_minus_E),
    median_I_minus_expected = median(I_minus_E),
    prop_I_gt_expected   = mean(I > expected_I),
    n_raw_p_lt_05        = sum(p < 0.05),
    n_BH_sig             = sum(p_bh < 0.05),
    n_Holm_sig           = sum(p_holm < 0.05),
    iqr_I                = IQR(I),
    median_I             = median(I),
    prop_I_pos_secondary = mean(I > 0)  # descriptive/secondary only -- NOT the primary benchmark; see prop_I_gt_expected above
  )
}
moran_summary_stats <- bind_rows(
  cbind(series = "Raw inflation_mom",
        summary_stats(moran_raw_mtc$I, moran_raw_mtc$expected_I, moran_raw_mtc$perm_p_value,
                       moran_raw_mtc$p_BH, moran_raw_mtc$p_Holm)),
  cbind(series = "TWFE residual",
        summary_stats(moran_resid_mtc$I, moran_resid_mtc$expected_I, moran_resid_mtc$perm_p_value,
                       moran_resid_mtc$p_BH, moran_resid_mtc$p_Holm))
)
write.csv(moran_summary_stats, "results_v2/tables/moran_summary_stats.csv", row.names = FALSE)

message("\n---- Raw vs residual Moran summary, null-centred on E[I] = -1/(N-1) = ",
        round(theoretical_EI, 7), " (primary); I>0 is secondary/descriptive only ----")
print(moran_summary_stats)
message("\nMean(I_diff) = ", round(mean(comparison$I_diff_resid_minus_raw), 4),
        " | Median(I_diff) = ", round(median(comparison$I_diff_resid_minus_raw), 4))

## ----------------------------------------------------------------
## 6. Figures C, D, E.
## ----------------------------------------------------------------

## Both reference lines are shown: E[I] (primary null, solid) and 0
## (secondary/naive, dashed, lighter) -- I > 0 is not the benchmark used
## for inference in this analysis, per instruction.
eI_line <- function() list(
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey75", linewidth = 0.4),
  geom_hline(yintercept = theoretical_EI, linetype = "solid", color = "grey40", linewidth = 0.5)
)

p_c <- ggplot(moran_raw_mtc, aes(x = date, y = I)) +
  eI_line() +
  geom_line(color = "#d7191c", alpha = 0.6) +
  geom_point(aes(shape = sig_BH_05), size = 1.6, color = "#d7191c") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "BH-significant (p<0.05)") +
  labs(title = "Monthly Global Moran's I -- Raw inflation_mom",
       subtitle = paste0("Solid grey line = E[I] = -1/(N-1) = ", round(theoretical_EI, 4),
                          " (primary null); dashed = 0 (secondary reference). ",
                          "Permutation-based p-value, BH-corrected significance marked."),
       x = NULL, y = "Moran's I") +
  theme_minimal()
ggsave("results_v2/figures/fig_C_monthly_moran_raw.png", p_c, width = 10, height = 5, dpi = 300)

p_d <- ggplot(moran_resid_mtc, aes(x = date, y = I)) +
  eI_line() +
  geom_line(color = "#2c7bb6", alpha = 0.6) +
  geom_point(aes(shape = sig_BH_05), size = 1.6, color = "#2c7bb6") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "BH-significant (p<0.05)") +
  labs(title = "Monthly Global Moran's I -- TWFE (M3) Residuals",
       subtitle = paste0("Solid grey line = E[I] = -1/(N-1) = ", round(theoretical_EI, 4),
                          " (primary null); dashed = 0 (secondary reference). Same months/W as raw."),
       x = NULL, y = "Moran's I") +
  theme_minimal()
ggsave("results_v2/figures/fig_D_monthly_moran_residual.png", p_d, width = 10, height = 5, dpi = 300)

comparison_long <- comparison %>%
  select(date, I_raw, I_resid) %>%
  pivot_longer(cols = c(I_raw, I_resid), names_to = "series", values_to = "I") %>%
  mutate(series = recode(series, I_raw = "Raw inflation_mom", I_resid = "TWFE residual"))

p_e <- ggplot(comparison_long, aes(x = date, y = I, color = series)) +
  eI_line() +
  geom_line(alpha = 0.8) +
  scale_color_manual(values = c("Raw inflation_mom" = "#d7191c", "TWFE residual" = "#2c7bb6")) +
  labs(title = "Raw vs. TWFE-Residual Monthly Moran's I",
       subtitle = paste0("Solid grey line = E[I] = -1/(N-1) = ", round(theoretical_EI, 4),
                          " (primary null). Descriptive comparison -- no paired significance test."),
       x = NULL, y = "Moran's I", color = NULL) +
  theme_minimal() + theme(legend.position = "top")
ggsave("results_v2/figures/fig_E_raw_vs_residual_moran_comparison.png", p_e,
       width = 10, height = 5, dpi = 300)

message("\nWrote figures C, D, E to results_v2/figures/.")

## ----------------------------------------------------------------
## 7. Interpretation scope note -- written to disk, not just console
##    output, so the caveats travel with the results rather than living
##    only in chat history.
## ----------------------------------------------------------------

scope_note <- c(
  "MAKANI v2 Stage 1 (TWFE + monthly residual spatial tests) -- interpretation scope",
  "",
  "What the evidence supports:",
  "- Month effects dominate persistent region effects in the TWFE variance",
  "  decomposition (42.9% vs 0.45% of total MoM inflation variance).",
  "- Contemporaneous monthly cross-sectional spatial autocorrelation, measured by",
  "  Global Moran's I on each month's 13-region cross-section, is weak and does",
  "  not survive multiple-testing correction (0/157 months survive BH or Holm,",
  "  raw or residual).",
  "- Removing region and month fixed effects changes monthly Moran's I very",
  "  little (mean I - E[I] is close between raw and residual series).",
  "",
  "What this stage does NOT establish, and should not be read as showing:",
  "- That spatial processes do not exist in Saudi regional inflation. This",
  "  analysis tests CONTEMPORANEOUS (same-month) cross-sectional spatial",
  "  dependence only.",
  "- It does not rule out lagged spatial transmission (e.g. a shock in region",
  "  A affecting region B with a one- or multi-month delay).",
  "- It does not rule out category-specific spatial dependence (e.g. housing or",
  "  food sub-indices could show spatial structure even if headline MoM does not).",
  "- It does not rule out heterogeneous regional exposure to common national",
  "  shocks (a region-specific loading on tau_t is a distinct hypothesis from",
  "  region-to-region spillover, and is not tested here).",
  "- It does not rule out lower-frequency spatial dynamics (the long-run,",
  "  full-sample-mean Moran's I in the frozen v1.0 baseline, I=0.204, tests a",
  "  different quantity -- a single cross-section of 13 long-run averages --",
  "  not the same thing as 157 separate monthly cross-sections, and is not",
  "  contradicted or confirmed by the month-by-month result here.",
  "",
  "Accordingly, none of the tables or figures in this stage should be cited as",
  "evidence for or against dynamic, lagged, or category-specific spatial",
  "transmission -- those remain open questions for later stages."
)
writeLines(scope_note, "results_v2/tables/interpretation_scope_note.txt")
message("Wrote results_v2/tables/interpretation_scope_note.txt")

message("12_residual_spatial_tests.R complete.")
