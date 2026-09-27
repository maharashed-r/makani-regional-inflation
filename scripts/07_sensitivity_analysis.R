# ================================
# Sensitivity & Robustness Analysis (reviewer-requested)
#
# Addresses reviewer points on the spatial analysis:
#   1. Are the correlation / Moran's I results driven by the two VAT policy
#      shock months (Jan 2018 introduction, Jul 2020 rate tripling)? ->
#      Sensitivity Analysis table, now testing not just the single shock
#      month but +/-1 month, the calendar quarter, and a six-month window
#      starting at each shock (a reviewer specifically asked why only one
#      month was tested, since VAT effects can persist longer).
#   2. Is the Global Moran's I result sensitive to the choice of spatial
#      weights matrix (Queen was used throughout)? -> Weights Sensitivity
#      table: Queen / Rook / KNN (k=4) / Distance band.
#   3. Global Moran's I only tells us whether spatial clustering exists
#      overall, not WHERE. -> Local Moran's I (LISA): classifies each
#      region into High-High / Low-Low / High-Low / Low-High / Not
#      significant, with a cluster map.
#   4. Correlation significance across 78 pairs inflates the false-positive
#      risk (multiple testing). -> Benjamini-Hochberg (FDR) correction
#      applied in addition to the raw p-values.
#   5. Global/Local Moran's I are both computed on each region's mean
#      inflation across the whole sample, collapsing the time dimension.
#      -> exploratory year-by-year Moran's I to show how spatial clustering
#      has evolved over time.
#
# Depends on: data_clean (01_load_data.R), inflation_sf, lw
# (04_spatial_analysis.R, 06_spatial_autocorrelation.R).
# ================================

library(dplyr)
library(tidyr)
library(spdep)
library(sf)
library(ggplot2)

dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Same Region -> map-name recoding used in 04_spatial_analysis.R, kept
# identical here so scenario region means join onto inflation_sf$name
# exactly as the main pipeline does.
recode_region_names <- function(region) {
  recode(region,
    "Aseer"             = "`Asir",
    "Jazan"             = "Jizan",
    "Riyadh"            = "Ar Riyad",
    "Madinah"           = "Al Madinah",
    "Al Qassim"         = "Al Quassim",
    "Eastern Province"  = "Ash Sharqiyah",
    "Jouf"              = "Al Jawf",
    "Al Jouf"           = "Al Jawf",
    "Al Baha"           = "Al Bahah",
    "Baha"              = "Al Bahah",
    "Tabouk"            = "Tabuk",
    "Northern Borders"  = "Al Hudud ash Shamaliyah"
  )
}

run_scenario <- function(scenario_name, excl_dates) {
  scenario_data <- data_clean %>% filter(!date %in% excl_dates)

  # --- correlation side ---
  wide_mat <- scenario_data %>%
    select(date, Region, inflation) %>%
    pivot_wider(names_from = Region, values_from = inflation) %>%
    select(-date) %>%
    as.matrix()

  cor_mat <- cor(wide_mat, use = "pairwise.complete.obs")
  diag(cor_mat) <- NA
  avg_cor <- mean(cor_mat, na.rm = TRUE)
  max_cor <- max(cor_mat, na.rm = TRUE)

  # --- Moran's I side: recompute regional means under this scenario and
  #     rejoin to the (fixed) map geometry, using the same weights `lw` ---
  region_means <- scenario_data %>%
    group_by(Region) %>%
    summarise(mean_inflation = mean(inflation, na.rm = TRUE), .groups = "drop") %>%
    mutate(Region = recode_region_names(Region))

  sf_scn <- inflation_sf %>%
    select(-mean_inflation) %>%
    left_join(region_means, by = c("name" = "Region"))

  moran_scn <- moran.test(
    sf_scn$mean_inflation, lw,
    zero.policy = TRUE, na.action = na.exclude
  )

  data.frame(
    Scenario             = scenario_name,
    N_months_excluded    = length(excl_dates),
    Average_Correlation  = round(avg_cor, 3),
    Maximum_Correlation  = round(max_cor, 3),
    Moran_I              = round(unname(moran_scn$estimate["Moran I statistic"]), 3),
    Moran_p_value        = signif(moran_scn$p.value, 3)
  )
}

# ============================================================
# 1. Sensitivity Analysis: VAT shock windows of increasing width
# ============================================================
jan2018 <- as.Date("2018-01-01")
jul2020 <- as.Date("2020-07-01")

# `before`/`after` are in whole months relative to the shock month.
# quarter  = the 3-month calendar quarter starting at the shock month
# 6 months = the 6-month window starting at the shock month
shift_month <- function(date, n) seq(date, by = paste(n, "months"), length.out = 2)[2]
window_dates <- function(center, before, after) {
  seq(shift_month(center, -before), shift_month(center, after), by = "month")
}

vat_scenarios <- list(
  "Full Sample"                       = as.Date(character(0)),
  "Excl. Jan 2018 (month)"            = window_dates(jan2018, 0, 0),
  "Excl. Jul 2020 (month)"            = window_dates(jul2020, 0, 0),
  "Excl. Both (month)"                = c(window_dates(jan2018, 0, 0), window_dates(jul2020, 0, 0)),
  "Excl. Jan 2018 (+/-1 month)"       = window_dates(jan2018, 1, 1),
  "Excl. Jul 2020 (+/-1 month)"       = window_dates(jul2020, 1, 1),
  "Excl. Both (+/-1 month)"           = c(window_dates(jan2018, 1, 1), window_dates(jul2020, 1, 1)),
  "Excl. Jan 2018 (quarter)"          = window_dates(jan2018, 0, 2),
  "Excl. Jul 2020 (quarter)"          = window_dates(jul2020, 0, 2),
  "Excl. Both (quarter)"              = c(window_dates(jan2018, 0, 2), window_dates(jul2020, 0, 2)),
  "Excl. Jan 2018 (6 months)"         = window_dates(jan2018, 0, 5),
  "Excl. Jul 2020 (6 months)"         = window_dates(jul2020, 0, 5),
  "Excl. Both (6 months)"             = c(window_dates(jan2018, 0, 5), window_dates(jul2020, 0, 5))
)

sensitivity_analysis <- bind_rows(
  lapply(names(vat_scenarios), function(nm) run_scenario(nm, vat_scenarios[[nm]]))
)
write.csv(sensitivity_analysis, "results/tables/sensitivity_analysis.csv", row.names = FALSE)
message("Sensitivity analysis (VAT shock windows) written.")
print(sensitivity_analysis)

# ============================================================
# 2. Spatial weights matrix sensitivity: Queen vs Rook vs KNN vs Distance band
# ============================================================
coords <- st_coordinates(st_centroid(st_geometry(inflation_sf)))

nb_queen <- poly2nb(inflation_sf, queen = TRUE)
nb_rook  <- poly2nb(inflation_sf, queen = FALSE)
nb_knn4  <- knn2nb(knearneigh(coords, k = 4))

# Distance band: smallest threshold that still gives every region >= 1
# neighbour (based on each region's nearest-neighbour distance), plus a
# small buffer.
nb_knn1   <- knn2nb(knearneigh(coords, k = 1))
band_dist <- max(unlist(nbdists(nb_knn1, coords))) * 1.05
nb_dist   <- dnearneigh(coords, 0, band_dist)

weights_list <- list(
  "Queen contiguity"        = nb_queen,
  "Rook contiguity"         = nb_rook,
  "K-nearest neighbours (k=4)" = nb_knn4,
  "Distance band"           = nb_dist
)

weights_rows <- lapply(names(weights_list), function(nm) {
  nb_scn <- weights_list[[nm]]
  lw_scn <- nb2listw(nb_scn, style = "W", zero.policy = TRUE)
  mt <- moran.test(
    inflation_sf$mean_inflation, lw_scn,
    zero.policy = TRUE, na.action = na.exclude
  )
  data.frame(
    Weight_Matrix  = nm,
    Avg_Neighbors  = round(mean(card(nb_scn)), 2),
    Moran_I        = round(unname(mt$estimate["Moran I statistic"]), 3),
    Z_value        = round(unname(mt$statistic), 2),
    P_value        = signif(mt$p.value, 3)
  )
})

weights_sensitivity <- bind_rows(weights_rows)
write.csv(weights_sensitivity, "results/tables/weights_sensitivity.csv", row.names = FALSE)
message("Spatial weights sensitivity analysis written.")
print(weights_sensitivity)

# ============================================================
# 3. Local Moran's I (LISA): where is the clustering, not just whether it exists
# ============================================================
z_val   <- as.numeric(scale(inflation_sf$mean_inflation))
lag_z   <- lag.listw(lw, z_val, zero.policy = TRUE)

local_moran <- localmoran(
  inflation_sf$mean_inflation, lw,
  zero.policy = TRUE, na.action = na.exclude
)

p_col <- grep("^Pr", colnames(local_moran), value = TRUE)[1]
lisa_p <- local_moran[, p_col]

quadrant <- dplyr::case_when(
  is.na(lisa_p) | lisa_p >= 0.05 ~ "Not significant",
  z_val > 0 & lag_z > 0 ~ "High-High",
  z_val < 0 & lag_z < 0 ~ "Low-Low",
  z_val > 0 & lag_z < 0 ~ "High-Low",
  z_val < 0 & lag_z > 0 ~ "Low-High",
  TRUE ~ "Not significant"
)

lisa_results <- data.frame(
  Region        = inflation_sf$name,
  Mean_Inflation = round(inflation_sf$mean_inflation, 4),
  Local_I       = round(local_moran[, "Ii"], 3),
  Z_value       = round(local_moran[, "Z.Ii"], 2),
  P_value       = signif(lisa_p, 3),
  Cluster       = quadrant
)
write.csv(lisa_results, "results/tables/lisa_results.csv", row.names = FALSE)
message("LISA (Local Moran's I) results written.")
print(lisa_results)

# LISA cluster map
inflation_sf_lisa <- inflation_sf %>% mutate(lisa_cluster = quadrant)

lisa_colors <- c(
  "High-High"       = "#d7191c",
  "Low-Low"         = "#2c7bb6",
  "High-Low"        = "#fdae61",
  "Low-High"        = "#abd9e9",
  "Not significant" = "grey90"
)

p_lisa <- ggplot(inflation_sf_lisa) +
  geom_sf(aes(fill = lisa_cluster), color = "white", linewidth = 0.3) +
  scale_fill_manual(values = lisa_colors, name = "LISA cluster (p < 0.05)") +
  labs(
    title = "Local Moran's I (LISA) Cluster Map: Regional Mean Inflation",
    subtitle = "Queen contiguity weights"
  ) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave("results/figures/lisa_cluster_map.png", p_lisa, width = 7, height = 6, dpi = 300)
message("LISA cluster map saved.")

# ============================================================
# 4. Correlation significance testing (+ Benjamini-Hochberg FDR correction)
#    Which of the 78 regional pairs are NOT significantly correlated, and
#    how many survive correction for testing 78 pairs simultaneously?
# ============================================================
full_wide <- data_clean %>%
  select(date, Region, inflation) %>%
  pivot_wider(names_from = Region, values_from = inflation)

region_names_all <- setdiff(names(full_wide), "date")
region_pairs <- combn(region_names_all, 2, simplify = FALSE)

corr_sig_rows <- lapply(region_pairs, function(p) {
  x  <- full_wide[[p[1]]]
  y  <- full_wide[[p[2]]]
  ok <- stats::complete.cases(x, y)
  test <- stats::cor.test(x[ok], y[ok])
  data.frame(
    Region_1             = p[1],
    Region_2             = p[2],
    Correlation          = round(unname(test$estimate), 3),
    P_value              = signif(test$p.value, 3),
    Significant_at_5pct  = test$p.value < 0.05
  )
})

correlation_significance <- bind_rows(corr_sig_rows) %>% arrange(P_value)

# Benjamini-Hochberg (FDR) correction across all 78 simultaneous tests
correlation_significance$P_value_BH <- signif(
  p.adjust(correlation_significance$P_value, method = "BH"), 3
)
correlation_significance$Significant_BH_5pct <- correlation_significance$P_value_BH < 0.05

write.csv(correlation_significance, "results/tables/correlation_significance.csv", row.names = FALSE)

n_pairs_total    <- nrow(correlation_significance)
n_pairs_sig      <- sum(correlation_significance$Significant_at_5pct)
n_pairs_sig_bh   <- sum(correlation_significance$Significant_BH_5pct)
message(
  "Correlation significance: ", n_pairs_sig, "/", n_pairs_total,
  " region pairs significant at raw 5%; ", n_pairs_sig_bh, "/", n_pairs_total,
  " significant after Benjamini-Hochberg correction."
)
print(correlation_significance %>% filter(!Significant_at_5pct))

# ============================================================
# 5. Exploratory: Global Moran's I year by year
#    Global/Local Moran's I above collapse the whole sample period into a
#    single mean per region. This section computes Moran's I separately
#    for each calendar year's mean inflation, to show (exploratorily,
#    without claiming a formal time-varying model) how spatial clustering
#    has moved over time.
# ============================================================
data_clean_yr <- data_clean %>% mutate(year = as.integer(format(date, "%Y")))
years <- sort(unique(data_clean_yr$year))

moran_by_year <- bind_rows(lapply(years, function(yr) {
  yr_means <- data_clean_yr %>%
    filter(year == yr) %>%
    group_by(Region) %>%
    summarise(mean_inflation = mean(inflation, na.rm = TRUE), .groups = "drop") %>%
    mutate(Region = recode_region_names(Region))

  sf_yr <- inflation_sf %>%
    select(-mean_inflation) %>%
    left_join(yr_means, by = c("name" = "Region"))

  n_valid <- sum(!is.na(sf_yr$mean_inflation))
  if (n_valid < 5) {
    return(data.frame(Year = yr, N_regions = n_valid, Moran_I = NA_real_, P_value = NA_real_))
  }

  mt <- tryCatch(
    moran.test(sf_yr$mean_inflation, lw, zero.policy = TRUE, na.action = na.exclude),
    error = function(e) NULL
  )
  if (is.null(mt)) {
    return(data.frame(Year = yr, N_regions = n_valid, Moran_I = NA_real_, P_value = NA_real_))
  }
  data.frame(
    Year      = yr,
    N_regions = n_valid,
    Moran_I   = round(unname(mt$estimate["Moran I statistic"]), 3),
    P_value   = signif(mt$p.value, 3)
  )
}))

write.csv(moran_by_year, "results/tables/moran_by_year.csv", row.names = FALSE)
message("Year-by-year Moran's I (exploratory) written.")
print(moran_by_year)

p_moran_time <- ggplot(moran_by_year, aes(x = Year, y = Moran_I)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(color = "steelblue") +
  geom_point(aes(color = P_value < 0.05), size = 2.5) +
  scale_color_manual(
    values = c(`TRUE` = "#d7191c", `FALSE` = "grey40"),
    name = "p < 0.05", na.value = "grey80"
  ) +
  labs(
    title = "Global Moran's I by Year (Exploratory)",
    subtitle = "Each point uses that calendar year's regional mean inflation only",
    x = "Year", y = "Moran's I"
  ) +
  theme_minimal()

ggsave("results/figures/moran_over_time.png", p_moran_time, width = 8, height = 5, dpi = 300)
message("Moran's I over time plot saved.")
