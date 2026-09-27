# ================================
# Spatial autocorrelation (Global Moran's I)
#
# Purpose: test formally whether regional mean inflation is spatially
# clustered (neighbouring regions have similar inflation) rather than
# randomly distributed across space. A simple correlation matrix between
# regions (03_exploratory_analysis.R) cannot answer this question because
# it ignores geography; Moran's I explicitly uses the map's adjacency
# structure.
#
# Depends on `inflation_sf` (built in 04_spatial_analysis.R): regional
# boundaries joined with each region's mean inflation.
# ================================

library(spdep)
library(sf)

dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# =========================
# 1. مصفوفة الجوار المكاني (Queen contiguity: مشاركة حدّ أو حتى نقطة واحدة)
# =========================
nb <- poly2nb(inflation_sf, queen = TRUE)

n_no_neighbors <- sum(card(nb) == 0)
if (n_no_neighbors > 0) {
  warning(n_no_neighbors, " region(s) have no neighbours under queen contiguity; ",
          "they are excluded via zero.policy = TRUE.")
}

# =========================
# 2. أوزان مكانية موحّدة الصفوف (row-standardized)
# =========================
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# =========================
# 3. اختبار Global Moran's I على متوسط التضخم لكل منطقة
# =========================
moran_result <- moran.test(
  inflation_sf$mean_inflation,
  lw,
  zero.policy   = TRUE,
  na.action     = na.exclude
)

moran_summary <- data.frame(
  statistic = "Global Moran's I",
  I         = unname(moran_result$estimate["Moran I statistic"]),
  expected  = unname(moran_result$estimate["Expectation"]),
  variance  = unname(moran_result$estimate["Variance"]),
  z_value   = unname(moran_result$statistic),
  p_value   = moran_result$p.value
)
write.csv(moran_summary, "results/tables/moran_test.csv", row.names = FALSE)

message(
  "Moran's I = ", round(moran_summary$I, 3),
  " (p = ", signif(moran_summary$p_value, 3), ")"
)

# =========================
# 4. Moran scatter plot (القيمة مقابل متوسط جيرانها المكاني)
# =========================
png("results/figures/moran_scatter_plot.png", width = 1800, height = 1500, res = 300)
moran.plot(
  inflation_sf$mean_inflation,
  lw,
  zero.policy = TRUE,
  xlab        = "Mean Inflation",
  ylab        = "Spatially Lagged Mean Inflation",
  main        = "Moran Scatter Plot: Regional Inflation",
  labels      = inflation_sf$name
)
dev.off()

# =========================
# 5. اختبار محاكاة مونت كارلو (Monte Carlo permutation test)
#
# مع عيّنة صغيرة جداً (n = 13 منطقة)، قيمة p الناتجة من moran.test() تعتمد
# على تقارب توزيعي (asymptotic) قد لا يكون دقيقاً بهذا الحجم. moran.mc()
# تعطي قيمة p قائمة على المحاكاة (9999 إعادة ترتيب عشوائي) وهي أكثر
# موثوقية هنا، وتُستخدم كفحص متانة (robustness check) لنتيجة moran.test().
# =========================
set.seed(123)
moran_mc_result <- moran.mc(
  inflation_sf$mean_inflation,
  lw,
  nsim        = 9999,
  zero.policy = TRUE
)

moran_mc_summary <- data.frame(
  statistic     = "Global Moran's I (Monte Carlo, 9999 permutations)",
  I             = unname(moran_mc_result$statistic),
  n_simulations = 9999,
  p_value       = moran_mc_result$p.value
)
write.csv(moran_mc_summary, "results/tables/moran_mc_test.csv", row.names = FALSE)

message(
  "Moran's I (Monte Carlo) = ", round(moran_mc_summary$I, 3),
  " (p = ", signif(moran_mc_summary$p_value, 3), ")"
)
