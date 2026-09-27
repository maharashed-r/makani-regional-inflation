# ================================
# Final plots for the paper
# Single source of truth for all three figures (histogram, correlation
# heatmap, map). Depends on objects built earlier in the pipeline:
#   - data_clean   (01_load_data.R)
#   - cor_melt     (03_exploratory_analysis.R)
#   - inflation_sf (04_spatial_analysis.R)
# ================================

library(ggplot2)
library(sf)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Histogram --------------------------------------------------------------
p_hist <- ggplot(data_clean, aes(x = inflation)) +
  geom_histogram(fill = "steelblue", bins = 30) +
  labs(title = "Distribution of Inflation")

# Map ----------------------------------------------------------------------
p_map <- ggplot(data = inflation_sf) +
  geom_sf(aes(fill = mean_inflation), color = "white", size = 0.3) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "Inflation",
    na.value = "grey90"
  ) +
  labs(
    title = "Regional Inflation in Saudi Arabia",
    subtitle = "Average Inflation by Region"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

# Correlation heatmap --------------------------------------------------------
# midpoint = 0 (not an arbitrary value) so the diverging scale is centered on
# "no correlation"; colors above/below 0 are then directly comparable.
p_corr <- ggplot(cor_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation"
  ) +
  theme_minimal() +
  labs(title = "Correlation Between Regions") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 14, face = "bold")
  )

# حفظ الرسوم
ggsave("results/figures/histogram.png",   p_hist, width = 6, height = 4)
ggsave("results/figures/map.png",         p_map,  width = 7, height = 5)
ggsave("results/figures/correlation.png", p_corr, width = 6, height = 5)
