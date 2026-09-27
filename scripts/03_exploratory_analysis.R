# ================================
# Exploratory analysis
# Regional inflation statistics + correlation structure
# (Runs after 01_load_data.R / 02_clean_data.R, which prepare `data_clean`)
#
# NOTE: this script only computes data/tables. Plotting is centralized in
# 05_plots.R so each figure has a single source of truth.
# ================================

library(dplyr)
library(tidyr)
library(reshape2)

# تحويل البيانات إلى wide لحساب الارتباط بين المناطق
inflation_wide <- data_clean %>%
  select(date, Region, inflation) %>%
  pivot_wider(names_from = Region, values_from = inflation)

inflation_matrix <- inflation_wide %>%
  select(-date)

# مصفوفة الارتباط (pairwise.complete.obs: أي شهر ناقص في منطقة لا يُسقط بقية
# الأزواج، فقط يُستبعد من ذلك الزوج بالتحديد)
cor_matrix <- cor(inflation_matrix, use = "pairwise.complete.obs")
cor_melt   <- melt(cor_matrix)

# إحصاءات وصفية لكل منطقة (المصدر الوحيد لهذا الجدول في المشروع)
region_stats <- data_clean %>%
  group_by(Region) %>%
  summarise(
    mean_inflation = mean(inflation, na.rm = TRUE),
    sd_inflation   = sd(inflation, na.rm = TRUE),
    min_inflation  = min(inflation, na.rm = TRUE),
    max_inflation  = max(inflation, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_inflation))

# حفظ النتائج
write.csv(region_stats, "results/tables/region_statistics.csv", row.names = FALSE)
write.csv(cor_matrix, "results/tables/regional_inflation_correlation.csv")
