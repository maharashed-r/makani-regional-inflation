# ================================
# Data cleaning / quality checks
# (Runs after 01_load_data.R, which loads `data_clean`)
# ================================

library(dplyr)

# 1. إزالة الصفوف المكررة تماماً (إن وجدت)
n_before <- nrow(data_clean)
data_clean <- data_clean %>% distinct()
n_dupes <- n_before - nrow(data_clean)
if (n_dupes > 0) {
  message(n_dupes, " duplicate row(s) removed.")
}

# 2. تحويل عمود التاريخ إلى Date فعلي (بدل نص)
data_clean$date <- as.Date(data_clean$date)

# 3. التحقق أن كل منطقة لديها نفس عدد الأشهر (لوحة بيانات متوازنة)
panel_check <- data_clean %>%
  count(Region, name = "n_months")

if (length(unique(panel_check$n_months)) > 1) {
  warning(
    "Unbalanced panel: regions do not all have the same number of observations:\n",
    paste(capture.output(print(panel_check)), collapse = "\n")
  )
}

# 4. التحقق من عدم وجود قيم تضخم غير منطقية (NA أو Inf)
n_invalid <- sum(!is.finite(data_clean$inflation))
if (n_invalid > 0) {
  warning(n_invalid, " non-finite inflation value(s) found (NA/Inf).")
}

# ملاحظة: حساب الإحصاءات الوصفية (region_statistics.csv) ومصفوفة الارتباط
# يتم في 03_exploratory_analysis.R فقط، تجنباً للتكرار.
