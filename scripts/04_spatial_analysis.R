# ================================
# Spatial analysis
# Builds `inflation_sf`: regional boundaries joined with mean inflation.
# (Runs after 01-03, which prepare `data_clean`.)
#
# NOTE: this script only builds the spatial data object. The actual map is
# drawn once, in 05_plots.R, so there is a single source of truth for the
# figure (previously this script AND 05_plots.R both saved map.png).
# ================================

library(dplyr)
library(sf)
library(rnaturalearth)

# =========================
# 1. متوسط التضخم لكل منطقة
# =========================
inflation_mean <- data_clean %>%
  group_by(Region) %>%
  summarise(mean_inflation = mean(inflation, na.rm = TRUE))

# =========================
# 2. ضبط أسماء المناطق لتطابق أسماء الخريطة (قبل الدمج)
# =========================
inflation_mean$Region <- recode(inflation_mean$Region,
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

# =========================
# 3. جلب حدود المناطق (مع تخزين محلي لتفادي الاعتماد على الإنترنت في كل تشغيل)
# =========================
shapefile_cache <- "data_raw/saudi_states_sf.rds"

if (file.exists(shapefile_cache)) {
  inflation_sf <- readRDS(shapefile_cache)
} else {
  inflation_sf <- ne_states(country = "Saudi Arabia", returnclass = "sf")
  saveRDS(inflation_sf, shapefile_cache)
}

# تعديل اسم حائل داخل الخريطة
inflation_sf$name <- recode(inflation_sf$name, "Ha'il" = "Hail")

# =========================
# 4. الدمج
# =========================
inflation_sf <- inflation_sf %>%
  left_join(inflation_mean, by = c("name" = "Region"))

# =========================
# 5. تحقق صريح من نجاح المطابقة (بدل الطباعة الصامتة فقط)
# =========================
unmatched <- inflation_sf$name[is.na(inflation_sf$mean_inflation)]
if (length(unmatched) > 0) {
  warning(
    "The following regions on the map did not match any region in the data ",
    "and will show as NA: ", paste(unmatched, collapse = ", ")
  )
}

unused_regions <- setdiff(inflation_mean$Region, inflation_sf$name)
if (length(unused_regions) > 0) {
  warning(
    "The following regions in the data did not match any region on the map: ",
    paste(unused_regions, collapse = ", ")
  )
}
