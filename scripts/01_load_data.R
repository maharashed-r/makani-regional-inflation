library(dplyr)
library(tidyr)

# التأكد من وجود مجلدات النتائج قبل أي كتابة لاحقة في السكربتات التالية
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# ================================================================
# بناء data_clean مباشرة من data_raw في كل تشغيل (بدل قراءة نسخة
# مخزّنة قديمة).
#
# السبب: النسخة المخزّنة سابقاً في data_clean/saudi_cpi_clean.csv كان
# عمود "inflation" فيها محسوباً بترتيب عمود الشهر (M) كنص (نصي/
# lexicographic) بدل رقمي. بما أن data_raw يحتوي أيضاً على صف "متوسط
# سنوي" لكل سنة برمز M == "-"، أصبح ترتيب كل سنة فعلياً:
#   "-", "1", "10", "11", "12", "2", "3", "4", "5", "6", "7", "8", "9"
# بدل الترتيب الزمني الصحيح (1..12). هذا جعل حساب الفرق اللوغاريتمي
# (log-diff) بين "الشهر الحالي" و"الشهر السابق حسب هذا الترتيب الخاطئ"
# ينتج قيم تضخم غير صحيحة لثلاثة أشهر من كل سنة تحديداً:
#   - يناير: يقارن بمتوسط السنة نفسها بدل ديسمبر من السنة السابقة
#   - فبراير: يقارن بديسمبر من نفس السنة (شهر لم يأتِ بعد) بدل يناير
#   - أكتوبر: يقارن بيناير من نفس السنة بدل سبتمبر
# تم التحقق من هذا رقمياً: القيم المخزّنة سابقاً طابقت هذه الحسابات
# الخاطئة بدقة أقل من 0.0000005، عبر عدة مناطق وسنوات مختلفة، ما يؤكد
# السبب بشكل قاطع. الأثر: 494 من أصل 2028 قيمة تضخم (24%) كانت خاطئة.
#
# الحل: إعادة بناء data_clean من data_raw مباشرة مع:
#   1. استبعاد صف المتوسط السنوي (M == "-") قبل أي ترتيب
#   2. تحويل M و Y إلى أرقام صحيحة (integer) بدل نص
#   3. الترتيب الزمني الصريح (arrange) قبل حساب log-diff لكل منطقة
# ================================================================

raw <- read.csv("data_raw/saudi_cpi_regions.csv", check.names = FALSE)

raw <- raw[raw$M != "-", ]          # استبعاد صف المتوسط السنوي
raw$Y <- as.integer(raw$Y)
raw$M <- as.integer(raw$M)          # تحويل الشهر إلى رقم صحيح (لا نص)

region_cols <- setdiff(names(raw), c("Y", "M"))

data_clean <- raw %>%
  pivot_longer(cols = all_of(region_cols), names_to = "Region", values_to = "CPI") %>%
  mutate(date = as.Date(sprintf("%04d-%02d-01", Y, M))) %>%
  arrange(Region, date) %>%         # ترتيب زمني صريح وصحيح لكل منطقة
  group_by(Region) %>%
  mutate(inflation = log(CPI) - log(lag(CPI))) %>%
  ungroup() %>%
  filter(!is.na(inflation)) %>%     # إسقاط أول شهر لكل منطقة (لا يوجد شهر سابق)
  select(date, Region, CPI, inflation) %>%
  arrange(Region, date)

# حفظ نسخة محدّثة وصحيحة من data_clean (تحل محل النسخة القديمة الخاطئة)
write.csv(data_clean, "data_clean/saudi_cpi_clean.csv", row.names = FALSE)

if (any(is.na(data_clean))) {
  message("Missing values detected in data_clean: ",
          sum(is.na(data_clean)), " cell(s). See 02_clean_data.R for the check.")
}