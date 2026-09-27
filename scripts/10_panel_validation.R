# ================================================================
# MAKANI v2 -- Panel validation report
#
# NEW v2 FILE. Reads only:
#   - data_v2/region_crosswalk.csv, data_v2/saudi_cpi_panel_v2.csv
#     (built by scripts/09_panel_build.R)
#   - data_raw/saudi_cpi_regions.csv, data_raw/saudi_states_sf.rds
#     (re-read independently, as ground truth, not assumed identical to
#     what 09 produced)
#
# Does not modify any file under data_clean/, results/, baseline/, or
# scripts/01-08. Writes only data_v2/panel_validation_report.txt.
#
# This script produces a human-readable validation report. The codified,
# CI-style pass/fail assertions live separately in
# tests/testthat/test-panel-v2.R (run via testthat, not this script).
#
# Run manually: Rscript scripts/10_panel_validation.R
# ================================================================

dir.create("data_v2", recursive = TRUE, showWarnings = FALSE)

panel     <- read.csv("data_v2/saudi_cpi_panel_v2.csv", stringsAsFactors = FALSE)
crosswalk <- read.csv("data_v2/region_crosswalk.csv",   stringsAsFactors = FALSE)
panel$date <- as.Date(panel$date)

raw_header       <- read.csv("data_raw/saudi_cpi_regions.csv", check.names = FALSE, nrows = 0)
raw_region_names <- setdiff(names(raw_header), c("Y", "M"))
geom_raw         <- readRDS("data_raw/saudi_states_sf.rds")
geom_df          <- as.data.frame(geom_raw)[, c("name", "name_ar", "iso_3166_2")]
registry         <- read.csv("data_v2/region_id_registry.csv", stringsAsFactors = FALSE, na.strings = "")

lines <- character(0)
add   <- function(...) lines <<- c(lines, sprintf(...))
ok_fail <- function(cond) if (isTRUE(cond)) "PASS" else "FAIL"

add("MAKANI v2 Panel Validation Report")
add("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
add(strrep("=", 70))

## ---- A. Panel structure -----------------------------------------------
add("")
add("A. PANEL STRUCTURE")

n_regions_panel <- length(unique(panel$region_id))
chk_a1 <- n_regions_panel == 13
add("  [%s] exactly 13 regions in panel (found %d)", ok_fail(chk_a1), n_regions_panel)

chk_a2 <- length(unique(crosswalk$region_id)) == 13 && nrow(crosswalk) == 13
add("  [%s] exactly 13 regions in crosswalk (found %d)", ok_fail(chk_a2), nrow(crosswalk))

dup_key <- paste(panel$region_id, panel$date)
n_dupe_keys <- sum(duplicated(dup_key))
chk_a3 <- n_dupe_keys == 0
add("  [%s] unique region_id x date (duplicate rows found: %d)", ok_fail(chk_a3), n_dupe_keys)

balance <- aggregate(date ~ region_id, data = panel, FUN = length)
names(balance)[2] <- "n_months"
first_last <- aggregate(date ~ region_id, data = panel,
                         FUN = function(x) c(min = format(min(x)), max = format(max(x))))
chk_a4 <- length(unique(balance$n_months)) == 1
add("  [%s] balanced panel (all regions have the same month count: %s)",
    ok_fail(chk_a4), paste(unique(balance$n_months), collapse = ","))

# Expected monthly date sequence: every region's dates must be a
# contiguous run of calendar months with no gaps.
expected_seq <- seq(min(panel$date), max(panel$date), by = "month")
gap_check <- sapply(split(panel$date, panel$region_id), function(d) {
  identical(sort(d), expected_seq)
})
chk_a5 <- all(gap_check)
add("  [%s] every region's date sequence matches the expected contiguous monthly run (%d/%d regions)",
    ok_fail(chk_a5), sum(gap_check), length(gap_check))

add("")
add("  Region balance table (months per region_id):")
for (i in seq_len(nrow(balance))) {
  add("    %-6s %d months", balance$region_id[i], balance$n_months[i])
}

add("")
add("  First / last date by region:")
fl <- do.call(rbind, lapply(split(panel$date, panel$region_id), function(d) {
  data.frame(first = format(min(d)), last = format(max(d)))
}))
fl$region_id <- rownames(fl)
for (i in seq_len(nrow(fl))) {
  add("    %-6s first=%s  last=%s", fl$region_id[i], fl$first[i], fl$last[i])
}

## ---- B. CPI integrity ---------------------------------------------------
add("")
add("B. CPI INTEGRITY")

chk_b1 <- is.numeric(panel$CPI)
add("  [%s] CPI column is numeric", ok_fail(chk_b1))

chk_b2 <- all(is.finite(panel$CPI)) && all(panel$CPI > 0)
add("  [%s] CPI strictly positive and finite for all %d rows (violations: %d)",
    ok_fail(chk_b2), nrow(panel), sum(!is.finite(panel$CPI) | panel$CPI <= 0))

chk_b3 <- all(panel$month >= 1 & panel$month <= 12) && all(!is.na(panel$date))
add("  [%s] no impossible dates (month in 1-12, all dates parse)", ok_fail(chk_b3))

raw_full <- read.csv("data_raw/saudi_cpi_regions.csv", check.names = FALSE)
n_annual_rows_in_source <- sum(raw_full$M == "-")
raw_full_no_annual <- raw_full[raw_full$M != "-", ]
chk_b4 <- nrow(panel) == 13 * nrow(raw_full_no_annual)
add("  [%s] no accidental annual-average rows in panel (%d annual rows existed in source and were excluded; panel row count %d == 13 x %d months)",
    ok_fail(chk_b4), n_annual_rows_in_source, nrow(panel), nrow(raw_full_no_annual))

# No text-sort chronology problem: within each region, the date column
# must already be monotonically increasing as stored (i.e. numeric Y/M
# sort, not lexicographic "1,10,11,12,2,...").
chrono_ok <- sapply(split(panel$date, panel$region_id), function(d) !is.unsorted(d))
chk_b5 <- all(chrono_ok)
add("  [%s] no text-sort chronology problem (dates monotonically increasing within every region)",
    ok_fail(chk_b5))

chk_b6 <- !any(is.na(panel$year)) && !any(is.na(panel$month)) && !any(is.na(panel$CPI))
add("  [%s] no silent coercion failures (no NA in year/month/CPI)", ok_fail(chk_b6))

## ---- C. Inflation integrity ---------------------------------------------
add("")
add("C. INFLATION INTEGRITY")

set.seed(42)
sample_idx <- sample(which(!is.na(panel$inflation_mom)), 5)
mom_ok <- logical(length(sample_idx))
for (i in seq_along(sample_idx)) {
  r <- panel[sample_idx[i], ]
  prev <- panel[panel$region_id == r$region_id &
                 panel$date == (seq(r$date, by = "-1 month", length.out = 2)[2]), "CPI"]
  recomputed <- 100 * (log(r$CPI) - log(prev))
  mom_ok[i] <- isTRUE(all.equal(recomputed, r$inflation_mom, tolerance = 1e-9))
  add("    sample %d: region=%s date=%s recomputed=%.10f stored=%.10f match=%s",
      i, r$region_id, r$date, recomputed, r$inflation_mom, mom_ok[i])
}
chk_c1 <- all(mom_ok)
add("  [%s] manual MoM recompute matches stored value on %d/%d sampled rows",
    ok_fail(chk_c1), sum(mom_ok), length(mom_ok))

sample_idx_yoy <- sample(which(!is.na(panel$inflation_yoy)), 5)
yoy_ok <- logical(length(sample_idx_yoy))
for (i in seq_along(sample_idx_yoy)) {
  r <- panel[sample_idx_yoy[i], ]
  prev <- panel[panel$region_id == r$region_id &
                 panel$date == (seq(r$date, by = "-12 month", length.out = 2)[2]), "CPI"]
  recomputed <- 100 * (log(r$CPI) - log(prev))
  yoy_ok[i] <- isTRUE(all.equal(recomputed, r$inflation_yoy, tolerance = 1e-9))
  add("    sample %d: region=%s date=%s recomputed=%.10f stored=%.10f match=%s",
      i, r$region_id, r$date, recomputed, r$inflation_yoy, yoy_ok[i])
}
chk_c2 <- all(yoy_ok)
add("  [%s] manual YoY recompute matches stored value on %d/%d sampled rows",
    ok_fail(chk_c2), sum(yoy_ok), length(yoy_ok))

first_mom_na <- sapply(split(panel, panel$region_id), function(d) {
  d <- d[order(d$date), ]
  is.na(d$inflation_mom[1]) && all(!is.na(d$inflation_mom[-1]))
})
chk_c3 <- all(first_mom_na)
add("  [%s] first MoM observation is structurally NA in every region, and no other MoM value is NA (%d/%d regions)",
    ok_fail(chk_c3), sum(first_mom_na), length(first_mom_na))

first12_yoy_na <- sapply(split(panel, panel$region_id), function(d) {
  d <- d[order(d$date), ]
  all(is.na(d$inflation_yoy[1:12])) && all(!is.na(d$inflation_yoy[-(1:12)]))
})
chk_c4 <- all(first12_yoy_na)
add("  [%s] first 12 YoY observations are structurally NA in every region, and no other YoY value is NA (%d/%d regions)",
    ok_fail(chk_c4), sum(first12_yoy_na), length(first12_yoy_na))

## ---- D. Spatial key integrity --------------------------------------------
add("")
add("D. SPATIAL KEY INTEGRITY")

chk_d1 <- all(crosswalk$region_name_gastat_en %in% raw_region_names) &&
          all(raw_region_names %in% crosswalk$region_name_gastat_en)
add("  [%s] every crosswalk region_name_gastat_en matches a raw CPI column, and vice versa",
    ok_fail(chk_d1))

geom_match <- match(crosswalk$geometry_spelling, geom_df$name)
chk_d2 <- !any(is.na(geom_match))
add("  [%s] every panel region maps to exactly one geometry (unmatched: %d)",
    ok_fail(chk_d2), sum(is.na(geom_match)))

chk_d3 <- !anyDuplicated(geom_match[!is.na(geom_match)])
add("  [%s] no duplicate geometry mappings (each geometry claimed at most once)", ok_fail(chk_d3))

chk_d4 <- length(unique(geom_match[!is.na(geom_match)])) == nrow(geom_df)
add("  [%s] all 13 geometries claimed by exactly one region (claimed: %d/%d)",
    ok_fail(chk_d4), length(unique(geom_match[!is.na(geom_match)])), nrow(geom_df))

add("")
add("  Geometry match table:")
for (i in seq_len(nrow(crosswalk))) {
  gname <- crosswalk$geometry_spelling[i]
  matched <- gname %in% geom_df$name
  add("    %-6s %-20s -> %-28s matched=%s", crosswalk$region_id[i],
      crosswalk$region_name_gastat_en[i], gname, matched)
}

## ---- E. Permanent region_id (immutable-ID hardening) ---------------------
add("")
add("E. PERMANENT region_id MAPPING")

expected_permanent_mapping <- c(
  "Al Baha" = "R01", "Al Jouf" = "R02", "Al Qassim" = "R03", "Aseer" = "R04",
  "Eastern Province" = "R05", "Hail" = "R06", "Jazan" = "R07", "Madinah" = "R08",
  "Makkah" = "R09", "Najran" = "R10", "Northern Borders" = "R11",
  "Riyadh" = "R12", "Tabouk" = "R13"
)

chk_e1 <- nrow(registry) == 13 && length(unique(registry$region_id)) == 13
add("  [%s] registry has exactly 13 unique region_id values", ok_fail(chk_e1))

registry_lookup <- setNames(registry$region_id, registry$canonical_region_name)
mapping_matches <- sapply(names(expected_permanent_mapping), function(rn) {
  identical(unname(registry_lookup[rn]), unname(expected_permanent_mapping[rn]))
})
chk_e2 <- all(mapping_matches)
add("  [%s] every region matches its predefined permanent region_id in the registry (%d/%d)",
    ok_fail(chk_e2), sum(mapping_matches), length(mapping_matches))

panel_lookup <- unique(panel[, c("Region", "region_id")])
panel_matches <- sapply(names(expected_permanent_mapping), function(rn) {
  identical(panel_lookup$region_id[panel_lookup$Region == rn], unname(expected_permanent_mapping[rn]))
})
chk_e3 <- all(panel_matches)
add("  [%s] built panel's region_id matches the same permanent mapping (%d/%d)",
    ok_fail(chk_e3), sum(panel_matches), length(panel_matches))

build_script_src <- readLines("scripts/09_panel_build.R")
chk_e4 <- !any(grepl("sort(raw_region_names)", build_script_src, fixed = TRUE)) &&
          any(grepl("region_id_registry.csv", build_script_src, fixed = TRUE))
add("  [%s] scripts/09_panel_build.R no longer derives region_id from alphabetical sort at runtime",
    ok_fail(chk_e4))

add("")
add("  Permanent region_id mapping:")
for (rn in names(expected_permanent_mapping)) {
  add("    %-6s %s", expected_permanent_mapping[rn], rn)
}

add("")
add("  Provenance check: Arabic-name / ISO columns are not mislabelled as GASTAT-sourced:")
chk_e5 <- all(c("region_name_ar_naturalearth", "iso_3166_2_naturalearth") %in% names(crosswalk)) &&
          !any(c("region_name_ar", "iso_3166_2") %in% names(crosswalk))
add("  [%s] crosswalk exposes region_name_ar_naturalearth / iso_3166_2_naturalearth, not a GASTAT-implying name",
    ok_fail(chk_e5))

## ---- F. Missing-value summary --------------------------------------------
add("")
add("F. MISSING-VALUE SUMMARY")
for (col in c("region_id", "Region", "date", "year", "month", "CPI", "inflation_mom", "inflation_yoy")) {
  n_na <- sum(is.na(panel[[col]]))
  add("    %-15s NA count = %d", col, n_na)
}
add("    Expected inflation_mom NA count = 13 (one structural first-obs per region): %s",
    ok_fail(sum(is.na(panel$inflation_mom)) == 13))
add("    Expected inflation_yoy NA count = 156 (12 structural first-obs per region x 13 regions): %s",
    ok_fail(sum(is.na(panel$inflation_yoy)) == 156))

## ---- G. National CPI differential ----------------------------------------
add("")
add("G. NATIONAL CPI DIFFERENTIAL")
add("  No official national Saudi CPI series exists anywhere in this repository")
add("  (data_raw/ contains only the 13-region file and the cached shapefile).")
add("  inflation_diff_mom / inflation_diff_yoy were NOT constructed, per instruction")
add("  not to fabricate a national series or silently use an unweighted 13-region")
add("  mean as an unlabelled proxy. See data_v2/README.md.")

## ---- Overall summary -------------------------------------------------
add("")
add(strrep("=", 70))
all_checks <- c(chk_a1, chk_a2, chk_a3, chk_a4, chk_a5,
                 chk_b1, chk_b2, chk_b3, chk_b4, chk_b5, chk_b6,
                 chk_c1, chk_c2, chk_c3, chk_c4,
                 chk_d1, chk_d2, chk_d3, chk_d4,
                 chk_e1, chk_e2, chk_e3, chk_e4, chk_e5)
add("OVERALL: %d/%d checks passed", sum(all_checks), length(all_checks))
if (!all(all_checks)) add("*** AT LEAST ONE CHECK FAILED -- DO NOT USE THIS PANEL UNTIL RESOLVED ***")

writeLines(lines, "data_v2/panel_validation_report.txt")
cat(paste(lines, collapse = "\n"), "\n")
