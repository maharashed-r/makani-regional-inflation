# ================================================================
# MAKANI v2 -- Canonical region x month panel builder
#
# NEW v2 FILE. Does not source, call, or modify any of scripts 01-08.
# Does not read or write data_clean/saudi_cpi_clean.csv or anything under
# results/ or baseline/. Reads only the authoritative raw inputs:
#   - data_raw/saudi_cpi_regions.csv   (GASTAT regional CPI, wide format)
#   - data_raw/saudi_states_sf.rds     (Natural Earth admin-1 geometry)
#
# Written in base R only (no dplyr/tidyr/sf) so it runs in any R
# installation without an renv restore -- this repo's own README already
# documents that the sf/units stack can fail to load on Windows under a
# non-ASCII project path, and the geometry file is only used here for its
# plain-text attribute columns (name/name_ar/iso_3166_2), not for any
# spatial operation, so no spatial package is actually needed.
#
# Outputs (all new, under data_v2/ -- nothing in data_clean/ or results/
# is touched):
#   data_v2/region_crosswalk.csv   Canonical region_id <-> every spelling
#                                   used anywhere in this repository.
#   data_v2/saudi_cpi_panel_v2.csv One row per region_id x month, CPI +
#                                   inflation_mom + inflation_yoy, with
#                                   structural lag NAs preserved (not
#                                   dropped, unlike data_clean).
#
# region_id is READ from data_v2/region_id_registry.csv, a permanent,
# version-controlled, hand-frozen lookup table -- it is never generated
# from alphabetical sorting (or any other runtime derivation) of the raw
# column names. This is deliberate: alphabetical order is a function of
# the CURRENT spelling and CURRENT column ordering in the raw source
# file, so an id derived that way could silently shift if either ever
# changes. The registry's ids (R01-R13) were frozen from the alphabetical
# assignment used in the first version of this script, so no existing
# panel value or id changes as a result of this hardening -- only the
# assignment *mechanism* changed, from "derive" to "look up."
#
# Run manually: Rscript scripts/09_panel_build.R
# ================================================================

dir.create("data_v2", recursive = TRUE, showWarnings = FALSE)

## ----------------------------------------------------------------
## 1. Region crosswalk
##
## The one and only place a raw-CPI-spelling <-> geometry-spelling alias
## map is defined for v2. Every future v2 script should read
## data_v2/region_crosswalk.csv instead of hard-coding its own recode()
## table (unlike scripts/04_spatial_analysis.R and
## scripts/07_sensitivity_analysis.R in the frozen v1.0 baseline, which
## each carry their own copy -- those are left untouched).
## ----------------------------------------------------------------

raw_header <- read.csv("data_raw/saudi_cpi_regions.csv", check.names = FALSE, nrows = 0)
raw_region_names <- setdiff(names(raw_header), c("Y", "M"))

if (length(raw_region_names) != 13) {
  stop("Expected exactly 13 region columns in data_raw/saudi_cpi_regions.csv, found ",
       length(raw_region_names), ": ", paste(raw_region_names, collapse = ", "))
}

# region_id assignment: looked up from the permanent registry, matched by
# exact name first, then by known alias. A raw column name that matches
# neither halts the build loudly (stop()) rather than silently minting a
# new id or falling back to any derived scheme -- a real future rename
# needs a human to extend data_v2/region_id_registry.csv, not code that
# guesses.
registry <- read.csv("data_v2/region_id_registry.csv", stringsAsFactors = FALSE, na.strings = "")

if (anyDuplicated(registry$region_id) || anyDuplicated(registry$canonical_region_name)) {
  stop("data_v2/region_id_registry.csv contains duplicate region_id or canonical_region_name values.")
}

lookup_region_id <- function(raw_name) {
  exact <- registry$region_id[registry$canonical_region_name == raw_name]
  if (length(exact) == 1) return(exact)

  alias_hit <- vapply(registry$known_aliases, function(a) {
    if (is.na(a) || !nzchar(a)) return(FALSE)
    raw_name %in% trimws(strsplit(a, ";")[[1]])
  }, logical(1))
  if (sum(alias_hit) == 1) return(registry$region_id[alias_hit])

  NA_character_
}

region_id <- vapply(raw_region_names, lookup_region_id, character(1))
names(region_id) <- raw_region_names

n_unmatched_id <- sum(is.na(region_id))
if (n_unmatched_id > 0) {
  stop(n_unmatched_id, " raw region name(s) have no permanent region_id in ",
       "data_v2/region_id_registry.csv (no exact or alias match): ",
       paste(raw_region_names[is.na(region_id)], collapse = ", "),
       ". Extend the registry -- do not derive an id at runtime.")
}
if (anyDuplicated(region_id)) {
  stop("Two or more raw region names resolved to the same region_id -- ambiguous registry match.")
}
if (!setequal(region_id, registry$region_id)) {
  stop("Not every permanent region_id in the registry was claimed by a raw region name -- ",
       "registry and raw source have diverged.")
}

# The one alias map this project needs (raw GASTAT spelling -> Natural
# Earth geometry spelling). Cross-checked for consistency against the
# recode() tables already present in the frozen baseline
# (scripts/04_spatial_analysis.R lines 25-38, scripts/07 lines 42-57) --
# reproduced here as the same correspondence, not altered, but defined
# once instead of twice. Correctness is verified programmatically below
# (join must hit all 13 geometries with none left over / none missing),
# not just asserted.
geometry_alias <- c(
  "Riyadh"            = "Ar Riyad",
  "Makkah"            = "Makkah",
  "Madinah"           = "Al Madinah",
  "Al Qassim"         = "Al Quassim",
  "Eastern Province"  = "Ash Sharqiyah",
  "Aseer"             = "`Asir",
  "Tabouk"            = "Tabuk",
  "Hail"              = "Ha'il",
  "Northern Borders"  = "Al Hudud ash Shamaliyah",
  "Jazan"             = "Jizan",
  "Najran"            = "Najran",
  "Al Baha"           = "Al Bahah",
  "Al Jouf"           = "Al Jawf"
)

missing_alias <- setdiff(raw_region_names, names(geometry_alias))
if (length(missing_alias) > 0) {
  stop("No geometry alias defined for raw region(s): ", paste(missing_alias, collapse = ", "))
}

# Geometry attribute table, read with base R only (readRDS() does not
# require the `sf` package to access a plain character column of an sf
# object -- verified: an sf object is a data.frame subclass).
geom_raw <- readRDS("data_raw/saudi_states_sf.rds")
geom_df  <- as.data.frame(geom_raw)[, c("name", "name_ar", "iso_3166_2")]
# Note: matched on the geometry's RAW name (e.g. "Ha'il"), unlike the
# baseline's 04_spatial_analysis.R which renames "Ha'il" -> "Hail" in
# place before joining. geometry_alias above maps "Hail" -> "Ha'il" to
# match this raw spelling directly.

if (any(duplicated(geom_df$name))) {
  stop("Duplicate geometry names found in data_raw/saudi_states_sf.rds -- cannot build a 1:1 crosswalk.")
}

crosswalk <- data.frame(
  region_id                      = region_id[raw_region_names],
  region_name_gastat_en          = raw_region_names,
  raw_data_spelling               = raw_region_names,
  geometry_spelling               = geometry_alias[raw_region_names],
  # Provenance made explicit in the column name itself, not just in a
  # note: these two columns are NOT from GASTAT. GASTAT publishes no
  # Arabic-language region file and no ISO code anywhere in this
  # repository -- both values are attributes of the cached Natural Earth
  # geometry file (data_raw/saudi_states_sf.rds), sourced there, not here.
  region_name_ar_naturalearth      = NA_character_,
  iso_3166_2_naturalearth          = NA_character_,
  legacy_aliases                  = NA_character_,
  notes                           = "",
  stringsAsFactors = FALSE
)

geom_match_idx <- match(crosswalk$geometry_spelling, geom_df$name)
n_unmatched_to_geometry <- sum(is.na(geom_match_idx))
if (n_unmatched_to_geometry > 0) {
  stop(n_unmatched_to_geometry, " crosswalk row(s) did not match any geometry name: ",
       paste(crosswalk$region_name_gastat_en[is.na(geom_match_idx)], collapse = ", "))
}
crosswalk$region_name_ar_naturalearth <- geom_df$name_ar[geom_match_idx]
crosswalk$iso_3166_2_naturalearth      <- geom_df$iso_3166_2[geom_match_idx]

n_unused_geometries <- sum(!seq_len(nrow(geom_df)) %in% geom_match_idx)
if (n_unused_geometries > 0) {
  stop(n_unused_geometries, " geometry row(s) were not claimed by any crosswalk row -- ",
       "expected exactly 13-to-13.")
}

# Known legacy/dead aliases: values that appear as recode() *targets* in
# the baseline scripts but never match any current raw column header --
# recorded for audit traceability, not because they are used by anything.
legacy_map <- list(
  "Al Jouf" = "Jouf",
  "Al Baha" = "Baha"
)
for (nm in names(legacy_map)) {
  crosswalk$legacy_aliases[crosswalk$region_name_gastat_en == nm] <- legacy_map[[nm]]
}
crosswalk$notes[!is.na(crosswalk$region_name_ar_naturalearth) & crosswalk$region_name_ar_naturalearth != ""] <-
  "region_name_ar_naturalearth / iso_3166_2_naturalearth sourced from data_raw/saudi_states_sf.rds (Natural Earth admin-1), NOT from GASTAT -- GASTAT publishes no Arabic-language region file or ISO code anywhere in this repository."

crosswalk <- crosswalk[order(crosswalk$region_id), ]
write.csv(crosswalk, "data_v2/region_crosswalk.csv", row.names = FALSE, na = "")
message("Wrote data_v2/region_crosswalk.csv (", nrow(crosswalk), " regions).")

## ----------------------------------------------------------------
## 2. Long-format region x month panel (base R pivot, explicit chronological
##    ordering by integer year/month -- same defensive ordering rationale
##    as the baseline's 01_load_data.R, re-derived independently here
##    rather than by sourcing that script).
## ----------------------------------------------------------------

raw <- read.csv("data_raw/saudi_cpi_regions.csv", check.names = FALSE)
n_raw_rows_incl_annual <- nrow(raw)

raw <- raw[raw$M != "-", ]                 # drop the annual-average row per year
n_annual_rows_dropped <- n_raw_rows_incl_annual - nrow(raw)
raw$Y <- as.integer(raw$Y)
raw$M <- as.integer(raw$M)

if (any(is.na(raw$Y)) || any(is.na(raw$M))) {
  stop("Non-numeric Y or M value survived coercion to integer -- silent coercion failure.")
}
if (any(raw$M < 1 | raw$M > 12)) {
  stop("Month value(s) outside 1-12 found after coercion -- impossible date.")
}

panel_list <- lapply(raw_region_names, function(rn) {
  data.frame(
    region_id = unname(region_id[rn]),
    Region    = rn,
    year      = raw$Y,
    month     = raw$M,
    CPI       = raw[[rn]],
    stringsAsFactors = FALSE
  )
})
panel <- do.call(rbind, panel_list)
panel$date <- as.Date(sprintf("%04d-%02d-01", panel$year, panel$month))

if (any(!is.finite(panel$CPI))) {
  stop(sum(!is.finite(panel$CPI)), " non-finite CPI value(s) found -- cannot proceed.")
}
if (any(panel$CPI <= 0)) {
  stop(sum(panel$CPI <= 0), " non-positive CPI value(s) found -- CPI must be strictly positive.")
}

# Explicit chronological sort -- the exact class of bug documented in the
# frozen baseline's 01_load_data.R header comment (text-sorted month order)
# is re-derived independently here rather than assumed away.
panel <- panel[order(panel$region_id, panel$year, panel$month), ]

if (anyDuplicated(paste(panel$region_id, panel$date))) {
  stop("Duplicate region_id x date rows found after pivot.")
}

## ----------------------------------------------------------------
## 3. Inflation: month-over-month and year-over-year, matched on an
##    explicit linear month index (year*12+month) per region_id rather
##    than on row position, so a hypothetical gap in the monthly sequence
##    could never silently misalign a lag (row-position lag would not
##    catch that; date-indexed lookup does, by construction).
## ----------------------------------------------------------------

panel$month_index <- panel$year * 12L + panel$month
panel$lookup_key   <- paste(panel$region_id, panel$month_index)

cpi_by_key <- setNames(panel$CPI, panel$lookup_key)

lag_cpi <- function(n) {
  key <- paste(panel$region_id, panel$month_index - n)
  unname(cpi_by_key[key])
}

cpi_lag1  <- lag_cpi(1)
cpi_lag12 <- lag_cpi(12)

panel$inflation_mom <- 100 * (log(panel$CPI) - log(cpi_lag1))
panel$inflation_yoy <- 100 * (log(panel$CPI) - log(cpi_lag12))

panel$month_index <- NULL
panel$lookup_key   <- NULL
panel <- panel[, c("region_id", "Region", "date", "year", "month", "CPI",
                    "inflation_mom", "inflation_yoy")]

write.csv(panel, "data_v2/saudi_cpi_panel_v2.csv", row.names = FALSE, na = "NA")
message("Wrote data_v2/saudi_cpi_panel_v2.csv (", nrow(panel), " region-month rows).")

## ----------------------------------------------------------------
## 4. National CPI differential -- checked, not fabricated.
## ----------------------------------------------------------------

national_source_found <- FALSE  # no national/aggregate CPI column or file exists anywhere in this repository (verified by inspection of data_raw/ and the full repo tree; see data_v2/README.md)
if (national_source_found) {
  stop("national_source_found is TRUE but no code path to build inflation_diff_* was written -- ",
       "this branch is intentionally unimplemented pending an actual authoritative national series.")
} else {
  message("No official national Saudi CPI series found in this repository or its raw inputs. ",
          "inflation_diff_mom / inflation_diff_yoy were NOT constructed (per instruction: no ",
          "unweighted-13-region-mean proxy without explicit labelling). See data_v2/README.md.")
}

## ----------------------------------------------------------------
## 5. Console summary (also captured by scripts/10_panel_validation.R)
## ----------------------------------------------------------------

message("---- Build summary ----")
message("Raw source rows (incl. annual-average rows): ", n_raw_rows_incl_annual)
message("Annual-average rows dropped: ", n_annual_rows_dropped)
message("Regions: ", length(raw_region_names))
message("Panel rows (region x month, all preserved): ", nrow(panel))
message("Date range: ", format(min(panel$date)), " to ", format(max(panel$date)))
message("inflation_mom non-NA: ", sum(!is.na(panel$inflation_mom)), " / ", nrow(panel))
message("inflation_yoy non-NA: ", sum(!is.na(panel$inflation_yoy)), " / ", nrow(panel))
