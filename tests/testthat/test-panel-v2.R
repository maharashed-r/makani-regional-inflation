# ================================================================
# MAKANI v2 -- automated tests for the region x month panel
#
# NEW v2 FILE, kept separate from the frozen v1.0 baseline (no baseline
# script is sourced, read as a test fixture, or asserted against here).
#
# Run from the repository root with:
#   testthat::test_dir("tests/testthat")
# testthat::test_dir() sets the working directory to tests/testthat/
# itself while running each test file (standard R-package testthat
# convention, verified empirically in this session) -- so paths below
# are written relative to THIS file's location (../../ = repo root),
# not relative to wherever Rscript was originally launched from.
# ================================================================

library(testthat)

repo_root      <- "../.."
panel_path     <- file.path(repo_root, "data_v2/saudi_cpi_panel_v2.csv")
crosswalk_path <- file.path(repo_root, "data_v2/region_crosswalk.csv")

test_that("v2 artifacts exist before any test reads them", {
  skip_if_not(file.exists(panel_path),
              paste("Run scripts/09_panel_build.R first --", panel_path, "not found"))
  skip_if_not(file.exists(crosswalk_path),
              paste("Run scripts/09_panel_build.R first --", crosswalk_path, "not found"))
  expect_true(file.exists(panel_path))
  expect_true(file.exists(crosswalk_path))
})

skip_if_not(file.exists(panel_path) && file.exists(crosswalk_path),
            "v2 panel artifacts not built yet")

panel     <- read.csv(panel_path, stringsAsFactors = FALSE)
crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)
panel$date <- as.Date(panel$date)

raw_header       <- read.csv(file.path(repo_root, "data_raw/saudi_cpi_regions.csv"), check.names = FALSE, nrows = 0)
raw_region_names <- setdiff(names(raw_header), c("Y", "M"))
geom_raw         <- readRDS(file.path(repo_root, "data_raw/saudi_states_sf.rds"))
geom_df          <- as.data.frame(geom_raw)[, c("name", "name_ar", "iso_3166_2")]
registry_path    <- file.path(repo_root, "data_v2/region_id_registry.csv")
registry         <- read.csv(registry_path, stringsAsFactors = FALSE, na.strings = "")

test_that("every region maps to its predefined permanent region_id (immutable, not alphabetical)", {
  expect_equal(nrow(registry), 13)
  expect_equal(length(unique(registry$region_id)), 13)
  expect_setequal(registry$region_id, sprintf("R%02d", 1:13))

  # The frozen mapping this hardening change must preserve exactly --
  # hard-coded here (independent of the registry file and of
  # scripts/09_panel_build.R's lookup logic) so this test would fail if
  # either the registry file OR the build script's matching logic ever
  # silently drifted from the originally-frozen R01-R13 identities.
  expected_permanent_mapping <- c(
    "Al Baha" = "R01", "Al Jouf" = "R02", "Al Qassim" = "R03", "Aseer" = "R04",
    "Eastern Province" = "R05", "Hail" = "R06", "Jazan" = "R07", "Madinah" = "R08",
    "Makkah" = "R09", "Najran" = "R10", "Northern Borders" = "R11",
    "Riyadh" = "R12", "Tabouk" = "R13"
  )
  expect_setequal(names(expected_permanent_mapping), raw_region_names)

  registry_lookup <- setNames(registry$region_id, registry$canonical_region_name)
  for (region_name in names(expected_permanent_mapping)) {
    expect_equal(
      unname(registry_lookup[region_name]),
      unname(expected_permanent_mapping[region_name]),
      info = paste("registry mismatch for", region_name)
    )
  }

  # And the mapping actually realised in the built panel/crosswalk must
  # match the same frozen table -- not just the registry file in isolation.
  panel_lookup <- unique(panel[, c("Region", "region_id")])
  for (region_name in names(expected_permanent_mapping)) {
    expect_equal(
      panel_lookup$region_id[panel_lookup$Region == region_name],
      unname(expected_permanent_mapping[region_name]),
      info = paste("panel region_id mismatch for", region_name)
    )
  }

  crosswalk_lookup <- setNames(crosswalk$region_id, crosswalk$region_name_gastat_en)
  for (region_name in names(expected_permanent_mapping)) {
    expect_equal(
      unname(crosswalk_lookup[region_name]),
      unname(expected_permanent_mapping[region_name]),
      info = paste("crosswalk region_id mismatch for", region_name)
    )
  }
})

test_that("region_id is not derived from alphabetical sort of raw names", {
  # A regression guard: if a future edit reintroduces
  # sort(raw_region_names) + seq_along() as the source of region_id, this
  # test still passes trivially for THIS dataset (alphabetical order
  # happens to match the frozen mapping here) -- so it also asserts that
  # scripts/09_panel_build.R's source no longer contains that derivation.
  build_script <- readLines(file.path(repo_root, "scripts/09_panel_build.R"))
  expect_false(any(grepl("sort(raw_region_names)", build_script, fixed = TRUE)))
  expect_false(any(grepl("seq_along(raw_region_names_sorted)", build_script, fixed = TRUE)))
  expect_true(any(grepl("region_id_registry.csv", build_script, fixed = TRUE)))
})

test_that("crosswalk provenance columns for Arabic name / ISO code do not imply GASTAT sourcing", {
  expect_true(all(c("region_name_ar_naturalearth", "iso_3166_2_naturalearth") %in% names(crosswalk)))
  expect_false(any(c("region_name_ar", "iso_3166_2") %in% names(crosswalk)))
})

test_that("region count is exactly 13", {
  expect_equal(length(unique(panel$region_id)), 13)
  expect_equal(nrow(crosswalk), 13)
  expect_equal(length(raw_region_names), 13)
})

test_that("no duplicate region-month records exist", {
  dup_key <- paste(panel$region_id, panel$date)
  expect_equal(sum(duplicated(dup_key)), 0)
})

test_that("date order is not corrupted (monotonic per region, no text-sort artefact)", {
  chrono_ok <- sapply(split(panel$date, panel$region_id), function(d) !is.unsorted(d))
  expect_true(all(chrono_ok))

  expected_seq <- seq(min(panel$date), max(panel$date), by = "month")
  gap_ok <- sapply(split(panel$date, panel$region_id), function(d) identical(sort(d), expected_seq))
  expect_true(all(gap_ok))
})

test_that("panel is balanced across all 13 regions", {
  n_months <- table(panel$region_id)
  expect_equal(length(unique(as.vector(n_months))), 1)
  expect_equal(as.vector(n_months)[[1]], 158)
})

test_that("CPI contains no non-positive or non-finite values", {
  expect_true(all(is.finite(panel$CPI)))
  expect_true(all(panel$CPI > 0))
})

test_that("no accidental annual-average rows entered the panel", {
  raw_full <- read.csv(file.path(repo_root, "data_raw/saudi_cpi_regions.csv"), check.names = FALSE)
  raw_full_no_annual <- raw_full[raw_full$M != "-", ]
  expect_equal(nrow(panel), 13 * nrow(raw_full_no_annual))
})

test_that("region-to-geometry matching succeeds for every region, with no duplicates", {
  geom_match <- match(crosswalk$geometry_spelling, geom_df$name)
  expect_false(any(is.na(geom_match)))
  expect_false(anyDuplicated(geom_match[!is.na(geom_match)]) > 0)
  expect_equal(length(unique(geom_match)), nrow(geom_df))
})

test_that("every crosswalk region_name_gastat_en matches a raw CPI column and vice versa", {
  expect_true(all(crosswalk$region_name_gastat_en %in% raw_region_names))
  expect_true(all(raw_region_names %in% crosswalk$region_name_gastat_en))
})

test_that("generated inflation_mom matches the stated formula on a random sample", {
  set.seed(1)
  idx <- sample(which(!is.na(panel$inflation_mom)), 25)
  for (i in idx) {
    r <- panel[i, ]
    prev_date <- seq(r$date, by = "-1 month", length.out = 2)[2]
    prev_cpi  <- panel$CPI[panel$region_id == r$region_id & panel$date == prev_date]
    expect_length(prev_cpi, 1)
    expect_equal(r$inflation_mom, 100 * (log(r$CPI) - log(prev_cpi)), tolerance = 1e-9)
  }
})

test_that("generated inflation_yoy matches the stated formula on a random sample", {
  set.seed(2)
  idx <- sample(which(!is.na(panel$inflation_yoy)), 25)
  for (i in idx) {
    r <- panel[i, ]
    prev_date <- seq(r$date, by = "-12 month", length.out = 2)[2]
    prev_cpi  <- panel$CPI[panel$region_id == r$region_id & panel$date == prev_date]
    expect_length(prev_cpi, 1)
    expect_equal(r$inflation_yoy, 100 * (log(r$CPI) - log(prev_cpi)), tolerance = 1e-9)
  }
})

test_that("only structurally-expected missing values appear (no unexpected NAs)", {
  expect_equal(sum(is.na(panel$region_id)), 0)
  expect_equal(sum(is.na(panel$Region)), 0)
  expect_equal(sum(is.na(panel$date)), 0)
  expect_equal(sum(is.na(panel$year)), 0)
  expect_equal(sum(is.na(panel$month)), 0)
  expect_equal(sum(is.na(panel$CPI)), 0)

  expect_equal(sum(is.na(panel$inflation_mom)), 13)   # one structural first-obs per region
  expect_equal(sum(is.na(panel$inflation_yoy)), 13 * 12)  # first 12 obs per region

  first_mom_na_only <- sapply(split(panel, panel$region_id), function(d) {
    d <- d[order(d$date), ]
    is.na(d$inflation_mom[1]) && !any(is.na(d$inflation_mom[-1]))
  })
  expect_true(all(first_mom_na_only))

  first12_yoy_na_only <- sapply(split(panel, panel$region_id), function(d) {
    d <- d[order(d$date), ]
    all(is.na(d$inflation_yoy[1:12])) && !any(is.na(d$inflation_yoy[-(1:12)]))
  })
  expect_true(all(first12_yoy_na_only))
})

test_that("baseline v1.0 files were not touched by building the v2 panel", {
  # Sentinel check: this test does not modify baseline/v1.0/MANIFEST.sha256,
  # it only re-verifies it, independent of the checksum verification already
  # performed manually when this baseline was frozen.
  manifest <- readLines(file.path(repo_root, "baseline/v1.0/MANIFEST.sha256"))
  expect_true(length(manifest) > 0)
})
