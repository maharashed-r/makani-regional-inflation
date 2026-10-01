# ================================================================
# MAKANI -- tests for the spatial-weight construction (post-v1.0
# great-circle distance correction, 2026-10).
#
# Rebuilds the weights exactly as scripts/12-16 and 18 do and checks:
# geographic coordinates, great-circle distances, link counts, the
# distance-band threshold, connectivity, Queen invariance and
# determinism. Paths are relative to this file (../../ = repo root),
# as in test-panel-v2.R.
# ================================================================

library(testthat)
suppressMessages({ library(sf); library(spdep); library(dplyr) })

repo_root <- "../.."
geom_path      <- file.path(repo_root, "data_raw/saudi_states_sf.rds")
crosswalk_path <- file.path(repo_root, "data_v2/region_crosswalk.csv")
registry_path  <- file.path(repo_root, "data_v2/region_id_registry.csv")

build_weights <- function() {
  crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)
  region_order <- sort(unique(crosswalk$region_id))
  geom_raw <- readRDS(geom_path)
  sf_frame <- crosswalk %>%
    select(region_id, region_name_gastat_en, geometry_spelling) %>%
    arrange(match(region_id, region_order)) %>%
    left_join(geom_raw %>% select(name), by = c("geometry_spelling" = "name"))
  sf_frame <- st_as_sf(sf_frame)
  coords <- st_coordinates(st_centroid(st_geometry(sf_frame)))
  nb_queen <- poly2nb(sf_frame, queen = TRUE)
  nb_knn4  <- knn2nb(knearneigh(coords, k = 4, longlat = TRUE))
  nb_knn1  <- knn2nb(knearneigh(coords, k = 1, longlat = TRUE))
  band_km  <- max(unlist(nbdists(nb_knn1, coords, longlat = TRUE))) * 1.05
  nb_dist  <- dnearneigh(coords, 0, band_km, longlat = TRUE)
  list(sf = sf_frame, coords = coords, region_order = region_order, nb_queen = nb_queen,
       nb_knn4 = nb_knn4, nb_knn1 = nb_knn1, nb_dist = nb_dist, band_km = band_km)
}

skip_if_not(file.exists(geom_path) && file.exists(crosswalk_path), "geometry or crosswalk not found")
w <- build_weights()

haversine_km <- function(a, b) {
  r <- pi / 180
  h <- sin((b[2] - a[2]) * r / 2)^2 + cos(a[2] * r) * cos(b[2] * r) * sin((b[1] - a[1]) * r / 2)^2
  2 * 6371.0088 * asin(sqrt(h))
}

test_that("centroids are WGS84 longitude/latitude", {
  expect_true(isTRUE(st_is_longlat(w$sf)))
  expect_equal(st_crs(w$sf)$epsg, 4326L)
  expect_true(all(w$coords[, 1] > 34 & w$coords[, 1] < 56))   # longitude, Saudi Arabia
  expect_true(all(w$coords[, 2] > 16 & w$coords[, 2] < 33))   # latitude
})

test_that("every active stage builds distance-based weights with longlat = TRUE", {
  scripts <- file.path(repo_root, "scripts", c(
    "12_residual_spatial_tests.R", "13_lagged_spatial_transmission.R",
    "14_distributed_spatial_lags.R", "15_nonlinear_spatial_dynamics.R",
    "16_structural_shocks.R", "18_final_robustness_model_lock.R"))
  for (f in scripts) {
    src <- readLines(f, warn = FALSE)
    calls <- grep("knearneigh\\(coords|nbdists\\(nb_knn1, coords|dnearneigh\\(coords", src, value = TRUE)
    expect_equal(length(calls), 4L, info = f)
    expect_true(all(grepl("longlat = TRUE", calls)), info = f)
  }
})

test_that("great-circle distances are active (km, not degrees)", {
  d <- unlist(nbdists(w$nb_knn1, w$coords, longlat = TRUE))
  hv <- unlist(lapply(seq_along(w$nb_knn1), function(i) haversine_km(w$coords[i, ], w$coords[w$nb_knn1[[i]], ])))
  expect_true(all(d > 100))                 # degree-scale values would be < 15
  expect_lt(max(abs(d - hv)), 2)            # spdep great circle vs haversine, km
})

test_that("KNN k=4 and k=1 graphs", {
  expect_equal(sum(card(w$nb_knn4)), 52L)
  expect_true(all(card(w$nb_knn4) == 4L))
  expect_equal(sum(card(w$nb_knn1)), 13L)
})

test_that("distance band: unchanged rule, corrected threshold and link count", {
  expect_equal(w$band_km, 500.9, tolerance = 0.5, scale = 1)
  expect_equal(sum(card(w$nb_dist)), 46L)   # directed links (23 undirected)
  expect_true(is.symmetric.nb(w$nb_dist, verbose = FALSE))
  expect_true(all(card(w$nb_dist) >= 1L))
})

test_that("all weight graphs are connected (one component, no islands)", {
  for (nb in list(w$nb_queen, w$nb_knn4, w$nb_dist)) {
    expect_equal(n.comp.nb(nb)$nc, 1L)
    expect_true(all(card(nb) > 0L))
  }
})

test_that("Queen contiguity does not depend on distance and is unchanged", {
  expect_equal(sum(card(w$nb_queen)), 52L)  # 26 undirected edges, as in spatial_weights_comparison.csv
  W <- listw2mat(nb2listw(w$nb_queen, style = "W", zero.policy = TRUE))
  expect_true(all(abs(rowSums(W) - 1) < 1e-8))
})

test_that("weight construction is deterministic", {
  w2 <- build_weights()
  expect_identical(w2$nb_knn4, w$nb_knn4)
  expect_identical(w2$nb_dist, w$nb_dist)
  expect_identical(w2$band_km, w$band_km)
})
