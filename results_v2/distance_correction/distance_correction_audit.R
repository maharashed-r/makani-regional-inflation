# ================================================================
# MAKANI -- post-v1.0 distance correction: BEFORE/AFTER audit.
#
# BEFORE = the committed state at main f88d776 (Euclidean distance on
# longitude/latitude degrees for KNN and distance-band weights).
# AFTER  = the working tree after scripts 12-16 and 18 were re-run with
# great-circle distances (spdep longlat = TRUE).
#
# Reads BEFORE files with `git show`, so it must be run from the
# repository root of a Git checkout:
#   Rscript results_v2/distance_correction/distance_correction_audit.R
# Writes only into results_v2/distance_correction/.
# ================================================================

suppressMessages({ library(sf); library(spdep); library(dplyr) })
BASE <- "f88d77661bffbc33460d0d404c24859bcbea457a"
out_dir <- "results_v2/distance_correction"

git_csv <- function(path) {
  txt <- system2("git", c("show", paste0(BASE, ":", path)), stdout = TRUE)
  read.csv(text = paste(txt, collapse = "\n"), stringsAsFactors = FALSE, check.names = FALSE)
}

## ---- 1. Weight graphs, rebuilt both ways -------------------------------------
crosswalk <- read.csv("data_v2/region_crosswalk.csv", stringsAsFactors = FALSE)
region_order <- sort(unique(crosswalk$region_id))
geom_raw <- readRDS("data_raw/saudi_states_sf.rds")
sf_frame <- crosswalk %>% select(region_id, region_name_gastat_en, geometry_spelling) %>%
  arrange(match(region_id, region_order)) %>%
  left_join(geom_raw %>% select(name), by = c("geometry_spelling" = "name")) %>% st_as_sf()
coords <- st_coordinates(st_centroid(st_geometry(sf_frame)))

build <- function(longlat) {
  knn4 <- knn2nb(knearneigh(coords, k = 4, longlat = longlat))
  knn1 <- knn2nb(knearneigh(coords, k = 1, longlat = longlat))
  band <- max(unlist(nbdists(knn1, coords, longlat = longlat))) * 1.05
  list(knn4 = knn4, knn1 = knn1, band = band, dist = dnearneigh(coords, 0, band, longlat = longlat))
}
before <- build(FALSE); after <- build(TRUE)
edges <- function(nb) unlist(lapply(seq_along(nb), function(i)
  if (length(nb[[i]]) && nb[[i]][1] != 0) paste0(region_order[i], "->", region_order[nb[[i]]])))
graph_row <- function(name, a, b, unit_a, unit_b, thr_a = NA, thr_b = NA) {
  ea <- edges(a); eb <- edges(b)
  data.frame(matrix = name,
             directed_links_before = length(ea), directed_links_after = length(eb),
             links_removed = paste(setdiff(ea, eb), collapse = "; "),
             links_added = paste(setdiff(eb, ea), collapse = "; "),
             components_before = n.comp.nb(a)$nc, components_after = n.comp.nb(b)$nc,
             min_card_before = min(card(a)), max_card_before = max(card(a)),
             min_card_after = min(card(b)), max_card_after = max(card(b)),
             threshold_before = thr_a, threshold_unit_before = unit_a,
             threshold_after = thr_b, threshold_unit_after = unit_b)
}
graphs <- rbind(
  graph_row("KNN k=4", before$knn4, after$knn4, NA, NA),
  graph_row("KNN k=1 (defines the band rule)", before$knn1, after$knn1, NA, NA),
  graph_row("Distance band (1.05 x max nearest-neighbour distance)", before$dist, after$dist,
            "Euclidean degrees (no physical unit)", "great-circle km", before$band, after$band))
write.csv(graphs, file.path(out_dir, "weights_graph_comparison.csv"), row.names = FALSE)

## ---- 2. Static Moran's I of mean monthly inflation by region (diagnostic) ---
panel <- read.csv("data_v2/saudi_cpi_panel_v2.csv", stringsAsFactors = FALSE)
mean_infl <- tapply(panel$inflation_mom, panel$region_id, mean, na.rm = TRUE)[region_order]
moran_row <- function(name, nb_b, nb_a) {
  mb <- moran.test(mean_infl, nb2listw(nb_b, style = "W"), randomisation = TRUE)
  ma <- moran.test(mean_infl, nb2listw(nb_a, style = "W"), randomisation = TRUE)
  data.frame(matrix = name, moran_I_before = unname(mb$estimate[1]), p_before = mb$p.value,
             moran_I_after = unname(ma$estimate[1]), p_after = ma$p.value)
}
moran <- rbind(moran_row("KNN k=4", before$knn4, after$knn4),
               moran_row("Distance band", before$dist, after$dist))
write.csv(moran, file.path(out_dir, "moran_static_comparison.csv"), row.names = FALSE)

## ---- 3. Every changed/unchanged result table of stages 12-18 ---------------
changed <- system2("git", c("diff", "--name-only", BASE, "--", "results_v2"), stdout = TRUE)
changed <- changed[grepl("\\.csv$", changed) & !grepl("^results_v2/distance_correction/", changed)]
cells <- list(); status <- list()
for (f in changed) {
  a <- git_csv(f); b <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(dim(a), dim(b)) || !identical(names(a), names(b))) {
    status[[f]] <- data.frame(file = f, cells_changed = NA, max_abs_numeric_change = NA, note = "shape changed"); next
  }
  n <- 0; mx <- 0
  for (j in seq_along(a)) for (i in seq_len(nrow(a))) {
    va <- a[i, j]; vb <- b[i, j]
    same <- if (is.numeric(va) && is.numeric(vb)) isTRUE(all.equal(va, vb, tolerance = 0)) || (is.na(va) && is.na(vb)) else identical(as.character(va), as.character(vb))
    if (!same) {
      n <- n + 1
      d <- if (is.numeric(va) && is.numeric(vb)) abs(vb - va) else NA
      if (!is.na(d)) mx <- max(mx, d)
      cells[[length(cells) + 1]] <- data.frame(file = f, row = i, row_label = as.character(a[i, 1]), column = names(a)[j],
                                               before = as.character(va), after = as.character(vb), abs_change = d)
    }
  }
  status[[f]] <- data.frame(file = f, cells_changed = n, max_abs_numeric_change = mx, note = "")
}
write.csv(do.call(rbind, cells), file.path(out_dir, "before_after_changed_cells.csv"), row.names = FALSE)
write.csv(do.call(rbind, status), file.path(out_dir, "before_after_changed_files.csv"), row.names = FALSE)

## ---- 4. Headline comparison (Stage 18, Section D; Stage 15 spline) ---------
fr_b <- git_csv("results_v2/final_robustness/tables/fr_D_spatial_weight_robustness.csv")
fr_a <- read.csv("results_v2/final_robustness/tables/fr_D_spatial_weight_robustness.csv", stringsAsFactors = FALSE)
nl_b <- git_csv("results_v2/nonlinear_spatial_dynamics/tables/nl_F_B_spline_robustness.csv")
nl_a <- read.csv("results_v2/nonlinear_spatial_dynamics/tables/nl_F_B_spline_robustness.csv", stringsAsFactors = FALSE)
headline <- rbind(
  data.frame(source = "fr_D (Stage 18, final model, DK)", spec = fr_b$weight_matrix,
             joint_p_before = fr_b$joint_p_value, joint_p_after = fr_a$joint_p_value,
             Theta3_before = fr_b$Theta_3, Theta3_after = fr_a$Theta_3,
             ci_before = sprintf("[%.4f, %.4f]", fr_b$ci_low, fr_b$ci_high),
             ci_after = sprintf("[%.4f, %.4f]", fr_a$ci_low, fr_a$ci_high),
             conclusion_before = ifelse(fr_b$joint_p_value > 0.05, "not significant", "significant"),
             conclusion_after = ifelse(fr_a$joint_p_value > 0.05, "not significant", "significant")),
  data.frame(source = "nl_F_B (Stage 15, spline robustness, DK)", spec = nl_b$spec,
             joint_p_before = nl_b$p_value_asymptotic, joint_p_after = nl_a$p_value_asymptotic,
             Theta3_before = NA, Theta3_after = NA, ci_before = NA, ci_after = NA,
             conclusion_before = ifelse(nl_b$p_value_asymptotic > 0.05, "not significant", "nominally significant"),
             conclusion_after = ifelse(nl_a$p_value_asymptotic > 0.05, "not significant", "nominally significant")))
ss_b <- git_csv("results_v2/structural_shocks/tables/ss_J_robustness_weights_by_event.csv")
ss_a <- read.csv("results_v2/structural_shocks/tables/ss_J_robustness_weights_by_event.csv", stringsAsFactors = FALSE)
headline <- rbind(headline, data.frame(
  source = "ss_J (Stage 16, event x weights, spatial-dynamics block, DK asymptotic)",
  spec = paste(ss_b$event_id, ss_b$weight_spec), joint_p_before = ss_b$p_value_asymptotic, joint_p_after = ss_a$p_value_asymptotic,
  Theta3_before = NA, Theta3_after = NA, ci_before = NA, ci_after = NA,
  conclusion_before = ifelse(ss_b$p_value_asymptotic > 0.05, "not significant", "nominally significant"),
  conclusion_after = ifelse(ss_a$p_value_asymptotic > 0.05, "not significant", "nominally significant")))
write.csv(headline, file.path(out_dir, "headline_before_after.csv"), row.names = FALSE)

## ---- 5. Queen / non-distance invariance: every other result table ----------
all_csv <- system2("git", c("ls-tree", "-r", "--name-only", BASE, "--", "results_v2"), stdout = TRUE)
all_csv <- all_csv[grepl("\\.(csv|txt)$", all_csv)]
inv <- data.frame(file = all_csv,
                  status = ifelse(all_csv %in% changed, "CHANGED: EXPECTED_DISTANCE_CORRECTION", "IDENTICAL to f88d776"))
write.csv(inv, file.path(out_dir, "invariance_all_result_tables.csv"), row.names = FALSE)
cat("Result tables identical:", sum(!all_csv %in% changed), " changed:", sum(all_csv %in% changed), "\n")
print(graphs[, 1:7]); print(moran); print(do.call(rbind, status)); print(headline)
