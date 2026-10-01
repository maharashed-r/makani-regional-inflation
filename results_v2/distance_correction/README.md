# Distance correction (post-v1.0, October 2026)

## What changed

The KNN (k=4) and distance-band spatial-weight matrices are built from the
centroids of the 13 regions. The centroids are WGS84 longitude/latitude
(EPSG:4326).

- **Up to v1.0:** Euclidean distance was applied directly to degrees.
- **Now:** distances are great-circle distances in km (`spdep`: `longlat = TRUE`
  in `knearneigh()`, `nbdists()` and `dnearneigh()`).

The change was made in scripts 12, 13, 14, 15, 16 and 18. Nothing else
changed:
- data, sample, regions, k = 4;
- the distance-band rule (1.05 x the largest nearest-neighbour distance);
- row standardisation, seeds and bootstrap designs;
- the BH families and estimators;
- the primary Queen matrix, which uses no distances.

The decision was made on geographic grounds before any result was seen. A
degree of longitude is about 94–107 km at these latitudes, a degree of
latitude about 111 km. The threshold of 4.88 "degrees" had no physical unit.

## Weight graphs (`weights_graph_comparison.csv`)

| Matrix | Directed links before → after | Change | Connected |
|---|---|---|---|
| KNN k=4 | 52 → 52 | R09→R07 replaced by R09→R12 | yes / yes |
| KNN k=1 (defines the band) | 13 → 13 | none | — |
| Distance band | 54 → 46 | removed R01↔R10, R02↔R08, R03↔R09, R10↔R12 | yes / yes |

Distance-band threshold: 4.884 degrees → **500.9 km**. This is the same rule,
not a tuned value.

## Results (`headline_before_after.csv`, `before_after_changed_cells.csv`)

| Result | Before | After |
|---|---|---|
| Final model (Stage 18), joint spatial test, KNN | p = 0.9979, Θ3 = 0.0001 | p = 0.9849, Θ3 = 0.0004 |
| Final model (Stage 18), joint spatial test, distance band | p = 0.8317, Θ3 = −0.0605 | p = 0.8949, Θ3 = −0.0235 |
| Final model (Stage 18), Queen (primary) | p = 0.3156, Θ3 = −0.0098 | unchanged |
| Spline robustness (Stage 15), KNN | p = 0.092 | p = 0.132 |
| Spline robustness (Stage 15), distance band | p = 0.022 | p = 0.006 |
| Event × weights (Stage 16), E3 VAT 2020, KNN | p = 0.075 | **p = 0.011** |
| Event × weights (Stage 16), E3 VAT 2020, distance band | p = 0.085 | **p = 5.6e-09** |
| Static Moran's I of mean inflation, KNN | I = 0.071, p = 0.161 | I = 0.084, p = 0.141 |
| Static Moran's I of mean inflation, distance band | I = 0.014, p = 0.287 | I = 0.038, p = 0.264 |

All p-values above are Driscoll–Kraay asymptotic unless stated. Every changed
cell is listed in `before_after_changed_cells.csv`.

## Conclusions

- **Unchanged:** the primary conclusion (no significant joint spatial
  predictive dependence) holds under all three matrices, before and after.
- **Unchanged:** the Stage-15 spline result still does not replicate under KNN.
- **Unchanged:** the Stage-16 pre-specified family (Queen, wild-cluster
  bootstrap, BH) does not change.
- **New qualification:** in the Stage-16 weight-robustness check, the E3 (July
  2020 VAT increase, COVID-confounded) spatial-dynamics test becomes
  nominally significant under KNN and distance-band weights (asymptotic DK; no
  bootstrap was run for robustness variants). The Queen primary test is not
  significant (DK p = 0.144; wild-cluster bootstrap p = 0.646). The manuscript
  now reports this.

## Invariance (`invariance_all_result_tables.csv`)

- Of the 137 result tables and classification files under `results_v2/`, 128 are identical to `f88d776`.
- The 9 changed tables are all distance-based (KNN or distance band).
- Every Queen-based table is identical, including all wild-cluster bootstraps
  (B = 9999), LORO, BH families and the final model.

Two regenerated files differed only because of the audit environment and
were restored to their committed versions:
- `results_v2/tables/lag_dk_implementation_audit.csv`: the recorded
  `sandwich` version string;
- `results_v2/nonlinear_spatial_dynamics/tables/nl_B_marginal_effect_by_state.csv`:
  only the `support_density` column, from `stats::density()`, which changed
  in R 4.4.0. All estimates were identical.

Of the regenerated figures, only the four that plot KNN or distance-band
results are updated: `fig_F2`, `fig_G3`, and final-robustness `fig_2` and
`fig_4`. The others were restored.

## Environment and reproduction

Re-run with R 4.3.3, spdep 1.3.1, sf 1.0.15 and sandwich 3.1.0 (the lock file
pins R 4.4.1 and spdep 1.4-2). With the original code, this environment
reproduced the committed `fr_D` table exactly (|Δp| ≤ 3e-16).

To reproduce this audit:

```
Rscript results_v2/distance_correction/distance_correction_audit.R
```

Run it from a Git checkout. It reads the BEFORE state with `git show`.

The v1.0 baseline (`scripts/07`, `results/`) is frozen and unchanged; see
`baseline/v1.0/METHODOLOGY_NOTE.md`, section 5.
