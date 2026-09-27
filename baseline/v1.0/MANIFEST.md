# MAKANI Exploratory Baseline v1.0

**Freeze date:** 2026-09-21
**Git commit at freeze:** `ef2ee9473b33fd65c3ea7ef99e08aa2ee75b9a00`
**Commit subject:** "Add Project Architecture overview page to the report" (2026-08-05T23:13:53+03:00)
**Repository state at freeze:** clean except one pre-existing untracked directory (`.github/`, deliberately dropped from tracking per commit `d438f55` — unrelated to this baseline, not included in it)

This document freezes the current MAKANI pipeline and its outputs as a named,
checksummed reference point: **"MAKANI Exploratory Baseline v1.0."** Nothing
in `data_raw/`, `data_clean/`, `scripts/`, `results/`, `paper.qmd`, or
`renv.lock` was modified to produce this manifest — every hash below was
computed by reading the files as they already existed. See
`baseline/v1.0/MANIFEST.sha256` for a machine-checkable checksum file
(`sha256sum -c baseline/v1.0/MANIFEST.sha256` from the repo root; verified
OK for all 36 entries at freeze time).

See also `baseline/v1.0/METHODOLOGY_NOTE.md` for the interpretive caveats
that apply to every result listed below.

---

## 1. Input files

| File | SHA256 | Size |
|---|---|---|
| `data_raw/saudi_cpi_regions.csv` | `9d7b4ff3c7368fde6d0403397ad0cab8a855e283d99e71efb2177eeb29ffb1e` | 39,121 B |
| `data_raw/saudi_states_sf.rds` | `25267763331d0d81a17676e55a6b7c6f9dd99742cbddcf17853ee2e56486b40` | 56,405 B |

`saudi_cpi_regions.csv`: monthly regional CPI, base 2023=100, GASTAT, wide
format, Jan 2013–Feb 2026, 13 regions + an annual-average row per year.
`saudi_states_sf.rds`: cached `rnaturalearth` Saudi Arabia state-boundary
polygons (avoids a network fetch on every run).

## 2. Rebuilt cache (not a second source of truth)

| File | SHA256 | Size |
|---|---|---|
| `data_clean/saudi_cpi_clean.csv` | `b6deca38d6813416390826fca55dea50a89b951e22aeb6a338135ee867b3cc0` | 122,691 B |

Deterministically regenerated from the two input files above by
`scripts/01_load_data.R` on every pipeline run. 2,041 rows = 13 regions ×
157 months (Feb 2013–Feb 2026, first month per region dropped for the
log-diff lag). Verified at freeze time: 0 missing values, 0 duplicate
`(date, Region)` pairs, fully balanced panel (157 rows/region, confirmed
directly).

## 3. Scripts

| Script | SHA256 | Run by `run_project.R`? |
|---|---|---|
| `scripts/01_load_data.R` | `f29b5b2ef86ee85e2b109ff8e15854768f12fb5e1f8858df50ba7989a10d3b0` | yes (1st) |
| `scripts/02_clean_data.R` | `d4d84c96f0d17edbd30cd03b76f6fed8f7340447ac2c3a430ac772dc64478fc` | yes (2nd) |
| `scripts/03_exploratory_analysis.R` | `3cf8a2052e8cb3409af65e19609a4c5d7b721d6e37062ffbd1f31c8af6782fb` | yes (3rd) |
| `scripts/04_spatial_analysis.R` | `96bbfe2bf8ad6f48b43112e71ea84fe60163f7d17bdf61991444b6a2a84c397` | yes (4th) |
| `scripts/06_spatial_autocorrelation.R` | `1857ecfb624cd45c5c06543f094c095351313fef4bcc0dc4dca4b0938c9c968` | yes (5th) |
| `scripts/07_sensitivity_analysis.R` | `55414184bb3b50287669ceaec0b245d4c1badedf61261dfcbb3ab1607d06bff` | yes (6th) |
| `scripts/05_plots.R` | `7e330b93b7b8ed757b3b977936fe291889156d4af1f88862733ffc014423922` | yes (7th, last) |
| `scripts/08_spatial_regression.R` | `3b9772ac838cc505a5fb83f3712f70249b5d7e26ac6a61bde6954ffb0c01457` | **no** — manual only, intercept-only SAR/SEM, deliberately excluded (see §5 caveats) |
| `scripts/run_project.R` | `85157b286066d06e2ca9addc709c9cc4ffc46b5a1452e92e76d258a8864f903` | orchestrator itself |

Execution order confirmed identical in `run_project.R` and `paper.qmd`'s
setup chunk: `01 → 02 → 03 → 04 → 06 → 07 → 05`.

## 4. Report source

| File | SHA256 |
|---|---|
| `paper.qmd` | `0d60868a32b2d885c84244f076f5d03c35e56d1c278a0b523e98c89d3fb11d4` |

`paper.html`/`paper_files/` are build output, not tracked, regenerated with
`quarto render paper.qmd`; `paper.docx` is the tracked submission
deliverable and is not separately hashed here since it is a Quarto render
artifact of `paper.qmd`, not an independent source file.

## 5. Environment / package versions

| | |
|---|---|
| R version | 4.4.1 |
| `renv.lock` SHA256 | `9e206e58fa2ed071fa54248f89816fd97c5c62efe8c467700a14ecf183a676b` |
| Total locked packages | 116 |

Key packages actually exercised by the analysis scripts (full list of 116 is
in `renv.lock`; restore the exact environment with `renv::restore()`):

| Package | Version | Used by |
|---|---|---|
| `sf` | 1.1-0 | 04, 06, 07, 08 (geometry) |
| `spdep` | 1.4-2 | 06, 07, 08 (Moran's I, LISA, weights matrices) |
| `spatialreg` | 1.4-3 | 08 (SAR/SEM) |
| `dplyr` | 1.2.1 | all scripts |
| `tidyr` | 1.3.2 | 01, 03, 07 |
| `reshape2` | 1.4.5 | 03 (correlation melt) |
| `ggplot2` | 4.0.3 | 05, 07 (all plots) |
| `rnaturalearth` | 1.2.0 | 04 (map fallback if cache absent) |
| `classInt` | 0.4-11 | transitive (sf/spdep classification) |
| `units` | 1.0-1 | transitive (sf) |
| `Rcpp` | 1.1.1 | transitive (sf/spdep compiled backends) |
| `DT` | 0.34.0 | paper.qmd HTML tables |
| `flextable` | 0.9.11 | paper.qmd Word tables |
| `plotly` | 4.12.0 | paper.qmd HTML interactive figures |
| `knitr` | 1.51 | paper.qmd rendering |
| `rmarkdown` | 2.31 | paper.qmd rendering |

## 6. Output tables

### 6a. Current-pipeline-generated (reproduced by `run_project.R` + manual `08`)

| File | SHA256 | Written by |
|---|---|---|
| `results/tables/region_statistics.csv` | `ed798444099e728ff3ebb3db88f9e00f35ac55b794260d07a07a0d9bdb43fed` | `03` |
| `results/tables/regional_inflation_correlation.csv` | `aa8c9a80a26848ef16898ddcee28489bd6fcfed2cf0193b7ff6292ea436573c` | `03` |
| `results/tables/moran_test.csv` | `056fa25fd738c717772d5e47eabbb6142b47e2e5ea65f7a80f5b9cdb8a41545` | `06` (asymptotic Moran's I) |
| `results/tables/moran_mc_test.csv` | `2b13a75364f870c6dd9a133751dfcd859c9dee58fc7887662be4c867549474a` | `06` (permutation Moran's I) |
| `results/tables/sensitivity_analysis.csv` | `3c5cc6fb1b95846dcf7317f9ed3fd07988ddbd7546f2cb5099be6e6f4c7e01f` | `07` (VAT-window scenarios) |
| `results/tables/weights_sensitivity.csv` | `670f54a8b5c4539714bad60b79da45cca7942686e8ec5176e0b71fe03658b47` | `07` |
| `results/tables/lisa_results.csv` | `c75ccbfc61708f16c8c2168d9d86fbac50cf465c40527ad56e3b55ce007fa56` | `07` |
| `results/tables/correlation_significance.csv` | `6b82400c14cd6cf58b8274cb36d4132640d5af29b81e5efcd3340d81593466d` | `07` (+ BH-FDR) |
| `results/tables/moran_by_year.csv` | `9d6ab9d58ad0cd781511bcfb0278dbdc0a26ba7e51bb71bf59c3ce4faa221e9` | `07` |
| `results/tables/spatial_model_comparison.csv` | `3280a154c40317d949271f97ed103c7946e61119fb23642d95a1384b045584a` | `08` (manual — see caveats) |

### 6b. Orphaned / legacy — NOT written by any current script

| File | SHA256 | Note |
|---|---|---|
| `results/tables/region_inflation_stats.csv` | `c88d8ef2f8b5fcebafbef9711b77bd15432559107566936209180f33804e9c3` | Duplicate byte-identical copy also at `results/region_inflation_stats.csv`. **Values do not match `region_statistics.csv`** (e.g. Al Baha mean inflation 0.00104 here vs. 0.00218 in the current table) — consistent with the pre-fix text-sort bug documented in `scripts/01_load_data.R`'s header comment. Predates the current pipeline; kept on disk for audit traceability only, **not to be cited as a current result.** |
| `results/region_inflation_stats.csv` | `c88d8ef2f8b5fcebafbef9711b77bd15432559107566936209180f33804e9c3` | Same file as above, stray top-level copy. |
| `results/tables/summary_statistics.csv` | `705e1d915530b3e2974053aebe3a30b23320aa2cdfef46bb0703f84c53707a9` | Single-row national aggregate; not produced by any current script. Its min/max values match the stale `region_inflation_stats.csv`, confirming the same pre-fix vintage. Not to be cited as current. |

None of these three are read by `paper.qmd` (verified by search) — they are
inert on disk, not part of the reporting chain, but were **not deleted or
overwritten** per the freeze requirement to preserve everything exactly as
found.

## 7. Output figures

### 7a. Current-pipeline-generated

| File | SHA256 | Written by |
|---|---|---|
| `results/figures/histogram.png` | `62593b14968bc76ff6c99e62e8e3737dcb574c5da19045682109e70cf344244` | `05` |
| `results/figures/map.png` | `f3603abb06a66ef2a92d8a1d091fe23e841497fb1c2b188486b832eb397054f` | `05` |
| `results/figures/correlation.png` | `e2d7908367d59f3c650d30944a9c206c4a4d2103688d6ef44edf5324933300` | `05` |
| `results/figures/moran_scatter_plot.png` | `8c3e74aa696d8a33ae272cad6bdf90e602425b9f0d989f06daccfc5784985ed` | `06` |
| `results/figures/lisa_cluster_map.png` | `91df2f596277dc3bf144c593dc0d85167bc7b75ad5e30b058dd137c28e54205` | `07` |
| `results/figures/moran_over_time.png` | `8f8d54a851a7bd8ab67a12423ff3568e81b08c917123f8a23b6d8e7fb4903b4` | `07` |

### 7b. Orphaned / legacy — NOT written by any current script

| File | SHA256 |
|---|---|
| `results/figures/inflation_correlation_heatmap.png` | `51074a5b2577dc6a6a703199b162490be26c4219aaba2734501460874bee3e4` |
| `results/figures/regional_correlation_heatmap.png` | `3e11092cf555f023d1f0b5083325c1ed37ae5e23078f38f4c0a911fbbeb5d83` |
| `results/figures/regional_inflation_facets.png` | `64dbaa64f558ec2f4a88afec73d496141bdce8da5e282acccc0583d859fddd6` |
| `results/figures/regional_inflation_trends.png` | `c5b3b9074ddc84ed49e7c2f9a102e54bb48741aae27660f73293c53059d9821` |
| `results/figures/saudi_inflation_map.png` | `50f0efc666b56530e51e5a79d405e988045185845eb720ffb1f1d5b659170716` |

Same caveat as §6b: predate the current script set (all last touched in the
same commit, `3cc941f`, as everything else), not referenced by `paper.qmd`,
preserved unchanged rather than deleted.

## 8. Headline numeric results at freeze (for quick reference only — full values in the tables above)

- Global Moran's I (asymptotic): I=0.204, z=1.68, p=0.046
- Global Moran's I (Monte Carlo, 9999 perms): I=0.204, p=0.0575
- LISA: 2/13 regions significant (Makkah High-High, Al Jawf Low-Low, p≈0.031 each, uncorrected)
- SAR/SEM (intercept-only): rho=lambda=0.433, LR p=0.262, OLS AIC beats both

## 9. How to verify this baseline later

```bash
cd <repo root>
sha256sum -c baseline/v1.0/MANIFEST.sha256
```

Any line reporting `FAILED` means that file has changed since this freeze —
investigate before treating current `results/` as equivalent to v1.0.
