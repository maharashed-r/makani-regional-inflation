# Reproducibility Guide

This guide takes a researcher from a fresh clone of this repository to the
tables, figures, and rendered manuscript reported in `paper.qmd`. No step
requires data beyond what is already bundled (see `DATA_AVAILABILITY.md`).

## 1. Requirements

- **R**: developed and locked against R 4.4.1 (see `.Rprofile` / `renv.lock`
  for the exact version `renv` will install if it differs from your system R).
- **Quarto** (optional, only needed to re-render `paper.qmd`): any recent
  release from <https://quarto.org/>.
- **Disk**: the full repository (code + bundled data + all committed
  results) is on the order of tens of megabytes; `renv::restore()` will
  additionally download and cache R package sources/binaries (typically a
  few hundred MB, kept in your `renv` cache, not in this repository).
- **No LaTeX / PDF renderer is required** to reproduce any table, figure, or
  the HTML/Word manuscript. A LaTeX distribution (e.g. TinyTeX) was **not**
  available in the environment this project was developed in, so PDF
  rendering of `paper.qmd` was never exercised here; `quarto render
  paper.qmd --to pdf` may work if you separately install a LaTeX engine
  (`tinytex::install_tinytex()`), but this is untested by the authors.

## 2. Environment setup

```r
# from the repository root, in R or RStudio (open MAKANI.Rproj first)
install.packages("renv")   # if not already installed
renv::restore()            # installs the exact package versions in renv.lock
```

### Windows note: non-ASCII project paths

If your local path to this repository contains non-ASCII characters (e.g.
Arabic), the `units`/`sf` package stack can fail to load with `error: no
database found!` (the underlying `udunits2` C library cannot open its XML
unit database through a non-ASCII path on Windows). Fixes, in order of
preference:

1. Clone/move the repository to an ASCII-only path (e.g. `C:\projects\makani`).
2. Or, without moving it, create a local `.Renviron` file (already
   `.gitignore`d — never commit it) in the repository root pointing `renv`'s
   package library at an ASCII-only location:
   ```
   RENV_PATHS_LIBRARY=C:/Users/<you>/AppData/Local/renv/library/<project-name>
   ```
   then run `renv::restore()` again.

## 3. Directory structure this project expects

Everything is relative to the repository root — no external directories or
environment-specific absolute paths are required. Scripts read from
`data_raw/`, `data_v2/`, and (for the v1 baseline pipeline) `data_clean/`,
and write to `results/` or `results_v2/<stage-name>/{tables,figures}/`,
all already present as of this commit and safe to regenerate in place.

## 4. Two pipelines, run in this order if reproducing from scratch

This repository contains **two** linked analyses that were developed
sequentially. Both are fully reproducible; the second builds on the first's
canonical panel but does not depend on the first's *results*.

### 4a. v1 — exploratory baseline (`scripts/01`–`08`)

```bash
Rscript scripts/run_project.R
```

Runs `01_load_data.R → 02_clean_data.R → 03_exploratory_analysis.R →
04_spatial_analysis.R → 06_spatial_autocorrelation.R →
07_sensitivity_analysis.R → 05_plots.R` in sequence, writing to
`data_clean/` and `results/`. `08_spatial_regression.R` (intercept-only
SAR/SEM diagnostic, not part of the main narrative) is run separately:
`Rscript scripts/08_spatial_regression.R`.

### 4b. v2 — confirmatory dynamic-panel sequence (`scripts/09`–`20`)

Run in numeric order; each stage is independent enough to inspect on its
own, but later stages read earlier stages' committed outputs as a
read-only reference for sample/consistency checks, so numeric order is
recommended for a from-scratch reproduction:

```bash
Rscript scripts/09_panel_build.R
Rscript scripts/10_panel_validation.R
Rscript scripts/11_twfe_common_shocks.R
Rscript scripts/12_residual_spatial_tests.R
Rscript scripts/13_lagged_spatial_transmission.R
Rscript scripts/14_distributed_spatial_lags.R
Rscript scripts/15_nonlinear_spatial_dynamics.R
Rscript scripts/16_structural_shocks.R
Rscript scripts/17_stationarity_temporal_dynamics.R
Rscript scripts/18_final_robustness_model_lock.R
Rscript scripts/20_power_mde_presubmission.R
```

Each writes to its own `results_v2/<stage-name>/{tables,figures}/`
directory (see the top-level `README.md` for the exact mapping) plus a
stage-level `*_REPORT.md`.

**Runtime warning:** several stages (15, 16, 18, 20) run Monte Carlo /
wild-cluster-bootstrap procedures with thousands to tens of thousands of
model refits and can take from tens of minutes to several hours on a
typical laptop, dominated by Stage 20's power simulation. Every such
script prints progress per grid point / section as it runs, and every
random element uses the fixed seed `20260921` for exact reproducibility of
the reported numbers (subject to normal cross-platform floating-point and
BLAS variation).

### Validation tests

```bash
Rscript -e 'testthat::test_dir("tests/testthat", package = "MAKANI")'
```

or open `MAKANI.Rproj` in RStudio and use **Build → Test**. The suite
validates the v2 canonical panel (`data_v2/saudi_cpi_panel_v2.csv`)
structurally and statistically (177 expectations across 16 test blocks).

## 5. Rendering the manuscript

```bash
quarto render paper.qmd --to html
quarto render paper.qmd --to docx
```

`paper.qmd` sources the v1 pipeline scripts directly in its setup chunk
(for the Appendix A exploratory material) and reads every v2 (Stages
9–20) number programmatically from the already-committed
`results_v2/*/tables/*.csv` files — it does **not** re-run any Stage
9–20 script itself. Run the v2 pipeline (§4b) first if those files are
not already present. `paper.html` and `paper_files/` are build output
(`.gitignore`d) and are not part of this repository; regenerate them with
the command above.

## 6. What "reproduce the core results" means here

Re-running §4b end to end and re-rendering `paper.qmd` should reproduce,
to within Monte Carlo/floating-point tolerance:

- the TWFE variance decomposition (Table 2);
- the final K=3 dynamic model's coefficients and inference (Table 3);
- the distributed-lag and robustness tables (Tables 4–5);
- the final spatial classification **B**, temporal classification **T2**,
  and power classification **P3** (`results_v2/final_robustness/` and
  `results_v2/presubmission_refinement/power_mde/`).

These three classifications are locked conclusions of this research
program as released; a faithful re-run should reproduce them exactly (the
underlying joint-test p-values and power estimates may vary by a small
Monte Carlo margin, not the qualitative classification).
