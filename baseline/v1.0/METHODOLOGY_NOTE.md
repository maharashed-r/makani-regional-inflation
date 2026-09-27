# Methodology Note — MAKANI Exploratory Baseline v1.0

This note states the interpretive limits that apply to every result frozen
in `baseline/v1.0/MANIFEST.md`. It does not change any script or any
number — it documents how the existing numbers should and should not be
read. See the accompanying full code audit (prior session) for the
line-by-line basis of each point below.

## 1. The current analysis primarily collapses the time dimension

`data_clean` holds 2,041 region-month observations (13 regions × 157
months, Feb 2013–Feb 2026), but every spatial-statistical step — the Global
and Local Moran's I in `06`/`07`, the weights-matrix sensitivity checks, and
the SAR/SEM models in `08` — operates on `inflation_sf$mean_inflation`: one
number per region, produced by averaging each region's inflation series
across the *entire* sample (`04_spatial_analysis.R`, lines 18–20). The
VAT-window sensitivity scenarios and the year-by-year Moran's I in `07`
re-collapse to a regional mean under a filtered sub-sample or within a
single calendar year respectively — still a fresh cross-sectional collapse
each time, never a joint model over both dimensions at once. No script in
this baseline fits a space-time (panel) spatial model. Consequently, the
richness of the monthly panel is used for the *descriptive* and
*correlation* results (`03`'s region statistics and correlation matrix,
which do operate row-by-row on the full panel) but not for any of the
formally spatial statistics.

## 2. n=13 spatial inference is exploratory

Every Moran's I, LISA, and SAR/SEM result in this baseline is computed on a
cross-section of 13 regions. At this sample size:

- Asymptotic p-values (`moran.test()`'s default) rely on a
  normal/randomisation approximation whose accuracy at n=13 is not
  guaranteed — this is why `06_spatial_autocorrelation.R` also runs a
  permutation test (see §3).
- The year-by-year Moran's I (`07`, §5) reruns the same n=13 test 14 times
  independently; 0 of 14 years reach conventional significance, which is
  expected given how little power a 13-unit test has, not evidence that no
  yearly pattern exists.
- LISA classifies each of the 13 regions with its own local p-value with
  **no multiple-testing correction applied** (unlike the 78-pair regional
  correlation test, which does get Benjamini-Hochberg correction in `07`,
  §4). At n=13 simultaneous tests, a Bonferroni-type threshold would
  require p≈0.0038 to survive; the two "significant" LISA clusters in this
  baseline (Makkah High-High, Al Jawf Low-Low) sit at p≈0.031 each — well
  short of that bar. Read LISA output in this baseline as
  **hypothesis-generating, not confirmatory**, regardless of the p-value
  column's face value.

## 3. Permutation Moran's I is preferred over asymptotic inference

`06_spatial_autocorrelation.R` computes Global Moran's I twice: once via
`moran.test()` (asymptotic, p=0.046) and once via `moran.mc()` (9,999
permutations, p=0.0575). The permutation test makes no distributional
assumption beyond exchangeability under the null and is the more defensible
figure at this sample size. **The asymptotic p-value should not be quoted
on its own** as evidence of significant spatial clustering in this
baseline — it crosses the conventional 0.05 threshold while its permutation
counterpart does not, and the permutation result is the one that should
govern the headline conclusion ("suggestive, not conclusive" spatial
clustering, as `paper.qmd`'s own Conclusion already states).

## 4. SAR/SEM intercept-only models are not economic models

`08_spatial_regression.R` fits `mean_inflation ~ 1` (no covariates) under
both a Spatial Lag (SAR) and Spatial Error (SEM) specification. With no
regional covariates, this is a mathematically degenerate special case: for
an intercept-only mean structure, the SAR and SEM formulations reduce to
observationally equivalent estimates (`rho` and `lambda` come out identical
at 0.433 in this baseline — not a coincidence, but a property of having no
X beyond the constant). Neither model's likelihood-ratio test against OLS
is significant here (p=0.262), and OLS has the lower (better) AIC than
either spatial specification. **This script does not estimate what drives
regional inflation** — it only re-expresses the same spatial dependence
already summarised by Moran's I in regression form, with no explanatory
covariates behind it. This is why `run_project.R` excludes it from the main
pipeline and why it must not be cited as evidence of a working spatial
econometric model of inflation. It is retained in this baseline exactly as
a documented negative/placeholder result, pending future covariate work
(explicitly out of scope for this baseline).

---

**Summary for any reader of `baseline/v1.0/`:** this baseline detects a
modest, sample-size-limited signal of spatial structure in regional
inflation (via the more consistently robust correlation and permutation-Moran
results) but does not yet explain it. Treat every number in
`MANIFEST.md` §8 as descriptive/exploratory evidence, not as a settled
finding — that distinction is the entire point of labelling this
"Exploratory Baseline v1.0" rather than a final result.
