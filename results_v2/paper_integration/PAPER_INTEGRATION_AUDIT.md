# MAKANI v2 — Stage 19: Paper Integration Audit

This document records the audit of `paper.qmd` performed before editing
(per the task's explicit instruction: audit first, identify obsolete
claims, replace only what v2 evidence supersedes, preserve valid prose),
and the final consistency check performed after integration and rendering.

## 1. Audit of the pre-existing `paper.qmd`

The document as found was entirely a **v1 (Stages 1-8) exploratory
analysis**: it sourced `scripts/01`-`07` directly, and its content was
limited to descriptive statistics, a regional correlation matrix, a
single full-sample Global Moran's I test (plus a Monte Carlo variant),
Local Indicators of Spatial Association (LISA), a VAT-window sensitivity
check, and a year-by-year exploratory Moran's I breakdown. It contained
**no reference to the MAKANI v2 confirmatory sequence** (Stages 9-18:
panel construction, TWFE decomposition, monthly permutation Moran's I,
lag-1 and distributed spatial lags, nonlinear/state-dependent tests,
structural-shock tests, stationarity, seasonal/temporal tests, or final
robustness). It also contained one internal defect independent of v1/v2
status: **its table numbering skipped from Table 4 to Table 6** (no Table
5 existed anywhere in the document).

## 2. Obsolete / superseded claims identified

| # | v1 claim (as originally written) | Why it is superseded by v2 evidence |
|---|---|---|
| 1 | Full-sample Global Moran's I (I ≈ real, "clears conventional significance... though close to the threshold") presented as the paper's primary spatial-dependence evidence | Stage 12's monthly permutation Moran's I (157 months, raw + TWFE-residual) finds 0/157 months survive BH or Holm correction — the single full-sample snapshot does not represent robust, multiply-tested spatial dependence |
| 2 | Economic Interpretation: "regions actively transmitting inflation to one another" implicitly entertained as a live possibility alongside common shocks | Stages 13-18 (lag-1, distributed K=3/K=6, final robustness) find no robust spatial *predictive* dependence at any tested horizon, weight matrix, or inference method |
| 3 | Conclusion: "The evidence for spatial clustering... is real but modest... best read as suggestive rather than conclusive" presented as the paper's final word | Superseded by the final v2 classification (Stage 18): **B — Weak/suggestive but non-robust spatial dependence**, arrived at through a full confirmatory sequence, not a single test |
| 4 | Limitations: "does not yet include... theoretical framework... covariate-based spatial econometric model... is the most valuable next step" | Stages 14-18 *did* build the dynamic panel / distributed-lag / robustness framework the original limitations section called for (though region-level covariates, e.g. population/trade/sectoral data, genuinely remain unavailable — this specific sub-point is preserved, updated) |
| 5 | Reproducibility section referencing only `scripts/01`-`07` | Superseded by the full v2 pipeline (Stages 9-18) and its commit history |
| 6 | Table numbering skip (Table 4 → Table 6, no Table 5) | Structural defect, fixed independent of v1/v2 content (renumbered A1-A5 in the relocated appendix) |

## 3. What was preserved verbatim (or near-verbatim)

All v1 descriptive/statistical prose — the correlation-significance
discussion, the LISA regional interpretation (Makkah/Al Jawf), the VAT
sensitivity-window analysis, the alternative-spatial-weights check, and
the year-by-year Moran's I breakdown — was preserved in full in
**Appendix A: Exploratory Baseline Analysis (v1, Preliminary)**, since
these remain valid *descriptive* facts about the data (a real correlation
structure exists; LISA clusters exist as computed). Only their
*interpretation* as evidence of spatial interaction was reframed, via a
short prefatory note and targeted forward-references to the main text,
rather than deleted.

## 4. What was added

- Abstract, Literature Review (6 verified citations — see below),
  Empirical Strategy (final K=3 TWFE specification + 10-stage pipeline
  table), a new Results section built entirely from programmatic
  `read.csv()` calls against Stage 12-17 output, a Robustness and
  Additional Analyses section built from Stage 18 output, Discussion,
  an updated Limitations section, an updated Conclusion stating the
  final spatial (B) and temporal (T2) classifications, References, and
  an expanded Reproducibility section listing every v2 stage's commit
  hash.
- Five new main-text figures (Figures 1-5) and five main-text tables
  (Tables 1-5), all sourced from already-committed `results_v2/*`
  outputs — none newly computed by this stage.

## 5. Literature sources (verified, not invented)

1. Ciccarelli, M., & Mojon, B. (2010). Global Inflation. *The Review of
   Economics and Statistics*, 92(3), 524-535.
2. Marques, H., Pino, G., & Tena, J. D. Regional Inflation Dynamics Using
   Space-Time Models. *Empirical Economics*.
3. Elhorst, J. P. (2017). Spatial Panel Data Analysis. In *Encyclopedia
   of GIS* (2nd ed.). Springer.
4. Auer, R., Levchenko, A. A., & Sauré, P. (2017). International
   Inflation Spillovers through Input Linkages. BIS Working Paper No. 623.
5. Hasan, M. M., & Alogeel, H. (2008). Understanding the Inflationary
   Process in the GCC Region: The Case of Saudi Arabia and Kuwait. IMF
   Working Paper WP/08/193.
6. Fareed, F., Rezghi, A., & Sandoz, C. (2023). Inflation Dynamics in the
   Gulf Cooperation Council (GCC): What is the Role of External Factors?
   IMF Working Paper WP/23/263.

Each was located and its exact title/authors/venue/year verified via web
search before citing (see conversation record); none was cited from
memory alone without verification.

## 6. Final consistency audit (post-integration)

| Check | Result |
|---|---|
| Abstract matches final results | PASS — all abstract numbers are `r`-computed from `results_v2/*` CSVs at render time |
| Introduction does not overclaim | PASS — explicitly states "We do not begin from the premise that regional spatial spillovers exist" |
| Methods match actual code | PASS — Empirical Strategy section's formula matches `fr_ref`/`f_ref` exactly as estimated in Stage 18 |
| All reported numbers match Stage 14-18 outputs | PASS — every numeric result in Results/Robustness is read programmatically, not retyped |
| No outdated v1 Moran-only conclusion remains as the paper's final word | PASS — v1 conclusion relocated to Appendix A with explicit superseding note; main Conclusion states the v2 classification |
| Nonlinear findings described as suggestive/non-robust | PASS — "suggestive, non-actionable evidence, not... an established nonlinear relationship" |
| Structural-shock results described as null | PASS — "No robust event-conditioned parameter instability is found" |
| Stationarity findings stated correctly | PASS — inflation stationary (13/13 regions), CPI level non-stationary, matches Stage 17 exactly |
| Final model is Stage-14 K=3 TWFE | PASS — stated explicitly in Empirical Strategy and Table 3 |
| Spatial classification is B | PASS — rendered inline from `fr_J_final_classifications.txt` |
| Temporal classification is T2 | PASS — rendered inline from the same file |
| No causal wording appears (asserting an established causal claim) | PASS — every occurrence of "causal", "transmission", "spillover", "contagion" is either (a) negated/hedged, (b) an explicit limitation statement, (c) a literal citation title, or (d) a script filename |
| Every citation exists | PASS — 6/6 verified via web search against publisher/institutional sources |
| Every table/figure referenced in text exists | PASS — Tables 1-5 and Figures 1-5 (main text) and Tables A1-A5 and Figures A0-A6 (appendix, renumbered to fix the original skip) all confirmed present in rendered output |
| Render status | HTML: success (exit 0, zero errors). DOCX: success (exit 0, zero errors). PDF: **not rendered** — no LaTeX engine (`tinytex::is_tinytex()` = FALSE) is installed in this environment, and per instruction no new software was installed to enable it; this is documented as a limitation rather than worked around |

## 7. Files not modified, and one deliberate, instructed baseline exception

No file under `baseline/`, `scripts/01`-`18`, `data_v2/`, `tests/`,
`research_notes/`, or any `results_v2/*` subdirectory other than
`results_v2/paper_integration/` (this directory) was created or modified
by Stage 19. `results/` (the v1 outputs read by the Appendix) was read
but not modified.

**One expected, deliberate exception:** `baseline/v1.0/MANIFEST.sha256`
checksums the *live* repo-root `paper.qmd` (and `paper.html`/`paper.docx`)
as part of the frozen v1.0 snapshot — `sha256sum -c` on the manifest now
correctly reports `paper.qmd: FAILED` (37/38 entries OK, not 38/38). This
is the **intended, explicitly instructed** consequence of Stage 19's own
task specification ("Update paper.qmd... Render the manuscript"), not an
accidental baseline modification: `paper.qmd` is the project's living
manuscript, and integrating the completed v2 evidence into it is this
stage's entire purpose. `baseline/v1.0/MANIFEST.sha256`,
`baseline/v1.0/MANIFEST.md`, and `baseline/v1.0/METHODOLOGY_NOTE.md`
themselves are **not** edited — the historical v1.0 freeze record remains
exactly as written at freeze time (2026-09-21, commit `72147ab`/`ef2ee94`);
only the live document it once also checksummed has now advanced, as
instructed. Every other one of the manifest's 37 remaining entries
(`scripts/`, `results/`, `data_raw/`, `data_clean/`, `renv.lock`) verifies
unchanged.
