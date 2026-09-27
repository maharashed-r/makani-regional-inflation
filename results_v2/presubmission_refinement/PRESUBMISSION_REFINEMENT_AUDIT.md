# MAKANI v2 — Pre-Submission Scientific Refinement Audit

This stage performed ONLY the four tasks authorized: (1) tighten framing/
contribution language, (2) strengthen the N=13 limitation discussion, (3)
add a simulation-based power/MDE analysis, (4) rewrite Future Research as a
disciplined program. **No new econometric model was added anywhere, and no
Stage 14–18 empirical value was changed.**

## 1. Manuscript phrases materially changed, and why

| Location | Before | After | Why |
|---|---|---|---|
| Title | "...Common Shocks, Local Persistence..." | "...Common Temporal Variation, Local Persistence..." | "Common shocks" implies an identified mechanism; month fixed effects identify only shared temporal variation, not which economic force produces it |
| Abstract, opening question | "common national temporal shocks" | "common temporal variation" | Same reason, applied throughout |
| Abstract, closing sentence | "...common national shocks and local region-specific dynamics... dominate the empirical pattern..." | Rewritten to state the exact preferred central conclusion, end with "no robust spatial predictive dependence is detected in the available data," and add one sentence on the power-analysis implication | Matches the task's explicitly preferred phrasing; avoids the prohibited "no spatial dependence exists" framing; discloses the power-analysis caveat at the abstract level since the P3 classification justifies it |
| Literature Review, "How this study differs" | "those common national forces (month fixed effects)" | "common temporal variation (captured statistically via month fixed effects, without asserting a specific identified mechanism)" | Prevents month FE being read as an identified causal mechanism |
| Results, §5.1 heading | "Common month effects dominate the explainable variance" | "Common temporal variation dominates the explainable variance" | Terminology consistency; "dominates" is retained here because it describes a directly *statistically estimated* quantity (variance share), which is legitimate under category A |
| Results, §5.1 prose | plain "Month fixed effects alone account for X%..." | Added an explicit gloss that month FE absorb shared variation "without identifying which specific mechanism... drives it" and a closing sentence distinguishing statistical decomposition from causal identification | Same reason |
| Discussion, intro sentence | "a large, statistically dominant common monthly component" | Rewritten to explicitly label the A/B/C distinction (statistically estimated / economically interpreted / possible explanation) required by Task 1 | Task 1's explicit instruction to distinguish these three categories |
| Discussion, bullet 3 | "...so the *residual* spatial test correctly finds little left to explain once the common shock is removed" | "...has, by construction, absorbed whatever common temporal variation exists that month, from whatever source... regardless of which specific national force... generated that common variation" | Avoids asserting month FE = a specific identified "shock" |
| Conclusion, opening | Generic "Across ten confirmatory stages..." opening | Opens with the exact task-specified sentence: *"Common temporal variation accounts for a substantial share of regional inflation dynamics, while the available 13-region panel provides no robust evidence of spatial predictive dependence after accounting for regional dynamics."* | Direct instruction to use this as the preferred central conclusion |
| Conclusion, "no spatial dependence exists" paragraph | "...after accounting for common national shocks..." / closing "...are what the data support as the dominant structure" | "...after accounting for common temporal variation..." / added explicit power-analysis caveat sentence | Same terminology fix; adds the absence-of-evidence-vs-evidence-of-absence point required by Task 2 |

No number, coefficient, p-value, classification (B/T2), or Stage 14–18 result was altered — only the prose framing around already-locked numbers.

## 2. N=13 limitation language added

Added in three places, each required by Task 2:

- **Robustness and Additional Analyses** (new subsection): a concise, main-text-appropriate summary of the power analysis and its headline implication.
- **Limitations**: the N=13 bullet was substantially rewritten to state explicitly that T is reasonably long (154 months), that N=13 limits *power* (not validity), that this is an absence-of-evidence/evidence-of-absence distinction, and that the observed CI remains wide enough to admit small effects.
- **Conclusion**: one sentence added tying the final classification to the power-analysis caveat.
- **Abstract**: one sentence added (justified because the power analysis returned classification P3), per the task's explicit conditional instruction ("if justified by the final power results, add one concise sentence to the Abstract").

## 3. Simulation assumptions (full detail in `scripts/20_power_mde_presubmission.R` header and `results_v2/presubmission_refinement/power_mde/simulation_design.csv`)

- **Estimand:** cumulative Theta_3 = theta_1+theta_2+theta_3.
- **Effect-size grid (§5):** {0, ±0.025, ±0.05, ±0.075, ±0.10, ±0.15, ±0.20} — pre-specified before any simulation, not adjusted afterward (see §11 below for why it was not widened despite MDE_80 falling outside it).
- **Lag-shape (§6):** primary = equal allocation (theta_k = Theta_3/3); secondary sensitivity = front-loaded (0.6, 0.3, 0.1) weights on Theta_3.
- **Residual-resampling method (§7):** moving block bootstrap, block length 6 months primary (≈ T^(1/3) = 5.36 for T=154, a standard rule-of-thumb — not tuned), 12 months as a documented sensitivity check on a reduced grid; blocks are resampled as full 13-region rows (preserving observed cross-sectional dependence within a block) and concatenated across time (preserving short-range temporal dependence).
- **Random seed (§8):** 20260921 (the seed used consistently across every MAKANI stage), set once at the top of the script — no per-replication reseeding.
- **Monte Carlo repetitions (§9):** 2,000 per point (primary grid, 13 points), 1,000 per point (front-loaded sensitivity, 13 points), 500 per point (block-length sensitivity, 5-point reduced grid), 100 outer × B=199 inner (DK-vs-WCB validation, 3 points).
- **MDE calculation (§10):** linear interpolation of power-vs-|Theta_3| (folding the signed grid into magnitudes by averaging the ± pair at each point) to find the smallest |Theta_3| reaching 80%/90% power; an analytic single-parameter z-approximation using the locked observed Driscoll-Kraay SE is reported alongside as a sanity check only, never substituted for the simulation.
- **Sensitivity analyses (§11):** front-loaded lag-shape (different profile, same grid); L=12 block length (reduced grid); DK-vs-WCB validation (inference-rule check).

## 4. Limitations of the power analysis (as stated in the paper and the script)

Fixed-regressor design (not a fully recursive dynamic re-simulation); moving-block-bootstrap residuals approximate but do not perfectly reproduce the true dependence structure; the DK-vs-WCB validation used a modest sample (100×3); MDE values are interpolated, not separately simulated at every possible magnitude; only the locked K=3/Queen design was assessed, not K=6 or alternative weights. All stated explicitly in `POWER_MDE_REPORT.md` and Appendix B of `paper.qmd`.

## 5. Future Research changes

The paper previously had **no** Future Research section. One was added between Limitations and Conclusion, structured as the four sequential steps specified in the task: (1) finer geography, (2) category-level inflation, (3) economic connectivity beyond geographic contiguity, (4) heterogeneous interactions (Z_i-conditioned spatial effects). Each step is explicitly framed as prospective and hypothesis-generating, not demonstrated; Step 1 explicitly does not claim sub-regional data currently exist; Step 2 explicitly references the already-documented region×category CPI data limitation; Step 4 explicitly states the interaction is NOT added to the current empirical model.

## 6. Confirmation: no new empirical model added to the main study

`scripts/20_power_mde_presubmission.R` estimates only variants of the ALREADY-LOCKED Stage-14 K=3 model (the restricted own-lag-only model, used solely to generate fitted values and residuals for the simulation DGP; and the same full K=3 spatial model, re-fit only on *simulated* data for power estimation). No interaction term, no nonlinear term, no SAR/SEM/SDM/threshold/regime-switching specification, and no new model family was introduced. `paper.qmd` was edited for prose/framing and to add the power-analysis appendix and Future Research section; no table, figure, or number from Stages 14–18 was regenerated or altered by this stage.

## 7. Scientific integrity check (pre-commit)

| Check | Result |
|---|---|
| No existing empirical estimate changed | PASS — `git diff` against all Stage 9–19 output directories is empty (verified below) |
| Final Stage-14 K=3 model remains locked | PASS — unchanged formula, unchanged coefficients (Table 3 reads from the same `fr_K_final_result_table.csv`) |
| Spatial classification remains B | PASS — read from the same `fr_J_final_classifications.txt`, unmodified |
| Temporal classification remains T2 | PASS — same file, unmodified |
| No causal language introduced | PASS — power-analysis prose describes detectability, not causal effects; "no causal identification" limitation retained verbatim |
| Month FE not called causal shocks | PASS — every instance rewritten to "common temporal variation," with explicit disclaimers that no specific mechanism is identified |
| Absence of evidence not presented as proof of zero effect | PASS — explicit distinction added in Discussion/Limitations/Conclusion/Abstract |
| Power analysis assumptions explicit | PASS — full pre-specification in script header, `simulation_design.csv`, and Appendix B |
| Simulated results clearly distinguished from observed results | PASS — every mention of simulated power/MDE is labelled as simulated; the observed Theta_3 estimate and CI (Stage 18, locked) are labelled "observed" throughout, including in Figure B3 |
| No cherry-picking of grid/lag-shape/block-length/spec | PASS — grid, shapes, and block lengths were fixed in the script header before the first replication ran; MDE_80 falling outside the grid was reported as such, not used to justify widening the grid post hoc |
| Future Research claims clearly prospective | PASS — see §5 above |
| All tables/figures referenced exist | PASS — verified in rendered `paper.html` (Tables B1–B4, Figures B1–B3, plus all main-text Tables 1–5/Figures 1–5 and Appendix A tables/figures) |
| Paper renders without errors | PASS — HTML and DOCX both exit 0 with zero errors (see §9) |

## 8. Baseline / prior-stage integrity

- `baseline/v1.0/MANIFEST.sha256`: 37/38 entries verify OK, with `paper.qmd` FAILED — this is the **same documented, deliberate exception first recorded in Stage 19** (`paper.qmd` is checksummed as part of the frozen v1.0 snapshot but is also the project's living manuscript, which Stage 19 and this stage are both explicitly instructed to update). `baseline/v1.0/MANIFEST.sha256`, `MANIFEST.md`, and `METHODOLOGY_NOTE.md` themselves were **not** touched. `paper.docx` is likewise expected to differ from its frozen-snapshot hash for the same reason.
- `git diff --stat` against every prior stage's script/data/results directory (scripts/01–20 other than 20 itself, `data_v2/`, `tests/`, `research_notes/`, `baseline/`, and every `results_v2/*` subdirectory other than `results_v2/presubmission_refinement/`) is empty.

## 9. Render status

- HTML: success, exit 0, zero errors.
- DOCX: success, exit 0, zero errors.
- PDF: not rendered — no LaTeX engine installed in this environment (unchanged from Stage 19); per instruction, no new software was installed to enable it.

## 10. Files created/modified by this stage

- Created: `scripts/20_power_mde_presubmission.R`, `results_v2/presubmission_refinement/PRESUBMISSION_REFINEMENT_AUDIT.md` (this file), `results_v2/presubmission_refinement/power_mde/` (POWER_MDE_REPORT.md, power_curve.csv, mde_summary.csv, simulation_design.csv, simulation_validation.csv, power_classification.txt, 3 figures).
- Modified: `paper.qmd` (framing edits, N=13 limitation strengthening, new Future Research section, new power-analysis subsection and Appendix B), `paper.docx` (re-rendered).
- Not modified: any Stage 9–18 script, any `results_v2/*` output directory other than `results_v2/presubmission_refinement/`, `data_v2/`, `tests/`, `research_notes/`, `baseline/v1.0/*`.
