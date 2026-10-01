# MAKANI v2 -- Final Robustness Report

Stage 7 (final) of the MAKANI v2 econometric sequence. Consolidates
Stages 14-17 (`fd7ef10`, `c63dba5`, `bc4fde0`, `e0566f5`) into one final,
robustness-checked specification. No new significance search, no new
model family.

## Section A: final reference specification

Stage-14 K=3 model, reproduced on N=13, T=154, obs=2002. Joint spatial test p=0.316; cumulative Theta_3=-0.0098 (95% CI [-0.1778, 0.1582]).

## Section B: leave-one-region-out

13/13 omissions leave the conclusion unchanged. Theta_3 range [-0.0670, 0.0347]; joint p range [0.077, 0.797].
The result is not driven by any single region.

## Section C: influential-month sensitivity

All exclusion specs (single most influential month, top 3 months, top 1% observations) leave spatial joint inference non-significant -- see `fr_C_influential_month_sensitivity.csv`.

## Section D: spatial-weight robustness

Queen p=0.316, KNN p=0.985, Distance p=0.895. Conclusion does not depend on W.

## Section E: dynamic-horizon robustness

K=1 p=0.911, K=3 p=0.320, K=6 p=0.408 (common sample). No materially different conclusion across horizons.

## Section F: temporal-specification robustness

T0 p=0.320, T1 (+lag 6) p=0.271 -- confirms Stage 17's finding: the spatial conclusion is unchanged.

## Section G: inference robustness

Driscoll-Kraay p=0.316, region-clustered p=0.517, classical (diagnostic) p=0.842, wild-cluster-bootstrap (B=9999) p=0.5754. All four agree.

## Section H: consolidated hypothesis-family summary

See `fr_H_consolidated_hypothesis_families.csv` -- ten major test families across Stages 12-17, none surviving multiple-testing correction except the trivially-expected month-FE-absorbed contemporaneous dependence check.

## Sections J-K: final classifications and result table

**Spatial: B: Weak/suggestive but non-robust spatial dependence**
**Temporal: T2: Moderate but non-decisive temporal persistence**
See `fr_K_final_result_table.csv` and `fr_K_robustness_summary.csv`.

## Files

All tables: `results_v2/final_robustness/tables/`. All figures: `results_v2/final_robustness/figures/`.
Script: `scripts/18_final_robustness_model_lock.R`.
