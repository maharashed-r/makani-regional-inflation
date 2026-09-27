# MAKANI v2 -- Final Model Decision

## Decision

The **Stage-14 K=3 dynamic two-way-fixed-effects specification** is retained
as MAKANI's final empirical model:

```
inflation_mom ~ region_id_f + date_f + own_lag_1 + own_lag_2 + own_lag_3 + w_lag_1 + w_lag_2 + w_lag_3
```

## Rationale (pre-specified hierarchy)

1. **Theoretical coherence**: a parsimonious distributed-lag TWFE model
   with region and month fixed effects is the natural, defensible
   starting point for a regional inflation panel with an unknown, likely
   weak, spatial structure -- no stage in this sequence produced a
   theoretically motivated reason to prefer a more complex alternative.
2. **Robustness**: Section B-G above show the model's central conclusion
   (non-significant joint spatial predictive dependence) is unchanged
   across every region omission, every influential-month exclusion,
   all three spatial weight matrices, all three dynamic horizons, the
   T0/T1 temporal extension, and four different inference methods
   (including a B=9999 wild-cluster bootstrap).
3. **Parsimony**: Stage 15's nonlinear extensions, Stage 16's
   event-interaction terms, and Stage 17's seasonal/extended-lag terms
   each add parameters without a correspondingly robust gain -- Stage
   15's one nominal result fails multiple-testing correction and
   cross-specification replication; Stage 16 finds 0/9 tests surviving
   correction; Stage 17 finds no material AIC/BIC improvement (best
   delta-BIC = 2.15, versus a conventional 'very strong evidence'
   threshold of 10).
4. **Residual adequacy**: no extension tested materially reduces
   residual serial dependence relative to the K=3 reference (Stage 17's
   own preferred extension actually shows a slightly WORSE raw Ljung-Box
   count, 7/13 vs. 6/13) -- residual structure is acknowledged (Section
   J's T2 classification) but is not resolved by any tested extension,
   so adding one would not be justified by this criterion either.
5. **AIC/BIC**: used only as supporting evidence throughout, never as
   the basis for a decision on their own, per the pre-specified
   hierarchy above.

No model was chosen because it produced the smallest p-value -- every
tested extension across Stages 15-18 returned a null or non-robust
result, and the final decision reflects that consistently.

## Final classifications

- **Spatial evidence: B: Weak/suggestive but non-robust spatial dependence**
- **Temporal dynamics: T2: Moderate but non-decisive temporal persistence**

These are reported separately, per instruction, and are not combined
into a single summary label.

## What remains open (explicitly, not resolved by this stage)

- Moderate residual temporal persistence (Section J, T2) is real but not
  materially reduced by any tested extension -- future work outside
  this MAKANI v2 sequence's scope (e.g. a genuinely different dynamic
  structure, not explored here) would be needed to resolve it.
- Stage 15's one nominal (raw p=0.0168, not surviving correction) spline
  result remains a documented, non-actionable curiosity, not evidence.

## Files

`results_v2/final_robustness/FINAL_ROBUSTNESS_REPORT.md`,
`results_v2/final_robustness/tables/`, `results_v2/final_robustness/figures/`.
