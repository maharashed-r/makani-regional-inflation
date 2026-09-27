# MAKANI v2 -- Simulation-Based Power / Minimum-Detectable-Effect (MDE) Analysis

**This is a design-sensitivity analysis, not a new empirical finding.** It
does not re-estimate, reinterpret, or alter any Stage 14-18 result. See
`scripts/20_power_mde_presubmission.R` for the full pre-specified design.

## Question answered

Given N=13, T=154 (the locked Stage-14 K=3 estimation sample), the
observed Queen spatial-weight structure, the locked dynamic
specification, and the observed residual dependence, what magnitude of
cumulative 3-month spatial predictive effect (Theta_3) could this design
plausibly detect at conventional power (80%/90%)?

## Design summary

- Effect-size grid (pre-specified, NOT adjusted after seeing results): -0.2, -0.15, -0.1, -0.075, -0.05, -0.025, 0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2
- Primary lag-shape: equal allocation across lags 1-3. Secondary (sensitivity): front-loaded (0.6/0.3/0.1).
- Fixed-regressor DGP with moving-block-bootstrapped residuals (block length 6 months primary, 12 months sensitivity).
- Primary inference: Driscoll-Kraay joint Wald test (same as Stage 14/18). Validated against a B=199 wild-cluster bootstrap at 3 effect sizes.
- Replications: 2000 (primary), 1000 (front-loaded sensitivity), 500 (block-length sensitivity), 100 x WCB(B=199) (validation). Seed=20260921, set once.

## Validation: DK asymptotic vs. nested wild-cluster bootstrap

See `simulation_validation.csv`. The DK-based rejection rate is used as
the primary simulation inference rule throughout the full grid (running
a full B=9999 WCB inside thousands of replications is not computationally
tractable); this validation checks that DK-based power tracks WCB-based
power reasonably at 3 selected effect sizes, on a reduced replication
count -- it is an approximation check, not a claim that the full-grid
power curve IS wild-bootstrap power.

## Power curve

See `power_curve.csv` and Figures 1-2.
Maximum power achieved anywhere within the pre-specified grid (|Theta_3| <= 0.20): 0.577.

## MDE

See `mde_summary.csv`.
- MDE_80 (simulated): NOT reached within the pre-specified grid (|Theta_3| <= 0.20)
- MDE_90 (simulated): NOT reached within the pre-specified grid
- MDE_80 (analytic sanity check, single-parameter z-approximation from the locked observed DK SE = 0.0857): 0.2402
- MDE_90 (analytic sanity check): 0.2779
The analytic approximation is a sanity check only and is not substituted for the simulation result.

## Power classification

**P3: Low power even for moderate effects**

This design cannot reliably distinguish a true zero cumulative spatial effect from a small-to-moderate one within the pre-specified grid. Small spatial effects (well within the range of what would be economically unsurprising) CANNOT be ruled out by this study's null/weak finding -- absence of evidence is not evidence of absence here. Only effects at or above the simulated MDE (or, if MDE_80 was not reached within the grid, effects larger than the largest grid point tested, |Theta_3|=0.20) would have been reliably detected.

## Limitations of this power analysis

- Fixed-regressor design: own-lag and spatial-lag regressor VALUES are held at their observed values rather than being regenerated recursively from a fully simulated dynamic panel (which would itself require assuming a data-generating process for how a true spatial effect would alter neighbours' inflation period by period). This is a standard simplification for power analysis, stated explicitly.
- Moving-block-bootstrap residuals approximate, but do not perfectly reproduce, the true (unknown) joint temporal/cross-sectional dependence structure of the actual error process.
- The DK-based rejection rule is validated against WCB at only 3 effect sizes and 100 replications each -- a modest validation sample, not an exhaustive one; Monte Carlo standard errors for that validation are reported in `simulation_validation.csv`.
- Power is estimated only at the pre-specified grid points; values between them are linearly interpolated for MDE_80/MDE_90, not separately simulated.
- This analysis assumes the SAME lag horizon (K=3) and spatial-weight structure (Queen) as the locked design; it does not assess power for K=6, alternative weights, or nonlinear specifications, which were not reopened per this stage's guardrails.

## Files

`power_curve.csv`, `mde_summary.csv`, `simulation_design.csv`, `simulation_validation.csv`, `power_classification.txt`, `fig_1_power_curve.png`, `fig_2_power_by_lag_shape.png`, `fig_3_observed_vs_detectable_scale.png`.
Script: `scripts/20_power_mde_presubmission.R`.
