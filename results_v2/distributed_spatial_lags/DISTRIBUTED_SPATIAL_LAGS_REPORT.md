# MAKANI v2 -- Distributed Temporal and Spatial Lags

Stage 3 of the MAKANI v2 econometric sequence, building on the hardened
lag-1 stage (`scripts/13_lagged_spatial_transmission.R`, commit `f7f5d5d`).
Tests whether a **distributed** block of lagged neighbour inflation
(spatial lags 1..K) carries additional predictive content beyond own
inflation dynamics, region fixed effects, and common month effects.

**Language discipline (per instruction): results below are described as
"evidence consistent with lagged spatial dependence" or its absence --
never as causal transmission, contagion, or structural spillovers.**

## Estimation samples

| Sample | N | T | Obs | Date range |
|---|---|---|---|---|
| K=3 common (M0-M4) | 13 | 154 | 2002 | 2013-05-01 to 2026-02-01 |
| K=6 robustness | 13 | 151 | 1963 | 2013-08-01 to 2026-02-01 |

All 13 construction-validation checks passed (no look-ahead, correct
region ordering in W, exactly one lag value per region-month per order,
no region-ID mixing, structural missingness exactly as expected for
k=1, 3, 6). See `dsl_construction_validation.csv`.

## M0-M4 model comparison (common K=3 sample)

See `dsl_model_comparison_M0_M4.csv` for the full table (N, params, R2,
adjusted R2, residual SE, AIC, BIC, within-R2 of the cumulative lag
block for each model). Incremental fit from adding the K=3 spatial-lag
block (M3 -> M4): see `dsl_incremental_spatial_block.csv` --
delta R2 = 0.000248, delta AIC = 5.09 (higher/worse), delta BIC = 21.90 (higher/worse), F-test p = 0.842.

## Primary K=3 coefficients (Driscoll-Kraay, M4)

See `dsl_k3_coefficients_dk.csv` for phi_1..3 and theta_1..3 with DK SE,
p-value, and 95% CI.

**Joint tests (DK asymptotic Wald):**
- H0: theta_1=theta_2=theta_3=0 -- stat=3.540, df=3, p=0.3156
- H0: phi_1=phi_2=phi_3=0 -- stat=5.096, df=3, p=0.1649

**Wild cluster bootstrap (B=9999, seed=20260921), joint spatial-lag test: p = 0.5754.**

## Cumulative effects

- Theta_3 = sum(theta_k) = -0.0098, DK SE = 0.0857, DK p = 0.9087, 95% CI [-0.1778, 0.1582]; WCB p = 0.8979
- Phi_3 = sum(phi_k) = -0.1323, DK SE = 0.1267, DK p = 0.2964, 95% CI [-0.3807, 0.1160]
- Theta_6 (K=6 robustness) = -0.0378, DK p = 0.7284; WCB p (B=1999) = 0.8259

## Spatial-weight robustness (K=3)

See `dsl_spatial_weight_robustness.csv` and `dsl_cumulative_Theta3_by_weight.csv`.
Sign of cumulative Theta_3 stable across Queen/KNN/distance: FALSE.
Rook not treated as independent robustness (already shown identical to
Queen in this sample, scripts/12).

## Residual diagnostics

See `dsl_residual_serial_correlation_comparison.csv` and
`dsl_k3_residual_ljungbox_by_region.csv` (with BH/Holm correction).
K=3 model (M4): 6/13 regions significant at raw 5%, 5/13 survive BH, 3/13 survive Holm.

## Robustness interpretation

See `dsl_robustness_summary.csv`, `dsl_robustness_k3_vs_k6_theta1.csv`,
and `dsl_robustness_cumulative_k3_vs_k6.csv` for the full numeric basis.
No lag length was selected based on statistical significance -- K=3 was
specified as primary and K=6 as robustness before either was estimated.

## Scientific guardrails

No claim of causality, contagion, established transmission, or
structural spillovers is made anywhere in this stage's outputs. Where
evidence is discussed, the language used is "evidence consistent with
lagged spatial dependence" (or its absence), never a causal claim.

## Files

All tables: `results_v2/distributed_spatial_lags/tables/`.
All figures: `results_v2/distributed_spatial_lags/figures/`.
Script: `scripts/14_distributed_spatial_lags.R`.
