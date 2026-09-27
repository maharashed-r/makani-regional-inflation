# MAKANI v2 -- Nonlinear and State-Dependent Regional Inflation Dynamics

Stage 4 of the MAKANI v2 econometric sequence, building on the
distributed-lag stage (`scripts/14_distributed_spatial_lags.R`, commit
`fd7ef10`) as the authoritative linear reference (K=3 own and Queen
spatial lags, region FE, month FE). That stage is not rebuilt or
reinterpreted here -- its K=3 sample is reproduced independently for
an exactly matched common estimation sample, and cross-checked against
its own saved output (match confirmed: TRUE).

**Nonlinearity was not assumed.** Every design choice below (state
variable, regime cutoff, spline df, threshold trimming/grid, composite
spatial variable) was pre-specified in the script's header before any
model was estimated -- see that header for the full rationale.

## Section A: high-inflation regime interaction

delta_1=0.1904 (p=0.098), delta_2=-0.0477 (p=0.661), delta_3=-0.0315 (p=0.720) [Driscoll-Kraay].
Joint Wald (DK), H0: delta_1=delta_2=delta_3=0: p=0.258. Wild cluster bootstrap (B=9999): p=0.770.
Cumulative spatial effect: normal-regime and high-regime estimates both
small and statistically indistinguishable from zero (see
`nl_A_cumulative_effect_by_regime.csv`).

## Section B: smooth nonlinearity (spline)

Joint Wald (DK) test of the 3 spline-interaction terms: stat=10.213, p=0.01684.
This is the only nominally significant result in this stage at the raw
5% level in its primary (Queen, K=3) specification. It does **not**
survive Section G's multiple-testing correction, does not replicate
under KNN weights, and does not replicate at K=6 (see Section F below)
-- see `nl_B_spline_coefficients.csv`, `nl_B_marginal_effect_by_state.csv`.

## Section C: quadratic diagnostic (secondary only)

Joint Wald (DK): p=0.076 (not significant at 5%). Implied turning point = 0.673 (original scale), 
within observed support: TRUE, with 143 of 2002 observations within 0.5 SD of it -- but since the 
underlying joint test is not significant, the turning point is reported for completeness, not interpreted
as an established economic extremum.

## Section D: threshold model

Estimated gamma_hat = 0.2833 (n_below=1703, n_above=299). **Caveat: gamma_hat landed at the upper edge
of the trimmed search range** -- a sign the SSR-minimizing point may be an artefact of the trim boundary
rather than a well-identified interior optimum; this is stated explicitly, not hidden.
Bootstrap existence test (B=499, 30 grid points): p=0.754.
**No statistically robust evidence of a threshold relationship.**

## Section E: asymmetry

coef_pos=0.0524, coef_neg=-0.1316, difference=0.1840, DK p=0.367, wild-cluster-bootstrap (B=999) p=0.616.
No evidence of asymmetric positive/negative spatial predictive content.

## Section F: robustness

See `nl_F_A_regime_weight_robustness.csv`, `nl_F_A_regime_definition_robustness.csv`,
`nl_F_B_spline_robustness.csv`, `nl_F_B_spline_extreme_obs_sensitivity.csv`.
Section B's nominal significance (Queen K=3, p=0.017) does NOT replicate under KNN weights (p=0.092)
or at K=6 (p=0.471); it IS present under distance-band weights (p=0.022) and is not driven by extreme
observations (excluding the top 1% Cook's-distance points leaves it similarly nominal). This is a mixed,
not-fully-robust picture.

## Section G: multiple-testing correction

Pre-specified 5-test family (A, B, C, D, E). See `nl_G_multiple_testing_family.csv` for the full table.
1/5 significant at raw 5%; 0/5 survive BH; 0/5 survive Holm.

## Section H: residual diagnostics

Stage-14 linear K=3 reference: 6/13 raw, 5/13 BH, 3/13 Holm. Stage-15 regime-interaction model: 6/13 raw, 4/13 BH, 2/13 Holm.

## Section I: evidence classification

**B: Suggestive but not robust nonlinear evidence**

Statistical nonlinearity, economic magnitude, predictive content, and
causal interpretation are explicitly distinguished: a nominal p<0.05 in
one specification (Section B, Queen K=3) is NOT robust across spatial
weights or dynamic horizon, does NOT survive multiple-testing
correction, and implies only small cumulative effect sizes even where
nominal. No causal, transmission, or contagion claim is made anywhere
in this stage.

## Files

All tables: `results_v2/nonlinear_spatial_dynamics/tables/`.
All figures: `results_v2/nonlinear_spatial_dynamics/figures/`.
Script: `scripts/15_nonlinear_spatial_dynamics.R`.
