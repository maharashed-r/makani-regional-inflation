# MAKANI v2 -- Stationarity, Seasonality, and Remaining Temporal Dynamics

Stage 6 of the MAKANI v2 econometric sequence. Builds on Stage 14
(`fd7ef10`, linear K=3 reference), Stage 15 (`c63dba5`), and Stage 16
(`bc4fde0`); none is rebuilt or reinterpreted here.

**No package was installed.** tseries/urca/plm/punitroots/CADFtest are
absent from the locked renv library; ADF, KPSS, Zivot-Andrews, and
Pesaran's CIPS are implemented directly in base R with standard
published critical values -- see the script header for full rationale
and the documented precision limitation on exact p-values.

## Section A: stationarity audit

Inflation (`inflation_mom`), per region (ADF+KPSS, drift spec): 13/13 stationary (both tests agree), 0/13 inconclusive, 0/13 non-stationary (both agree).
CPI level (both drift and trend specs) and Stage-14 reference-model residuals were also tested per region -- see `st_A_stationarity_tests_by_region.csv`.
Panel-level CIPS (cross-sectionally robust): CPI=-1.344, inflation=-2.715, residuals=-2.771 (see `st_A_panel_CIPS_summary.csv`; exact Pesaran (2007) critical values for N=13 not transcribed -- documented limitation).

## Section B: break-aware stationarity sensitivity

Zivot-Andrews (Model C) applied to inflation_mom per region: 1/13 regions' stationarity conclusion CHANGES once a single break is allowed for.
This is a descriptive sensitivity check only -- Stage 16's own break-diagnostic finding (BIC-selected m=0 breaks in the national series) is not reinterpreted or treated as causal evidence here.

## Section C: seasonal lag-12

own_lag_12 coefficient (Driscoll-Kraay): estimate=-0.0116, p=0.688. Incremental R2 over S0: 0.0001. AIC/BIC: S0=4224.0/5132.8, S1=4225.7/5140.1.

## Section D: extended temporal memory

BIC-preferred specification: T1 (own lag 1-3 + 6). See `st_D_extended_temporal_model_comparison.csv` for the full T0-T3 comparison.

## Section E: residual diagnostics

Stage 17 preferred model: 7/13 raw, 5/13 BH, 2/13 Holm (vs. Stage 14: 6/13, 5/13, 3/13; Stage 15: 6/13, 4/13, 2/13; Stage 16: 5/13, 3/13, 3/13). A one-region change is not treated as substantive on its own.

## Section F: spatial-conclusion check

Joint spatial-lag test under the preferred Stage-17 temporal specification: p=0.369 (vs. p=0.427 under the Stage-14-equivalent model on the same common sample). Conclusion DOES NOT CHANGE -- Stage 14's weak/no robust spatial predictive dependence finding still holds.

## Section G: robustness

Alternative (shorter) ADF lag-selection rule changes the stationarity conclusion for 0/13 regions. Extreme-observation exclusion leaves the spatial joint test's conclusion materially unchanged (see `st_G_robustness_extreme_obs.csv`). K=6 spatial horizon was not re-run -- not required for comparability since the preferred model retains Stage 14's original K=3 spatial block.

## Section H: evidence classification

**B: Inflation broadly stationary but modest residual temporal structure remains**

## Files

All tables: `results_v2/stationarity_temporal_dynamics/tables/`.
All figures: `results_v2/stationarity_temporal_dynamics/figures/`.
Script: `scripts/17_stationarity_temporal_dynamics.R`.
