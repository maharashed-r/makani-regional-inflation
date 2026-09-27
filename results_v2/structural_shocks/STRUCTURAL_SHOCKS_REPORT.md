# MAKANI v2 -- Structural Stability and Pre-Specified National Shock Episodes

Stage 5 of the MAKANI v2 econometric sequence, building on the
distributed-lag stage (`scripts/14_distributed_spatial_lags.R`, commit
`fd7ef10`) as the linear reference and the nonlinear stage
(`scripts/15_nonlinear_spatial_dynamics.R`, commit `c63dba5`). Neither is
rebuilt or reinterpreted here. **A structural break was not assumed** --
every test below can, and mostly does, return a null result.

## Section A: verified event registry

Full sourcing in `research_notes/structural_shock_registry.md`. Three
events verified against official Saudi government sources (or reporting
that itself cites one): energy price reform wave 1 (effective
2015-12-29, Council of Ministers Resolution No. 95); the Q1-2018 fiscal
reform package (VAT introduction + energy price reform wave 2, BOTH
effective 2018-01-01 -- combined into one event window since the two are
not separately identifiable in monthly data); and the VAT rate increase
from 5% to 15% (effective 2020-07-01, overlapping the COVID-19 pandemic
period).

## Sections D/E/F: per-event stability tests (primary window, Queen, K=3)

See `ss_D_wald_wcb_by_event.csv`, `ss_E_wald_wcb_by_event.csv`, `ss_F_combined_wald_wcb_by_event.csv`.
- E1 (E1: Energy Price Reform (Wave 1)): D (spatial) WCB p=0.976; E (own) WCB p=0.307; F (combined) WCB p=0.093; descriptive source label: neither.
- E2 (E2: Fiscal Reform Package (VAT Introduction + Energy Price Reform Wave 2)): D (spatial) WCB p=0.970; E (own) WCB p=0.404; F (combined) WCB p=0.569; descriptive source label: neither.
- E3 (E3: VAT Rate Increase (5% -> 15%)): D (spatial) WCB p=0.646; E (own) WCB p=0.953; F (combined) WCB p=0.965; descriptive source label: neither.

## Section G: data-driven break diagnostic (descriptive only)

Base-R dynamic-programming replica of the Bai-Perron global L2-optimal partition (min segment=24 months, max breaks=5, BIC selection) applied to the national inflation series -- BIC-selected number of breaks: m=0.
Data-detected breaks are reported purely descriptively and are NOT automatically attributed to any named policy; see `ss_G_detected_breaks.csv` for descriptive proximity (in months) to the nearest verified registry event.

## Section H: COVID-19 period

The VAT-rate-increase event (E3, Jul 2020) falls entirely within the pre-specified broad pandemic-disruption window (2020-03-01 to 2021-12-31). No attempt is made to separate a VAT-increase-specific effect from general pandemic-era disruption in this stage.

## Section I: multiple-testing correction

Pre-specified 9-test family (3 events x D/E/F, each using its designated wild-cluster-bootstrap primary p-value). See `ss_I_multiple_testing_family.csv`.
0/9 significant at raw 5%; 0/9 survive BH; 0/9 survive Holm.

## Section J: robustness

See `ss_J_robustness_weights_by_event.csv`, `ss_J_robustness_k6_by_event.csv`, `ss_J_robustness_window_length_by_event.csv`.

## Section K: residual diagnostics

Stage 14 reference: 6/13 raw, 5/13 BH, 3/13 Holm. Stage 15 reference: 6/13 raw, 4/13 BH, 2/13 Holm. Stage 16 all-events combined model: 5/13 raw, 3/13 BH, 3/13 Holm.

## Section L: evidence classification

**A: No evidence of parameter instability around verified episodes**

Parameter instability, temporal persistence, spatial predictive dependence, and causal policy effects are explicitly distinct. No causal, transmission,
contagion, or policy-impact claim is made anywhere in this stage. Two mechanism confounds are documented and left unresolved by design: the Jan-2018
VAT-introduction/energy-price-reform-wave-2 confound (identical effective date), and the Jul-2020 VAT-rate-increase/COVID-19 confound (overlapping window).

## Files

All tables: `results_v2/structural_shocks/tables/`.
All figures: `results_v2/structural_shocks/figures/`.
Registry: `research_notes/structural_shock_registry.md`.
Script: `scripts/16_structural_shocks.R`.
