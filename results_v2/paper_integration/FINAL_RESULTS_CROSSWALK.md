# MAKANI v2 — Stage 19: Final Results Crosswalk

Every quantitative claim in `paper.qmd`'s main text (Abstract through
Reproducibility) is read programmatically from one of the files below at
render time — none is hand-typed. This table lets any reader verify a
specific number in the manuscript against its exact source file and, for
each stage, its commit hash.

| Manuscript location | Number(s) | Source file | Column(s) | Stage / commit |
|---|---|---|---|---|
| Abstract; Table 2; §5.1 | Month FE share = 42.9%; Region FE share = 0.45% | `results_v2/tables/twfe_variance_decomposition_mom.csv` | `share_of_total` | Stage 11 / `e626dee` |
| Abstract; Table 3; §6 | Final K=3 coefficients, CIs, cumulative effects | `results_v2/final_robustness/tables/fr_K_final_result_table.csv` | all | Stage 18 / `5565b1a` |
| Abstract; Table 3; §6 | N=13, T=154, obs=2002, R², AIC, BIC, WCB p=0.575 | `results_v2/final_robustness/tables/fr_K_final_model_metadata.csv` | all | Stage 18 / `5565b1a` |
| §5.2 | Monthly Moran's I: 11/157 raw-sig, 0 BH, 0 Holm | `results_v2/tables/moran_multiple_testing_summary.csv` | `n_raw_p_lt_05`, `n_surviving_BH`, `n_surviving_Holm` | Stage 12 / `e626dee` |
| §5.2 | Lag-1 spatial coefficient, DK SE/p | `results_v2/tables/lag_inference_phi_theta.csv` | row: `method`="Driscoll-Kraay panel-corrected (PRIMARY)", `coef`="w_lag_1" | Stage 13 / `f7f5d5d` |
| §5.2 | Lag-1 wild-cluster-bootstrap p | `results_v2/tables/lag_wild_cluster_bootstrap.csv` | `p_value_wcb`, `n_bootstrap_reps` | Stage 13 / `f7f5d5d` |
| §5.3; Table 4 | K=3 joint theta test and cumulative Theta_3 | `results_v2/distributed_spatial_lags/tables/dsl_k3_wild_bootstrap_joint_and_cumulative.csv` | `statistic`, `p_value_wcb` | Stage 14 / `fd7ef10` |
| §5.3; Table 4 | K=1/K=3/K=6 common-sample comparison | `results_v2/final_robustness/tables/fr_E_dynamic_horizon_robustness.csv` | all | Stage 18 / `5565b1a` |
| §5.4 | Nonlinear 5-test family (raw/BH/Holm p) | `results_v2/nonlinear_spatial_dynamics/tables/nl_G_multiple_testing_family.csv` | `raw_p_value`, `p_BH`, `p_Holm` | Stage 15 / `c63dba5` |
| §5.5 | Structural-shock 9-test family (raw/BH sig counts) | `results_v2/structural_shocks/tables/ss_I_multiple_testing_family.csv` | `sig_raw_05`, `sig_BH_05` | Stage 16 / `bc4fde0` |
| §5.6 | Inflation/CPI stationarity region counts | `results_v2/stationarity_temporal_dynamics/tables/st_A_stationarity_summary_counts.csv` | `n_regions`, `conclusion` | Stage 17 / `e0566f5` |
| §5.6 | Seasonal lag-12 coefficient/p | `results_v2/stationarity_temporal_dynamics/tables/st_C_seasonal_lag12_coefficient.csv` | `estimate`, `p_value` | Stage 17 / `e0566f5` |
| §6; Table 5 | Leave-one-region-out p range | `results_v2/final_robustness/tables/fr_B_leave_one_region_out_summary.csv` | `min_p`, `max_p` | Stage 18 / `5565b1a` |
| §6; Table 5 | Spatial-weight robustness (Queen/KNN/distance) | `results_v2/final_robustness/tables/fr_D_spatial_weight_robustness.csv` | all | Stage 18 / `5565b1a` |
| §6; Table 5 | Final robustness summary (all 5 checks) | `results_v2/final_robustness/tables/fr_K_robustness_summary.csv` | all | Stage 18 / `5565b1a` |
| §6 | Inference-method robustness (DK/cluster/classical/WCB) | `results_v2/final_robustness/tables/fr_G_inference_robustness.csv` | `p_value` | Stage 18 / `5565b1a` |
| §9 Conclusion | Final spatial classification (B); final temporal classification (T2) | `results_v2/final_robustness/tables/fr_J_final_classifications.txt` | line 1; line 2 | Stage 18 / `5565b1a` |
| Table 1 | Panel structure (N, T, obs, date range, NA counts) | `data_v2/saudi_cpi_panel_v2.csv`, `data_v2/region_crosswalk.csv` | computed live from raw panel | Stage 9 / `9b4bc89` |
| Figure 1 | Regional inflation time series | `results_v2/figures/fig_A_regional_inflation_mom_over_time.png` | — | Stage 12 / `e626dee` |
| Figure 2 | Month fixed effects over time | `results_v2/figures/fig_B_month_fe_tau_over_time.png` | — | Stage 11 / `e626dee` |
| Figure 3 | Monthly Moran's I, raw vs. residual | `results_v2/figures/fig_E_raw_vs_residual_moran_comparison.png` | — | Stage 12 / `e626dee` |
| Figure 4 | Distributed spatial-lag coefficients (K=3) | `results_v2/distributed_spatial_lags/figures/fig_G1_k3_theta_by_lag.png` | — | Stage 14 / `fd7ef10` |
| Figure 5 | Leave-one-region-out forest plot | `results_v2/final_robustness/figures/fig_1_loro_theta3_forest_plot.png` | — | Stage 18 / `5565b1a` |

## Appendix A (v1) crosswalk

All Appendix A numbers are read from `results/tables/*.csv` and
`results/figures/*.png`, produced by the frozen v1 pipeline
(`scripts/01`-`07`), unchanged since the `baseline-v1.0` tag. These are
descriptive facts about the data, not re-verified against v2 output —
their *interpretation* in the main text is what changed, not their
values.

## Verification note

This crosswalk was built by re-reading each source file listed above
during Stage 19 (commit hashes confirmed via `git log`), not from
memory of earlier stages' console output. Any future stage that
modifies one of these files without a corresponding re-render of
`paper.qmd`/`paper.html`/`paper.docx` will make the manuscript stale;
re-run `quarto render paper.qmd --to html` and `--to docx` after any
such change.
