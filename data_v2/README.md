# MAKANI v2 Panel Documentation

Built by `scripts/09_panel_build.R`, validated by `scripts/10_panel_validation.R`
and `tests/testthat/test-panel-v2.R`. This is new MAKANI v2 development,
separate from the frozen `MAKANI Exploratory Baseline v1.0`
(`baseline/v1.0/`, git tag `baseline-v1.0`) — nothing here reads, sources,
or writes any v1.0 baseline file.

## Observation unit

One row per `region_id x month` — a canonical region identifier (see
below) crossed with a calendar month, distinct from the baseline's
`(date, Region)` key which uses the raw GASTAT spelling directly.

## Date coverage

2013-01-01 to 2026-02-01 (158 calendar months), identical source coverage
to the v1.0 baseline — same raw input file, `data_raw/saudi_cpi_regions.csv`.

## Number of regions

13, verified three independent ways: (a) 13 region columns in the raw CPI
file header, (b) 13 rows in the cached Natural Earth geometry file, (c) a
successful, unambiguous 1:1 match between the two (§ Geometry matching).

## Number of raw panel rows

**2,054** = 13 regions × 158 months. This is *larger* than the v1.0
baseline's `data_clean/saudi_cpi_clean.csv` (2,041 rows) because v2
**keeps** the structurally-missing first-lag row per region instead of
dropping it (see "Treatment of structural lag missingness" below).

## CPI definition

Unchanged from the source: GASTAT Consumer Price Index, base year
2023 = 100, one value per region per month, taken directly from
`data_raw/saudi_cpi_regions.csv` after excluding that file's per-year
annual-average row (`M == "-"`).

## MoM inflation formula

```
inflation_mom_it = 100 * (log(CPI_it) - log(CPI_i,t-1))
```

Matched on an explicit linear month index (`year*12 + month`) per
`region_id`, not on row position — so a hypothetical gap in the monthly
sequence could never silently misalign a lag. (No such gap exists in this
data; verified in `data_v2/panel_validation_report.txt` §A.) Expressed in
percent, consistent with the requested formula (the v1.0 baseline's
`inflation` column uses the same log-difference but left as a fraction,
not multiplied by 100 — v2's `inflation_mom`/`inflation_yoy` are
percentages; do not average the two scales together without adjusting).

## YoY inflation formula

```
inflation_yoy_it = 100 * (log(CPI_it) - log(CPI_i,t-12))
```

Same date-indexed (not row-indexed) lag-matching approach, 12 months back.

## Treatment of structural lag missingness

**Preserved, not dropped.** Every region's first month (Jan 2013) has
`inflation_mom = NA` (no prior month exists) and every region's first 12
months have `inflation_yoy = NA` (no observation 12 months earlier
exists). This differs deliberately from
`scripts/01_load_data.R`'s baseline behaviour, which filters out the
first NA row per region before writing `data_clean/saudi_cpi_clean.csv`.
Any future v2 model script must explicitly decide how to handle these NAs
(e.g. `na.omit()`, a balanced-subsample restriction, or an unbalanced
panel estimator) rather than inheriting a silent drop.

Verified counts: `inflation_mom` NA = 13 (one per region);
`inflation_yoy` NA = 156 (twelve per region × 13 regions). No other NA
appears anywhere in the panel.

## Region identifier design

`region_id` = `R01`–`R13`, **permanent and immutable**, read from
`data_v2/region_id_registry.csv` — a small, hand-frozen,
version-controlled lookup table (`region_id, canonical_region_name,
known_aliases, frozen_date, notes`). `scripts/09_panel_build.R` matches
each raw CPI column name against the registry's `canonical_region_name`
(exact match) or `known_aliases` (fallback match) and **halts the build**
if any raw name matches neither — it never derives, re-sorts, or
re-numbers an id at runtime. This is deliberate: the original version of
this script assigned `region_id` by alphabetically sorting the raw column
names, which meant the id was a function of the *current* spelling and
*current* column order in the raw source file and could in principle
shift if either changed. The registry freezes the ids that alphabetical
assignment originally produced (R01=Al Baha … R13=Tabouk — unchanged, see
`data_v2/region_id_registry.csv`), so no existing panel value or id
changed when this hardening was applied; only the assignment *mechanism*
did, from "derive" to "look up." A future spelling change, alias, or
reordering in the raw source is handled by a human extending the
registry's `known_aliases` column, not by code silently reassigning ids.

Do not confuse `data_v2/region_id_registry.csv` (the permanent identity
source of truth — minimal, stable, rarely touched) with
`data_v2/region_crosswalk.csv` (richer reference data built *from* the
registry on every run — geometry spelling, Arabic name, ISO code, legacy
aliases). The crosswalk is regenerated each build; the registry is not.

An `iso_3166_2_naturalearth` column is included in the crosswalk as
supplementary reference only (sourced directly from the cached Natural
Earth geometry file's own `iso_3166_2` attribute, not typed in by hand,
and — as the column name states — **not from GASTAT**), and should be
independently checked against an authoritative ISO source before being
relied on externally.

The full mapping — `region_id`, the canonical GASTAT English spelling,
the raw-data spelling (identical to the GASTAT spelling in this dataset),
the geometry spelling, the Arabic name, the ISO code, and any known legacy
aliases — lives in `data_v2/region_crosswalk.csv`. This is the **one**
place this correspondence is defined for v2 onward; it replaces the need
for any future v2 script to carry its own copy of a `recode()` table like
the ones in the frozen baseline's `scripts/04_spatial_analysis.R` and
`scripts/07_sensitivity_analysis.R` (both left untouched).

The Arabic region names (`region_name_ar_naturalearth`) are **not**
sourced from GASTAT — no Arabic-language column or file exists anywhere
in this repository. The column name says so explicitly, not just a note:
these values are taken from the `name_ar` attribute already present in
the cached `data_raw/saudi_states_sf.rds` (Natural Earth admin-1
dataset), which was already part of the frozen baseline's own raw inputs.
This is additionally disclosed per-row in the crosswalk's `notes` column.

## Geometry matching method

Base R only — no `sf`/`spdep` package required, since only the plain
character attribute columns (`name`, `name_ar`, `iso_3166_2`) of the
cached shapefile are read, not any spatial operation — these are the
*source* file's own column names (unchanged by v2); the crosswalk that v2
*produces* renames them to `region_name_ar_naturalearth` and
`iso_3166_2_naturalearth` for unambiguous provenance, per above.
`readRDS()` on an `sf` object does not require the `sf` package to be
loaded to access its non-geometry columns (verified in this session: an
`sf` object is a
`data.frame` subclass).

Matching is a plain 1:1 `match()` between the crosswalk's
`geometry_spelling` column (the alias map defined once in
`scripts/09_panel_build.R`, cross-checked for consistency against the
baseline's existing `recode()` tables) and the geometry file's raw `name`
column — verified programmatically, not just asserted: the build script
fails loudly (`stop()`) if any region is unmatched, if any geometry is
claimed twice, or if any geometry is left unclaimed. Result: all 13
regions matched to all 13 geometries, 1:1, with zero ambiguity (see
`data_v2/panel_validation_report.txt` §D for the full table).

## National CPI differential

**Not available, not constructed.** This repository's only raw input is
the 13-region file; no official national/aggregate Saudi CPI series exists
anywhere in `data_raw/`, `data_clean/`, `results/`, or `baseline/`.
Per instruction, an unweighted mean of the 13 regions was **not** silently
substituted as an unlabelled national proxy — `inflation_diff_mom` and
`inflation_diff_yoy` do not exist in this panel. If a genuine GASTAT
national CPI series is obtained later, it should be added as a new raw
input file (not fabricated from the regional data) before these columns
are constructed.

## Known raw-source asymmetry (informational only, not imputed)

`data_raw/saudi_cpi_regions.csv` is missing 2013's annual-average row.
Every other complete calendar year (2014–2025) has 13 rows (12 months + 1
GASTAT-provided annual average, `M=="-"`); 2013 has only 12 (no average
row). 2026 correctly has 2 rows (partial year, no average expected yet).
This is **documented here only** — the missing 2013 average row is not
manufactured, imputed, or otherwise reconstructed anywhere in this
pipeline, and it has **no effect** on `saudi_cpi_panel_v2.csv`: every
`M=="-"` row is excluded from the monthly panel regardless of which years
have one. Relevant only if someone later tries to reconcile this repo's
data against GASTAT's own published annual figures for 2013.
