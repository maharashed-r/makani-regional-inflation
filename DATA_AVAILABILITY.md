# Data Availability

This repository **bundles all data it uses**. Nothing is excluded or
withheld — both raw inputs are small, official, and openly licensed.

## Raw inputs (`data_raw/`)

| File | Content | Source | License / terms |
|---|---|---|---|
| `saudi_cpi_regions.csv` | Monthly Consumer Price Index by region, wide format, base year 2023 = 100, January 2013 – February 2026 | General Authority for Statistics (GASTAT), Kingdom of Saudi Arabia (<https://www.stats.gov.sa/en/>) | Published under Saudi Arabia's national Open Data Policy (Saudi Open Data Portal, <https://od.data.gov.sa/en/policies/open-data>), which permits sharing, reuse, and adaptation with attribution and without distorting the data or its source. This repository attributes GASTAT as the source throughout (`paper.qmd`, `README.md`, script headers). |
| `saudi_states_sf.rds` | Cached region-boundary polygons for Saudi Arabia's 13 administrative regions | [Natural Earth](https://www.naturalearthdata.com/) (via the R `rnaturalearth` package), cached locally to avoid a network fetch on every pipeline run | Natural Earth data is in the public domain — no attribution legally required, though Natural Earth is credited here as a courtesy. |

## Derived data

Everything else under `data_clean/`, `data_v2/`, and `results_v2/*/tables/`
is **derived deterministically from the two files above** by the scripts in
`scripts/` — none of it is a second, independent data source, and all of it
is regenerable by re-running the pipeline (see `REPRODUCIBILITY.md`).

## What is not restricted, and why this note exists anyway

Neither input file contains personal data, confidential business
information, or any indicator of a restrictive license. This note exists to
make the provenance and license basis explicit for readers and reviewers,
consistent with standard reproducibility-repository practice, not because
any material was withheld.

## If you update the CPI series

GASTAT periodically republishes and revises regional CPI series (see the
methodology note in `data_v2/README.md` regarding the August 2025
methodology update that this project's data already reflects). If you
substitute a newer `data_raw/saudi_cpi_regions.csv`, re-run the full
pipeline (`REPRODUCIBILITY.md`) rather than mixing outputs from different
vintages of the raw file.
