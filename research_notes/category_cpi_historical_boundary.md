# Regional × Category CPI: Historical Boundary Research Note

**Status:** Documentation only. No panel built, no model run.
**Branch:** `makani-v2-panel`
**Scope:** Determines the earliest date from which an official
`region × expenditure_category × month` CPI panel can be defensibly
constructed for MAKANI, and records why that date is later than the
existing headline (`region × month`) panel's start.

This note documents a *research decision*, not an analysis result. It
supersedes nothing already committed — `data_v2/saudi_cpi_panel_v2.csv`
(region × month, headline CPI, 2013–2026) is unaffected and remains the
primary MAKANI v2 panel. No downloaded GASTAT file is part of this
repository; all source files referenced below live only in the session
scratchpad that produced this note.

---

## 1. The verified July 2025 → August 2025 transition

Established by downloading and directly inspecting the two adjacent
monthly GASTAT bulletins that bracket the boundary — not by
interpolation between distant sample points.

| | July 2025 (last city-level) | August 2025 (first region-level) |
|---|---|---|
| Table title | "Index Numbers **by City** and Expenditure Category, July 2025" | "Indices **By Region** and Expenditure Category, August 2025" |
| Geographic units | 16 cities + "All Cities" (Riyadh, Makkah, Jeddah, Dammam, Medina, Taif, Alhofof, Abha, Buraydah, Tabuk, Hail, Jazzan, Najran, Baha, Skaka, Arar) | **13 administrative regions** + "All Regions" — exact string match to MAKANI's immutable `region_id` registry, zero aliasing needed |
| Divisions (COICOP level 1) | 12 | **13** |
| Groups (COICOP level 2) | 43 | 48 |
| Base year | 2018 = 100 | **2023 = 100** |
| Sheet language / format | Arabic only | Bilingual AR-EN, new numbered table scheme (1.1–5.3), matching the current (Aug 2026) format |

Three changes — base-year rebase, COICOP classification revision, and
city→region geography — landed in the **same single monthly release**.
This was not a gradual drift; it was one coordinated methodology
overhaul.

The new division introduced in the 13-division scheme is **"Insurance
and Financial Services,"** split out of the prior "Miscellaneous Goods
and Services." All other divisions are name-stable or minor renames of
an apparently unchanged concept (e.g. "Communication" → "Information and
Communication"; "Restaurants and Hotels" → "Restaurants and
Accommodation Services"), though these renames have not themselves been
individually verified as definitionally identical.

## 2. Exact conclusion

**A clean `region × category × month` panel is directly available only
from August 2025.** No official file, API, or dashboard was found that
cross-tabulates GASTAT's 13 administrative regions by expenditure
category for any date before August 2025.

## 3. Strict clean window

| | |
|---|---|
| Start | **August 2025** |
| End | **August 2026** (latest bulletin available at time of research) |
| T | **13 months** |
| N | **13 regions** |
| Eligible categories | All 13 divisions (zero harmonisation needed — single vintage throughout the window) |
| Expected observations (division level) | 13 × 13 × 13 = **2,197** |
| Expected observations (group level) | 13 × 48 × 13 = **8,112** |

## 4. Why this window is currently too short for the main MAKANI longitudinal/spatial-panel analysis

The MoM headline pipeline already built for MAKANI v2
(`scripts/11_twfe_common_shocks.R`, `scripts/12_residual_spatial_tests.R`)
runs on **157 months** of region-level data — enough to fit two-way
fixed effects with meaningful residual degrees of freedom, and enough to
run 157 independent monthly permutation Moran's I tests, apply
Benjamini-Hochberg/Holm correction across them, and still say something
about the distribution of monthly spatial dependence. A 13-month window:

- covers barely one seasonal cycle plus one extra month — no ability to
  separate seasonal effects from trend or from genuine shocks;
- cannot support a two-way fixed-effects decomposition with any useful
  residual power (13 region dummies + 13 month dummies against 13×13=169
  observations per category leaves almost no residual variation to work
  with);
- would yield at most ~13 monthly cross-sectional Moran's I tests per
  category — too few for the multiple-testing-corrected, permutation-based
  inference standard already established for the headline series;
- contains no known VAT-style national shock or other identifiable event
  to anchor a sensitivity analysis against, unlike the 2013–2026 headline
  window which spans two VAT changes and the 2015/16 energy-price reform.

Using this window for anything beyond a structural feasibility check
would risk exactly the kind of underpowered, easily-overinterpreted
result the v1.0 baseline and v2 TWFE stage have both been careful to
avoid.

## 5. Why older city-level CPI must not be self-aggregated into administrative regions

Before August 2025, GASTAT's region×category-equivalent table used
**16 sample cities**, not the 13 administrative regions, as its
geographic unit. Aggregating those cities into MAKANI's 13 regions
ourselves (e.g. by simple or population-weighted averaging of the
cities nominally located within each region) would silently manufacture
a "regional" series GASTAT never validated or published. This would:

- reintroduce, by a different route, exactly the kind of
  analyst-constructed proxy the v2 panel-build stage already refused to
  do for a national CPI series ("do not construct ... unless explicitly
  labelled as an analytical proxy");
- conflate a city's price level with its administrative region's
  price level without GASTAT's own weighting methodology, which is not
  published at the city→region level;
- produce a series that looks region-shaped but carries an undocumented,
  self-introduced aggregation error, compounding across 12+ years if
  pushed back to 2013.

Per the explicit scope of this research stage, **no such aggregation was
performed.** If a longer city-level series is ever useful on its own
terms (not as a region proxy), that would be a separate, explicitly
labelled workstream — not a silent substitute for regional data.

## 6. Evidence on chain-linked historical reconstruction, by series (not inferred from each other)

From GASTAT's own Methodology and Quality Report for Consumer Price
Index Statistics, v4.0 (§3.1, Data description):

| Series | GASTAT's own "since 2013" language | Status |
|---|---|---|
| National headline CPI | "Time series of indices by administrative region since 2013" *(this line covers regional headline, see below; national headline covered by the General Index row of the category series)* | **Explicitly documented** |
| National category CPI | "Time series of indices by expenditure category since 2013" / "...monthly rates of change... since 2013" / "...annual rates of change... since 2014" | **Explicitly documented** |
| Regional headline CPI | "Time series of indices by administrative region since 2013" / corresponding monthly (2013) and annual (2014) rate-of-change series | **Explicitly documented** |
| **Regional × category CPI** | Described only as "Indices by region and expenditure category **(for the reference month)**" — no "since 2013" or any continuous-series language anywhere in the document | **Not documented as a reconstructed historical series** |

This is deliberately **not** an inference drawn from the first three
rows applying to the fourth. Two independent lines of evidence converge
on the fourth row separately:

1. **Textual**: the methodology document's own language for this
   specific table omits the continuity claim it makes explicitly for the
   other three.
2. **Empirical**: direct inspection of the actual bulletin files (§1
   above) shows the region-level geography for this cross-tabulation
   simply did not exist in any published file before August 2025 — a
   series cannot have been chain-linked back to 2013 at a geographic
   grain that was never published before August 2025.

## 7. Inspected vintage table and SHA256 hashes

All files below were downloaded into the session scratchpad only, under
explicit authorization scoped to transition-date verification. **None
are part of this repository.**

| Reference month | Geography | Geographic units | Categories | Classification | Base year | Source file | SHA256 |
|---|---|---|---|---|---|---|---|
| Feb 2019 | City | 16 cities | 12 divisions | Pre-COICOP-2018 (12-division) | 2013=100 | `consumer_price_index_february_2019_en.pdf` | `34117b04cdf29a82bec61a3ceb6de17a6bdc150c4869308377032f83b65fa615` |
| Jul 2025 | City | 16 cities + All | 12 div / 43 grp | Pre-COICOP-2018 (12-division) | 2018=100 | `CPI_Tables-Jul_2025-AR-EN.xlsx` | `e13f5aadabd519150605057ba6edb2c983cb367aafbfbbb22a7d0106a1271153` |
| **Aug 2025** | **Region** | **13 regions + All** | **13 div / 48 grp** | **COICOP 2018 (13-division)** | **2023=100** | `CPI_Tables-Aug_2025-AR-EN.xlsx` | `c58de026aeb74800798121319781c0cf13993be9f7b154d4a9a3833d391a862c` |
| Aug 2026 | Region | 13 regions + All | 13 div / 48 grp | COICOP 2018 (13-division) | 2023=100 | `CPI Tables-Aug 2026-AR-EN.xlsx` | `d7366297d310a0be815eef427df894968076569de6d36d4eff35b35b4a671d4a` |
| n/a (reference doc) | — | — | — | — | — | `Methodology_and_Quality_Report_of_Consumer_Price_Index_EN_v4.0.pdf` | `c8c7cc89721bf2f5ea01ec68dd208aa244ee827b0dce8db204990ae59375a36a` |

Source for all files: General Authority for Statistics (GASTAT),
`stats.gov.sa`, Consumer Prices Statistics → Publications.

## 8. Implication for MAKANI

**Category-specific regional analysis is deferred from the primary
longitudinal design, not discarded permanently.** The existing
`region × month` headline panel (`data_v2/saudi_cpi_panel_v2.csv`,
2013–2026, 157 usable MoM months) remains MAKANI v2's primary dataset
and is unaffected by this finding. A `region × category × month` panel
is real and buildable starting August 2025, but is not yet long enough
to support the TWFE / monthly-Moran / multiple-testing-corrected
pipeline already validated on the headline series. Revisit once enough
additional months have accumulated to give the category panel
comparable statistical power — no fixed re-check date is set here, since
that depends on how many months are judged sufficient when the question
is next revisited, not on a calendar deadline.

## 9. Theoretically relevant candidates for future category analysis

Once the category panel is long enough to analyze, these five remain
the strongest first candidates — chosen on data-quality and economic-
mechanism grounds in the prior feasibility audit, before any spatial
result was inspected, and unaffected by the boundary finding here since
all five are name-stable divisions present on both sides of the
Jul/Aug 2025 transition:

- **Food and Beverages**
- **Clothing and Footwear**
- **Housing, Water, Electricity, Gas and Other Fuels**
- **Transport**
- **Restaurants and Accommodation Services**

Within the August-2025-onward window these five (like all 13 divisions)
require zero classification harmonisation — the constraint on using them
sooner is statistical power (short T), not category comparability.
