# MAKANI v2 — Structural Shock / National Episode Registry

Stage 16 input document. Verified **before** any Stage-16 model was estimated,
per the task's requirement that no event date be hard-coded from memory or
assumed. Every date below traces to an official Saudi government source (or,
where the primary government portal was not directly retrievable in this
session, to independent corroborating reporting that itself cites the same
official government action — noted explicitly per event). No causal
inference is drawn from this registry; it documents *when* verified national
policy/administrative-price episodes occurred, nothing more.

All three verified events fall within the MAKANI v2 panel's coverage window
(January 2013 – February 2026, monthly).

---

## Monthly coding rule

GASTAT's regional CPI series (and therefore `inflation_mom`) is published at
calendar-month frequency with no intra-month timing information. An event
whose effective date falls within calendar month *M* is coded to event month
*M* in its entirety — finer-than-month timing cannot be distinguished in a
monthly index. All three verified events below have effective dates either
at the very start of a calendar month (2018-01-01, 2020-07-01) or within the
final three days of one (2015-12-29), so this rule is unambiguous in every
case; there is no genuinely mid-month effective date requiring a judgment
call.

---

## E1: Energy / Administered-Price Reform — Wave 1

| Field | Value |
|---|---|
| Event name | Energy and administered-price reform, first wave (electricity, fuel, gas, water) |
| Announcement date | 2015-12-28 |
| Effective date | 2015-12-29 (electricity, per Council of Ministers Resolution No. 95); other administered prices (gasoline, natural gas/ethane, water) changed within the same turn-of-year window, 2015-12-29 to 2016-01-01 |
| Event month (primary) | December 2015 (2015-12) |
| Economic relevance | Council of Ministers Resolution No. 95 raised electricity tariffs and, in the same policy action, raised administered fuel, gas, and water prices — the opening move of the Kingdom's five-year subsidy-reform program, undertaken in response to the 2014–2015 oil-revenue decline. Directly affects the "Housing, Water, Electricity, Gas" and "Transport" CPI divisions. |
| Source | Council of Ministers Resolution No. 95 (issued 2015-12-28), announced via the Saudi Press Agency (SPA), the Kingdom's official state news agency. Corroborated by contemporaneous international reporting (Reuters, Al Jazeera) that itself cites the SPA/Council-of-Ministers announcement. |
| Verification status | **VERIFIED** — resolution number, announcement date, and effective date are consistently corroborated across independent sources citing the same official government action. The primary Umm Al-Qura Gazette publication text was not independently retrieved in this session; this is noted as a minor sourcing limitation, not a basis for doubting the date, which is corroborated by multiple independent citations of the same Council of Ministers resolution. |
| Primary analysis window | Event month + following 2 months → Dec 2015 – Feb 2016 |
| Sensitivity windows | (1) Event month only: Dec 2015. (2) Event month + following 5 months: Dec 2015 – May 2016 |

---

## E2: Fiscal Reform Package — VAT Introduction + Energy Price Reform Wave 2 (COMBINED — see confound note)

| Field | Value |
|---|---|
| Event name | Q1-2018 fiscal reform package (VAT introduction and the second wave of energy-price reform) |
| Sub-component 1 | **VAT introduction.** Effective date: **2018-01-01**. ZATCA's official VAT Law page states: *"It came into effect on January 1, 2018, with its original version published on 4 Dhul-Qi'dah 1438 AH."* Standard rate 5%, introduced under the GCC Unified VAT Agreement. Source: ZATCA (Zakat, Tax and Customs Authority) official website, VAT Law page (`zatca.gov.sa`). |
| Sub-component 2 | **Energy price reform, second wave.** Effective date: **2018-01-01**, coinciding with the fiscal reform package. Gasoline and residential electricity prices rose sharply (premium gasoline +127%, regular gasoline +83%, per contemporaneous reporting on the Fiscal Balance Program); the Citizen's Account cash-transfer program (established Dec 2017) was introduced alongside it to offset the impact on lower-income households. Source: multiple independent reports of the Fiscal Balance Program's implementation, consistent with the government's own published Fiscal Balance Program documents. |
| Event month (primary) | January 2018 (2018-01) |
| Economic relevance | Directly affects nearly all CPI divisions (VAT is broad-based); energy-price component directly affects Transport and Housing/Utilities divisions. |
| **CONFOUND NOTE — read before interpreting any E2 result** | Both sub-components share the **identical** effective date (2018-01-01). With monthly regional CPI data, these two policy actions are **not separately identifiable** — any estimated change in regional inflation dynamics around January 2018 cannot be attributed to VAT introduction *or* energy-price reform individually. The registry documents them as two distinct verified policy actions for completeness, but the Stage-16 regression analysis treats them as **one combined event window** ("Jan-2018 Fiscal Reform Package"), consistent with the task's own guardrail against combining distinct mechanisms into one interpretation merely because dates coincide — here the dates do not merely coincide, they are literally the same, so no attempt is made to disentangle them at all. |
| Verification status | **VERIFIED** for both sub-components' effective dates (ZATCA official source for VAT; multiple independently corroborating sources for the energy-price component, consistent with the government's own Fiscal Balance Program). |
| Primary analysis window | Event month + following 2 months → Jan 2018 – Mar 2018 |
| Sensitivity windows | (1) Event month only: Jan 2018. (2) Event month + following 5 months: Jan 2018 – Jun 2018 |

---

## E3: VAT Rate Increase (5% → 15%)

| Field | Value |
|---|---|
| Event name | VAT rate increase from 5% to 15% |
| Announcement date | 2020-05-11 (Ministry of Finance press release; Minister of Finance H.E. Mohammed Al-Jadaan) |
| Effective date | **2020-07-01** |
| Event month (primary) | July 2020 (2020-07) |
| Economic relevance | Tripling of the standard VAT rate — the single largest one-time indirect-tax change in the panel's coverage window. Broad-based across nearly all CPI divisions. |
| Source | Saudi Ministry of Finance official press release (11–12 May 2020) announcing the increase effective 1 July 2020, as one of several fiscal measures addressing the COVID-19-era fiscal shortfall; ZATCA (then GAZT) transitional-provisions guidelines published 2020-06-07 confirming the 1 July 2020 effective date and specifying that the 15% rate applies to supplies made on or after that date. |
| Verification status | **VERIFIED** — announcement date, effective date, and rate change are consistently confirmed by the Ministry of Finance announcement and the tax authority's own transitional guidance. |
| **CONFOUND NOTE — read before interpreting any E3 result** | The effective date (July 2020) falls squarely within the broad COVID-19 pandemic macro-disruption period (see Section H of the Stage-16 analysis). No attempt is made in this stage to separate a VAT-increase-specific effect from general pandemic-era demand, supply, and mobility disruption — with a single 2020 national episode and no comparable non-pandemic tax-only counterfactual in the sample, this stage's data cannot support that separation, and the report states this explicitly rather than implying one mechanism over the other. |
| Primary analysis window | Event month + following 2 months → Jul 2020 – Sep 2020 |
| Sensitivity windows | (1) Event month only: Jul 2020. (2) Event month + following 5 months: Jul 2020 – Dec 2020 |

---

## Candidate events considered and excluded

No candidate event from the task's suggested list was excluded for lack of
verification — all three (energy/administered-price reform, VAT
introduction, VAT rate increase) were successfully verified against official
or officially-sourced material, as documented above. No additional
candidate events were added beyond the three suggested in the task
specification; this registry does not attempt an exhaustive search for
every possible national policy episode in 2013–2026, only the pre-specified
candidates named in the task.

## Explicit non-causal scope note

This registry records *dates of verified national policy/administrative-price
actions*. It does not, by itself, imply that regional inflation dynamics
changed around any of these dates — that is an empirical question tested
(without presupposing a positive answer) in
`scripts/16_structural_shocks.R` and reported in
`results_v2/structural_shocks/STRUCTURAL_SHOCKS_REPORT.md`.
