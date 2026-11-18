# GEMS: The German Epidemic Microsimulation System — Paper Examples

This repository contains the code and model configurations behind the example applications in

> Ponge, J., Patzner, J., Suer, J., et al. *GEMS: The German Epidemic Microsimulation System.*

All simulations use [GEMS](https://github.com/IMMIDD/GEMS) **v1.0.2** (pinned in [Manifest.toml](Manifest.toml)).
This README documents the parameterization of every example: population, contacts, transmission, disease progression, and interventions.
All disease examples share one **reference model**, [models/DE_model_first_wave.toml](models/DE_model_first_wave.toml) ([Section 4](#4-reference-model-covid-19-wild-type-first-wave-2020)). Each example sets its own initial infections, horizon, and interventions in code, via the `Simulation` constructor or the scenario function's keyword arguments; [Section 5](#5-example-specific-parameterization) lists only these deviations.

## Contents

1. [Repository layout](#1-repository-layout)
2. [Reproducing the results](#2-reproducing-the-results)
3. [Overview: examples, code, and configuration](#3-overview-examples-code-and-configuration)
4. [Reference model (COVID-19, wild type, first wave 2020)](#4-reference-model-covid-19-wild-type-first-wave-2020)
   - [4.1 Simulation and initial infections](#41-simulation-and-initial-infections)
   - [4.2 Contacts](#42-contacts)
   - [4.3 Transmission](#43-transmission)
   - [4.4 Progression assignment](#44-progression-assignment)
   - [4.5 Transition durations](#45-transition-durations)
5. [Example-specific parameterization](#5-example-specific-parameterization)
   - [5.1 First wave 2020 (Sec. 3.1.1, Fig 10–11)](#51-first-wave-2020-sec-311-fig-1011)
   - [5.2 Scenario modeling (Sec. 3.2.1, Fig 14–15)](#52-scenario-modeling-sec-321-fig-1415)
   - [5.3 Parameter sweep (Sec. 3.2.2, Fig 16)](#53-parameter-sweep-sec-322-fig-16)
   - [5.4 Multi-criteria strategy evaluation (Sec. 3.2.3, Fig 17, Tab 2)](#54-multi-criteria-strategy-evaluation-sec-323-fig-17-tab-2)
   - [5.5 Performance benchmark (Sec. 3.3, Fig 18)](#55-performance-benchmark-sec-33-fig-18)
6. [Population analyses (Fig 4–7)](#6-population-analyses-fig-47)
7. [References](#7-references)

---

## 1. Repository layout

| Path | Content |
|------|---------|
| [models/DE_model_first_wave.toml](models/DE_model_first_wave.toml) | Reference model ([Section 4](#4-reference-model-covid-19-wild-type-first-wave-2020)) |
| [src/GEMSPaper.jl](src/GEMSPaper.jl) | Package entry point and `runall()` |
| [src/infectiousness_profile.jl](src/infectiousness_profile.jl) | `InfectiousnessProfile` transmission function |
| `src/sc_*.jl` | One file per example scenario |
| [src/first_wave/](src/first_wave/) | RKI data, imported and introduced infections, intervention phases of March 2020 (fitted contact reductions), calibration, and figures of the first wave example |
| [src/demographics.jl](src/demographics.jl), [src/settings.jl](src/settings.jl), [src/mobility.jl](src/mobility.jl), [src/contacts.jl](src/contacts.jl) | Population analyses (Fig 4–7) |
| [src/benchmarking.jl](src/benchmarking.jl) | Performance benchmark (Fig 18) |
| [calibrate_first_wave.jl](calibrate_first_wave.jl) | Calibration of the reference model (transmission scale, generation time, contact factors); not part of `runall()` |
| [inputs/](inputs/) | Calibration results read by the scenarios, county income data |
| `data/` | RKI case data, downloaded on the first run |

## 2. Reproducing the results

Requires Julia 1.12. In the repository root:

```julia
pkg> activate .
pkg> instantiate
```

```julia
using GEMSPaper
runall()                  # all analyses and scenarios, written to results/<timestamp>/
```

Each example can also be run on its own: `scenario_first_wave()`, `scenario_household_isolation()`, `scenario_sweeping()`, `scenario_testing()`, `benchmark()`, `contacts()`, `demographics()`, `mobility()`, `settinganalyses()`.

The first wave example downloads the RKI case data (about 400 MB, checked against a SHA-256 checksum) to `data/` on its first run. Its calibrated contact factors are read from `inputs/`; to recalibrate, run [calibrate_first_wave.jl](calibrate_first_wave.jl).

The parameter sweep stores its results in `sc_sweep_data.csv`; `replot_sweeps(folder)` redraws the heatmap and the variability summary from that file without re-running the 8,405 simulations.

All simulations use fixed random seeds ([Section 3](#3-overview-examples-code-and-configuration)). Bit-identical results require the same number of Julia threads.

**Machines used for the paper**

| Machine | Hardware | Julia | Threads (`JULIA_NUM_THREADS`) | Used for |
|---|---|---|---|---|
| Server | AMD EPYC 9575F (64 cores, 128 threads), 755 GB DDR5 RAM (5200 MT/s, CL52), NVMe SSD | TODO (1.12?) | TODO | TODO: all examples? performance benchmark (Sec. 3.3) |
| Desktop | Intel i9-14900K, 128 GB DDR5 RAM, 4 TB NVMe SSD | TODO (1.12?) | TODO | desktop runtime of the national model (Sec. 3.3) |

*Terminology:* a *random seed* initializes the random number generator of a simulation; infections that start an epidemic are called *initial infections* (at tick 0) or *introduced infections* (added during the run).

## 3. Overview: examples, code, and configuration

| Paper element | Script → function | Config | Population | Horizon | Runs / random seeds |
|---|---|---|---|---|---|
| Fig 4–6 (population, demographics, commuting) | `demographics.jl`, `settings.jl`, `mobility.jl` | none | DE | – | – |
| Fig 7 (contact patterns) | [contacts.jl](src/contacts.jl) → `contacts()` | reference model | DE | – | 0 |
| Sec. 3.1.1 First wave 2020 (Fig 10–11) | [sc_first_wave.jl](src/sc_first_wave.jl) → `scenario_first_wave()` | reference model | DE | 10.02.–31.05.2020 | see [5.1](#51-first-wave-2020-sec-311-fig-1011) |
| Sec. 3.1.2 RSV immunization 2025/26 (Fig 12–13) | separate repository: [IMMIDD/2025-RSV-scenario](https://github.com/IMMIDD/2025-RSV-scenario) | – | – | – | – |
| Sec. 3.2.1 Scenario modeling (Fig 14A, baseline readouts Fig 15) | [sc_household_isolation.jl](src/sc_household_isolation.jl) → `scenario_household_isolation()` | reference model | DE | 180 d | 1 per scenario / 0 |
| Fig 14B (code listing) | [sc_household_isolation_screenshot_file.jl](src/sc_household_isolation_screenshot_file.jl) | GEMS defaults | on-the-fly, 100,000 | 365 d | 1 per scenario / none (random) |
| Sec. 3.2.2 Parameter sweep (Fig 16) | [sc_sweeping.jl](src/sc_sweeping.jl) → `scenario_sweeping()` | reference model | SL | 365 d | 41 × 41 grid × 5 / 1–8405 |
| Sec. 3.2.3 Multi-criteria evaluation (Fig 17, Tab 2) | [sc_testing.jl](src/sc_testing.jl) → `scenario_testing()` | reference model | DE | 365 d | 1 per scenario / 0, 1, 2 |
| Sec. 3.3 Performance (Fig 18) | [benchmarking.jl](src/benchmarking.jl) → `benchmark()` | GEMS defaults | 16 states + DE; on-the-fly 0.1–50 M | 365 d | 11 per population (median) / 1–11 |

## 4. Reference model (COVID-19, wild type, first wave 2020)

Config: [models/DE_model_first_wave.toml](models/DE_model_first_wave.toml). The config comments contain the same sources and derivations.

**Source tags:** `[L]` literature · `[C]` calibrated · `[A]` assumption · `[D]` GEMS default · `[I]` illustrative / scenario setting

### 4.1 Simulation and initial infections

| Parameter | Value | Source |
|---|---|---|
| Time step | 1 day | [D] |
| Start / end date | 2020-02-10 / 2020-05-31 | [I] before the first imported infections of the first wave (reported from week 9 on); end of the forecast |
| Initial infections | none (`InfectedFraction`, fraction 0.0) | [I] the scenarios set their own initial infections |
| Duration | 112 days (`TimesUp`) | [I] the scenarios set their own horizon |
| Global setting | none | [D] |
| `[Population]` | n = 100,000, household size 3, office size 5, school size 100 | [D] only used for generated populations; all examples use the GEMS population models |

### 4.2 Contacts

Daily contacts ~ Poisson(λ), drawn uniformly from the members of the setting (`ContactparameterSampling`).

| Setting type | Level | λ | Source |
|---|---|---|---|
| Household | – | 3.3 | [L] Mossong et al. 2008 |
| Workplace | Office / Department / Workplace / WorkplaceSite | 1.4 / 0.7 / 0.35 / 0.35 | [L] Mossong et al. 2008 (total 2.8); [A] split 1/2, 1/4, 1/8, 1/8 |
| School | SchoolClass / SchoolYear / School / SchoolComplex | 0.9 / 0.45 / 0.225 / 0.225 | [L] Mossong et al. 2008 (total 1.8); [A] split 1/2, 1/4, 1/8, 1/8 |
| Municipality | – | 4.0 | [L] Mossong et al. 2008 |

The split weights contacts towards the most immediate level (own office or class): each level gets half of the one below, the two outermost levels the same.

### 4.3 Transmission

The `InfectiousnessProfile` transmission function ([src/infectiousness_profile.jl](src/infectiousness_profile.jl)) scales the per-contact transmission probability by an infectiousness profile over the days τ since the infector's infectiousness onset:

p(τ) = `transmission_scale` · w(τ), with w(τ) the profile mass in [τ, τ + 1). The profile mass before τ = 0 is dropped and the rest normalized to 1, so `transmission_scale` is the sum of the daily per-contact probabilities over the infection. Recovered individuals are immune.

| Parameter | Value | Source |
|---|---|---|
| Profile | shifted gamma relative to symptom onset: shape 97.1875, rate 3.71875, shift 25.625 d | [L] He et al. 2020, corrected by Ashcroft et al. 2020, Table 1 |
| Config parameters | `distribution = "Gamma"`, `parameters = [97.1875, 0.26891]` (shape, scale = 1 / rate), `offset = 21.625` | [L] adapted: offset = shift − pre-symptomatic period (25.625 − 4.0); profile mass before infectiousness onset (4 %) dropped |
| `transmission_scale` | 0.50486 | [C] fitted to R0 = 3.28 with [calibrate_first_wave.jl](calibrate_first_wave.jl) |
| Target R0 | 3.28 | [L] Liu et al. 2020 |

**Fit of the transmission scale:** population Rhineland-Palatinate (contact structure closest to Germany), 0.1 % random initial infections, 50 days per run, 10 runs per candidate; R0 = mean R_eff of ticks 1–10; linear search (candidate × target / result, at most 15 iterations) from the value in the config until |R0 − 3.28| ≤ 0.02.

**Transmission limits of GEMS:** agents transmit from infectiousness onset until recovery (= RNA clearance, [Section 4.5](#45-transition-durations)), from severeness onset on only in their household, and not while hospitalized. Isolated and quarantined individuals only have household contacts.

### 4.4 Progression assignment

Each infected individual's progression category is drawn from Categorical(p), where p is the matrix row of their age group (`AgeBasedProgressionAssignment` [D], groups `-39`, `40-59`, `60-79`, `80-`). Only `Critical` cases can die, after ICU admission: death ~ Bernoulli(p = `death_probability`).

| Input | Value | Source |
|---|---|---|
| Asymptomatic share | 0.17 | [L] Byambasuren et al. 2020 |
| IHR and IFR by age | see below | [L] adapted: Verity et al. 2020 (IHR: Table 3, IFR: Table 1), means of the 10-year bins weighted by the GEMS DE population |
| `Critical.death_probability` | 0.53 | [L] adapted: mortality of invasively ventilated patients (Karagiannidis et al. 2020) |
| `Severe` share | = IHR | [A] as many severely ill at home as in hospital |

**IHR and IFR per age group** (weight = n / sum of n in the group):

| Group | Bin | IHR % | IFR % | n (GEMS DE) | Weight |
|---|---|---|---|---|---|
| -39 | 0–9 | 0.00 | 0.00161 | 8,615,732 | 0.2289 |
| | 10–19 | 0.0408 | 0.00695 | 7,215,600 | 0.1917 |
| | 20–29 | 1.04 | 0.0309 | 10,607,181 | 0.2819 |
| | 30–39 | 3.43 | 0.0844 | 11,194,255 | 0.2975 |
| | **weighted** | **1.3212** | **0.0355** | 37,632,768 | 1 |
| 40-59 | 40–49 | 4.25 | 0.161 | 10,601,426 | 0.4434 |
| | 50–59 | 8.16 | 0.595 | 13,308,631 | 0.5566 |
| | **weighted** | **6.4264** | **0.4026** | 23,910,057 | 1 |
| 60-79 | 60–69 | 11.8 | 1.93 | 10,119,814 | 0.5716 |
| | 70–79 | 16.6 | 4.28 | 7,584,813 | 0.4284 |
| | **weighted** | **13.8564** | **2.9368** | 17,704,627 | 1 |
| 80- | 80+ | 18.4 | 7.80 | 5,316,274 | 1 |

Matrix columns: `Critical = IFR / 0.53`, `Hospitalized = IHR − Critical`, `Severe = IHR`, `Symptomatic = 1 − Asymptomatic − Severe − Hospitalized − Critical` (rounded to 6 decimals).

| Age group | Asymptomatic | Symptomatic | Severe | Hospitalized | Critical |
|---|---|---|---|---|---|
| -39 | 0.17 | 0.803576 | 0.013212 | 0.012542 | 0.000670 |
| 40-59 | 0.17 | 0.701472 | 0.064264 | 0.056668 | 0.007596 |
| 60-79 | 0.17 | 0.552872 | 0.138564 | 0.083153 | 0.055411 |
| 80- | 0.17 | 0.462000 | 0.184000 | 0.036830 | 0.147170 |

### 4.5 Transition durations

All durations are Poisson(λ) distributed, in days. GEMS adds 1 day to every exposure → infectiousness onset draw.

| Category | Transition | λ | Source |
|---|---|---|---|
| all | exposure → infectiousness onset | 1 + Poisson(0.8) | [L] adapted: incubation 5.8 (McAloon et al. 2020, mean) − pre-symptomatic period 4.0 − GEMS offset 1.0 |
| Asymptomatic | infectiousness onset → recovery | 9.5 | [L] Hu et al. 2020: communicable period of asymptomatic cases, median 9.5 d |
| Symptomatic, Severe, Hospitalized, Critical | infectiousness onset → symptom onset | 4.0 | [L] adapted: 5 % quantile of the infectiousness profile, 3.7 d before symptom onset (He et al. 2020, corrected by Ashcroft et al. 2020), rounded up |
| Symptomatic | symptom onset → recovery | 13.4 | [L] Byrne et al. 2020: symptom onset → first of two negative RT-PCR tests, mean 13.4 d |
| Severe, Hospitalized, Critical | symptom onset → severeness onset | 2.0 | [A] bedridden 2 d after symptom onset, half of the time to hospital admission |
| Severe | severeness onset → severeness offset | 8.0 | [A] same as hospital stay without ICU |
| Hospitalized | severeness onset → hospital admission | 2.0 | [L] adapted: symptom onset → hospital admission, median 4 d (Schilling et al. 2020), minus 2 d at home |
| Hospitalized | hospital admission → hospital discharge | 8.0 | [L] Schilling et al. 2020: hospital stay without ICU, median 8 d |
| Critical | severeness onset → hospital admission | 3.0 | [L] adapted: symptom onset → hospital admission of ICU patients, median 5 d (Schilling et al. 2020), minus 2 d at home |
| Critical | hospital admission → ICU admission | 0.0 | [L] Schilling et al. 2020: median 0 d |
| Critical (survives) | ICU admission → ICU discharge | 11.0 | [L] Schilling et al. 2020: ICU stay of survivors, median 11 d |
| Critical (survives) | ICU discharge → hospital discharge | 3.0 | [L] adapted: hospital stay of ICU patients 14 d − ICU stay 11 d (Schilling et al. 2020) |
| Critical (dies) | ICU admission → death | 8.0 | [L] Schilling et al. 2020: ICU stay of deceased, median 8 d |
| Hospitalized, Critical | hospital discharge → severeness offset | 3.0 | [A] discharged before the severe phase is over, severely ill at home for 3 days |
| Severe, Hospitalized, Critical | severeness offset → recovery | 4.0 | [A] 4-day convalescence between the end of severe illness and recovery |

## 5. Example-specific parameterization

Interventions use the Trigger–Strategy–Measure (TriSM) framework ([Ponge et al. 2024](https://doi.org/10.1109/WSC63780.2024.10838778)), listed as trigger → strategy (condition) → measures (offset).

### 5.1 First wave 2020 (Sec. 3.1.1, Fig 10–11)

Code: [sc_first_wave.jl](src/sc_first_wave.jl), [src/first_wave/](src/first_wave/). Config: reference model. Calibrated on the RKI cases by symptom onset; compared with the RKI cases by reporting date (also by symptom onset with `diagnostics = true`).

**Data**

| Item | Value | Source |
|---|---|---|
| RKI cases | county, reporting date, symptom onset; version of 14.10.2026 | Robert Koch-Institut 2026 (Zenodo, CC BY 4.0) |
| Unknown onsets | onsets > 21 d before the report, or later after the report than the 99 % quantile of the model's incubation period, are treated as unknown | [A] data errors / a case reported before onset is at most one incubation period before it |
| Counties | Berlin districts merged into 11000; Eisenach (16056) merged into the Wartburgkreis (16063) | data structure |
| Reporting delays | RKI onset → report delays per onset week (all of Germany, onsets 01.02.–31.05.; weeks with < 200 cases take the nearest week) | [L] adapted: Robert Koch-Institut 2026 |
| Median reporting delay | 5 d (`MEDIAN_REPORT_DELAY`, fixed in the code) | [L] adapted: median onset → report delay of the RKI cases with a known onset 01.02.–31.05.2020 (Robert Koch-Institut 2026) |
| Imports, countries of infection | RKI situation reports of 04.03., 06.03., 08.03., 22.03., 24.03., 26.03., 27.04.2020 | Robert Koch-Institut 2020 |
| Tyrol overnight stays | guests from Germany, winter 2018/19, 7 origin regions | [L] Amt der Tiroler Landesregierung 2019, Tab. 6, p. 11 |
| Origin regions of the Tyrol data | states mapped to whole regions: Nord = SH, HH, NI, HB; MD = HE, RP, SL; Ost = BB, MV, SN, ST, TH; BY, BW, NW, BE | [A] approximation: Statistik Austria assigns guests by postal code (e.g. 68/69 Mannheim/Heidelberg to MD, 89 Ulm to BY, 38 incl. parts of ST to Nord; Statistik Austria n.d.) |
| County geography | land areas and population (`EWZ`) of the GEMS county shapes; straight-line distance from the area-weighted centroid in Gauss–Krüger zone 3 (EPSG:31467) | [D] GEMS county shapes; [A] straight-line distance from the centroid |
| County income | disposable income of private households (incl. non-profit institutions) per inhabitant 2019 ([inputs/county_income_2019.csv](inputs/county_income_2019.csv)); Eisenach (16056) takes the value of the Wartburgkreis; used relative to the population-weighted national mean | [L] Statistische Ämter des Bundes und der Länder 2026, Tab. 2.4 |
| School holidays | KMK holiday calendar 2019/20 | [L] Sekretariat der Kultusministerkonferenz 2018 |

**Introduced infections** (added via `stepmod`, as `RegionalSeeds` in GEMS v1.0.2 only infects at tick 0)

| Item | Value | Source |
|---|---|---|
| Imported cases per reporting week (weeks 10–17) | 262, 2,781, 6,454, 3,814, 1,095, 267, 113, 31 | [L] RKI situation report of 27.04.2020, p. 3, Fig. 3 (exposure "Ausland"; n = 97,150 cases with a known place of exposure) |
| Imported cases, week 9 | 35 | [L] adapted: log-linear fit of the cases with exposure abroad (sum of the listed countries) on 04., 06., 08.03. (64, 130, 223; Tab. 2), evaluated at 02.03.; [A] stands for all imports reported up to the end of week 9 (fit not in the code) |
| Which cases are imports | cases with exposure abroad among the cases with a known place of exposure (63 % of the cases up to 27.04.) | [A] the cases with an unknown place of exposure (37 %) are not imputed |
| Share of Austria and Italy | 0.90 (weeks 9–12), 0.84 (weeks 13–17) | [L] adapted: weeks 9–10: Italy 200 of 223 foreign exposures on 08.03. (Tab. 2; no Austria); weeks 11–12 and 13: share in the increase of the country mentions 08.→22.03. and 22.→26.03. (Tab. 3 and text); [A] report intervals taken as reporting weeks, mentions as cases; 0.84 carried over to weeks 14–17 (the report of 27.04. bounds the share for 22.03.→27.04. at ≤ 0.814) |
| Share of Italy in the Alpine imports | 1.0 (weeks 9–10), 0.21 (11–12), 0.07 (13–17) | [L] adapted: no Austria on 08.03.; Italy's share in the increase of the Alpine mentions 08.→22.03. (995 of 4,849) and 22.03.→27.04. (474 of 6,443; Tab. 2); [A] reference date 27.04. here, 26.03. for the share of Austria and Italy |
| Import curve | infections per day (10.02.–26.04.) back-calculated from the weekly counts (reporting week = ISO week) over the time to onset (1 + Poisson(4.8) d; asymptomatic: infection + 6 d) and the RKI onset → report delays of the onset week, separately for Austria, Italy, and other countries; EMS with smoothing weights 1/4, 1/2, 1/4, flat start, 200 iterations; from 01.02. for the R(t) back-calculation | time to onset: reference model ([Section 4.5](#45-transition-durations)); delays: RKI (Data); [L] Becker et al. 1991 (method); [A] 200 iterations without a convergence criterion; no imports before 10.02. (config start; 01.02. for R(t)) or after 26.04. (end of week 17); asymptomatic imports reported with onset at infection + 6 d |
| Infections per reported import | 1 / 0.22 ≈ 4.5 | [A] imports detected like all infections; returning travellers from risk areas were probably tested more often, so the imports may be overestimated |
| Seeding of the imports | per day and source, expected count × 4.5 rounded stochastically (floor + Bernoulli draw); counties drawn with replacement by the allocation weights (computed with the county populations of the model); a random susceptible resident of the county (any age) infected; imports into counties outside the population model dropped | [A] |
| Reference points | Ischgl (47.0128 N, 10.2919 E) for Austria, Bozen (46.498 N, 11.354 E) for Italy | [A] Tyrol accounts for 2,737 of the 2,995 mentions with an Austrian region in the RKI report of 26.03. (Tab. 3), with Ischgl as its main cluster; South Tyrol/Trentino is the Italian region named most often (reports of 08.–26.03.) |
| Allocation to counties: Austria | Tyrol overnight stays per German origin region, distributed within the region by population × relative income^1.0 × (distance to Ischgl / 100 km)^−0.4, × (1 + 1.0 × school holidays of the state on the infection date) | see the following rows |
| Distance decay | 0.4 | [L] adapted: slope of ln(Tyrol overnight stays per inhabitant) on ln(population-weighted mean distance to Ischgl) over the 7 origin regions (−0.39; fit not in the code), applied within the regions; [A] not adjusted for income (−0.25 with income^1.0 removed) |
| Income elasticity | 1.0 | [A] trips proportional to income (also used for the other imports) |
| Holiday factor | 1.0 | [A] weight doubled during school holidays |
| Allocation: Italy | 70 % like Austria with Bozen as reference point (same Tyrol origin regions and holiday factor), 30 % by population × income | [A] ski share 0.7; share of South Tyrol/Trentino in the RKI reports: 58 % of the Italian exposures (08.03., Tab. 2), 85–86 % of the Italian regions named (24./26.03., Tab. 3), 21 % of all Italy mentions (22.03.); plausible range 0.2–0.86 |
| Allocation: other countries | population × relative income^1.0 | [A] |
| School holidays (mid-February–March 2020) | BY, BW 22.02.–01.03.; HH 29.02.–15.03.; RP 15.–23.02.; SL 15.–25.02.; SN, MV 08.–23.02.; BW: 50 % of the schools | [L] Sekretariat der Kultusministerkonferenz 2018 (all states except BW); [A] adjacent weekends added; [A] BW: carnival week as in BY (no state-wide winter holidays, carnival holidays set by each school), coverage 0.5 (no source, range 0–1); not included: ST, TH (10.–14.02.), movable days in other states (e.g. carnival in NRW), Easter holidays |
| Heinsberg cluster | RKI cases of Kreis Heinsberg with onset (or, if unknown, reporting date) up to 29.02. and reported by 07.03., infected 1:1 at onset − 6 d (rounded mean incubation 1 + 0.8 + 4.0), or report − 11 d (+ median reporting delay 5 d) if the onset is unknown; 180 cases seeded (170 with onset); 132 cases with onset up to 29.02. reported later are not seeded | [A] without a direct source: cut-off 29.02. (two weeks after the carnival session of 15.02., Streeck et al. 2020), reporting limit 07.03., 1:1 (not divided by the ascertainment); imports are allocated to Heinsberg in addition |

**Contact phases and testing**

| Item | Value | Source |
|---|---|---|
| Phase start dates | phase 2: 09.03. (large events cancelled), phase 3: 16.03. (schools closed), phase 4: 23.03. (contact restrictions) | [L] dates of the measures |
| Transitions between phases | linear over 3, 3, 4 days | [L] Dehning et al. 2020 |
| Contact factors (phases 1–4) | 0.570 / 0.163 / 0.165 / 0.119, multiplied with the contact parameters of all settings except households | [C] [inputs/first_wave_factors_DE_model_first_wave.csv](inputs/first_wave_factors_DE_model_first_wave.csv) |
| Care homes | households with at least 11 members, scaled like the other settings | [A] |
| School closure | all school settings closed from 16.03.2020 | [L] |
| Ascertainment | 0.22 of all infections detected (`ASCERTAINMENT_2020`); also used to scale the reported imports, not the Heinsberg cluster | [L] adapted: under-ascertainment factor of 4 to 5 in the first half of 2020 from seroprevalence studies (Robert Koch-Institut 2021, p. 1), i.e. 20–25 % of infections detected; 1 / 4.5 = 0.22; [A] the same factor for imports |
| Testing | perfect PCR test at symptom onset with probability p = (0.22 − hospitalized share) / (symptomatic share − hospitalized share), and at hospital admission (if not detected before); positives are not isolated | [A] |
| Reported cases by reporting date | detected infections, reported at symptom onset (asymptomatic: infection + 6 d) plus a delay drawn from the RKI delays of the onset week | [A] |
| Reported cases by onset | detected infections at symptom onset (asymptomatic: infection + 6 d), each kept with the share of RKI cases with a known onset (per onset week); compared with the RKI cases with a known onset | [A] |

**Calibration** ([calibrate_first_wave.jl](calibrate_first_wave.jl) → `calibrate_first_wave()`, not part of `runall()`)

| Step | Setup |
|---|---|
| 1. Transmission scale | [Section 4.3](#43-transmission) |
| 2a. Generation time of the model | Saarland, 70 days, scale 0.35, 20 extra infections per day in the Regionalverband Saarbrücken (ticks 7–27), ascertainment 1.0, random seed 2; cohort of ticks 10–22; mean 5.68 d ([inputs/generation_time_DE_model_first_wave.csv](inputs/generation_time_DE_model_first_wave.csv)) |
| 2b. Contact factors | fitted to the weekly RKI cases by onset, weeks 10–16 (up to 19.04., weight 0.5 for week 10), with the R per phase of GEMS kept inside the 95 % interval of the RKI R per phase (R(t) back-calculated from 01.02., infections 01.03.–12.04., block bootstrap with 1,000 samples, block length 7 d); surrogate fit, then damped Gauss–Newton steps, at most 7 evaluations with 2 runs each (random seeds 0, 1); runs until 17.05. |
| 2c. Band | 20 factor samples (log-normal) from the approximate covariance of the fit, σ²(AᵀA)⁻¹ with A the weighted Jacobian of the residuals at the optimum and σ² = residual sum of squares / (7 weeks − 4 factors) (random seeds 1001–1020) ([inputs/first_wave_band_factors_DE_model_first_wave.csv](inputs/first_wave_band_factors_DE_model_first_wave.csv)) |

**Evaluation** (`rki_county_metrics` in [src/first_wave/rki_comparison.jl](src/first_wave/rki_comparison.jl))

| Item | Setup |
|---|---|
| Windows | fitting period (weeks 10–16, up to 19.04.) and forecast (weeks 17–22, up to 31.05.), cumulative; single weeks 10, 11, 13, 15, 17 |
| County agreement | Spearman correlation of the RKI and GEMS incidence per county (cases per 100,000; RKI with the official, GEMS with the model population; mean of the runs) |
| Uncertainty | bootstrap over states: 2,000 samples of the 16 states drawn with replacement, each with all its counties (random seed 0); SD, 95 % interval, and two-sided p-value |
| State agreement | Spearman correlation of the incidence of the 16 states |
| Baselines | correlation of the RKI incidence with the population, the introduced infections per inhabitant, the closeness to Ischgl, and the earlier RKI incidence (forecast window: incidence of the fitting period; fitting window and single weeks: RKI incidence by reporting date up to 01.03.) |
| Gain over the introductions | Spearman of GEMS minus Spearman of the introduced infections per inhabitant, with interval and p-value from the same bootstrap samples (paired); for `:random_seeds`, this is the contribution of the simulated spread (population and contact model) |

**Parts of the scenario** (`FIRST_WAVE_PARTS`)

By default, `scenario_first_wave()` runs only the parts of the paper (`FIRST_WAVE_PAPER_PARTS`, marked *paper*) and writes the outputs of the paper (cases by reporting date, the paper figures `sc_first_wave_paper_figure_1_reported.png` and `…_2_reported.png`, the seeding summary, and the county metrics) and the files they are built from. The other parts can be selected with `parts`, e.g. `scenario_first_wave(parts = (:projections, :figures), diagnostics = true)`; `diagnostics = true` also writes the comparisons by symptom onset and the figures of all parts.

| Part | Setup | Runs / random seeds | Paper |
|---|---|---|---|
| `:final` | calibrated factors until 31.05. (six weeks of forecast after 19.04.) | 3 / 0–2 | paper |
| `:projections` | from 20.04. on, the contacts of phase 2, 3, or 4, or of phase 1 (until 17.05.); schools reopen in the projections to phase 1 and 2 | 2 per variant / 0, 1 | – |
| `:counterfactuals` | contacts of phase k kept from the start of phase k + 1 on (k = 1–3) | 1 per variant / 0 | – |
| `:band` | factor samples of the calibration, 5–95 % band, national and per state | 20 / 1001–1020 | paper |
| `:seeds_robustness` | all RKI cases infected before 09.03. as introduced infections (divided by the ascertainment, Heinsberg 1:1), no transmission before 09.03. | 3 / 0–2 | – |
| `:random_seeds` | the seeds of `:final` (imports and Heinsberg cluster), the same number per day, in counties drawn in proportion to their population (drawn anew for each random seed) | 3 / 0–2 | paper |
| `:figures` | paper figures, seeding summary, county metrics; with `diagnostics`, the figures of all parts whose results exist | – | paper |

### 5.2 Scenario modeling (Sec. 3.2.1, Fig 14–15)

Code: [sc_household_isolation.jl](src/sc_household_isolation.jl). Config: reference model.

| Parameter | Reference model | This example |
|---|---|---|
| Initial infections | none | 0.01 % of the population, random (`infected_fraction = 0.0001`) |
| Duration | 112 days | 180 days |
| Random seed | – | 0 (both scenarios) |

| Scenario | Trigger | Strategy (condition) | Measures (offset) → follow-up |
|---|---|---|---|
| Baseline | – | – | – |
| Scenario | `SymptomTrigger` | Find Household Members | `FindSettingMembers` in Household, incl. the symptomatic individual (0) → *Isolation* |
| | | Isolation | `SelfIsolation(14)` (0) |

All household members, including the symptomatic individual, isolate at home for 14 days. Each new symptomatic household member triggers the strategy again and restarts the 14 days for the whole household.

The baseline run also produces the GEMS readout figure (Fig 15, `sc_household_isolation_baseline_combined.png`): daily infections and deaths, effective reproduction number, household attack rate, infection duration, cumulative disease progressions, progression categories, and maps of the initial reproduction number and the infection fatality per county. Its ResultData object uses the default style, since the progression categories and the fatality map need the raw infections.

The Fig 14B listing ([sc_household_isolation_screenshot_file.jl](src/sc_household_isolation_screenshot_file.jl)) uses the same intervention with GEMS defaults ([Section 5.5](#55-performance-benchmark-sec-33-fig-18)) on 100,000 agents, without a fixed random seed.

### 5.3 Parameter sweep (Sec. 3.2.2, Fig 16)

Code: [sc_sweeping.jl](src/sc_sweeping.jl). Config: reference model.

To focus on the regime where self-isolation can be decisive, we assume that other measures already prevent 67 % of all transmissions, such that without isolation R_init lies slightly above 1.

| Parameter | Reference model | This example |
|---|---|---|
| Population | – | Saarland |
| Initial infections | none | 0.1 % of the population, random |
| Duration | 112 days | 365 days |
| Transmission scale | 0.50486 | × (1 − 0.67): 67 % of transmissions prevented by other measures (`transmission_reduction = 0.67`) [I] |
| Asymptomatic share (x-axis) | 0.17 | 0 to 1, step 0.025, set in all age groups; the other categories rescaled proportionally |
| Adherence (y-axis) | – | 0 to 1, step 0.025 |
| Repetitions | – | 5 per combination; random seeds 1–8405 (one per run) |
| Outcome | – | R_init per run (`r0(PostProcessor(sim))`: infections caused by the first 5 % of all infections, divided by their number) |

| Trigger | Strategy (condition) | Measures |
|---|---|---|
| `SymptomTrigger` | Self Isolation (adheres ~ Bernoulli(p = adherence), drawn from the simulation's RNG) | `SelfIsolation(14)` (0) |

The heatmap shows the mean R_init per combination. Alongside the data (`sc_sweep_data.csv`) and the heatmap (`sc_sweep_heatmap.png`), the scenario writes `sc_sweep_variability.txt`, which summarizes the variability between replicates: the per-cell sample SD across replicates, aggregated over the grid as the pooled SD (root of the mean variance), mean, 95th percentile, and maximum; the per-cell range (maximum − minimum) across replicates, aggregated as its minimum, maximum, mean, and median; the standard error of the cell means; and the number of cells whose replicates lie on both sides of R = 1.

### 5.4 Multi-criteria strategy evaluation (Sec. 3.2.3, Fig 17, Tab 2)

Code: [sc_testing.jl](src/sc_testing.jl). Config: reference model.

To focus on the regime where tracing can be decisive, we assume that other measures (for example, limits on gatherings or closed restaurants and leisure venues) already prevent half of the community contacts outside households, workplaces, and schools.

| Parameter | Reference model | This example |
|---|---|---|
| Initial infections | none | 50 infections in Gangelt (municipality key 05370008, county Heinsberg) |
| Duration | 112 days | 365 days |
| Municipality contact parameter | 4.0 | 2.0: 50 % of community contacts prevented (`municipality_contact_reduction = 0.5`) [I] |
| Random seeds | – | 0 (Isolation Only), 1 (Isolation & Testing), 2 (Isolation, Testing & Tracing) |

**Test type** "Antigen Test" [I]: positive ~ Bernoulli(p = 0.80) if infected (sensitivity), Bernoulli(p = 0.01) if not (specificity 0.99).

| Scenario | Trigger | Strategy (condition) | Measures (offset) → follow-up |
|---|---|---|---|
| Isolation Only | `SymptomTrigger` | Self Isolation | `SelfIsolation(14)` (0) |
| Isolation & Testing | `SymptomTrigger` | Self Isolation | `SelfIsolation(14)` (0); `Test(Antigen)` (0) |
| Isolation, Testing & Tracing | `SymptomTrigger` | Isolation & Testing (not yet detected) | `SelfIsolation(14)` (0); `Test(Antigen)` (0) → positive: *Contact Identification*; negative: *Cancel Isolation* |
| | | Contact Identification | `FindSettingMembers` in Household, Office (if working), SchoolClass (if student) (1 d) → *Isolate & Test Contacts* |
| | | Isolate & Test Contacts (not yet detected) | `SelfIsolation(14)` (0); `Test(Antigen)` (0) → negative: *Cancel Isolation* |
| | | Cancel Isolation | `CancelSelfIsolation` (0) |

The tests in Isolation & Testing have no follow-up: the scenario detects cases without changing the spread, so it differs from Isolation Only only by the random seed. The tracing scenario corresponds to Fig 9 of the paper.

### 5.5 Performance benchmark (Sec. 3.3, Fig 18)

GEMS default parameterization ([data/DefaultConf.toml @ v1.0.2](https://github.com/IMMIDD/GEMS/blob/v1.0.2/data/DefaultConf.toml)), also used for Fig 14B.

| Parameter | Value |
|---|---|
| Initial infections | 0.1 % random |
| Duration | 365 days |
| Transmission probability | 0.2 |
| Daily contacts | Poisson(λ = 1.0) on every setting level |
| Progression category (ages ≤14 / 15–65 / 66+) | Categorical(p): Asymptomatic 0.40 / 0.25 / 0.15; Symptomatic 0.58 / 0.60 / 0.40; Severe 0.01 / 0.11 / 0.25; Hospitalized 0.007 / 0.03 / 0.12; Critical 0.003 / 0.01 / 0.08 |
| Durations (days) | exposure → infectiousness 1 + Poisson(λ = 1); infectiousness → symptoms Poisson(λ = 1); symptom → recovery Poisson(λ = 7) (asymptomatic: infectiousness → recovery Poisson(λ = 8)); symptom → severeness onset Poisson(λ = 1); severeness duration (Severe) Poisson(λ = 7); to hospital Poisson(λ = 2); hospital stay Poisson(λ = 7); to ICU Poisson(λ = 2); ICU stay Poisson(λ = 7); ICU discharge → hospital discharge Poisson(λ = 7); ICU → death Poisson(λ = 10); hospital discharge → severeness offset Poisson(λ = 3); severeness offset → recovery Poisson(λ = 4) |
| Death (Critical) | Bernoulli(p = 0.3) |
| On-the-fly population | 100,000 agents, mean household size 3, office size 5, school size 100 |

Protocol: one warm-up run, then 11 runs per population (random seeds 1–11), medians reported. Populations: 16 German states and Germany, plus on-the-fly populations of 0.1, 1, 5, 10, 20, and 50 million agents.

## 6. Population analyses (Fig 4–7)

No disease parameters are involved.

- **Fig 4–6:** Germany (84,563,726 individuals).
- **Fig 7:** Germany, contacts of the reference model.

## 7. References

- Amt der Tiroler Landesregierung, Abteilung Raumordnung und Statistik (2019). Der Tourismus im Winter 2018/2019. Innsbruck, August 2019. Tab. 6 (p. 11). https://www.tirol.gv.at/fileadmin/themen/statistik-budget/statistik/downloads/Tourismus/Der_Tourismus_im_Winter_18_19.pdf
- Ashcroft, P., et al. (2020). COVID-19 infectivity profile correction. *Swiss Med Wkly* 150:w20336. https://doi.org/10.4414/smw.2020.20336
- Becker, N. G., Watson, L. F., Carlin, J. B. (1991). A method of non-parametric back-projection and its application to AIDS data. *Stat Med* 10(10):1527–1542. https://doi.org/10.1002/sim.4780101005
- Byambasuren, O., et al. (2020). Estimating the extent of asymptomatic COVID-19 and its potential for community transmission: Systematic review and meta-analysis. *JAMMI* 5(4):223–234. https://doi.org/10.3138/jammi-2020-0030
- Byrne, A. W., et al. (2020). Inferred duration of infectious period of SARS-CoV-2: rapid scoping review and analysis of available evidence for asymptomatic and symptomatic COVID-19 cases. *BMJ Open* 10(8):e039856. https://doi.org/10.1136/bmjopen-2020-039856
- Dehning, J., et al. (2020). Inferring change points in the spread of COVID-19 reveals the effectiveness of interventions. *Science* 369(6500):eabb9789. https://doi.org/10.1126/science.abb9789
- He, X., et al. (2020). Temporal dynamics in viral shedding and transmissibility of COVID-19. *Nat Med* 26:672–675. https://doi.org/10.1038/s41591-020-0869-5; author correction *Nat Med* 26:1491–1493. https://doi.org/10.1038/s41591-020-1016-z
- Hu, Z., et al. (2020). Clinical characteristics of 24 asymptomatic infections with COVID-19 screened among close contacts in Nanjing, China. *Sci China Life Sci* 63(5):706–711. https://doi.org/10.1007/s11427-020-1661-4
- Karagiannidis, C., et al. (2020). Case characteristics, resource use, and outcomes of 10 021 patients with COVID-19 admitted to 920 German hospitals: an observational study. *Lancet Respir Med* 8(9):853–862. https://doi.org/10.1016/S2213-2600(20)30316-7
- Liu, Y., et al. (2020). The reproductive number of COVID-19 is higher compared to SARS coronavirus. *J Travel Med* 27(2):taaa021. https://doi.org/10.1093/jtm/taaa021
- McAloon, C., et al. (2020). Incubation period of COVID-19: a rapid systematic review and meta-analysis of observational research. *BMJ Open* 10(8):e039652. https://doi.org/10.1136/bmjopen-2020-039652
- Mossong, J., et al. (2008). Social contacts and mixing patterns relevant to the spread of infectious diseases. *PLoS Med* 5(3):e74. https://doi.org/10.1371/journal.pmed.0050074
- Ponge, J., et al. (2024). A standardized framework for modeling non-pharmaceutical interventions in individual-based infectious disease simulations. *WSC 2024*, 1106–1117. https://doi.org/10.1109/WSC63780.2024.10838778
- Robert Koch-Institut (2020). Täglicher Lagebericht des RKI zur Coronavirus-Krankheit-2019 (COVID-19), reports of 04.03., 06.03., 08.03. (Tab. 2), 22.03., 24.03., 26.03. (Tab. 3), and 27.04.2020 (Fig. 3, Tab. 2). https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/2020-MM-DD-de.pdf?__blob=publicationFile (no longer online; archived copies: https://web.archive.org/web/2020/ followed by that URL)
- Robert Koch-Institut (2021). SARS-CoV-2-Seroprävalenz in der Allgemeinbevölkerung in Deutschland [Factsheet], as of 23.11.2021, p. 1. https://doi.org/10.25646/9288
- Robert Koch-Institut (2026). SARS-CoV-2 Infektionen in Deutschland [Data set], version of 14.10.2026. Zenodo. https://doi.org/10.5281/zenodo.23357370. License: CC BY 4.0
- Schilling, J., Lehfeld, A.-S., Schumacher, D., Ullrich, A., Diercke, M., Buda, S., Haas, W., RKI COVID-19 Study Group (2020). Krankheitsschwere der ersten COVID-19-Welle in Deutschland basierend auf den Meldungen gemäß Infektionsschutzgesetz. *Journal of Health Monitoring* 5(S11):2–20. https://doi.org/10.25646/7169
- Sekretariat der Kultusministerkonferenz (2018). Ferien im Schuljahr 2019/2020 in den Ländern in der Bundesrepublik Deutschland, as of 11.07.2018. https://www.kmk.org/fileadmin/Dateien/pdf/Ferienkalender/FER19_20.pdf
- Statistik Austria (n.d.). Zuordnungsverzeichnis Deutschland (German postal codes → origin regions of the Austrian tourism statistics), edited by Landesstatistik Steiermark. https://www.landesentwicklung.steiermark.at/cms/dokumente/12656884_141979459/b19d32c7/Zuordnungsverzeichnis%20Oesterreich%20Deutschland-PLZ-PDF.pdf
- Statistische Ämter des Bundes und der Länder (2026). Einkommen der privaten Haushalte in den kreisfreien Städten und Landkreisen der Bundesrepublik Deutschland 1995 bis 2023. Volkswirtschaftliche Gesamtrechnungen der Länder, Reihe 2, Kreisergebnisse Band 3, Tab. 2.4 (Berechnungsstand Februar 2025, published January 2026). https://www.statistikportal.de/de/vgrdl/publikationenbsalt. License: Datenlizenz Deutschland – Namensnennung – Version 2.0.
- Streeck, H., et al. (2020). Infection fatality rate of SARS-CoV2 in a super-spreading event in Germany. *Nat Commun* 11:5829. https://doi.org/10.1038/s41467-020-19509-y
- Verity, R., et al. (2020). Estimates of the severity of coronavirus disease 2019: a model-based analysis. *Lancet Infect Dis* 20(6):669–677. https://doi.org/10.1016/S1473-3099(20)30243-7
