# GEMS: The German Epidemic Microsimulation System — Paper Examples

This repository contains the code and model configurations behind the example applications in

> Ponge, J., Patzner, J., Suer, J., et al. *GEMS: The German Epidemic Microsimulation System.*

All simulations use [GEMS](https://github.com/IMMIDD/GEMS) **v1.0.1** (pinned in [Manifest.toml](Manifest.toml)).
This README documents the parameterization of every example: population, contacts, transmission, disease progression, and interventions.
All disease examples share one **reference model**, [models/DE_model_first_wave.toml](models/DE_model_first_wave.toml) ([Section 4](#4-reference-model-covid-19-wild-type-first-wave-2020)). Each example sets its own initial infections, horizon, and interventions in code, via the `Simulation` constructor or the scenario function's keyword arguments; [Section 5](#5-example-specific-parameterization) lists only these deviations.

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
| [src/benchmarking.jl](src/benchmarking.jl) | Performance benchmark (Fig 15) |
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

All simulations use fixed random seeds ([Section 3](#3-overview-examples-code-and-configuration)). Bit-identical results require the same number of Julia threads.

*Terminology:* a *random seed* initializes the random number generator of a simulation; infections that start an epidemic are called *initial infections* (at tick 0) or *introduced infections* (added during the run).

## 3. Overview: examples, code, and configuration

| Paper element | Script → function | Config | Population | Horizon | Runs / random seeds |
|---|---|---|---|---|---|
| Fig 4–6 (population, demographics, commuting) | `demographics.jl`, `settings.jl`, `mobility.jl` | none | DE | – | – |
| Fig 7 (contact patterns) | [contacts.jl](src/contacts.jl) → `contacts()` | GEMS defaults | SL | – | 0 |
| Sec. 3.1.1 First wave 2020 | [sc_first_wave.jl](src/sc_first_wave.jl) → `scenario_first_wave()` | reference model | DE | 10.02.–31.05.2020 | see [5.1](#51-first-wave-2020-sec-311) |
| Sec. 3.1.2 Scenario modeling (Fig 12A) | [sc_household_isolation.jl](src/sc_household_isolation.jl) → `scenario_household_isolation()` | reference model | DE | 100 d | 1 per scenario / 0 |
| Fig 12B (code listing) | [sc_household_isolation_screenshot_file.jl](src/sc_household_isolation_screenshot_file.jl) | GEMS defaults | on-the-fly, 100,000 | 365 d | 1 per scenario / none (random) |
| Sec. 3.1.3 Parameter sweep (Fig 13) | [sc_sweeping.jl](src/sc_sweeping.jl) → `scenario_sweeping()` | reference model | SL | 365 d | 41 × 41 grid × 3 / 1–5043 |
| Sec. 3.1.4 Multi-criteria evaluation (Fig 14, Tab 2) | [sc_testing.jl](src/sc_testing.jl) → `scenario_testing()` | reference model | DE | 1095 d | 1 per scenario / 0, 1, 2 |
| Sec. 3.2 Performance (Fig 15) | [benchmarking.jl](src/benchmarking.jl) → `benchmark()` | GEMS defaults | 16 states + DE; on-the-fly 0.1–50 M | 365 d | 11 per population (median) / 1–11 |

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

### 5.1 First wave 2020 (Sec. 3.1.1)

Code: [sc_first_wave.jl](src/sc_first_wave.jl), [src/first_wave/](src/first_wave/). Config: reference model. Compared with the RKI cases by reporting date and by symptom onset.

**Data**

| Item | Value | Source |
|---|---|---|
| RKI cases | county, reporting date, symptom onset; version of 14.10.2026 | Robert Koch-Institut 2026 (Zenodo, CC BY 4.0) |
| Unknown onsets | onsets > 21 d before the report, or later after the report than the 99 % quantile of the model's incubation period, are treated as unknown | [A] data errors / a case reported before onset is at most one incubation period before it |
| Counties | Berlin districts merged into 11000; Eisenach (16056) merged into the Wartburgkreis (16063) | data structure |
| Reporting delays | RKI onset → report delays per onset week (all of Germany) | RKI cases |

**Introduced infections** (added via `stepmod`, as `RegionalSeeds` in GEMS v1.0.1 only infects at tick 0)

| Item | Value | Source |
|---|---|---|
| Imported cases per reporting week (weeks 9–17) | 35, 262, 2,781, 6,454, 3,814, 1,095, 267, 113, 31 | [L] RKI situation report of 27.04.2020, Fig. 3 (weeks 10–17); week 9 extrapolated from the reports of 04.–08.03. |
| Share of Austria and Italy | 0.90 (weeks 9–12), 0.84 (weeks 13–17) | [L] country tables of the RKI situation reports up to 08.03., 08.–22.03., and 22.–26.03. |
| Share of Italy in the Alpine imports | 1.0 (weeks 9–10), 0.21 (11–12), 0.07 (13–17) | [L] RKI situation reports of 08.03., 22.03. (995 of 4,849), and 27.04. (474 of 6,443) |
| Import curve | back-calculated to infection dates over the model's time to onset and the RKI reporting delays (200 iterations), divided by the ascertainment | [A] |
| Allocation to counties: Austria | Tyrol overnight stays per German origin region (winter 2018/19), distributed within the region by population × relative income^1.0 × (distance to Ischgl / 100 km)^−0.4, × (1 + 1.0 × school holidays of the state) | [L] Amt der Tiroler Landesregierung 2019, Tab. 6; [L] income elasticity 1.0 (T-MONA guest survey); [L] adapted: distance decay 0.4 (slope of the Tyrol overnight stays across the regions); [A] holiday factor 1.0 |
| Allocation: Italy | 70 % like Austria with Bozen as reference point, 30 % by population × income | [L] ski share 0.7 (RKI situation reports) |
| Allocation: other countries | population × income | [A] |
| School holidays (mid-February–March 2020) | BY, BW 22.02.–01.03.; HH 29.02.–15.03.; RP 15.–23.02.; SL 15.–25.02.; SN, MV 08.–23.02.; BW: 50 % of the schools | [L] school holiday calendars; [A] BW coverage 0.5 |
| County income | disposable income of private households per inhabitant 2019 | [L] VGR der Länder, Reihe 2 Band 3, Tab. 2.4 ([inputs/county_income_2019.csv](inputs/county_income_2019.csv)) |
| Heinsberg cluster | RKI cases of Kreis Heinsberg with onset (or, if unknown, reporting date) up to 29.02. and reported by 07.03., infected 1:1 at onset − 6 d (rounded mean incubation 1 + 0.8 + 4.0), or report − 11 d (+ median reporting delay 5 d) if the onset is unknown | [A] |

**Contact phases and testing**

| Item | Value | Source |
|---|---|---|
| Phase start dates | phase 2: 09.03. (large events cancelled), phase 3: 16.03. (schools closed), phase 4: 23.03. (contact restrictions) | [L] dates of the measures |
| Transitions between phases | linear over 3, 3, 4 days | [L] Dehning et al. 2020 |
| Contact factors (phases 1–4) | 0.570 / 0.164 / 0.162 / 0.120, multiplied with the contact parameters of all settings except households | [C] [inputs/first_wave_factors_DE_model_first_wave.csv](inputs/first_wave_factors_DE_model_first_wave.csv) |
| Care homes | households with at least 11 members, scaled like the other settings | [A] |
| School closure | all school settings closed from 16.03.2020 | [L] |
| Ascertainment | 0.22 of all infections detected | – |
| Testing | perfect PCR test at symptom onset with probability p = (0.22 − hospitalized share) / (symptomatic share − hospitalized share), and at hospital admission (if not detected before); positives are not isolated | [A] |
| Reported cases by reporting date | detected infections, reported at symptom onset (asymptomatic: infection + 6 d) plus a delay drawn from the RKI delays of the onset week | [A] |
| Reported cases by onset | detected infections at symptom onset (asymptomatic: infection + 6 d), each kept with the share of RKI cases with a known onset (per onset week); compared with the RKI cases with a known onset | [A] |

**Calibration** ([calibrate_first_wave.jl](calibrate_first_wave.jl) → `calibrate_first_wave()`, not part of `runall()`)

| Step | Setup |
|---|---|
| 1. Transmission scale | [Section 4.3](#43-transmission) |
| 2a. Generation time of the model | Saarland, 70 days, scale 0.35, 20 extra infections per day in the Regionalverband Saarbrücken (ticks 7–27), ascertainment 1.0, random seed 2; cohort of ticks 10–22; mean 5.77 d ([inputs/generation_time_DE_model_first_wave.csv](inputs/generation_time_DE_model_first_wave.csv)) |
| 2b. Contact factors | fitted to the weekly RKI cases by onset, weeks 10–16 (up to 19.04., weight 0.5 for week 10), with the R per phase of GEMS kept inside the 95 % interval of the RKI R per phase (R(t) back-calculated from 01.02., infections 01.03.–12.04., block bootstrap with 1,000 samples, block length 7 d); surrogate fit, then damped Gauss–Newton steps, at most 7 evaluations with 2 runs each (random seeds 0, 1); runs until 17.05. |
| 2c. Band | 20 factor samples from the calibration uncertainty (random seeds 1001–1020) ([inputs/first_wave_band_factors_DE_model_first_wave.csv](inputs/first_wave_band_factors_DE_model_first_wave.csv)) |

**Parts of the scenario** (`FIRST_WAVE_PARTS`)

| Part | Setup | Runs / random seeds |
|---|---|---|
| `:final` | calibrated factors until 31.05. (six weeks of forecast after 19.04.) | 3 / 0–2 |
| `:projections` | from 20.04. on, the contacts of phase 2, 3, or 4, or of phase 1 (until 17.05.); schools reopen in the projections to phase 1 and 2 | 2 per variant / 0, 1 |
| `:counterfactuals` | contacts of phase k kept from the start of phase k + 1 on (k = 1–3) | 1 per variant / 0 |
| `:band` | factor samples of the calibration, 5–95 % band | 20 / 1001–1020 |
| `:seeds_robustness` | all RKI cases infected before 09.03. as introduced infections (divided by the ascertainment, Heinsberg 1:1), no transmission before 09.03. | 3 / 0–2 |
| `:figures` | figures of all parts whose results exist | – |

### 5.2 Scenario modeling (Sec. 3.1.2, Fig 12)

Code: [sc_household_isolation.jl](src/sc_household_isolation.jl). Config: reference model.

| Parameter | Reference model | This example |
|---|---|---|
| Initial infections | none | 0.01 % of the population, random (`infected_fraction = 0.0001`) |
| Duration | 112 days | 100 days |
| Random seed | – | 0 (both scenarios) |

| Scenario | Trigger | Strategy (condition) | Measures (offset) → follow-up |
|---|---|---|---|
| Baseline | – | – | – |
| Scenario | `SymptomTrigger` | Find Household Members | `FindSettingMembers` in Household, incl. the symptomatic individual (0) → *Isolation* |
| | | Isolation | `SelfIsolation(14)` (0) |

All household members, including the symptomatic individual, isolate at home for 14 days. Each new symptomatic household member triggers the strategy again and restarts the 14 days for the whole household.

The Fig 12B listing ([sc_household_isolation_screenshot_file.jl](src/sc_household_isolation_screenshot_file.jl)) uses the same intervention with GEMS defaults ([Section 5.5](#55-performance-benchmark-sec-32-fig-15)) on 100,000 agents, without a fixed random seed.

### 5.3 Parameter sweep (Sec. 3.1.3, Fig 13)

Code: [sc_sweeping.jl](src/sc_sweeping.jl). Config: reference model.

| Parameter | Reference model | This example |
|---|---|---|
| Population | – | Saarland |
| Initial infections | none | 0.1 % of the population, random |
| Duration | 112 days | 365 days |
| Transmission scale | 0.50486 | × (1 − 0.6): 60 % of transmissions prevented by other measures [I] |
| Asymptomatic share (x-axis) | 0.17 | 0 to 1, step 0.025, set in all age groups; the other categories rescaled proportionally |
| Adherence (y-axis) | – | 0 to 1, step 0.025 |
| Repetitions | – | 3 per combination; random seeds 1–5043 (one per run) |
| Outcome | – | R_init per run (`r0(PostProcessor(sim))`) |

| Trigger | Strategy (condition) | Measures |
|---|---|---|
| `SymptomTrigger` | Self Isolation (adheres ~ Bernoulli(p = adherence), drawn from the simulation's RNG) | `SelfIsolation(14)` (0) |

### 5.4 Multi-criteria strategy evaluation (Sec. 3.1.4, Fig 14, Tab 2)

Code: [sc_testing.jl](src/sc_testing.jl). Config: reference model.

| Parameter | Reference model | This example |
|---|---|---|
| Initial infections | none | 10 infections in Gangelt (municipality key 05370008, county Heinsberg) |
| Duration | 112 days | 1095 days |
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

The tracing scenario corresponds to Fig 9 of the paper.

### 5.5 Performance benchmark (Sec. 3.2, Fig 15)

GEMS default parameterization ([data/DefaultConf.toml @ v1.0.1](https://github.com/IMMIDD/GEMS/blob/v1.0.1/data/DefaultConf.toml)), also used for Fig 7 and Fig 12B.

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
- **Fig 7:** Saarland, GEMS defaults.

## 7. References

- Amt der Tiroler Landesregierung (2019). Der Tourismus im Winter 2018/19, Tab. 6.
- Ashcroft, P., et al. (2020). COVID-19 infectivity profile correction. *Swiss Med Wkly* 150:w20336. https://doi.org/10.4414/smw.2020.20336
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
- Robert Koch-Institut (2020). Täglicher Lagebericht des RKI zur Coronavirus-Krankheit-2019 (COVID-19), reports of 04.03.–27.04.2020.
- Robert Koch-Institut (2026). SARS-CoV-2 Infektionen in Deutschland [Data set], version of 14.10.2026. Zenodo. https://doi.org/10.5281/zenodo.23357370. License: CC BY 4.0
- Schilling, J., Lehfeld, A.-S., Schumacher, D., Diercke, M., Buda, S., Haas, W., RKI COVID-19 Study Group (2020). Krankheitsschwere der ersten COVID-19-Welle in Deutschland basierend auf den Meldungen gemäß Infektionsschutzgesetz. *Journal of Health Monitoring* 5(S11). https://doi.org/10.25646/7169
- Statistische Ämter des Bundes und der Länder (2025). Volkswirtschaftliche Gesamtrechnungen der Länder, Reihe 2 Band 3, Tab. 2.4 (Berechnungsstand Februar 2025). statistikportal.de. License: Datenlizenz Deutschland – Namensnennung – Version 2.0.
- T-MONA guest survey (income elasticity of ski trips).
- Verity, R., et al. (2020). Estimates of the severity of coronavirus disease 2019: a model-based analysis. *Lancet Infect Dis* 20(6):669–677. https://doi.org/10.1016/S1473-3099(20)30243-7
