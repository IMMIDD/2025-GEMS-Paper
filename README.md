# GEMS: The German Epidemic Microsimulation System — Paper Examples

This repository contains the code and model configurations behind the example applications in

> Ponge, J., Patzner, J., Suer, J., et al. *GEMS: The German Epidemic Microsimulation System.*

All simulations use [GEMS](https://github.com/IMMIDD/GEMS) **v1.0.1** (pinned in [Manifest.toml](Manifest.toml)).
This README documents the parameterization of every example: population, contacts, transmission, disease progression, and interventions.
Most examples share one **reference parameterization** ([Section 4](#4-reference-parameterization-covid-19-wild-type-2020)); each example in [Section 5](#5-example-specific-parameterization) lists only what it changes.

---

## 0. 🚧 Open items (remove before release)

This README describes the **target** parameterization. Markers:

- **🔄 CHANGE**: differs from the current config or code (current value given as `config: …`)
- **🔍 VERIFY**: value or source not yet checked against the primary source
- **🚧 DECISION** (not decided; proposed value shown) · **🚧 UPDATE** (decided, not yet applied) · **🚧 INPUT** (information missing)

| ID | Type | Item | Section |
|----|------|------|---------|
| [P1](#p1) | DECISION | Death pathway. **B1** (shown): built-in categories, every death goes through ICU. **B2:** custom category "dies in hospital without ICU" | 4.4 |
| [P2](#p2) | UPDATE | `Severe` = one severe case per hospitalized case | 4.4 |
| [P3](#p3) | UPDATE | Progression matrix from Verity et al. 2020 🔍, `death_probability` from Karagiannidis et al. 2020 🔍 | 4.4 |
| [P4](#p4) | DECISION | Latent period 3.55 (incubation 6.55 d incl. GEMS' +1 day) instead of 4.55 (7.55 d) | 4.5 |
| [P5](#p5) | UPDATE | Severity, hospital, and ICU durations from Schilling et al. 2020 🔍 (C9); check that GEMS accepts a fixed 0 for same-day transitions (paper line 430–431 says durations can be set to zero) | 4.5 |
| [P6](#p6) | DECISION | Asymptomatic infectious period 15.4 d (as symptomatic) instead of 13.4 d | 4.5 |
| [P7](#p7) | INPUT | Reason for 0.5 municipality contacts in the multi-criteria example | 5.4 |
| [P8](#p8) | DECISION | Fit transmission rates to R_init instead of mean R_eff over days 1–10 | 4.3 |
| [P9](#p9) | DECISION | Include the RKI validation? | 5.7 |
| [P10](#p10) | INPUT | Beryl (RSV): code repository and parameter table | 5.5 |
| [P12](#p12) | UPDATE | `runall()` builds Fig 7 from DE, the figure shows Saarland | 2 |
| [P13](#p13) | UPDATE | Delete unused [DE_model_R0_3.28_Beckingen_Seed.toml](models/DE_model_R0_3.28_Beckingen_Seed.toml) | 7 |
| P14 | UPDATE | Apply all 🔄 changes, refit transmission rates, rerun Fig 10, 11, 12A, 14, Tab 2, RKI validation | all |

<a id="config-changes"></a>**🔄 Config changes:**

| # | Parameter | Current | Target | Status |
|---|---|---|---|---|
| C1 | Progression matrix | WSC 2025 fatality values split 70/21/9 % into Severe/Hospitalized/Critical | from Verity et al. 2020 🔍 (Sec. 4.4) | [P3](#p3) |
| C2 | `Severe` share | 0.0014 / 0.0063 / 0.0350 / 0.1036 | = IHR per age group | [P2](#p2) |
| C3 | `Hospitalized` / `Critical` split | 70 : 30 | `Critical = IFR / 0.53` | [P1](#p1) |
| C4 | `Critical.death_probability` | 1.0 | 0.53 | [P3](#p3) |
| C5 | `exposure_to_infectiousness_onset` | Poisson(λ = 4.55) | Poisson(λ = 3.55) | [P4](#p4) |
| C6 | `Asymptomatic.infectiousness_onset_to_recovery` | Poisson(λ = 13.4) | Poisson(λ = 15.4) | [P6](#p6) |
| C7 | Fit target in [fit_r0.jl](fit_r0.jl) | mean R_eff, days 1–10 | R_init | [P8](#p8) |
| C8 | Transmission rates | 0.03142, 0.06888, 0.077, RKI grid | refit after C1–C7, C9 | P14 |
| C9 | Severity, hospital, and ICU durations | see Sec. 4.5 (`config: …`) | from Schilling et al. 2020 🔍 | [P5](#p5) |

<a id="verify"></a>**🔍 To verify:**

| ID | What |
|----|------|
| V1 | Population weights (currently approximate Destatis 31.12.2019 counts; table ID unverified). Use exact Destatis counts or the GEMS DE population |
| V2 | Karagiannidis et al. 2020: 53 % mortality of invasively ventilated patients (so far only from abstract/press release) |
| V3 | Verity et al. 2020: IFR (Table 1) and IHR (Table 3) by 10-year age bin |
| V4 | Schilling et al. 2020, Table 5: durations used in Sec. 4.5 |
| V5 | Salje et al. 2020: 19 % ICU admission among hospitalized (only needed for B2, [P1](#p1)) |
| V6 | DOIs, volumes, and pages in [Section 8](#8-references) |

GEMS defaults and the +1 offset are taken from the GEMS v1.0.1 source code. Adopted from the WSC 2025 model with its sources: Mossong, Liu, Byambasuren, Wu, Byrne.

**Manuscript fixes:**

| ID | Item |
|----|------|
| M1 | Line 365 cites COVIMOD for contact numbers; configs use POLYMOD (Mossong), and Fig 7 does not depend on contact numbers. |
| M2 | Sec. 3.1.2: baseline is seeded with 0.1 % random infections, not in Gangelt. Fig 12B listing runs GEMS defaults, not the Fig 12A setup. |
| M3 | Sec. 3.1.4 does not mention its deviations (municipality contacts, transmission rate, 3 years, one seed per scenario). |
| M4 | Sec. 3.1.5 placeholder "(add table)" (line 854). |
| M5 | Fig 8B shows ICU survival, but currently every ICU patient dies (`death_probability = 1.0`). Resolved by C4. |
| M6 | Check stated incubation and infectious periods if P4 and P6 are adopted. |
| M7 | Pathogen label. Line 747–748 (Sec. 3.1.1): "parameterize GEMS for a COVID-19-like pathogen" → "for COVID-19 (wild type, 2020; full parameterization in the paper repository)". Line 808–809 (Sec. 3.1.4): "a COVID-19-like outbreak with illustrative parameter values that are not calibrated to observed data" → "the COVID-19 parameterization of the previous examples; results are not calibrated to observed data". |
| M8 | Numbers in the text that change with the rerun (P14): line 750–751 ("approximately eight months"), line 834–837 (attack rates 66.1 % / 38.5 %, 398 million tests, 649 million / 4.9 billion quarantine days), Fig 10, 11 (incl. regional R0 and fatality maps), 12A, 14, Tab 2. |
| M9 | Data Availability (line 1040–1042): mention that the repository README contains the full parameterization of all examples. |
| M10 | If the RKI validation is included ([P9](#p9)): line 1008–1012 states that the paper contains no validation against surveillance data. |

---

## 1. Repository layout

| Path | Content |
|------|---------|
| [models/](models/) | TOML configuration files ([Section 7](#7-configuration-files)) |
| [src/GEMSPaper.jl](src/GEMSPaper.jl) | Package entry point and `runall()` |
| `src/sc_*.jl` | One file per example scenario |
| [src/demographics.jl](src/demographics.jl), [src/settings.jl](src/settings.jl), [src/mobility.jl](src/mobility.jl), [src/contacts.jl](src/contacts.jl) | Population analyses (Fig 4–7) |
| [src/benchmarking.jl](src/benchmarking.jl) | Performance benchmark (Fig 15) |
| `src/rki_*.jl` | RKI data and seeding for the validation scenario |
| [fit_r0.jl](fit_r0.jl) | Transmission-rate fitting |

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

Each example can also be run on its own: `scenario_outbreak()`, `scenario_self_isolation()`, `scenario_sweeping()`, `scenario_testing()`, `benchmark()`, `contacts()`, `demographics()`, `mobility()`, `settinganalyses()`.

All simulations use fixed seeds ([Section 3](#3-overview-examples-code-and-configuration)). Bit-identical results require the same number of Julia threads.

> <a id="p12"></a>🚧 **UPDATE (P12):** Fig 7 was created with `contacts()` on Saarland; `runall()` passes the DE simulation.

## 3. Overview: examples, code, and configuration

| Paper element | Script → function | Config | Population | Horizon | Runs / seeds |
|---|---|---|---|---|---|
| Fig 4–6 (population, demographics, commuting) | `demographics.jl`, `settings.jl`, `mobility.jl` | none | DE | – | – |
| Fig 7 (contact patterns) | [contacts.jl](src/contacts.jl) → `contacts()` | GEMS defaults | SL | – | seed 0 |
| Sec. 3.1.1 Outbreak projection (Fig 10, 11) | [sc_outbreak.jl](src/sc_outbreak.jl) → `scenario_outbreak()` | [DE_model_R0_3.28_Heinsberg_Seed.toml](models/DE_model_R0_3.28_Heinsberg_Seed.toml) 🔄 | DE | 365 d | 1 / seed 0 |
| Sec. 3.1.2 Scenario modeling (Fig 12A) | [sc_self_isolation.jl](src/sc_self_isolation.jl) → `scenario_self_isolation()` | [DE_model_R0_3.28.toml](models/DE_model_R0_3.28.toml) 🔄 | DE | 365 d | 1 per scenario / seed 0 |
| Fig 12B (code listing) | [sc_self_isolation_screenshot_file.jl](src/sc_self_isolation_screenshot_file.jl) | GEMS defaults | on-the-fly, 100,000 | 365 d | seed 0 |
| Sec. 3.1.3 Parameter sweep (Fig 13) | [sc_sweeping.jl](src/sc_sweeping.jl) → `scenario_sweeping()` | defined in script + GEMS defaults | SL | 365 d | 41 × 41 grid × 3 / seeds 1–3 |
| Sec. 3.1.4 Multi-criteria evaluation (Fig 14, Tab 2) | [sc_testing.jl](src/sc_testing.jl) → `scenario_testing()` | [DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml](models/DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml) 🔄 | DE | 1095 d | 1 per scenario / seeds 0, 1, 2 |
| Sec. 3.1.5 Beryl (RSV) | 🚧 [P10](#p10) | – | – | – | – |
| Sec. 3.2 Performance (Fig 15) | [benchmarking.jl](src/benchmarking.jl) → `benchmark()` | GEMS defaults | 16 states + DE; on-the-fly 0.1–50 M | 365 d | 11 per population (median) / seeds 1–11 |
| 🚧 RKI validation ([P9](#p9)) | [sc_rki_validation.jl](src/sc_rki_validation.jl) | [DE_model_RKI_2020.toml](models/DE_model_RKI_2020.toml) 🔄, [DE_model_RKI_2020_short_gt.toml](models/DE_model_RKI_2020_short_gt.toml) 🔄 | DE | 42 d | seed 0 |

## 4. Reference parameterization (COVID-19, wild type, 2020)

Config: [models/DE_model_R0_3.28.toml](models/DE_model_R0_3.28.toml).

**Source tags:** `[L]` literature · `[C]` calibrated · `[A]` assumption · `[D]` GEMS default · `[I]` illustrative

### 4.1 Simulation and seeding

| Parameter | Value | Source |
|---|---|---|
| Time step | 1 day | [D] |
| Duration | 365 days | [I] |
| Initial infections | 0.1 % of the population, random | [I] |

### 4.2 Contacts

Contact partners are drawn uniformly from the setting.

| Setting type | Level | Daily contacts | Source |
|---|---|---|---|
| Household | – | Poisson(λ = 3.3) | [L] Mossong et al. 2008 |
| Workplace | Office / Department / Workplace / WorkplaceSite | Poisson(λ = 1.4 / 0.7 / 0.35 / 0.35) | [L] Mossong et al. 2008 (total 2.8); [A] split 1/2, 1/4, 1/8, 1/8 |
| School | SchoolClass / SchoolYear / School / SchoolComplex | Poisson(λ = 0.9 / 0.45 / 0.225 / 0.225) | [L] Mossong et al. 2008 (total 1.8); [A] split 1/2, 1/4, 1/8, 1/8 |
| Municipality | – | Poisson(λ = 4.0) | [L] Mossong et al. 2008 |

### 4.3 Transmission

| Parameter | Value | Source |
|---|---|---|
| Transmission function | constant per-contact probability (`ConstantTransmissionRate`) | [D] |
| Transmission probability | 🚧 refit 🔄 (`config: 0.03142`) | [C] fitted to R0 = 3.28 with [fit_r0.jl](fit_r0.jl) (10 runs per candidate) |
| Target R0 | 3.28 | [L] Liu et al. 2020 |

> <a id="p8"></a>🚧 **DECISION (P8):** Fit target R_init 🔄 (`code: mean R_eff, days 1–10`).

### 4.4 Progression assignment

Each infected individual's progression category is drawn from Categorical(p), where p is the matrix row of their age group (`-39`, `40-59`, `60-79`, `80-`). In GEMS, only `Critical` cases can die, after ICU admission: death ~ Bernoulli(p = `death_probability`).

<a id="p3"></a>

| Input | Value | Source |
|---|---|---|
| Asymptomatic share | 0.17 | [L] Byambasuren et al. 2020 |
| IFR and IHR by age | 🔄 | [L] Verity et al. 2020, Tables 1 and 3 🔍 [V3](#verify); adapted: population-weighted means of the 10-year bins per age group, weights 🔍 [V1](#verify) from Destatis 2019 |
| `Critical.death_probability` | 0.53 🔄 (`config: 1.0`) | [L] Karagiannidis et al. 2020 🔍 [V2](#verify); adapted: mortality of invasively ventilated patients |
| `Severe` share | = IHR 🔄 | <a id="p2"></a>[A] one severe case per hospitalized case |

Matrix: `Critical = IFR / 0.53`, `Hospitalized = IHR − Critical`, `Severe = IHR`, `Symptomatic = remainder`.

**Progression matrix** 🔄 (🔍 preliminary, approximate weights [V1](#verify)):

| Age group | Asymptomatic | Symptomatic | Severe | Hospitalized 🚧 | Critical 🚧 |
|---|---|---|---|---|---|
| 0–39 | 0.17 | 0.803818 | 0.013091 | 0.012426 | 0.000665 |
| 40–59 | 0.17 | 0.700496 | 0.064752 | 0.057054 | 0.007698 |
| 60–79 | 0.17 | 0.554284 | 0.137858 | 0.083099 | 0.054759 |
| 80+ | 0.17 | 0.462000 | 0.184000 | 0.036830 | 0.147170 |

> <a id="p1"></a>🚧 **DECISION (P1):** B1 (shown) or B2 (custom category "dies in hospital without ICU"). Under B1, 40 % (60–79) and 80 % (80+) of hospitalized patients go to ICU.

### 4.5 Transition durations

In days. GEMS adds 1 day to every exposure → infectiousness-onset draw.

| Category | Transition | Distribution | Basis |
|---|---|---|---|
| all | exposure → infectiousness onset | <a id="p4"></a>🚧 1 + Poisson(λ = 3.55) 🔄 (`config: 1 + Poisson(λ = 4.55)`) | [L] adapted: incubation 6.55 (Wu et al. 2022) minus pre-symptomatic period 2.0 (Byrne et al. 2020) minus 1 |
| Asymptomatic | infectiousness onset → recovery | <a id="p6"></a>🚧 Poisson(λ = 15.4) 🔄 (`config: Poisson(λ = 13.4)`) | [A] same as symptomatic (2.0 + 13.4) |
| Symptomatic, Severe, Hospitalized, Critical | infectiousness onset → symptom onset | Poisson(λ = 2.0) | [L] Byrne et al. 2020 |
| Symptomatic | symptom onset → recovery | Poisson(λ = 13.4) | [L] Byrne et al. 2020 |
| Severe, Hospitalized | symptom onset → severeness onset | Poisson(λ = 4) 🔄 (`config: Poisson(λ = 3.0)`) | <a id="p5"></a>[L] Schilling et al. 2020 🔍 [V4](#verify): symptom onset → hospital admission, median 4 d |
| Critical | symptom onset → severeness onset | Poisson(λ = 5) 🔄 (`config: Poisson(λ = 3.0)`) | [L] Schilling et al. 2020 🔍 [V4](#verify): symptom onset → hospital admission of ICU patients, median 5 d |
| Hospitalized, Critical | severeness onset → hospital admission | 🚧 0 🔄 (`config: Poisson(λ = 4.0)`) | [A] admission at the onset of severe illness |
| Critical | hospital admission → ICU admission | 🚧 0 🔄 (`config: Poisson(λ = 2.0)`) | [L] Schilling et al. 2020 🔍 [V4](#verify): median 0 d |
| Hospitalized | hospital admission → hospital discharge | Poisson(λ = 8) 🔄 (`config: Poisson(λ = 7.0)`) | [L] Schilling et al. 2020 🔍 [V4](#verify): hospital stay without ICU, median 8 d |
| Critical (survives) | ICU admission → ICU discharge | Poisson(λ = 11) 🔄 (`config: Poisson(λ = 7.0)`) | [L] Schilling et al. 2020 🔍 [V4](#verify): ICU stay of survivors, median 11 d |
| Critical (dies) | ICU admission → death | Poisson(λ = 8) 🔄 (`config: Poisson(λ = 10)`) | [L] Schilling et al. 2020 🔍 [V4](#verify): ICU stay of deceased, median 8 d |
| Critical (survives) | ICU discharge → hospital discharge | Poisson(λ = 3) 🔄 (`config: Poisson(λ = 7.0)`) | [L] adapted: Schilling et al. 2020 🔍 [V4](#verify) hospital stay of ICU patients (14 d) minus ICU stay of survivors (11 d) |
| Severe | severeness onset → severeness offset | Poisson(λ = 8) 🔄 (`config: Poisson(λ = 7.0)`) | [A] same as hospital stay without ICU |
| Hospitalized, Critical | hospital discharge → severeness offset | Poisson(λ = 3.0) | [A] patients are discharged before the severe phase is over and remain severely ill at home for 3 days |
| Severe, Hospitalized, Critical | severeness offset → recovery | Poisson(λ = 4.0) | [A] 4-day convalescence phase between the end of severe illness and recovery |

> 🚧 **UPDATE (P5):** Check that GEMS accepts a fixed 0 for the two same-day transitions.

## 5. Example-specific parameterization

Interventions use the Trigger–Strategy–Measure (TriSM) framework ([Ponge et al. 2024](https://doi.org/10.1109/WSC63780.2024.10838778)), listed as trigger → strategy (condition) → measures (offset).

### 5.1 Outbreak projection (Sec. 3.1.1, Fig 10, 11)

Config: [DE_model_R0_3.28_Heinsberg_Seed.toml](models/DE_model_R0_3.28_Heinsberg_Seed.toml) 🔄

| Parameter | Reference | This example |
|---|---|---|
| Initial infections | 0.1 % random | 10 infections in Gangelt (municipality key 05370008, county Heinsberg) |

Interventions: none.

### 5.2 Scenario modeling (Sec. 3.1.2, Fig 12)

Config: [DE_model_R0_3.28.toml](models/DE_model_R0_3.28.toml) 🔄 (reference)

| Scenario | Trigger | Strategy (condition) | Measures |
|---|---|---|---|
| Baseline | – | – | – |
| Scenario | `SymptomTrigger` | Self Isolation (none) | `SelfIsolation(14)` (0) |

The Fig 12B listing ([sc_self_isolation_screenshot_file.jl](src/sc_self_isolation_screenshot_file.jl)) uses the same intervention with GEMS defaults ([Section 5.6](#56-performance-benchmark-sec-32-fig-15)) on 100,000 agents.

### 5.3 Parameter sweep (Sec. 3.1.3, Fig 13)

Pathogen defined in [sc_sweeping.jl](src/sc_sweeping.jl); everything else GEMS defaults ([Section 5.6](#56-performance-benchmark-sec-32-fig-15)). Population: Saarland.

| Parameter | Value | Source |
|---|---|---|
| Transmission probability | 0.04 | [I] |
| Progression categories | Asymptomatic, Symptomatic (one age group) | [I] |
| Asymptomatic: exposure → infectiousness onset → recovery | 1 → 7 days (fixed, not drawn) | [I] |
| Symptomatic: exposure → infectiousness onset → symptom onset → recovery | 1 → 2 → 7 days (fixed, not drawn) | [I] |
| Symptomatic share (x-axis) | 0 to 1, step 0.025 | swept |
| Adherence (y-axis) | 0 to 1, step 0.025 | swept |
| Repetitions | 3 per combination (seeds 1–3) | |
| Outcome | R_init (mean of 3 runs) | |

| Trigger | Strategy (condition) | Measures |
|---|---|---|
| `SymptomTrigger` | Self Isolation (adheres ~ Bernoulli(p = adherence)) | `SelfIsolation(14)` (0) |

### 5.4 Multi-criteria strategy evaluation (Sec. 3.1.4, Fig 14, Tab 2)

Config: [DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml](models/DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml) 🔄

| Parameter | Reference | This example |
|---|---|---|
| Initial infections | 0.1 % random | 10 infections in Gangelt (05370008) |
| Municipality contacts | 4.0 | 0.5 <a id="p7"></a>🚧 INPUT (P7) |
| Transmission probability | – | 🚧 refit 🔄 (`config: 0.06888`), [C] fitted to R0 = 3.28 |
| Duration | 365 days | 1095 days |
| Seeds | – | 0 (Isolation Only), 1 (Isolation & Testing), 2 (Isolation, Testing & Tracing) |

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

### 5.5 Beryl — RSV (Sec. 3.1.5)

> <a id="p10"></a>🚧 **INPUT (P10):** Code repository and parameter table. Scenario definitions: [RespiCompass RSV Round 1 2025/26](https://github.com/european-modelling-hubs/RespiCompass/blob/main/Previous_Rounds/2025-2026_round_1/round1_2526_rsv.md).

### 5.6 Performance benchmark (Sec. 3.2, Fig 15)

GEMS default parameterization ([data/DefaultConf.toml @ v1.0.1](https://github.com/IMMIDD/GEMS/blob/v1.0.1/data/DefaultConf.toml)), also used for Fig 7, Fig 12B, and the parameter sweep.

| Parameter | Value |
|---|---|
| Initial infections | 0.1 % random |
| Duration | 365 days |
| Transmission probability | 0.2 |
| Daily contacts | Poisson(λ = 1.0) on every setting level |
| Progression category (ages ≤14 / 15–65 / 66+) | Categorical(p): Asymptomatic 0.40 / 0.25 / 0.15; Symptomatic 0.58 / 0.60 / 0.40; Severe 0.01 / 0.11 / 0.25; Hospitalized 0.007 / 0.03 / 0.12; Critical 0.003 / 0.01 / 0.08 |
| Durations (days) | exposure → infectiousness 1 + Poisson(λ = 1); infectiousness → symptoms Poisson(λ = 1); symptom → recovery Poisson(λ = 7) (asymptomatic: infectiousness → recovery Poisson(λ = 8)); severeness onset Poisson(λ = 1); severeness duration Poisson(λ = 7); to hospital Poisson(λ = 2); hospital stay Poisson(λ = 7); to ICU Poisson(λ = 2); ICU stay Poisson(λ = 7); ICU → death Poisson(λ = 10); post-hospital Poisson(λ = 3) + Poisson(λ = 4) |
| Death (Critical) | Bernoulli(p = 0.3) |
| On-the-fly population | 100,000 agents, mean household size 3, office size 5, school size 100 |

Protocol: one warm-up run, then 11 runs per population (seeds 1–11), medians reported. Populations: 16 German states and Germany, plus on-the-fly populations of 0.1, 1, 5, 10, 20, and 50 million agents.

### 5.7 🚧 RKI validation, first wave 2020

> <a id="p9"></a>🚧 **DECISION (P9):** Include?

Config: [DE_model_RKI_2020.toml](models/DE_model_RKI_2020.toml) 🔄

| Parameter | Value |
|---|---|
| Simulation period | 2020-02-10 to 2020-03-22 |
| Initial infections | none; imported cases from RKI data plus the Heinsberg cluster, back-dated by 7 days [A] |
| Transmission probability | 🚧 rerun 🔄; grid search against weekly RKI cases (weeks 10–11) [C] |

Variant [DE_model_RKI_2020_short_gt.toml](models/DE_model_RKI_2020_short_gt.toml) 🔄 (shorter generation time): exposure → infectiousness 1 + Poisson(λ = 3.0), asymptomatic infectiousness → recovery Poisson(λ = 6.0), symptom onset → recovery Poisson(λ = 4.0) [A]; transmission probability 🚧 refit 🔄 (`config: 0.077`, reference rate scaled by the ratio of mean infectious periods).

## 6. Population analyses (Fig 4–7)

No disease parameters are involved.

- **Fig 4–6:** Germany (84,563,726 individuals).
- **Fig 7:** Saarland, GEMS defaults.

## 7. Configuration files

| File | Used by | Differences from [DE_model_R0_3.28.toml](models/DE_model_R0_3.28.toml) |
|---|---|---|
| [DE_model_R0_3.28.toml](models/DE_model_R0_3.28.toml) 🔄 | 5.2 | reference |
| [DE_model_R0_3.28_Heinsberg_Seed.toml](models/DE_model_R0_3.28_Heinsberg_Seed.toml) 🔄 | 5.1 | seeding in Gangelt |
| [DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml](models/DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml) 🔄 | 5.4 | seeding in Gangelt; municipality contacts 0.5; transmission rate; 1095 days |
| [DE_model_RKI_2020.toml](models/DE_model_RKI_2020.toml) 🔄 | 5.7 | 2020 dates; no initial infections; 42 days |
| [DE_model_RKI_2020_short_gt.toml](models/DE_model_RKI_2020_short_gt.toml) 🔄 | 5.7 | as above; latent and infectious periods; transmission rate |
| <a id="p13"></a>🚧 [DE_model_R0_3.28_Beckingen_Seed.toml](models/DE_model_R0_3.28_Beckingen_Seed.toml) | unused | 🚧 delete |

| Parameter | TOML key |
|---|---|
| Seeding | `[Simulation.StartCondition]` |
| Duration | `[Simulation.StopCriterion]` |
| Transmission probability | `[Pathogens.Covid19.transmission_function.parameters] transmission_rate` |
| Durations | `[Pathogens.Covid19.progressions.<Category>.<from>_to_<to>]` |
| Death probability | `[Pathogens.Covid19.progressions.Critical] death_probability` |
| Age groups, progression matrix | `[Pathogens.Covid19.progression_assignment.parameters]` |
| Contacts | `[Settings.<SettingType>.contact_sampling_method.parameters] contactparameter` |

## 8. References

🚧 🔍 = to verify ([V1–V6](#verify)) · no marker = adopted from the WSC 2025 model.

- Byambasuren, O., et al. (2020). Estimating the extent of asymptomatic COVID-19 and its potential for community transmission. *JAMMI* 5(4):223–234. https://doi.org/10.3138/jammi-2020-0030
- Byrne, A. W., et al. (2020). Inferred duration of infectious period of SARS-CoV-2. *BMJ Open* 10(8):e039856. https://doi.org/10.1136/bmjopen-2020-039856
- 🔍 Destatis (Statistisches Bundesamt). Population by age, Germany, 31.12.2019. GENESIS-Online table 12411-0005. https://www-genesis.destatis.de
- 🔍 Karagiannidis, C., et al. (2020). Case characteristics, resource use, and outcomes of 10 021 patients with COVID-19 admitted to 920 German hospitals. *Lancet Respir Med* 8(9):853–862. https://doi.org/10.1016/S2213-2600(20)30316-7
- Liu, Y., et al. (2020). The reproductive number of COVID-19 is higher compared to SARS coronavirus. *J Travel Med* 27(2):taaa021. https://doi.org/10.1093/jtm/taaa021
- Mossong, J., et al. (2008). Social contacts and mixing patterns relevant to the spread of infectious diseases. *PLoS Med* 5(3):e74. https://doi.org/10.1371/journal.pmed.0050074
- Ponge, J., et al. (2024). A standardized framework for modeling non-pharmaceutical interventions in individual-based infectious disease simulations. *WSC 2024*, 1106–1117. https://doi.org/10.1109/WSC63780.2024.10838778
- 🔍 Salje, H., et al. (2020). Estimating the burden of SARS-CoV-2 in France. *Science* 369(6500):208–211. https://doi.org/10.1126/science.abc3517
- 🔍 Schilling, J., Lehfeld, A.-S., Schumacher, D., Ullrich, A., Diercke, M., Buda, S., Haas, W., RKI COVID-19 Study Group (2020). Krankheitsschwere der ersten COVID-19-Welle in Deutschland basierend auf den Meldungen gemäß Infektionsschutzgesetz. *Journal of Health Monitoring* 5(S11). https://doi.org/10.25646/7169
- 🔍 Verity, R., et al. (2020). Estimates of the severity of coronavirus disease 2019: a model-based analysis. *Lancet Infect Dis* 20(6):669–677. https://doi.org/10.1016/S1473-3099(20)30243-7
- Wu, Y., et al. (2022). Incubation period of COVID-19 caused by unique SARS-CoV-2 strains. *JAMA Netw Open* 5(8):e2228008. https://doi.org/10.1001/jamanetworkopen.2022.28008
