# GEMS: The German Epidemic Microsimulation System — Paper Examples

This repository contains the code and model configurations behind the example applications in

> Ponge, J., Patzner, J., Suer, J., et al. *GEMS: The German Epidemic Microsimulation System.*

All simulations use [GEMS](https://github.com/IMMIDD/GEMS) **v1.0.1** (pinned in [Manifest.toml](Manifest.toml)).
This README documents the full model parameterization of every example: population, contacts, transmission, natural disease history, and intervention setup.
Most examples share one **reference parameterization** ([Section 4](#4-reference-parameterization-covid-19-like-pathogen)); each example in [Section 5](#5-example-specific-parameterization) lists only what it changes.

---

## 0. 🚧 Open items (remove before release)

Everything marked 🚧 in this README is not final. Callouts use these labels:

- **DECISION**: we have not decided yet
- **UPDATE**: the decision is made, but configs, results, or text still need to change
- **INPUT**: needs information we have not written down yet

| ID | Type | Item | Affects |
|----|------|------|---------|
| [P1](#p1) | DECISION | Death pathway. **B1:** built-in categories only; every death goes through ICU, so ICU numbers for the elderly are overestimated. **B2:** custom "dies in hospital without ICU" progression category in this repo. | 4.4, 5.4 (Tab. 2) |
| [P2](#p2) | DECISION | Share of `Severe` progressions (bedridden, never hospitalized): 0, or a stated assumption | 4.4 |
| [P3](#p3) | UPDATE | Rebuild the progression matrix from age-specific IFR/IHR (Verity et al. 2020, 10-year bins) and set `death_probability` to an ICU fatality (Karagiannidis et al. 2020). Then refit transmission rates and rerun Fig 10/11, Tab 2, and RKI validation | 4.3, 4.4, 4.6, 5.x |
| [P4](#p4) | DECISION | GEMS adds +1 day to every latent period, so mean incubation is **7.55 d**, not the 6.55 d from Wu et al. 2022. Either set `exposure_to_infectiousness_onset = 3.55` or report 7.55 d | 4.5, 4.6 |
| [P5](#p5) | INPUT | Sources for severity, hospital, ICU, and death durations (currently GEMS defaults or assumptions without a source) | 4.5 |
| [P6](#p6) | DECISION | Asymptomatic infectious period: currently 13.4 d, i.e. 2 d shorter than symptomatic (15.4 d). The WSC 2025 model used the same length for both | 4.5 |
| [P7](#p7) | INPUT | Why the multi-criteria example uses 0.5 instead of 4.0 municipality contacts | 5.4 |
| [P8](#p8) | DECISION | Transmission-rate fit targets mean R_eff over days 1–10 (not R_init as reported in the paper). The Gangelt-seeded config is fitted from only 10 seed infections. Keep, or refit against R_init? | 4.3, 5.4 |
| [P9](#p9) | DECISION | Include the RKI validation (first wave 2020) in the paper/README? | 5.7 |
| [P10](#p10) | INPUT | Beryl (RSV): link to code repository and parameter table | 5.5 |
| [P11](#p11) | UPDATE | Fill in implied/realized quantities (IFR, IHR, ICU share, R_init per config) after the rerun | 4.6 |
| [P12](#p12) | UPDATE | `runall()` builds Fig 7 from the DE population, while the caption says Saarland. Fix `runall()` or document the direct `contacts()` call | 2, 6 |
| [P13](#p13) | UPDATE | Delete the unused config [DE_model_R0_3.28_Beckingen_Seed.toml](models/DE_model_R0_3.28_Beckingen_Seed.toml) | 7 |

**Manuscript fixes (not README):**

| ID | Item |
|----|------|
| M1 | Line 365 says contact numbers are "calibrated with empirical contact matrices from the COVIMOD study". Example configs use POLYMOD-based numbers, and Fig 7 does not depend on contact numbers at all. Suggested rewording: Fig 7 reflects the population structure under uniform within-setting mixing; comparison with survey data is in [46]; contact numbers are listed in this README. |
| M2 | Sec. 3.1.2 calls the baseline "similar to the preceding example", but it is seeded with 0.1 % random infections, not in Gangelt. The Fig 12B caption ("the full code") implies the listing reproduces Fig 12A, but the listing runs GEMS defaults on an on-the-fly population. |
| M3 | Sec. 3.1.4 does not mention the deviations from the reference parameterization (municipality contacts, refitted transmission rate, 3-year horizon, one seed per scenario). |
| M4 | Sec. 3.1.5 still contains the placeholder "(add table)" (line 854). |
| M5 | Fig 8B describes an ICU survival path. With the current `death_probability = 1.0`, every ICU patient dies, so ICU admissions equal deaths in Tab 2. Resolved by [P3](#p3). |

---

## 1. Repository layout

| Path | Content |
|------|---------|
| [models/](models/) | TOML configuration files (see [Section 7](#7-configuration-files)) |
| [src/GEMSPaper.jl](src/GEMSPaper.jl) | Package entry point and `runall()` |
| `src/sc_*.jl` | One file per example scenario |
| [src/demographics.jl](src/demographics.jl), [src/settings.jl](src/settings.jl), [src/mobility.jl](src/mobility.jl), [src/contacts.jl](src/contacts.jl) | Population and structural analyses (Fig 4–7) |
| [src/benchmarking.jl](src/benchmarking.jl) | Performance benchmark (Fig 15) |
| `src/rki_*.jl` | RKI data handling and seeding for the validation scenario |
| [fit_r0.jl](fit_r0.jl) | Transmission-rate fitting |

## 2. Reproducing the results

Requires Julia 1.12 (the Manifest was generated with 1.12.6). In the repository root:

```julia
pkg> activate .
pkg> instantiate
```

```julia
using GEMSPaper
runall()                  # all analyses and scenarios, written to results/<timestamp>/
```

Each example can also be run on its own: `scenario_outbreak()`, `scenario_self_isolation()`, `scenario_sweeping()`, `scenario_testing()`, `benchmark()`, `contacts()`, `demographics()`, `mobility()`, `settinganalyses()`.

All simulations use fixed seeds (see [Section 3](#3-overview-examples-code-and-configuration)). GEMS derives one random-number generator per thread from the simulation seed, so **bit-identical results require the same number of Julia threads**.

> [!WARNING]
> <a id="p12"></a>🚧 **UPDATE (P12):** Fig 7 was generated by calling `contacts()` directly (Saarland, GEMS defaults). `runall()` passes the DE simulation instead and would produce a Germany-wide version.

## 3. Overview: examples, code, and configuration

| Paper element | Script → function | Config | Population | Horizon | Runs / seeds |
|---|---|---|---|---|---|
| Fig 4–6 (population, demographics, commuting) | `demographics.jl`, `settings.jl`, `mobility.jl` | none (population only) | DE | – | – |
| Fig 7 (contact patterns) | [contacts.jl](src/contacts.jl) → `contacts()` | GEMS defaults (structure only, see [Sec. 6](#6-population-and-structural-analyses-fig-47)) | SL | – | seed 0 |
| Sec. 3.1.1 Outbreak projection (Fig 10, 11) | [sc_outbreak.jl](src/sc_outbreak.jl) → `scenario_outbreak()` | [DE_model_R0_3.28_Heinsberg_Seed.toml](models/DE_model_R0_3.28_Heinsberg_Seed.toml) | DE | 365 d | 1 / seed 0 |
| Sec. 3.1.2 Scenario modeling (Fig 12A) | [sc_self_isolation.jl](src/sc_self_isolation.jl) → `scenario_self_isolation()` | [DE_model_R0_3.28.toml](models/DE_model_R0_3.28.toml) | DE | 365 d | 1 per scenario / seed 0 |
| Fig 12B (code listing) | [sc_self_isolation_screenshot_file.jl](src/sc_self_isolation_screenshot_file.jl) | GEMS defaults | on-the-fly, 100,000 | 365 d | seed 0 |
| Sec. 3.1.3 Parameter sweep (Fig 13) | [sc_sweeping.jl](src/sc_sweeping.jl) → `scenario_sweeping()` | pathogen defined in script + GEMS defaults | SL | 365 d | 41 × 41 grid × 3 / seeds 1–3 |
| Sec. 3.1.4 Multi-criteria evaluation (Fig 14, Tab 2) | [sc_testing.jl](src/sc_testing.jl) → `scenario_testing()` | [DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml](models/DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml) | DE | 1095 d | 1 per scenario / seeds 0, 1, 2 |
| Sec. 3.1.5 Beryl (RSV) | 🚧 not in this repository ([P10](#p10)) | – | – | – | – |
| Sec. 3.2 Performance (Fig 15) | [benchmarking.jl](src/benchmarking.jl) → `benchmark()` | GEMS defaults | 16 states + DE; on-the-fly 0.1–50 M | 365 d | 11 per population (median) / seeds 1–11 |
| 🚧 RKI validation ([P9](#p9)) | [sc_rki_validation.jl](src/sc_rki_validation.jl) | [DE_model_RKI_2020.toml](models/DE_model_RKI_2020.toml), [DE_model_RKI_2020_short_gt.toml](models/DE_model_RKI_2020_short_gt.toml) | DE | 42 d | seed 0 |

## 4. Reference parameterization (COVID-19-like pathogen)

Defined in [models/DE_model_R0_3.28.toml](models/DE_model_R0_3.28.toml). Many values follow the earlier GEMS study [Ponge et al. 2025 (WSC)](https://doi.org/10.1109/WSC68292.2025.11339134); the tables cite the primary sources.

**Source tags:** `[L]` literature · `[C]` calibrated in this repository · `[A]` assumption · `[D]` GEMS default value · `[I]` illustrative choice for a methods demonstration

### 4.1 Simulation and seeding

| Parameter | Value | Source |
|---|---|---|
| Time step | 1 day | [D] |
| Duration | 365 days (`TimesUp`) | [I] |
| Initial infections | 0.1 % of the population, drawn uniformly at random (`InfectedFraction`) | [I] as in Ponge et al. 2025 |
| Global setting | disabled | [D] |

### 4.2 Contacts

Contacts are drawn for each infectious individual and day. The number of contacts per setting level follows a Poisson distribution with the mean given below (`ContactparameterSampling`). Contact partners are chosen uniformly among the members of that setting, including members of its subordinate settings. Quarantined, isolated, and severely ill individuals have household contacts only; hospitalized individuals have none.

| Setting type | Level | Mean daily contacts | Total per setting type | Source |
|---|---|---|---|---|
| Household | – | 3.3 | 3.3 | [L] Mossong et al. 2008 |
| Workplace | Office | 1.4 | 2.8 | [L] Mossong et al. 2008 (total); split across levels [A] |
| | Department | 0.7 | | |
| | Workplace | 0.35 | | |
| | WorkplaceSite | 0.35 | | |
| School | SchoolClass | 0.9 | 1.8 | [L] Mossong et al. 2008 (total); split across levels [A] |
| | SchoolYear | 0.45 | | |
| | School | 0.225 | | |
| | SchoolComplex | 0.225 | | |
| Municipality | – | 4.0 | 4.0 | [L] Mossong et al. 2008 |

### 4.3 Transmission

| Parameter | Value | Source |
|---|---|---|
| Transmission function | `ConstantTransmissionRate` (same per-contact probability for all settings, ages, and days) | [D] |
| Per-contact transmission probability | 0.03142 | [C] fitted to R0 = 3.28 |
| Target R0 | 3.28 | [L] Liu et al. 2020 |

**Fitting procedure** ([fit_r0.jl](fit_r0.jl)): linear search over the transmission rate. For each candidate rate, 10 simulations of 70 days are run. The fitted quantity is the mean daily effective reproduction number over days 1–10, averaged across runs; the search stops within a tolerance of 0.02.

> [!CAUTION]
> <a id="p8"></a>🚧 **DECISION (P8):** The paper reports R_init (secondary infections of the first 5 % of cases), but the fit targets R_eff over days 1–10. For Gangelt-seeded configs (10 seeds), days 1–10 contain very few infections. Keep this procedure or refit against R_init?
>
> 🚧 **UPDATE ([P3](#p3)):** Refit all transmission rates once the progression matrix changes.

### 4.4 Progression assignment

Every infected individual is assigned one progression category at exposure, with age-specific probabilities (`AgeBasedProgressionAssignment`). Only the `Critical` category can end in death, and only after ICU admission (GEMS v1.0.1).

> [!WARNING]
> <a id="p3"></a>🚧 **UPDATE (P3):** The matrix below is the **current** config. It was migrated from the GEMS v0.4.3 configuration of Ponge et al. (2025). That configuration's fatality values (0.2 / 0.9 / 5 / 14.8 %) were split into `Severe` (70 %), `Hospitalized` (21 %), and `Critical` (9 %). As a result, the implied IFR is about 10× lower than intended, and `death_probability = 1.0` means every ICU patient dies. **It will be replaced** by a matrix derived as follows:
>
> | Quantity | Source |
> |---|---|
> | IFR by age (10-year bins) | [L] Verity et al. 2020, Table 1 |
> | Infection–hospitalization ratio (IHR) by age | [L] Verity et al. 2020, Table 3 |
> | ICU fatality → `death_probability` | [L] Karagiannidis et al. 2020 (0.53, invasively ventilated patients, Germany, first wave) |
> | Asymptomatic share | [L] Byambasuren et al. 2020 (0.17) |
> | ICU admission among hospitalized (B2 only) | [L] Salje et al. 2020 |
>
> Derivation (B1): `Critical = IFR / d`, `Hospitalized = IHR − Critical`, `Symptomatic = 1 − 0.17 − IHR − Severe`.

> [!CAUTION]
> <a id="p1"></a>🚧 **DECISION (P1): B1 or B2?**
> - **B1:** Built-in categories only. Every death passes through ICU. The derivation above requires 49 % (age 70–79) and 80 % (age 80+) of hospitalized patients to go to ICU; Salje et al. 2020 report about 19 % overall. Deaths are correct, but ICU admissions and ICU load (Tab 2) are overestimated.
> - **B2:** Add a custom progression category "hospitalized, dies without ICU" in this repository. Deaths, hospitalizations, and ICU admissions can then all match their sources.
>
> <a id="p2"></a>🚧 **DECISION (P2):** `Severe` share (bedridden at home, never hospitalized): 0, or a stated assumption?

**Current values (to be replaced):**

| Age group | Asymptomatic | Symptomatic | Severe | Hospitalized | Critical |
|---|---|---|---|---|---|
| 0–39 | 0.17 | 0.828 | 0.0014 | 0.00042 | 0.00018 |
| 40–59 | 0.17 | 0.821 | 0.0063 | 0.00189 | 0.00081 |
| 60–79 | 0.17 | 0.780 | 0.0350 | 0.01050 | 0.00450 |
| 80+ | 0.17 | 0.682 | 0.1036 | 0.03108 | 0.01332 |

`Critical.death_probability = 1.0` 🚧

### 4.5 Transition durations

All durations are drawn from Poisson distributions with the given mean, in days. GEMS v1.0.1 **adds 1 day** to every exposure → infectiousness-onset draw so that no one becomes infectious on the day of exposure. The "effective mean" column includes this offset.

| Category | Transition | Config mean | Effective mean | Source |
|---|---|---|---|---|
| all | exposure → infectiousness onset | 4.55 | **5.55** 🚧 | [L] Wu et al. 2022 (incubation) minus [L] Byrne et al. 2020 (pre-symptomatic period) |
| Asymptomatic | infectiousness onset → recovery | 13.4 | 13.4 | 🚧 [A], see [P6](#p6) |
| Symptomatic, Severe, Hospitalized, Critical | infectiousness onset → symptom onset | 2.0 | 2.0 | [L] Byrne et al. 2020 |
| Symptomatic | symptom onset → recovery | 13.4 | 13.4 | [L] Byrne et al. 2020 |
| Severe, Hospitalized, Critical | symptom onset → severeness onset | 3.0 | 3.0 | 🚧 [D] ([P5](#p5)) |
| Severe | severeness onset → severeness offset | 7.0 | 7.0 | 🚧 [D] |
| Severe, Hospitalized, Critical | severeness offset → recovery | 4.0 | 4.0 | 🚧 [D] |
| Hospitalized, Critical | severeness onset → hospital admission | 4.0 | 4.0 | 🚧 [D] |
| Hospitalized | hospital admission → hospital discharge | 7.0 | 7.0 | 🚧 [D] |
| Hospitalized, Critical | hospital discharge → severeness offset | 3.0 | 3.0 | 🚧 [D] |
| Critical | hospital admission → ICU admission | 2.0 | 2.0 | 🚧 [D] |
| Critical (survives) | ICU admission → ICU discharge | 7.0 | 7.0 | 🚧 [D] |
| Critical (survives) | ICU discharge → hospital discharge | 7.0 | 7.0 | 🚧 [D] |
| Critical (dies) | ICU admission → death | 10 | 10 | 🚧 [D] |

> [!CAUTION]
> <a id="p4"></a>🚧 **DECISION (P4):** The resulting mean incubation period is 1 + 4.55 + 2.0 = **7.55 d**, not the 6.55 d reported by Wu et al. 2022 (the WSC 2025 model, GEMS v0.4.3, had no +1 offset). Change `exposure_to_infectiousness_onset` to 3.55 (requires refit and rerun) or report 7.55 d.
>
> <a id="p5"></a>🚧 **INPUT (P5):** Find literature sources for the severity, hospital, and ICU durations, or keep them as GEMS defaults and label them `[A]`.
>
> <a id="p6"></a>🚧 **DECISION (P6):** Asymptomatic cases are infectious for 13.4 d on average, symptomatic cases for 2.0 + 13.4 = 15.4 d. Intended?

### 4.6 Implied and realized quantities

Derived from the parameters above, so the inputs can be checked against the outputs.

| Quantity | Value | Note |
|---|---|---|
| Mean incubation period | 7.55 d 🚧 | [P4](#p4) |
| Mean infectious period, symptomatic | 15.4 d | |
| Mean infectious period, asymptomatic | 13.4 d 🚧 | [P6](#p6) |
| Implied IFR by age (current matrix) | 0.018 % / 0.081 % / 0.45 % / 1.33 % 🚧 | changes with [P3](#p3) |
| Realized IFR, IHR, ICU share by age | 🚧 | [P11](#p11), after rerun |
| Realized R_init per config | 🚧 | [P11](#p11), after rerun |

> [!WARNING]
> <a id="p11"></a>🚧 **UPDATE (P11):** Fill in the realized quantities from the rerun simulations.

## 5. Example-specific parameterization

Interventions use the Trigger–Strategy–Measure (TriSM) framework ([Ponge et al. 2024](https://doi.org/10.1109/WSC63780.2024.10838778)). The intervention tables list them as trigger → strategy (condition) → measures (offset).

### 5.1 Outbreak projection (Sec. 3.1.1, Fig 10, 11)

Config: [DE_model_R0_3.28_Heinsberg_Seed.toml](models/DE_model_R0_3.28_Heinsberg_Seed.toml)

| Parameter | Reference | This example |
|---|---|---|
| Initial infections | 0.1 % random | **10 infections in Gangelt** (municipality key 05370008, county Heinsberg), the first major German outbreak in 2020 (`RegionalSeeds`) |

Interventions: none.

### 5.2 Scenario modeling (Sec. 3.1.2, Fig 12)

Config: [DE_model_R0_3.28.toml](models/DE_model_R0_3.28.toml) (reference, unchanged; 0.1 % random initial infections)

| Scenario | Trigger | Strategy (condition) | Measures |
|---|---|---|---|
| Baseline | – | – | – |
| Scenario | `SymptomTrigger` | Self Isolation (none) | `SelfIsolation(14)` (immediately): no out-of-household contacts for 14 days |

The code listing in Fig 12B ([sc_self_isolation_screenshot_file.jl](src/sc_self_isolation_screenshot_file.jl)) shows the same intervention setup but runs with GEMS defaults ([Section 5.6](#56-performance-benchmark-sec-32-fig-15)) on a 100,000-agent on-the-fly population. It illustrates the API and does not reproduce Fig 12A.

### 5.3 Parameter sweep (Sec. 3.1.3, Fig 13)

No config file. The pathogen is defined in [sc_sweeping.jl](src/sc_sweeping.jl). Contacts, seeding, and run length come from the GEMS defaults ([Section 5.6](#56-performance-benchmark-sec-32-fig-15)). Population: Saarland (SL).

| Parameter | Value | Source |
|---|---|---|
| Transmission probability | 0.04 | [I] |
| Progression categories | Asymptomatic, Symptomatic (one age group, 0–120) | [I] |
| Asymptomatic: exposure → infectiousness onset → recovery | 0 (+1) → 7 days, fixed | [I] |
| Symptomatic: exposure → infectiousness onset → symptom onset → recovery | 0 (+1) → 2 → 7 days, fixed | [I] |
| Symptomatic share (x-axis) | 0 to 1, step 0.025 | swept |
| Adherence (y-axis) | 0 to 1, step 0.025 | swept |
| Repetitions | 3 per combination (seeds 1–3); 41 × 41 × 3 = 5,043 runs | |
| Outcome | R_init (mean of 3 runs) | |

Durations here are fixed values, not Poisson draws.

| Trigger | Strategy (condition) | Measures |
|---|---|---|
| `SymptomTrigger` | Self Isolation (adheres if `rand() ≤ adherence`, evaluated once per symptomatic individual) | `SelfIsolation(14)` (immediately) |

### 5.4 Multi-criteria strategy evaluation (Sec. 3.1.4, Fig 14, Tab 2)

Config: [DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml](models/DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml)

| Parameter | Reference | This example |
|---|---|---|
| Initial infections | 0.1 % random | 10 infections in Gangelt (05370008) |
| Municipality contacts | 4.0 | **0.5** 🚧 |
| Transmission probability | 0.03142 | **0.06888** ([C], refitted to R0 = 3.28 with [fit_r0.jl](fit_r0.jl)) |
| Duration | 365 days | **1095 days** |
| Seeds | – | 0 (Isolation Only), 1 (Isolation & Testing), 2 (Isolation, Testing & Tracing) |

> [!CAUTION]
> <a id="p7"></a>🚧 **INPUT (P7):** Why municipality contacts were reduced to 0.5 for this example.

**Test type** "Antigen Test": sensitivity 0.80, specificity 0.99 [I]; positive results are reported.

| Scenario | Trigger | Strategy (condition) | Measures (offset) → follow-up |
|---|---|---|---|
| Isolation Only | `SymptomTrigger` | Self Isolation | `SelfIsolation(14)` (0) |
| Isolation & Testing | `SymptomTrigger` | Self Isolation | `SelfIsolation(14)` (0); `Test(Antigen)` (0), no follow-up |
| Isolation, Testing & Tracing | `SymptomTrigger` | Isolation & Testing (not yet detected) | `SelfIsolation(14)` (0); `Test(Antigen)` (0) → positive: *Contact Identification*; negative: *Cancel Isolation* |
| | | Contact Identification | `FindSettingMembers(Household)` (1 d); `FindSettingMembers(Office)` (1 d, if working); `FindSettingMembers(SchoolClass)` (1 d, if student); each excluding the index case → *Isolate & Test Contacts* for every member found |
| | | Isolate & Test Contacts (not yet detected) | `SelfIsolation(14)` (0); `Test(Antigen)` (0) → negative: *Cancel Isolation*; positive: remains isolated, no further tracing |
| | | Cancel Isolation | `CancelSelfIsolation` (0) |

The tracing scenario corresponds to Fig 9 of the paper.

> [!WARNING]
> 🚧 **UPDATE ([P3](#p3)):** Tab 2 health-care outcomes (hospitalizations, ICU, deaths, loads) will change with the new progression matrix.

### 5.5 Beryl — RSV (Sec. 3.1.5)

> [!CAUTION]
> <a id="p10"></a>🚧 **INPUT (P10):** Link to the code repository and add the parameter table (or a link to it). Scenario definitions: [RespiCompass RSV Round 1 2025/26](https://github.com/european-modelling-hubs/RespiCompass/blob/main/Previous_Rounds/2025-2026_round_1/round1_2526_rsv.md).

### 5.6 Performance benchmark (Sec. 3.2, Fig 15)

GEMS default parameterization ([data/DefaultConf.toml @ v1.0.1](https://github.com/IMMIDD/GEMS/blob/v1.0.1/data/DefaultConf.toml)), also used by Fig 7, Fig 12B, and (for contacts, seeding, and duration) the parameter sweep:

| Parameter | Value |
|---|---|
| Initial infections | 0.1 % random |
| Duration | 365 days |
| Transmission probability | 0.2 |
| Contacts | 1.0 per day on every setting level |
| Progression matrix (ages ≤14 / 15–65 / 66+) | Asymptomatic 0.40 / 0.25 / 0.15; Symptomatic 0.58 / 0.60 / 0.40; Severe 0.01 / 0.11 / 0.25; Hospitalized 0.007 / 0.03 / 0.12; Critical 0.003 / 0.01 / 0.08 |
| Durations (Poisson means) | exposure → infectiousness 1 (+1); infectiousness → symptoms 1; symptom → recovery 7 (asymptomatic: infectiousness → recovery 8); severeness onset 1; severeness duration 7; to hospital 2; hospital stay 7; to ICU 2; ICU stay 7; ICU → death 10; post-hospital 3 + 4 |
| Critical death probability | 0.3 |
| On-the-fly default population | 100,000 agents, mean household size 3, office size 5, school size 100 |

Benchmark protocol: one warm-up run (excluded), then 11 runs per population (seeds 1–11), reported as medians. Populations: all 16 German states and Germany (Gesyland), plus on-the-fly populations of 0.1, 1, 5, 10, 20, and 50 million agents. Hardware: see Sec. 3.2 of the paper.

### 5.7 🚧 RKI validation, first wave 2020

> [!CAUTION]
> <a id="p9"></a>🚧 **DECISION (P9):** Include in the paper/README? If so, complete this section.

Config: [DE_model_RKI_2020.toml](models/DE_model_RKI_2020.toml), the reference parameterization with these changes:

| Parameter | Value |
|---|---|
| Simulation period | 2020-02-10 to 2020-03-22 (42 days) |
| Initial infections | none; seeding during the run from imported cases (RKI case data) plus the Heinsberg cluster, back-dated by an assumed 7-day incubation period |
| Transmission probability | grid search against national weekly RKI cases (weeks 10–11) |

Variant [DE_model_RKI_2020_short_gt.toml](models/DE_model_RKI_2020_short_gt.toml) shortens the generation time: exposure → infectiousness 3.0 (+1), asymptomatic infectious period 6.0, symptom onset → recovery 4.0; transmission probability 0.077.

## 6. Population and structural analyses (Fig 4–7)

These figures use the GEMS population models only; no disease parameters are involved.

- **Fig 4–6:** Germany (Gesyland synthetic population, 84,563,726 individuals).
- **Fig 7:** Saarland. Contact pairs are sampled from the setting structure under uniform within-setting mixing, and the figure shows shares of contacts. The number of contacts per setting does not affect the figure; it reflects only which individuals share households, schools, workplaces, and municipalities.

## 7. Configuration files

| File | Used by | Changes relative to [DE_model_R0_3.28.toml](models/DE_model_R0_3.28.toml) |
|---|---|---|
| [DE_model_R0_3.28.toml](models/DE_model_R0_3.28.toml) | 5.2 | reference |
| [DE_model_R0_3.28_Heinsberg_Seed.toml](models/DE_model_R0_3.28_Heinsberg_Seed.toml) | 5.1 | seeding in Gangelt |
| [DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml](models/DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml) | 5.4 | seeding in Gangelt; municipality contacts 0.5; transmission 0.06888; 1095 days |
| [DE_model_RKI_2020.toml](models/DE_model_RKI_2020.toml) | 5.7 | 2020 dates; no initial infections; 42 days |
| [DE_model_RKI_2020_short_gt.toml](models/DE_model_RKI_2020_short_gt.toml) | 5.7 | as above; shorter latent and infectious periods; transmission 0.077 |
| 🚧 [DE_model_R0_3.28_Beckingen_Seed.toml](models/DE_model_R0_3.28_Beckingen_Seed.toml) | unused | <a id="p13"></a>🚧 **UPDATE (P13):** delete |

**Where parameters live in the TOML files:**

| Parameter | TOML key |
|---|---|
| Seeding | `[Simulation.StartCondition]` |
| Duration | `[Simulation.StopCriterion]` |
| Transmission probability | `[Pathogens.Covid19.transmission_function.parameters] transmission_rate` |
| Durations | `[Pathogens.Covid19.progressions.<Category>.<from>_to_<to>]` |
| Death probability | `[Pathogens.Covid19.progressions.Critical] death_probability` |
| Progression matrix | `[Pathogens.Covid19.progression_assignment.parameters]` |
| Contacts | `[Settings.<SettingType>.contact_sampling_method.parameters] contactparameter` |

## 8. Known simplifications

- 🚧 Death pathway: depends on [P1](#p1).
- GEMS v1.0.1 adds 1 day to every exposure → infectiousness-onset draw.
- Transmission probability is constant across settings, ages, and time.
- Except for the RKI validation, the examples demonstrate GEMS' capabilities. Their parameter values are not calibrated to observed data, and differences between scenarios are within-model comparisons.

## 9. References

- Byambasuren, O., et al. (2020). Estimating the extent of asymptomatic COVID-19 and its potential for community transmission. *JAMMI* 5(4):223–234. https://doi.org/10.3138/jammi-2020-0030
- Byrne, A. W., et al. (2020). Inferred duration of infectious period of SARS-CoV-2. *BMJ Open* 10(8):e039856. https://doi.org/10.1136/bmjopen-2020-039856
- Karagiannidis, C., et al. (2020). Case characteristics, resource use, and outcomes of 10 021 patients with COVID-19 admitted to 920 German hospitals. *Lancet Respir Med* 8(9):853–862. https://doi.org/10.1016/S2213-2600(20)30316-7
- Liu, Y., et al. (2020). The reproductive number of COVID-19 is higher compared to SARS coronavirus. *J Travel Med* 27(2):taaa021. https://doi.org/10.1093/jtm/taaa021
- Mossong, J., et al. (2008). Social contacts and mixing patterns relevant to the spread of infectious diseases. *PLoS Med* 5(3):e74. https://doi.org/10.1371/journal.pmed.0050074
- Ponge, J., et al. (2024). A standardized framework for modeling non-pharmaceutical interventions in individual-based infectious disease simulations. *WSC 2024*, 1106–1117. https://doi.org/10.1109/WSC63780.2024.10838778
- Ponge, J., Patzner, J., Hellingrath, B., Karch, A. (2025). Targeted household quarantining: Enhancing the efficiency of epidemic responses. *WSC 2025*, 953–964. https://doi.org/10.1109/WSC68292.2025.11339134
- Salje, H., et al. (2020). Estimating the burden of SARS-CoV-2 in France. *Science* 369(6500):208–211. https://doi.org/10.1126/science.abc3517
- Verity, R., et al. (2020). Estimates of the severity of coronavirus disease 2019: a model-based analysis. *Lancet Infect Dis* 20(6):669–677. https://doi.org/10.1016/S1473-3099(20)30243-7
- Wu, Y., et al. (2022). Incubation period of COVID-19 caused by unique SARS-CoV-2 strains. *JAMA Netw Open* 5(8):e2228008. https://doi.org/10.1001/jamanetworkopen.2022.28008
