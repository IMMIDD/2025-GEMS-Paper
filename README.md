# GEMS: The German Epidemic Microsimulation System — Paper Examples

This repository contains the code and model configurations behind the example applications in

> Ponge, J., Patzner, J., Suer, J., et al. *GEMS: The German Epidemic Microsimulation System.*

All simulations use [GEMS](https://github.com/IMMIDD/GEMS) **v1.0.2**, the branch [`release_1.0.2`](https://github.com/IMMIDD/GEMS/tree/release_1.0.2) at commit [`7f3b79e`](https://github.com/IMMIDD/GEMS/commit/7f3b79e51752f1d316a836c2e9d78a76d73572fe) (pinned in [Manifest.toml](Manifest.toml); not yet tagged, the package version there is still 1.0.1). It gives every setting its own copy of its contact sampling method, which the age-based municipality contacts of the reference model need ([Section 4.2](#42-contacts)).
This README documents the parameterization of every example: population, contacts, transmission, disease progression, and interventions, as well as the validation of the population model against official statistics ([Section 6](#6-population-analyses-and-validation)).
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
   - [5.1 First wave 2020 (Sec. 3.1.1, Fig 11–12)](#51-first-wave-2020-sec-311-fig-1112)
   - [5.2 Scenario modeling (Sec. 3.2.1, Fig 15–16)](#52-scenario-modeling-sec-321-fig-1516)
   - [5.3 Parameter sweep (Sec. 3.2.2, Fig 17)](#53-parameter-sweep-sec-322-fig-17)
   - [5.4 Multi-criteria strategy evaluation (Sec. 3.2.3, Fig 18, Tab 2)](#54-multi-criteria-strategy-evaluation-sec-323-fig-18-tab-2)
   - [5.5 Performance benchmark (Sec. 3.3, Fig 19)](#55-performance-benchmark-sec-33-fig-19)
6. [Population analyses and validation](#6-population-analyses-and-validation)
   - [6.1 Population and contact figures (Fig 4–8)](#61-population-and-contact-figures-fig-48)
   - [6.2 Population validation](#62-population-validation)
   - [6.3 Commuting validation](#63-commuting-validation)
7. [References](#7-references)

---

## 1. Repository layout

| Path | Content |
|------|---------|
| [models/DE_model_first_wave.toml](models/DE_model_first_wave.toml) | Reference model ([Section 4](#4-reference-model-covid-19-wild-type-first-wave-2020)) |
| [src/GEMSPaper.jl](src/GEMSPaper.jl) | Package entry point and `runall()` |
| [src/infectiousness_profile.jl](src/infectiousness_profile.jl) | `InfectiousnessProfile` transmission function |
| `src/sc_*.jl` | One file per example scenario |
| [src/figure_style.jl](src/figure_style.jl) | Sizes, fonts, colors, and file names of the figures; titling, combining, and saving of panels |
| [src/first_wave/](src/first_wave/) | RKI data, imported and introduced infections, intervention phases of March 2020 (fitted contact reductions), calibration, comparison of municipality mixing assumptions, and figures of the first wave example |
| [src/contacts/](src/contacts/) | POLYMOD contact data: download, contacts per setting layer, household contact sampling (`HouseholdSizeSampling`), age-age matrix of the municipality contacts, derivation of the contact parameters of the config, and scaling of contact parameters |
| [src/demographics.jl](src/demographics.jl), [src/demographic_maps.jl](src/demographic_maps.jl), [src/settings_map.jl](src/settings_map.jl), [src/settings.jl](src/settings.jl), [src/mobility.jl](src/mobility.jl), [src/contacts.jl](src/contacts.jl) | Population analyses (Fig 4, 6, 8) |
| [src/population_validation.jl](src/population_validation.jl), [src/setting_validation.jl](src/setting_validation.jl) | Population validation against official statistics: age, sex, households, employment and education, setting sizes (Fig 5, [Section 6.2](#62-population-validation)) |
| [src/commuter_validation.jl](src/commuter_validation.jl) | Commuting balance and validation against the BA commuter statistics (Fig 7, [Section 6.3](#63-commuting-validation)) |
| [src/benchmarking.jl](src/benchmarking.jl) | Performance benchmark (Fig 19) |
| [calibrate_first_wave.jl](calibrate_first_wave.jl) | Calibration of the reference model (transmission scale, generation time, contact factors); not part of `runall()` |
| [derive_contact_parameters.jl](derive_contact_parameters.jl) | Derivation of the contact parameters of the reference model from POLYMOD; not part of `runall()` |
| [inputs/](inputs/) | Calibration results read by the scenarios, contact parameters derived from POLYMOD, reference data from official statistics (population, settings, commuting), county income; sources and licenses of every file in [inputs/README.md](inputs/README.md) |
| [figures/](figures/) | Drawn figures of the paper: PowerPoint deck ([drawn_figures.pptx](figures/drawn_figures.pptx)) and the panels that are combined with plots (`figures/panels/`), see [Section 2](#drawn-figures) |
| [tools/export_figures.ps1](tools/export_figures.ps1) | Exports the drawn figures from PowerPoint; not part of `runall()` |
| [tools/convert_to_tiff.ps1](tools/convert_to_tiff.ps1) | Converts the PNG figures of a folder to TIFF, see [Section 2](#tiff-files-for-the-submission); not part of `runall()` |
| `data/` | Downloaded on demand, not part of the repository: RKI case data, POLYMOD raw data, KMK class statistics ([inputs/README.md](inputs/README.md#downloaded-at-runtime-not-in-the-repository)) |

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

With a population of all 16 federal states, `demographics()` also runs the population validation (`population_validation()`) and `mobility()` the commuting validation (`commuter_validation()`). The settings map, the demographic maps, the contact patterns, and both validations save their data as CSV files, so the figures can be redrawn without loading the population: `settings_map(data = load_settings_per_pixel(folder))`, `contacts(panels = load_contact_matrices(folder))`, `demographic_maps(data = load_county_demographics(folder))`, `population_validation(data = load_validation_data(folder))` and `commuter_validation(data = load_commuter_data(folder))`.

The first wave example downloads the RKI case data (about 400 MB, checked against a SHA-256 checksum) to `data/` on its first run. Its calibrated contact factors are read from `inputs/`; to recalibrate, run [calibrate_first_wave.jl](calibrate_first_wave.jl). The population validation downloads the KMK class statistics to `data/` on its first run.

The contact parameters of the reference model are derived from the POLYMOD survey with [derive_contact_parameters.jl](derive_contact_parameters.jl). It downloads the POLYMOD raw data to `data/polymod/`, writes the parameters and the tables they come from to [inputs/contacts/](inputs/contacts/), and warns about config entries that differ from the derived values.

The parameter sweep stores its results in `sc_sweep_data.csv`; `replot_sweeps(folder)` redraws the figure and rewrites the variability summary from that file without re-running the 8,405 simulations. Likewise, `replot_first_wave(folder)` redraws the figures of the first wave, `replot_household_isolation(folder)` redraws the figures of the household isolation scenario, and `replot_testing(folder)` redraws the figure and rewrites the tables of the testing scenario from the stored results of its runs.

For a quick test, the testing scenario also runs with the population of a smaller state and an outbreak in one of its municipalities, e.g. `scenario_testing(population = "SL", municipality = 10042111)` (Beckingen, Saarland).

All simulations use fixed random seeds ([Section 3](#3-overview-examples-code-and-configuration)). Bit-identical results require the same number of Julia threads. Simulations with interventions are only exactly reproducible on a single thread in GEMS v1.0.2 (fixed in v1.1.1).

**Machines used for the paper**

| Machine | Hardware | Julia | Threads (`JULIA_NUM_THREADS`) | Used for |
|---|---|---|---|---|
| Server | AMD EPYC 9575F (64 cores, 128 threads), 755 GB DDR5 RAM (5200 MT/s, CL52), NVMe SSD | TODO (1.12?) | TODO | TODO: all examples? performance benchmark (Sec. 3.3) |
| Desktop | Intel i9-14900K, 128 GB DDR5 RAM, 4 TB NVMe SSD | TODO (1.12?) | TODO | desktop runtime of the national model (Sec. 3.3) |

### Drawn figures

The figures that are not generated from data (Fig 1–3, 9, 10, and the strategy of Fig 17) are drawn in [figures/drawn_figures.pptx](figures/drawn_figures.pptx), one figure per slide. The slides have the size of a full page of PLOS (7.5 × 8.75 in), so fonts and lines have their print sizes. On Windows with PowerPoint (closed), export them with

```powershell
powershell -ExecutionPolicy Bypass -File tools\export_figures.ps1
```

This writes a PNG of every figure (300 dpi, cropped to its width and height) to `results/<timestamp>_drawn_figures/`. The layout of a slide sets the width of its figure (full width 2250 px, text width 1560 px, panel width 750 px). The notes of a slide hold its file name, or `skip` if it is not exported (style sheet).

Slides with `panel` in the second line of their notes are panels that the Julia code combines with a plot. Julia reads them from `figures/panels/`, so PowerPoint is not needed to reproduce the results. After changing such a slide, update its panel with `-UpdatePanels`. The code listing of Fig 15 (`Fig_15B_code_listing.png`) is a screenshot of [sc_household_isolation_screenshot_file.jl](src/sc_household_isolation_screenshot_file.jl).

### TIFF files for the submission

PLOS requires the figures as TIFF files. On Windows, convert all PNG figures of a folder with

```powershell
powershell -ExecutionPolicy Bypass -File tools\convert_to_tiff.ps1 results\<folder>
```

This writes a TIFF of the same name next to every PNG (LZW compression, 300 dpi, white background, same pixels). Existing TIFF files are overwritten, subfolders are not searched. With `-Out <folder>`, the files are written to another folder. The script warns about figures above the limits of PLOS (2250 × 2625 px, 10 MB).

*Terminology:* a *random seed* initializes the random number generator of a simulation; infections that start an epidemic are called *initial infections* (at tick 0) or *introduced infections* (added during the run).

## 3. Overview: examples, code, and configuration

| Paper element | Script → function | Config | Population | Horizon | Runs / random seeds |
|---|---|---|---|---|---|
| Fig 4, 6 (population, demographics) | `demographics.jl`, `settings.jl` | none | DE | – | – |
| Fig 5 (population validation) | [population_validation.jl](src/population_validation.jl) → `population_validation()` (part of `demographics()`) | none | DE | – | – |
| Fig 7 (commuting validation) | [commuter_validation.jl](src/commuter_validation.jl) → `commuter_validation()` (part of `mobility()`) | none | DE | – | – |
| Fig 8 (contact patterns) | [contacts.jl](src/contacts.jl) → `contacts()` | reference model, plus uniform and strict mixing in the municipalities | DE | – | 1 day, all members of all settings |
| Contact parameters (Sec. 2.4, reference model) | [derive_contact_parameters.jl](derive_contact_parameters.jl) → `derive_contact_parameters()` | reference model | DE (age weights) | – | – |
| Sec. 3.1.1 First wave 2020 (Fig 11–12), incl. municipality mixing comparison | [sc_first_wave.jl](src/sc_first_wave.jl) → `scenario_first_wave()` | reference model | DE | 10.02.–31.05.2020 | see [5.1](#51-first-wave-2020-sec-311-fig-1112) |
| Sec. 3.1.2 RSV immunization 2025/26 (Fig 13–14) | separate repository: [IMMIDD/2025-RSV-scenario](https://github.com/IMMIDD/2025-RSV-scenario) | – | – | – | – |
| Sec. 3.2.1 Scenario modeling (Fig 15A, baseline readouts Fig 16) | [sc_household_isolation.jl](src/sc_household_isolation.jl) → `scenario_household_isolation()` | reference model | DE | 180 d | 1 per scenario / 0 |
| Fig 15B (code listing) | [sc_household_isolation_screenshot_file.jl](src/sc_household_isolation_screenshot_file.jl) | GEMS defaults | on-the-fly, 100,000 | 365 d | 1 per scenario / none (random) |
| Sec. 3.2.2 Parameter sweep (Fig 17) | [sc_sweeping.jl](src/sc_sweeping.jl) → `scenario_sweeping()` | reference model | SL | 365 d | 41 × 41 grid × 5 / 1–8405 |
| Sec. 3.2.3 Multi-criteria evaluation (Fig 18, Tab 2) | [sc_testing.jl](src/sc_testing.jl) → `scenario_testing()` | reference model | NRW | 365 d | 11 per scenario / 0–32 |
| Sec. 3.3 Performance (Fig 19) | [benchmarking.jl](src/benchmarking.jl) → `benchmark()` | reference model | 16 states + DE; on-the-fly 0.1–50 M | 365 d | 11 per population (median) / 1–11 |

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
| `[Population]` | n = 100,000, household size 3, office size 5, school size 100 | [D] only used for generated populations; all examples except the on-the-fly populations of the benchmark ([5.5](#55-performance-benchmark-sec-33-fig-19)) use the GEMS population models |

### 4.2 Contacts

All contact parameters are derived from the German participants of the POLYMOD survey (Mossong et al. 2008; data: Mossong et al. 2020) with [derive_contact_parameters.jl](derive_contact_parameters.jl) (code in [src/contacts/](src/contacts/)). [inputs/contacts/contact_parameters_DE.csv](inputs/contacts/contact_parameters_DE.csv) lists the value and derivation of every config entry. The config holds them rounded to 3 decimals.

**Derivation**

| Step | Rule |
|---|---|
| Layer of a contact | home → household; school → education, if the participant attends daycare, school, or university (children's diaries: `child_care`, without childminders; own diaries: full-time education); work → work, if the participant works (school contacts of working participants, e.g. teachers, also go to work); all other contacts (transport, leisure, other places, no location, school or work without membership) → municipality |
| Contacts per member and day | mean over the participants who are members of the layer (households and municipality: all participants), weighted by age (groups 0, 5, 10, 15, 20, 30, …, 70+) to the members of the layer in the GEMS DE population; daily average of weekdays (5/7) and weekend (2/7) |
| Visitors at home | home contacts beyond the other household members (household size − 1) count as municipality contacts |
| Split over setting levels | [A] contacts weighted towards the most immediate level (own office or class): each level gets half of the one below, the two outermost levels the same (1/2, 1/4, 1/8, 1/8) |

**Contact sampling per setting type**

| Setting type | Sampling method | Parameters | Source |
|---|---|---|---|
| Household (2–10 members) | `HouseholdSizeSampling` ([src/contacts/household_sampling.jl](src/contacts/household_sampling.jl)): each other present member is met once per day with probability q, by household size | q = 0.953 / 0.935 / 0.899 / 0.804 for 2 / 3 / 4 / 5+ members | [L] adapted: share of the other household members met per day by the German participants, by household size ([polymod_household_probabilities.csv](inputs/contacts/polymod_household_probabilities.csv)) |
| Household (care homes, 11+ members) | Poisson(λ) contacts drawn uniformly from the present members | λ = 2.706 | [L] adapted: mean daily home contacts of all German participants (POLYMOD has no care home residents); [A] care homes from 11 members on |
| Workplace | `ContactparameterSampling`: Poisson(λ) contacts drawn uniformly from the members | Office / Department / Workplace / WorkplaceSite: λ = 1.960 / 0.980 / 0.490 / 0.490 (total 3.920) | [L] adapted: work contacts per working participant; [A] split |
| Education (daycare to university) | `ContactparameterSampling` | SchoolClass / SchoolYear / School / SchoolComplex: λ = 1.337 / 0.668 / 0.334 / 0.334 (total 2.673) | [L] adapted: education contacts per participant in education; [A] split |
| Municipality | `AgeBasedContactSampling` (GEMS): Poisson(c) contacts per day, drawn from the members by the age-age matrix below | c = 4.155; matrix [polymod_municipality_matrix_DE.txt](inputs/contacts/polymod_municipality_matrix_DE.txt), 5-year age groups (0–4, …, 65–69, 70+) | [L] adapted: all other contacts, including visitors at home |

[polymod_contact_rates.csv](inputs/contacts/polymod_contact_rates.csv) compares the contacts per layer for Germany and all POLYMOD countries, with and without the cap on home contacts.

**Municipality contact matrix** ([src/contacts/polymod_matrix.jl](src/contacts/polymod_matrix.jl); steps in [polymod_municipality_matrix_DE.csv](inputs/contacts/polymod_municipality_matrix_DE.csv)): from the municipality contacts of the German participants, without the visitors at home (not identifiable in the data) and the 1.5 % of contacts without an age:

1. mean daily contacts c<sub>ij</sub> of a participant in age group i with persons in age group j (weekdays 5/7, weekend 2/7),
2. made reciprocal with the GEMS DE population N: (c<sub>ij</sub> N<sub>i</sub> + c<sub>ji</sub> N<sub>j</sub>) / (2 N<sub>i</sub>),
3. divided by N<sub>j</sub> (contacts per capita, independent of the age structure),
4. rows normalized to 1.

GEMS weights each row with the age structure of the municipality, computed from its members, so the age of the contacts reflects both the age preferences of POLYMOD and the local population. The expected number of contacts per member equals the contact parameter (0.9997 of c in the German first wave runs).

### 4.3 Transmission

The `InfectiousnessProfile` transmission function ([src/infectiousness_profile.jl](src/infectiousness_profile.jl)) scales the per-contact transmission probability by an infectiousness profile over the days τ since the infector's infectiousness onset:

p(τ) = `transmission_scale` · w(τ), with w(τ) the profile mass in [τ, τ + 1). The profile mass before τ = 0 is dropped and the rest normalized to 1, so `transmission_scale` is the sum of the daily per-contact probabilities over the infection. Recovered individuals are immune.

| Parameter | Value | Source |
|---|---|---|
| Profile | shifted gamma relative to symptom onset: shape 97.1875, rate 3.71875, shift 25.625 d | [L] He et al. 2020, corrected by Ashcroft et al. 2020, Table 1 |
| Config parameters | `distribution = "Gamma"`, `parameters = [97.1875, 0.26891]` (shape, scale = 1 / rate), `offset = 21.625` | [L] adapted: offset = shift − pre-symptomatic period (25.625 − 4.0); profile mass before infectiousness onset (4 %) dropped |
| `transmission_scale` | 0.45784 | [C] fitted to R0 = 3.28 with [calibrate_first_wave.jl](calibrate_first_wave.jl) |
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

### 5.1 First wave 2020 (Sec. 3.1.1, Fig 11–12)

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
| Imported cases, week 9 | 35 | [L] adapted: log-linear fit of the cases with exposure abroad (sum of the listed countries) on 04., 06., 08.03. (64, 130, 223; Tab. 2), evaluated at 02.03.; [A] stands for all imports reported up to the end of week 9 (fit: `GEMSPaper.imports_week9()` in [rki_imports.jl](src/first_wave/rki_imports.jl), 35.3) |
| Which cases are imports | cases with exposure abroad among the cases with a known place of exposure (63 % of the cases up to 27.04.) | [A] the cases with an unknown place of exposure (37 %) are not imputed |
| Share of Austria and Italy | 0.90 (weeks 9–12), 0.84 (weeks 13–17) | [L] adapted: weeks 9–10: Italy 200 of 223 foreign exposures on 08.03. (Tab. 2; no Austria); weeks 11–12 and 13: share in the increase of the country mentions 08.→22.03. and 22.→26.03. (Tab. 3 and text); [A] report intervals taken as reporting weeks, mentions as cases; 0.84 carried over to weeks 14–17 (the report of 27.04. bounds the share for 22.03.→27.04. at ≤ 0.814) |
| Share of Italy in the Alpine imports | 1.0 (weeks 9–10), 0.21 (11–12), 0.07 (13–17) | [L] adapted: no Austria on 08.03.; Italy's share in the increase of the Alpine mentions 08.→22.03. (995 of 4,849) and 22.03.→27.04. (474 of 6,443; Tab. 2); [A] reference date 27.04. here, 26.03. for the share of Austria and Italy |
| Import curve | infections per day (10.02.–26.04.) back-calculated from the weekly counts (reporting week = ISO week) over the time to onset (1 + Poisson(4.8) d; asymptomatic: infection + 6 d) and the RKI onset → report delays of the onset week, separately for Austria, Italy, and other countries; EMS with smoothing weights 1/4, 1/2, 1/4, flat start, 200 iterations; from 01.02. for the R(t) back-calculation | time to onset: reference model ([Section 4.5](#45-transition-durations)); delays: RKI (Data); [L] Becker et al. 1991 (method); [A] 200 iterations without a convergence criterion; no imports before 10.02. (config start; 01.02. for R(t)) or after 26.04. (end of week 17); asymptomatic imports reported with onset at infection + 6 d |
| Infections per reported import | 1 / 0.22 ≈ 4.5 | [A] imports detected like all infections; returning travellers from risk areas were probably tested more often, so the imports may be overestimated |
| Seeding of the imports | per day and source, expected count × 4.5 rounded stochastically (floor + Bernoulli draw); counties drawn with replacement by the allocation weights (computed with the county populations of the model); a random susceptible resident of the county (any age) infected; imports into counties outside the population model dropped | [A] |
| Reference points | Ischgl (47.0128 N, 10.2919 E) for Austria, Bozen (46.498 N, 11.354 E) for Italy | [A] Tyrol accounts for 2,737 of the 2,995 mentions with an Austrian region in the RKI report of 26.03. (Tab. 3), with Ischgl as its main cluster; South Tyrol/Trentino is the Italian region named most often (reports of 08.–26.03.) |
| Allocation to counties: Austria | Tyrol overnight stays per German origin region, distributed within the region by population × relative income^1.0 × (distance to Ischgl / 100 km)^−0.4, × (1 + 1.0 × school holidays of the state on the infection date) | see the following rows |
| Distance decay | 0.4 | [L] adapted: slope of ln(Tyrol overnight stays per inhabitant) on ln(population-weighted mean distance to Ischgl) over the 7 origin regions (slope −0.39; `GEMSPaper.fit_distance_decay()` in [rki_imports.jl](src/first_wave/rki_imports.jl) returns the decay 0.39), applied within the regions; [A] not adjusted for income (slope −0.25 with income^1.0 removed, `GEMSPaper.fit_distance_decay(income_elasticity = 1)`) |
| Income elasticity | 1.0 | [A] trips proportional to income (also used for the other imports) |
| Holiday factor | 1.0 | [A] weight doubled during school holidays |
| Allocation: Italy | 70 % like Austria with Bozen as reference point (same Tyrol origin regions and holiday factor), 30 % by population × income | [A] ski share 0.7; share of South Tyrol/Trentino in the RKI reports: 58 % of the Italian exposures (08.03., Tab. 2), 85–86 % of the Italian regions named (24./26.03., Tab. 3), 21 % of all Italy mentions (22.03.); plausible range 0.2–0.86 |
| Allocation: other countries | population × relative income^1.0 | [A] |
| School holidays (mid-February–March 2020) | BY, BW 22.02.–01.03.; HH 29.02.–15.03.; RP 15.–23.02.; SL 15.–25.02.; SN, MV 08.–23.02.; BW: 50 % of the schools | [L] Sekretariat der Kultusministerkonferenz 2018 (all states except BW); [A] adjacent weekends added; [A] BW: carnival week as in BY (no state-wide winter holidays, carnival holidays set by each school), coverage 0.5 (no source, range 0–1); not included: ST, TH (10.–14.02.), movable days in other states (e.g. carnival in NRW), Easter holidays |
| Heinsberg cluster | RKI cases of Kreis Heinsberg with onset (or, if unknown, reporting date) up to 29.02. and reported by 07.03., infected 1:1 at onset − 6 d (rounded mean incubation 1 + 0.8 + 4.0), or report − 11 d (+ median reporting delay 5 d) if the onset is unknown; 180 cases seeded (170 with onset); 132 cases with onset up to 29.02. reported later are not seeded | [A] without a direct source: cut-off 29.02. (two weeks after the carnival session of 15.02., Streeck et al. 2020), reporting limit 07.03., 1:1 (not divided by the ascertainment); imports are allocated to Heinsberg in addition |
| Plausibility check | share of introduced infections among all infections per phase, in the seeding summary (`sc_first_wave_final_seeding_summary.csv`, part `:figures`) | [L] compared with the share of reported cases with exposure abroad: up to 45.1 % in reporting week 11 (Frank et al. 2021, p. 5) |

**Contact phases and testing**

| Item | Value | Source |
|---|---|---|
| Phase start dates | phase 2: 09.03. (large events cancelled), phase 3: 16.03. (schools closed), phase 4: 23.03. (contact restrictions) | [L] dates of the measures to which Dehning et al. 2020 relate their three change points |
| Transitions between phases | linear over 3, 3, 4 days | [L] Dehning et al. 2020 |
| Contact factors (phases 1–4) | 0.641 / 0.207 / 0.182 / 0.139, multiplied with the contact parameters of all settings except households (municipalities keep their age-based sampling) | [C] [inputs/first_wave_factors_DE_model_first_wave.csv](inputs/first_wave_factors_DE_model_first_wave.csv) |
| Care homes | households with at least 11 members (`large_size` of `HouseholdSizeSampling`), λ = 2.706 scaled like the other settings | [A] |
| School closure | all school settings closed from 16.03.2020 | [L] date of the school closures (Dehning et al. 2020) |
| Ascertainment | 0.22 of all infections detected (`ASCERTAINMENT_2020`); also used to scale the reported imports, not the Heinsberg cluster | [L] adapted: under-ascertainment factor of 4 to 5 in the first half of 2020 from seroprevalence studies (Robert Koch-Institut 2021, p. 1), i.e. 20–25 % of infections detected; 1 / 4.5 = 0.22; [A] the same factor for imports |
| Testing | perfect PCR test at symptom onset with probability p = (0.22 − hospitalized share) / (symptomatic share − hospitalized share), and at hospital admission (if not detected before); positives are not isolated | [A] |
| Reported cases by reporting date | detected infections, reported at symptom onset (asymptomatic: infection + 6 d) plus a delay drawn from the RKI delays of the onset week | [A] |
| Reported cases by onset | detected infections at symptom onset (asymptomatic: infection + 6 d), each kept with the share of RKI cases with a known onset (per onset week); compared with the RKI cases with a known onset | [A] |

**Calibration** ([calibrate_first_wave.jl](calibrate_first_wave.jl) → `calibrate_first_wave()`, not part of `runall()`)

| Step | Setup |
|---|---|
| 1. Transmission scale | [Section 4.3](#43-transmission) |
| 2a. Generation time of the model | Saarland, 70 days, scale 0.35, 20 extra infections per day in the Regionalverband Saarbrücken (ticks 7–27), ascertainment 1.0, random seed 2; cohort of ticks 10–22; mean 5.76 d ([inputs/generation_time_DE_model_first_wave.csv](inputs/generation_time_DE_model_first_wave.csv)) |
| 2b. Contact factors | fitted to the weekly RKI cases by onset, weeks 10–16 (up to 19.04., weight 0.5 for week 10), with the R per phase of GEMS kept inside the 95 % interval of the RKI R per phase (renewal equation, Cori et al. 2013, with the generation time of the model from 2a, Gostic et al. 2020; R(t) back-calculated from 01.02., infections 01.03.–12.04., block bootstrap with 1,000 samples, block length 7 d); surrogate fit, then damped Gauss–Newton steps, at most 7 evaluations with 2 runs each (random seeds 0, 1); runs until 17.05. |
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

By default, `scenario_first_wave()` runs only the parts of the paper (`FIRST_WAVE_PAPER_PARTS`, marked *paper*) and writes the outputs of the paper (cases by reporting date, the paper figures `Fig_11_first_wave.png` and `Fig_12_first_wave_states.png`, the seeding summary, the county metrics, and the mixing comparison) and the files they are built from. The other parts can be selected with `parts`, e.g. `scenario_first_wave(parts = (:projections, :figures), diagnostics = true)`; `diagnostics = true` also writes the comparisons by symptom onset and the figures of all parts.

| Part | Setup | Runs / random seeds | Paper |
|---|---|---|---|
| `:final` | calibrated factors until 31.05. (six weeks of forecast after 19.04.) | 3 / 0–2 | paper |
| `:projections` | from 20.04. on, the contacts of phase 2, 3, or 4, or of phase 1 (until 17.05.); schools reopen in the projections to phase 1 and 2 | 2 per variant / 0, 1 | – |
| `:counterfactuals` | contacts of phase k kept from the start of phase k + 1 on (k = 1–3) | 1 per variant / 0 | – |
| `:band` | factor samples of the calibration, 5–95 % band, national and per state | 20 / 1001–1020 | paper |
| `:seeds_robustness` | all RKI cases infected before 09.03. as introduced infections (divided by the ascertainment, Heinsberg 1:1), no transmission before 09.03. | 3 / 0–2 | – |
| `:random_seeds` | the seeds of `:final` (imports and Heinsberg cluster), the same number per day, in counties drawn in proportion to their population (drawn anew for each random seed) | 3 / 0–2 | paper |
| `:mixing` | municipality contacts with uniform mixing and with same-age mixing instead of the POLYMOD matrix, each with the introductions of `:final` and of `:random_seeds`; compared with these two parts (run first if their outputs are not in the folder) | 3 per variant / 0–2 | paper |
| `:figures` | paper figures, seeding summary, county metrics; with `diagnostics`, the figures of all parts whose results exist | – | paper |

**Municipality mixing comparison** (`:mixing`, [src/first_wave/rki_mixing.jl](src/first_wave/rki_mixing.jl))

| Item | Setup |
|---|---|
| Mixing assumptions (`MUNICIPALITY_MIXING`) | reference: POLYMOD matrix ([Section 4.2](#42-contacts)); `:uniform`: `ContactparameterSampling`, contacts drawn uniformly from the members; `:same_age`: `AgeBasedContactSampling` with the identity matrix over the 15 age groups of the POLYMOD matrix (contacts only within the own age group) |
| Parameters | the contact parameter c = 4.155, the transmission scale, and the contact factors of the reference fit (no recalibration) |
| National fit | weekly cases GEMS / RKI in the fitted weeks (lowest, highest, RMS of the log ratio), reported cases of the whole run |
| Spatial pattern | county and state Spearman correlations with the RKI incidence as in the Evaluation above; difference alternative − reference with 95 % interval and p-value from the same bootstrap over states (paired); mean county Spearman of single runs between the random seeds of a model and with the reference model |
| Age and settings | age groups of cases and deaths, mean age of the infected, share of transmissions per setting type, share of municipality transmissions within the same 10-year age band, expected municipality contacts per member relative to c |
| Outputs | `sc_first_wave_mixing_summary.csv` (one column per model), `sc_first_wave_mixing_ages.csv`, and the outputs of each variant (`sc_first_wave_mixing_<variant>[_random_seeds]_*.csv`) |

### 5.2 Scenario modeling (Sec. 3.2.1, Fig 15–16)

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

The baseline run also produces the GEMS readout figure (Fig 16, `Fig_16_baseline_readouts.png`): daily infections, effective reproduction number, household attack rate (household sizes 2–6), infection duration, cumulative disease progressions, progression categories, and maps of the initial reproduction number and the infection fatality per county. The ResultData objects of both runs (style `HouseholdIsolationRD`) are stored in `household_isolation/` in the result folder. They hold the aggregated data that the figures need; of the raw infections, these are the infections per progression category and age group and the deaths per infection per county.

Outputs: `Fig_15_household_isolation.png` (the effect of household isolation above the code listing), `Fig_15A_household_isolation.png` (without the listing), and the readout figure.

The Fig 15B listing ([sc_household_isolation_screenshot_file.jl](src/sc_household_isolation_screenshot_file.jl)) uses the same intervention with the GEMS default parameterization ([data/DefaultConf.toml @ 7f3b79e](https://github.com/IMMIDD/GEMS/blob/7f3b79e51752f1d316a836c2e9d78a76d73572fe/data/DefaultConf.toml)) on 100,000 agents, without a fixed random seed.

### 5.3 Parameter sweep (Sec. 3.2.2, Fig 17)

Code: [sc_sweeping.jl](src/sc_sweeping.jl). Config: reference model.

To focus on the regime where self-isolation can be decisive, we assume that other measures already prevent 60 % of all transmissions, such that without isolation R_init lies above 1.

| Parameter | Reference model | This example |
|---|---|---|
| Population | – | Saarland |
| Initial infections | none | 0.1 % of the population, random |
| Duration | 112 days | 365 days |
| Transmission scale | 0.45784 | × (1 − 0.6): 60 % of transmissions prevented by other measures (`transmission_reduction = 0.6`) [I] |
| Asymptomatic share (x-axis) | 0.17 | 0 to 1, step 0.025, set in all age groups; the other categories rescaled proportionally |
| Adherence (y-axis) | – | 0 to 1, step 0.025 |
| Repetitions | – | 5 per combination; random seeds 1–8405 (one per run) |
| Outcome | – | R_init per run (`r0(PostProcessor(sim))`: infections caused by the first 5 % of all infections, divided by their number) |

| Trigger | Strategy (condition) | Measures |
|---|---|---|
| `SymptomTrigger` | Self Isolation (adheres ~ Bernoulli(p = adherence), drawn from the simulation's RNG) | `SelfIsolation(14)` (0) |

The heatmap shows the mean R_init per combination. Alongside the data (`sc_sweep_data.csv`) and the figure (`Fig_17_parameter_sweep.png`, text width, with the [drawn](#drawn-figures) strategy; the heatmap alone is `Fig_17B_sweep_heatmap.png`), the scenario writes `sc_sweep_variability.txt`, which summarizes the variability between replicates: the per-cell sample SD across replicates, aggregated over the grid as the pooled SD (root of the mean variance), mean, 95th percentile, and maximum; the per-cell range (maximum − minimum) across replicates, aggregated as its minimum, maximum, mean, and median; the standard error of the cell means; and the number of cells whose replicates lie on both sides of R = 1.

### 5.4 Multi-criteria strategy evaluation (Sec. 3.2.3, Fig 18, Tab 2)

Code: [sc_testing.jl](src/sc_testing.jl). Config: reference model.

To focus on the regime where tracing can be decisive, we assume that other measures (for example, limits on gatherings or closed restaurants and leisure venues) already prevent half of the community contacts outside households, workplaces, and schools.

| Parameter | Reference model | This example |
|---|---|---|
| Population | Germany | North Rhine-Westphalia |
| Initial infections | none | 50 infections in Gangelt (municipality key 05370008, county Heinsberg) |
| Duration | 112 days | 365 days |
| Municipality contact parameter | 4.155 | 2.078: 50 % of community contacts prevented (`municipality_contact_reduction = 0.5`), age-based sampling kept [I] |
| Runs / random seeds | – | 11 per scenario: 0, 3, …, 30 (Isolation Only), 1, 4, …, 31 (Isolation & Testing), 2, 5, …, 32 (Isolation, Testing & Tracing) |

**Test type** "Antigen Test" [I]: positive ~ Bernoulli(p = 0.80) if infected (sensitivity), Bernoulli(p = 0.01) if not (specificity 0.99).

| Scenario | Trigger | Strategy (condition) | Measures (offset) → follow-up |
|---|---|---|---|
| Isolation Only | `SymptomTrigger` | Self Isolation | `SelfIsolation(14)` (0) |
| Isolation & Testing | `SymptomTrigger` | Self Isolation | `SelfIsolation(14)` (0); `Test(Antigen)` (0) |
| Isolation, Testing & Tracing | `SymptomTrigger` | Isolation & Testing (not yet detected) | `SelfIsolation(14)` (0); `Test(Antigen)` (0) → positive: *Contact Identification*; negative: *Cancel Isolation* |
| | | Contact Identification | `FindSettingMembers` in Household, Office (if working), SchoolClass (if student) (1 d) → *Isolate & Test Contacts* |
| | | Isolate & Test Contacts (not yet detected) | `SelfIsolation(14)` (0); `Test(Antigen)` (0) → negative: *Cancel Isolation* |
| | | Cancel Isolation | `CancelSelfIsolation` (0) |

The tests in Isolation & Testing have no follow-up: the scenario detects cases without changing the spread, so it differs from Isolation Only only by the random seed. The tracing scenario corresponds to Fig 10 of the paper.

The top panels of the figure (`Fig_18_testing_scenarios.png`) show all runs, the bottom panels the run with the median peak day of infections of each scenario. The scenario stores the results of every run (`runs/seed_*/resultdata.jld2`, GEMS `ResultData` in the style `TestingRD`) and writes the outcomes of every run (`sc_testing_outcomes.csv`), their statistics per scenario (`sc_testing_summary.csv`: mean, standard deviation, minimum, maximum, and 95 % confidence interval of the mean), the comparison of tracing with the other scenarios (`sc_testing_comparison.csv`: averted infections, deaths, and workdays and schooldays of absence, additional tests, and additional tests per averted infection, death, workday of absence, and schoolday of absence; same statistics), and the daily data of every run (`sc_testing_daily.csv`). If it is interrupted, calling it again with the same `folder` and parameters skips the stored runs.

The outcomes are defined in `outcome_df()`. The duration of the wave is the number of days between 5 % and 95 % of all infections. Absent workers and students are isolated, severely ill, or in hospital. Their days are split into isolated and infected, isolated and not infected, and severely ill without isolation. The lost workdays and schooldays are the isolated days. Every day of the week counts, and students are the members of all educational settings. The isolations after false positive tests (individuals and days) are computed from the individual tests when a run is stored, and left empty for runs stored before this was added. The number of workers and students is stored with each run, too. If no run in the folder has it, `replot_testing(folder, population = "NRW")` counts them in the population.

### 5.5 Performance benchmark (Sec. 3.3, Fig 19)

Code: [benchmarking.jl](src/benchmarking.jl). Config: reference model. `benchmark(configfile = nothing)` runs the same protocol with the GEMS default parameterization instead, for comparison; `benchmark_machine.csv` records the config.

| Parameter | Reference model | This example |
|---|---|---|
| Initial infections | none | 0.1 % of the population, random (`infected_fraction = 0.001`) |
| Duration | 112 days | 365 days |
| Random seed | – | 1–11 (one per run) |
| On-the-fly populations | – | generated with the `[Population]` parameters of the reference model (household size 3, office size 5, school size 100) |

Protocol:

- **Runtime by population size:** 16 German states and Germany, plus on-the-fly populations of 0.1, 1, 5, 10, 20, and 50 million agents. The population files are downloaded and both population paths (population model and on-the-fly) are run once for compilation before anything is timed. Then 11 runs per population (random seeds 1–11), all with the threads of the running Julia process. Each run times three stages separately: initialization (reading population and settings files up to the ready simulation), simulation (`run!`), and post-processing (`ResultData` with all analyses). Post-processing with the `LightRD` style (fewer analyses) is timed separately and only reported in the summary. Each run also records the number of infections, the mean infection duration (exposure to removal), compile time, garbage collection time, allocated memory, and the live heap after a full garbage collection, once after the simulation and once after post-processing. Fig 19 shows medians with interquartile ranges.
- **Speedup by number of threads:** RP, 5 runs each with 1, 2, 4, 8, 16, 32, 64, and 128 threads, in one Julia process per thread count.
- **Peak memory:** maximum resident set size of one run of Germany (simulation and post-processing) in a fresh Julia process, i.e., the RAM a machine needs.

Outputs: `benchmark_machine.csv` (hardware, Julia version, GEMS version with repository, branch and tree from the Manifest, commit of this repository, threads, config), `benchmark_runs.csv` and `benchmark_summary.csv` (all runs; median and interquartile range per population), `benchmark_speedup_runs.csv` and `benchmark_speedup_summary.csv`, `benchmark_peak_memory.csv`, `Fig_19_benchmarks.png` (A–C: runtime of each stage by population size, D: speedup by number of threads; text width), and `benchmarks_per_agent.png` (runtime per million agents). `replot_benchmarks(folder)` redraws the figures from the CSV files; `replot_benchmarks(folder; on_the_fly = false)` leaves out the on-the-fly populations.

## 6. Population analyses and validation

No disease parameters are involved. All analyses use the GEMS population model of Germany (84,563,726 individuals, Gesyland).

### 6.1 Population and contact figures (Fig 4–8)

- **Fig 4:** settings map (`demographics()`, `Fig_04_settings_map.png` and `settings_map.csv.gz`). All households, school classes, and offices are counted per pixel of the figure (1 pixel is about 420 m), and the pixels are colored: households from dark gray to black (darker for more households), schools in blue and workplaces in green in pixels with fewer than 150 households. `coloring = household_opacity` draws the households black and partly opaque above schools and workplaces instead.
- **Fig 5:** population validation, see [Section 6.2](#62-population-validation).
- **Fig 6:** demographic maps per county (`demographics()`, `Fig_06_demographic_maps.png` and `demographic_maps.csv`). The values follow the definitions of the GEMS maps: individuals per km² (cadastral area, municipalities without an area in the GEMS shapefile left out), mean size of households with up to 10 members, share of single-person households, mean age, and shares of persons aged 0–14 and 65 or older. Eisenach (16056) is merged into the Wartburgkreis (16063), as in the GEMS county shapes.
- **Fig 7:** commuting. `mobility()` writes the shares of in- and out-commuters per county as CSV and draws the commuting figure ([Section 6.3](#63-commuting-validation)). The earlier maps of the commuter shares are drawn with `mobility(share_maps = true)`.
- **Fig 8:** contact patterns, as mean daily contacts per person by age (`contacts()`, `Fig_08_contact_matrices.png` and `contact_matrices.csv`). Every member of every setting draws contacts once with the GEMS sampling functions, and the sums are divided by all persons of the ego's age group. Panels A–C: households, workplaces, and educational settings (uniform mixing), summed over all levels of the setting hierarchy, with the contact parameters of the reference model ([Section 4.2](#42-contacts)). Panels D–F: municipalities with uniform mixing, with the POLYMOD matrix of the reference model, and with strict mixing (own age group only), i.e. the variants of the mixing comparison ([Section 5.1](#51-first-wave-2020-sec-311-fig-1112)). The age groups have five years (top group 85+). The contact matrices of E and F end with the group 70+.

`demographics()`, `settinganalyses()`, `mobility()`, and `contacts()` use the population of Germany by default, as in the paper. For quick tests, pass a smaller population, e.g. `sim = Simulation(population = "SL")`. The validations then are skipped, as they compare all 16 federal states. `contacts()` also needs the reference model, e.g. `sim = Simulation(configfile = "models/DE_model_first_wave.toml", population = "SL")`.

### 6.2 Population validation

`population_validation()` ([src/population_validation.jl](src/population_validation.jl), setting sizes in [src/setting_validation.jl](src/setting_validation.jl)) compares the population model with official statistics and writes `Fig_05_population_validation.png` and its data (`validation_*.csv`). The reference data are listed with table codes, reference dates, and licenses in [inputs/README.md](inputs/README.md#population-reference-data-for-the-validation-of-the-population-model).

| Panel | Comparison (Germany) | Reference |
|---|---|---|
| A) Age and sex | share of the population by single year of age (top group 90+) and sex | Destatis 12411-0013, 31.12.2019 |
| B) Household size | share of private households by size (1 to 5+); GEMS households with up to 9 members (Gesyland models care homes as households of 90+ residents) | Microcensus 12211-9034, 2019 |
| C) Family type | share of private households by type: couples without children, couples with children, single mothers, single fathers, other (single-person households and households without a core family) | Zensus 2022 5000H-1005, 15.05.2022 |
| D) Senior status | share of households by members aged 65+: only seniors (1, 2, 3+), seniors and younger persons, no seniors | Zensus 2022 5000H-1003, 15.05.2022 |
| E) Workplaces | share of workers by workplace size (1–9, 10–49, 50–249, 250–499, 500+) | BA establishments by size class, 30.06.2019 |
| F) Mean group size | members per daycare group and per school class (primary, secondary grades 5–10, special schools); secondary 11+ without reference (taught in courses) | Destatis 22541-0001, 01.03.2019; KMK class statistics 2019 |
| G) Universities | share of students by university size | Destatis 21311-0002, WS 2019/20 |
| H)–M) Federal states | GEMS vs. official share per state: population, age groups, employment and education (employed, daycare, pupils, not working 15–64, retired 65+), household sizes, family types, senior status | as above; employment: Zensus 2022 2000S-1002; pupils: Destatis 21111-0011, 21121-0010 (2019/20); daycare: Destatis 22541-0112 (01.03.2019) |

GEMS has no relationships between household members, so family types are classified by age: members at least 18 years younger than the oldest member count as children ([A]). Persons from 65 on without a workplace or education setting count as retired ([A]). Household types, senior status, and employment come from the Zensus 2022, which reports them for all federal states in categories that can be derived from the GEMS population.

### 6.3 Commuting validation

`commuter_validation()` ([src/commuter_validation.jl](src/commuter_validation.jl)) compares the commutes between counties in GEMS with the commuter statistics of the Bundesagentur für Arbeit (BA) of 30.06.2019 (employees subject to social insurance) and writes `Fig_07_commuting.png` and its data (`commuter_validation_*.csv`). Reference data: [inputs/README.md](inputs/README.md#commuters-reference-data-for-the-validation-of-commuting).

| Panel | Content |
|---|---|
| A) Commuting balance (GEMS) | net commuters per county (in- minus out-commuters) as peaks on the county centres (up: net in-commuting, down: net out-commuting), over the population density |
| B) Flows between counties | share of the out-commuters of the home county going to each other county, GEMS vs. BA, for the pairs with at least 10 commuters in both (smaller flows are only published in aggregates); agreement as Common Part of Commuters (CPC; Lenormand et al. 2012, 2016), with both flow sets scaled to the same total |
| C) Commuting distance | share of the workers who work in their home county and who commute over 0–25, 25–50, 50–100, 100–200, 200+ km (distance between the county centres), GEMS vs. BA |

Counties as in 2019 (the boundaries of the population and of the BA data), with Eisenach (16056) merged into the Wartburgkreis (16063) as in the GEMS county shapes. GEMS counts all workers, the BA only employees subject to social insurance (without civil servants and self-employed persons), so the panels compare shares.

## 7. References

- Amt der Tiroler Landesregierung, Abteilung Raumordnung und Statistik (2019). Der Tourismus im Winter 2018/2019. Innsbruck, August 2019. Tab. 6 (p. 11). https://www.tirol.gv.at/fileadmin/themen/statistik-budget/statistik/downloads/Tourismus/Der_Tourismus_im_Winter_18_19.pdf
- Ashcroft, P., et al. (2020). COVID-19 infectivity profile correction. *Swiss Med Wkly* 150:w20336. https://doi.org/10.4414/smw.2020.20336
- Becker, N. G., Watson, L. F., Carlin, J. B. (1991). A method of non-parametric back-projection and its application to AIDS data. *Stat Med* 10(10):1527–1542. https://doi.org/10.1002/sim.4780101005
- Byambasuren, O., et al. (2020). Estimating the extent of asymptomatic COVID-19 and its potential for community transmission: Systematic review and meta-analysis. *JAMMI* 5(4):223–234. https://doi.org/10.3138/jammi-2020-0030
- Byrne, A. W., et al. (2020). Inferred duration of infectious period of SARS-CoV-2: rapid scoping review and analysis of available evidence for asymptomatic and symptomatic COVID-19 cases. *BMJ Open* 10(8):e039856. https://doi.org/10.1136/bmjopen-2020-039856
- Cori, A., Ferguson, N. M., Fraser, C., Cauchemez, S. (2013). A new framework and software to estimate time-varying reproduction numbers during epidemics. *Am J Epidemiol* 178(9):1505–1512. https://doi.org/10.1093/aje/kwt133
- Dehning, J., et al. (2020). Inferring change points in the spread of COVID-19 reveals the effectiveness of interventions. *Science* 369(6500):eabb9789. https://doi.org/10.1126/science.abb9789
- Frank, C., Hellenbrand, W., Sievers, C. (2021). Betrachtung der reiseassoziierten COVID-19-Fälle im Sommer 2020 unter Berücksichtigung der Schulferien, Reisetätigkeit und Testkapazitäten. *Epidemiologisches Bulletin* 8/2021:3–16. https://doi.org/10.25646/7955
- Gostic, K. M., et al. (2020). Practical considerations for measuring the effective reproductive number, Rt. *PLoS Comput Biol* 16(12):e1008409. https://doi.org/10.1371/journal.pcbi.1008409
- He, X., et al. (2020). Temporal dynamics in viral shedding and transmissibility of COVID-19. *Nat Med* 26:672–675. https://doi.org/10.1038/s41591-020-0869-5; author correction *Nat Med* 26:1491–1493. https://doi.org/10.1038/s41591-020-1016-z
- Hu, Z., et al. (2020). Clinical characteristics of 24 asymptomatic infections with COVID-19 screened among close contacts in Nanjing, China. *Sci China Life Sci* 63(5):706–711. https://doi.org/10.1007/s11427-020-1661-4
- Karagiannidis, C., et al. (2020). Case characteristics, resource use, and outcomes of 10 021 patients with COVID-19 admitted to 920 German hospitals: an observational study. *Lancet Respir Med* 8(9):853–862. https://doi.org/10.1016/S2213-2600(20)30316-7
- Lenormand, M., Huet, S., Gargiulo, F., Deffuant, G. (2012). A universal model of commuting networks. *PLoS ONE* 7(10):e45985. https://doi.org/10.1371/journal.pone.0045985
- Lenormand, M., Bassolas, A., Ramasco, J. J. (2016). Systematic comparison of trip distribution laws and models. *J Transp Geogr* 51:158–169. https://doi.org/10.1016/j.jtrangeo.2015.12.008
- Liu, Y., et al. (2020). The reproductive number of COVID-19 is higher compared to SARS coronavirus. *J Travel Med* 27(2):taaa021. https://doi.org/10.1093/jtm/taaa021
- McAloon, C., et al. (2020). Incubation period of COVID-19: a rapid systematic review and meta-analysis of observational research. *BMJ Open* 10(8):e039652. https://doi.org/10.1136/bmjopen-2020-039652
- Mossong, J., et al. (2008). Social contacts and mixing patterns relevant to the spread of infectious diseases. *PLoS Med* 5(3):e74. https://doi.org/10.1371/journal.pmed.0050074
- Mossong, J., Hens, N., Jit, M., et al. (2020). POLYMOD social contact data (Version 2) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.3874557. License: CC BY-NC 4.0
- Ponge, J., et al. (2024). A standardized framework for modeling non-pharmaceutical interventions in individual-based infectious disease simulations. *WSC 2024*, 1106–1117. https://doi.org/10.1109/WSC63780.2024.10838778
- Robert Koch-Institut (2020). Täglicher Lagebericht des RKI zur Coronavirus-Krankheit-2019 (COVID-19), reports of 04.03., 06.03., 08.03. (Tab. 2), 22.03., 24.03., 26.03. (Tab. 3), and 27.04.2020 (Fig. 3, Tab. 2). https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Situationsberichte/2020-MM-DD-de.pdf?__blob=publicationFile (no longer online; archived copies: https://web.archive.org/web/2020/ followed by that URL)
- Robert Koch-Institut (2021). SARS-CoV-2-Seroprävalenz in der Allgemeinbevölkerung in Deutschland [Factsheet], as of 23.11.2021, p. 1. https://doi.org/10.25646/9288
- Robert Koch-Institut (2026). SARS-CoV-2 Infektionen in Deutschland [Data set], version of 14.10.2026. Zenodo. https://doi.org/10.5281/zenodo.23357370. License: CC BY 4.0
- Schilling, J., Lehfeld, A.-S., Schumacher, D., Ullrich, A., Diercke, M., Buda, S., Haas, W., RKI COVID-19 Study Group (2020). Krankheitsschwere der ersten COVID-19-Welle in Deutschland basierend auf den Meldungen gemäß Infektionsschutzgesetz. *Journal of Health Monitoring* 5(S11):2–20. https://doi.org/10.25646/7169
- Sekretariat der Kultusministerkonferenz (2018). Ferien im Schuljahr 2019/2020 in den Ländern in der Bundesrepublik Deutschland, as of 11.07.2018. https://www.kmk.org/fileadmin/Dateien/pdf/Ferienkalender/FER19_20.pdf
- Sekretariat der Kultusministerkonferenz (n.d.). Schüler/-innen, Klassen, Lehrkräfte und Absolvierende der Schulen 2010 bis 2019. Statistische Veröffentlichungen der KMK, Dokumentation 227, Teil A (allgemeinbildende und berufliche Schulen), retrieved 21.11.2026. https://www.kmk.org/ (downloaded at runtime, see [inputs/README.md](inputs/README.md))
- Statistik der Bundesagentur für Arbeit (2020). Tables of the employment statistics, reference date 30.06.2019: Betriebe und sozialversicherungspflichtig Beschäftigte; Gemeindedaten der sozialversicherungspflichtig Beschäftigten nach Wohn- und Arbeitsort; Pendlerverflechtungen der sozialversicherungspflichtig Beschäftigten nach Kreisen. Nürnberg. https://statistik.arbeitsagentur.de/ (files, URLs, and terms of use in [inputs/README.md](inputs/README.md))
- Statistik Austria (n.d.). Zuordnungsverzeichnis Deutschland (German postal codes → origin regions of the Austrian tourism statistics), edited by Landesstatistik Steiermark. https://www.landesentwicklung.steiermark.at/cms/dokumente/12656884_141979459/b19d32c7/Zuordnungsverzeichnis%20Oesterreich%20Deutschland-PLZ-PDF.pdf
- Statistische Ämter des Bundes und der Länder (2026). Einkommen der privaten Haushalte in den kreisfreien Städten und Landkreisen der Bundesrepublik Deutschland 1995 bis 2023. Volkswirtschaftliche Gesamtrechnungen der Länder, Reihe 2, Kreisergebnisse Band 3, Tab. 2.4 (Berechnungsstand Februar 2025, published January 2026). https://www.statistikportal.de/de/vgrdl/publikationenbsalt. License: Datenlizenz Deutschland – Namensnennung – Version 2.0.
- Statistisches Bundesamt (Destatis) (2026). GENESIS-Online, tables 12411-0013, 12211-9034, 21111-0011, 21121-0010, 21311-0002, 22541-0001, 22541-0112, retrieved 21.11.2026. https://www-genesis.destatis.de/. License: Datenlizenz Deutschland – Namensnennung – Version 2.0 (details per table in [inputs/README.md](inputs/README.md))
- Statistische Ämter des Bundes und der Länder (2026). Zensus 2022 database, tables 5000H-1003, 5000H-1005, 2000S-1002, retrieved 21.11.2026. https://ergebnisse.zensus2022.de/. License: Datenlizenz Deutschland – Namensnennung – Version 2.0
- Streeck, H., et al. (2020). Infection fatality rate of SARS-CoV2 in a super-spreading event in Germany. *Nat Commun* 11:5829. https://doi.org/10.1038/s41467-020-19509-y
- Verity, R., et al. (2020). Estimates of the severity of coronavirus disease 2019: a model-based analysis. *Lancet Infect Dis* 20(6):669–677. https://doi.org/10.1016/S1473-3099(20)30243-7
