# Inputs: sources, copyright and licenses

This folder contains the inputs read by the analyses: results of our own calibrations and data from official statistics. The repository's GPL-3.0 license applies to the code only. Data from third parties keeps its original license, listed below for every file.

## Rules for adding data

- Keep downloaded files **unchanged** (raw files). Tables derived from them are built by code and marked as own calculations.
- Add an entry to the file list below with: source and table code, reference date, data status ("Stand"), retrieval date, URL, license, required attribution, and changes.
- Files we extracted or reformatted ourselves carry a header that names the source, the license and the changes (see `county_income_2019.csv`).

## Licenses of the data providers

### Statistisches Bundesamt (Destatis), GENESIS-Online; Statistische Ämter des Bundes und der Länder

- License: [Datenlizenz Deutschland – Namensnennung – Version 2.0](https://www.govdata.de/dl-de/by-2-0) (dl-de/by-2-0). This also covers data retrieved through the GENESIS API.
- Copying and redistribution, including excerpts, are permitted with a source note.
- Attribution: "Statistisches Bundesamt (Destatis), Genesis-Online, <retrieval date>; Datenlizenz by-2-0", for derived tables with the addition "eigene Berechnung/eigene Darstellung". The note must refer to the license text and to the dataset (table code or URL).
- Terms: [Copyright for the GENESIS-Online database](https://www.destatis.de/EN/Service/Legal-Notice/CopyrightGENESISOnlineDatabase.html)

### Statistik der Bundesagentur für Arbeit (BA)

- Storing, copying and redistributing, including excerpts, are permitted with the source note "© Statistik der Bundesagentur für Arbeit". When published online, link to the [homepage of the BA statistics](https://statistik.arbeitsagentur.de/).
- The contents must **not be altered**. Own calculations are permitted but must be marked as such.
- Terms: [Impressum of the BA statistics](https://statistik.arbeitsagentur.de/DE/Home/Service/Impressum/impressum_node.html)

## Files

### Own results (repository license)

| File | Content |
|---|---|
| `first_wave_factors_DE_model_first_wave.csv` | Calibrated contact factors of the first wave model ([calibrate_first_wave.jl](../calibrate_first_wave.jl)) |
| `first_wave_band_factors_DE_model_first_wave.csv` | Factor samples for the uncertainty band of the first wave model |
| `generation_time_DE_model_first_wave.csv` | Generation time of the first wave model |

### Third-party data

| File | Source | License | Changes |
|---|---|---|---|
| `county_income_2019.csv` | Statistische Ämter des Bundes und der Länder, Volkswirtschaftliche Gesamtrechnungen der Länder, Reihe 2 Band 3, Tab. 2.4 (Berechnungsstand Februar 2025), [statistikportal.de](https://www.statistikportal.de) | dl-de/by-2-0 | Extract of the 2019 values per county; columns renamed (`ags`, `name`, `income_2019`) |

### `population/`: reference data for the validation of the population model

Raw downloads from [GENESIS-Online](https://www-genesis.destatis.de/datenbank/online/) and the [Zensus 2022 database](https://ergebnisse.zensus2022.de/datenbank/online/), unchanged, and the employment statistics of the BA (see below). Format of the GENESIS and Zensus files: semicolon-separated CSV export with decimal commas, UTF-8 with BOM, title rows above and notes and copyright rows below the table, and a value flag column ("e") after every value.

#### GENESIS-Online (Statistisches Bundesamt)

| File | Table | Content | Reference date | Retrieved |
|---|---|---|---|---|
| `12411-0013_de.csv` | [12411-0013](https://www-genesis.destatis.de/datenbank/online/statistic/12411/table/12411-0013) Bevölkerung: Bundesländer, Stichtag, Geschlecht, Altersjahre (Fortschreibung des Bevölkerungsstandes) | Persons by sex (male, female, total) and single year of age (0 to 90 and older), 16 federal states | 31.12.2019 and 31.12.2025 | 21.11.2026 |
| `12211-9034_de.csv` | [12211-9034](https://www-genesis.destatis.de/datenbank/online/statistic/12211/table/12211-9034) Privathaushalte: Bundesländer, Jahre (bis 2019), Haushaltsgröße (Mikrozensus) | Private households in thousands by size (1 to 5 and more persons), 16 federal states; households at main and secondary residence | 1991–2004 and 2019 | 21.11.2026 |
| `21111-0011_de.csv` | [21111-0011](https://www-genesis.destatis.de/datenbank/online/statistic/21111/table/21111-0011) Schüler: Bundesländer, Schuljahr, Geschlecht, Schulart, Jahrgangsstufen (Statistik der allgemeinbildenden Schulen) | Pupils at general schools by sex (male, female, total), school type and grade, 16 federal states | School year 2019/20 | 21.11.2026 |
| `21121-0010_de.csv` | [21121-0010](https://www-genesis.destatis.de/datenbank/online/statistic/21121/table/21121-0010) Schüler in beruflichen Schulen: Bundesländer, Schuljahr, Schulart (Statistik der beruflichen Schulen) | Pupils at vocational schools by school type, 16 federal states | School year 2019/20 | 21.11.2026 |
| `22541-0001_de.csv` | [22541-0001](https://www-genesis.destatis.de/datenbank/online/statistic/22541/table/22541-0001) Tageseinrichtungen, Pädagogisches Personal, Kinder, Genehmigte Plätze, Gruppen: Deutschland, Stichtag (Statistik der Kindertageseinrichtungen) | Daycare facilities, pedagogical staff, children, approved places and groups, Germany | 01.03.2019 | 21.11.2026 |
| `21311-0002_de.csv` | [21311-0002](https://www-genesis.destatis.de/datenbank/online/statistic/21311/table/21311-0002) Studierende: Deutschland, Semester, Nationalität, Geschlecht, Hochschulen (Statistik der Studierenden) | Students by nationality (German, foreign, total), sex (male, female, total) and university, Germany | WS 2019/20 | 21.11.2026 |
| `22541-0112_de.csv` | [22541-0112](https://www-genesis.destatis.de/datenbank/online/statistic/22541/table/22541-0112) Kinder in Tageseinrichtungen: Bundesländer, Stichtag, Alter der Nichtschul- und Schulkinder (Statistik der Kindertageseinrichtungen) | Children in daycare facilities by age, separately for children not yet attending school (under 1 to 7 and older) and school children, 16 federal states | 01.03.2019 | 21.11.2026 |

- Source note (as given in the files): "© Statistisches Bundesamt (Destatis), 2026".
- License: dl-de/by-2-0 (see above). Figures and tables derived from these files are own calculations ("eigene Berechnung").
- Changes: none.

#### Zensus 2022 database (Statistische Ämter des Bundes und der Länder)

| File | Table | Content | Reference date | Retrieved |
|---|---|---|---|---|
| `5000H-1003_de.csv` | [5000H-1003](https://ergebnisse.zensus2022.de/datenbank/online/statistic/5000H/table/5000H-1003) Haushalte: Seniorenstatus eines privaten Haushalts (ausführlich) | Private households by senior status (only seniors: single, two persons, more persons; seniors and younger persons; no seniors), 16 federal states (no Germany total); counts and shares | 15.05.2022 | 21.11.2026 |
| `5000H-1005_de.csv` | [5000H-1005](https://ergebnisse.zensus2022.de/datenbank/online/statistic/5000H/table/5000H-1005) Haushalte: Typ des privaten Haushalts | Private households by type (single-person households, couples without children, couples with children, single parents (mothers, fathers), multi-person households without a core family), 16 federal states (no Germany total); counts and shares | 15.05.2022 | 21.11.2026 |
| `2000S-1002_de.csv` | [2000S-1002](https://ergebnisse.zensus2022.de/datenbank/online/statistic/2000S/table/2000S-1002) Personen: Erwerbsstatus (ausführlich) | Persons by employment status (labour force: employed, unemployed; not in the labour force: below the minimum age of 15, recipients of pensions/capital income, pupils and students not economically active, others), 16 federal states (no Germany total); counts and shares | 15.05.2022 | 21.11.2026 |

- Source note (as given in the files): "© Statistische Ämter des Bundes und der Länder, Deutschland, 2026. Lizenziert unter der Datenlizenz Deutschland – Namensnennung – Version 2.0".
- License: dl-de/by-2-0 (see above). Figures and tables derived from these files are own calculations ("eigene Berechnung").
- Changes: none.

#### Statistik der Bundesagentur für Arbeit (BA)

| File | Table | Content | Reference date | Retrieved |
|---|---|---|---|---|
| `bstbetr-heft-d-0-201912-xls.xls` | [Betriebe und sozialversicherungspflichtig Beschäftigte, Deutschland (Jahreszahlen)](https://statistik.arbeitsagentur.de/DE/Navigation/Statistiken/Fachstatistiken/Beschaeftigung/Betriebe/Betriebe-Nav.html) (Beschäftigungsstatistik, Tabellen; created 30.03.2020) | Establishments (Betriebe, local units) and employees subject to social insurance by establishment size class, economic sector and county, Germany; sheet 2.3: time series 2009–2019 by size class (1–5 to 500 and more employees) | 30.06.2019 (sheet 2.3: 30.06.2009–30.06.2019) | 21.11.2026 |
| `ba_establishments_by_size_2019.csv` | Extract of `bstbetr-heft-d-0-201912-xls.xls`, sheet 2.3 | Establishments and employees subject to social insurance by size class, Germany | 30.06.2019 | – |

- Source note: "© Statistik der Bundesagentur für Arbeit"; citation as given in the file: "Statistik der Bundesagentur für Arbeit, Tabellen, Betriebe und sozialversicherungspflichtige Beschäftigung, Nürnberg, 30. Juni 2019".
- License: BA terms (see above), as stated in the file: storing, copying and redistributing, also in excerpts, with a source note; the contents must not be altered; own calculations must be marked; online, link to the [homepage of the BA statistics](https://statistik.arbeitsagentur.de/).
- Changes: the raw file is unchanged. The extract contains the row 30.06.2019 of sheet 2.3 with unchanged values (renamed columns, class bounds added), because the code cannot read the old Excel format. Its header names the source, the license and the changes.

### Downloaded at runtime (not in the repository)

These files are not part of the repository, because their license does not permit redistribution (KMK) or because they are too large (RKI). The code downloads them into the gitignored `data/` folder when they are missing and warns if their checksum differs from the version used in the paper.

| File | Source | Reference date | Version | Downloaded by |
|---|---|---|---|---|
| `Aktuell_Deutschland_SarsCov2_Infektionen.csv` | Robert Koch-Institut (2026). SARS-CoV-2 Infektionen in Deutschland [Data set]. Zenodo. [doi:10.5281/zenodo.23357370](https://doi.org/10.5281/zenodo.23357370) (about 400 MB) | 2020 onwards | version of 14.10.2026 | `obtain_rki_cases` ([rki_data.jl](../src/first_wave/rki_data.jl)) |
| `SKL_Teil_A_Schueler_2019.xlsx` | Sekretariat der Kultusministerkonferenz (KMK), [Dokumentation 227](https://www.kmk.org/downloads-dokumente/statistik/schulstatistik/schueler/-innen-klassen-lehrkraefte-und-absolvierende/archiv-schueler/-innen-klassen-lehrkraefte-und-absolvierende.html): Schüler/-innen, Klassen, Lehrkräfte und Absolventen/-innen der Schulen 2010–2019, Teil A (allgemeinbildende und berufliche Schulen) | 2010–2019 | retrieved 21.11.2026 | `obtain_kmk_classes` ([setting_validation.jl](../src/setting_validation.jl)) |

- License of the RKI data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Attribution as in the citation above.
- License of the KMK file: the [terms of use of the KMK](https://www.kmk.org/impressum.html) ("Nutzungsbedingungen") permit copies and downloads for personal, private and non-commercial use only. Redistribution requires the written consent of the KMK. Values derived from the file (mean class sizes) are own calculations.
