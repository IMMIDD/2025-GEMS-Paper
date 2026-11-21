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

Raw downloads from [GENESIS-Online](https://www-genesis.destatis.de/datenbank/online/) and the [Zensus 2022 database](https://ergebnisse.zensus2022.de/datenbank/online/), unchanged. Format: semicolon-separated CSV export with decimal commas, UTF-8 with BOM, title rows above and notes and copyright rows below the table, and a value flag column ("e") after every value.

#### GENESIS-Online (Statistisches Bundesamt)

| File | Table | Content | Reference date | Retrieved |
|---|---|---|---|---|
| `12411-0013_de.csv` | [12411-0013](https://www-genesis.destatis.de/datenbank/online/statistic/12411/table/12411-0013) Bevölkerung: Bundesländer, Stichtag, Geschlecht, Altersjahre (Fortschreibung des Bevölkerungsstandes) | Persons by sex (male, female, total) and single year of age (0 to 90 and older), 16 federal states | 31.12.2019 and 31.12.2025 | 21.11.2026 |
| `12211-9034_de.csv` | [12211-9034](https://www-genesis.destatis.de/datenbank/online/statistic/12211/table/12211-9034) Privathaushalte: Bundesländer, Jahre (bis 2019), Haushaltsgröße (Mikrozensus) | Private households in thousands by size (1 to 5 and more persons), 16 federal states; households at main and secondary residence | 1991–2004 and 2019 | 21.11.2026 |

- Source note (as given in the files): "© Statistisches Bundesamt (Destatis), 2026".
- License: dl-de/by-2-0 (see above). Figures and tables derived from these files are own calculations ("eigene Berechnung").
- Changes: none.

#### Zensus 2022 database (Statistische Ämter des Bundes und der Länder)

| File | Table | Content | Reference date | Retrieved |
|---|---|---|---|---|
| `1000A-2027_de.csv` | [1000A-2027](https://ergebnisse.zensus2022.de/datenbank/online/statistic/1000A/table/1000A-2027) Personen: Geschlecht – Alter (Jahresschritte) | Persons by sex (total, male, female) and single year of age (0 to 100 and older), Germany and the 16 federal states; counts and shares | 15.05.2022 | 20.11.2026 |
| `5000H-1001_de.csv` | [5000H-1001](https://ergebnisse.zensus2022.de/datenbank/online/statistic/5000H/table/5000H-1001) Haushalte: Größe des privaten Haushalts (bis 6 und mehr Personen) | Private households by size (1 to 6 and more persons), 16 federal states (no Germany total); counts and shares | 15.05.2022 | 20.11.2026 |
| `5000H-1003_de.csv` | [5000H-1003](https://ergebnisse.zensus2022.de/datenbank/online/statistic/5000H/table/5000H-1003) Haushalte: Seniorenstatus eines privaten Haushalts (ausführlich) | Private households by senior status (only seniors: single, two persons, more persons; seniors and younger persons; no seniors), 16 federal states (no Germany total); counts and shares | 15.05.2022 | 21.11.2026 |
| `5000H-1005_de.csv` | [5000H-1005](https://ergebnisse.zensus2022.de/datenbank/online/statistic/5000H/table/5000H-1005) Haushalte: Typ des privaten Haushalts | Private households by type (single-person households, couples without children, couples with children, single parents (mothers, fathers), multi-person households without a core family), 16 federal states (no Germany total); counts and shares | 15.05.2022 | 21.11.2026 |
| `2000S-1002_de.csv` | [2000S-1002](https://ergebnisse.zensus2022.de/datenbank/online/statistic/2000S/table/2000S-1002) Personen: Erwerbsstatus (ausführlich) | Persons by employment status (labour force: employed, unemployed; not in the labour force: below the minimum age of 15, recipients of pensions/capital income, pupils and students not economically active, others), 16 federal states (no Germany total); counts and shares | 15.05.2022 | 21.11.2026 |
| `1000A-2075_de.csv` | [1000A-2075](https://ergebnisse.zensus2022.de/datenbank/online/statistic/1000A/table/1000A-2075) Personen: Größe des privaten Haushalts – Seniorenstatus eines privaten Haushalts | Persons by size of their private household (1 to 6 and more persons) and its senior status (only seniors, seniors and younger persons, no seniors), 16 federal states (no Germany total); counts | 15.05.2022 | 21.11.2026 |

- Source note (as given in the files): "© Statistische Ämter des Bundes und der Länder, Deutschland, 2026. Lizenziert unter der Datenlizenz Deutschland – Namensnennung – Version 2.0".
- License: dl-de/by-2-0 (see above). Figures and tables derived from these files are own calculations ("eigene Berechnung").
- Changes: none.
