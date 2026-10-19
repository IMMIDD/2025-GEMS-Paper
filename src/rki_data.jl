###
### RKI DATA (FIRST WAVE 2020)
###

export RKI_CASES_FILE
export load_rki_cases
export rki_county
export rki_cases_by_county
export county_population
export rki_incidence
export rki_onset_delays


"""
    RKI_CASES_FILE

Default location of the RKI case data
(`Aktuell_Deutschland_SarsCov2_Infektionen.csv` from
https://github.com/robert-koch-institut/SARS-CoV-2-Infektionen_in_Deutschland).
"""
const RKI_CASES_FILE = joinpath("data", "Aktuell_Deutschland_SarsCov2_Infektionen.csv")


"""
    rki_county(id::Integer)

Maps an RKI `IdLandkreis` to a GEMS county id (Berlin districts are merged into 11000).
"""
rki_county(id::Integer) = 11001 <= id <= 11012 ? 11000 : Int(id)


"""
    load_rki_cases(; file = RKI_CASES_FILE, max_onset_lead = 21)

Loads the RKI case data. Columns:

| Name      | Description                                        |
| :-------- | :------------------------------------------------- |
| `county`  | County id (Berlin districts merged)                |
| `ags`     | County `AGS`                                       |
| `report`  | Reporting date (`Meldedatum`)                      |
| `onset`   | Symptom onset, `missing` if unknown                |
| `refdate` | Onset if known, otherwise reporting date           |
| `cases`   | Number of cases (`AnzahlFall`)                     |

Onsets more than `max_onset_lead` days before the reporting date are treated as unknown.
"""
function load_rki_cases(; file = RKI_CASES_FILE, max_onset_lead = 21)
    isfile(file) || throw(ArgumentError("RKI case data not found at $file."))

    df = CSV.read(file, DataFrame,
        select = [:IdLandkreis, :Meldedatum, :Refdatum, :IstErkrankungsbeginn, :NeuerFall, :AnzahlFall],
        dateformat = "yyyy-mm-dd")

    # NeuerFall == -1: removed since the previous day's version
    df = df[df.NeuerFall .!= -1, :]

    county = rki_county.(df.IdLandkreis)
    onset = [
        e == 1 && Dates.value(m - r) <= max_onset_lead ? r : missing
        for (e, m, r) in zip(df.IstErkrankungsbeginn, df.Meldedatum, df.Refdatum)]

    return DataFrame(
        county = county,
        ags = AGS.(county .* 1000),
        report = df.Meldedatum,
        onset = onset,
        refdate = coalesce.(onset, df.Meldedatum),
        cases = df.AnzahlFall)
end


"""
    rki_cases_by_county(df; date = :refdate, from, to)

RKI cases per county with `date` (`:refdate`, `:onset`, or `:report`) between
`from` and `to`. Counties without cases are not included.
"""
function rki_cases_by_county(df::DataFrame; date::Symbol = :refdate, from::Date, to::Date)
    return df[coalesce.(from .<= df[!, date] .<= to, false), :] |>
        d -> groupby(d, :ags) |>
        d -> combine(d, :cases => sum => :cases)
end


"""
    county_population(sim::Simulation)

Population per county, summed up from the municipalities of the population model.
"""
function county_population(sim::Simulation)
    return sim |> municipalities |>
        m -> DataFrame(ags = county.(ags.(m)), population = size.(m)) |>
        d -> groupby(d, :ags) |>
        d -> combine(d, :population => sum => :population)
end


"""
    rki_incidence(df, pop; date = :refdate, from, to)

RKI incidence per 100,000 per county between `from` and `to` (`pop` from
`county_population`). The result can be passed to `agsmap`.
"""
function rki_incidence(df::DataFrame, pop::DataFrame; date::Symbol = :refdate, from::Date, to::Date)
    cases = rki_cases_by_county(df, date = date, from = from, to = to)

    unmatched = setdiff(cases.ags, pop.ags)
    if !isempty(unmatched)
        @warn "RKI counties without a match in the population model (dropped): $unmatched"
    end

    return leftjoin(pop, cases, on = :ags) |>
        d -> transform(d, :cases => ByRow(c -> coalesce(c, 0)) => :cases) |>
        d -> transform(d, [:cases, :population] => ((c, p) -> c .* 100_000 ./ p) => :incidence) |>
        d -> sort(d, :ags, by = a -> a.id) |>
        d -> DataFrames.select(d, :ags, :incidence)
end


"""
    rki_onset_delays(df; from = Date(2020, 2, 1), to = Date(2020, 4, 30))

Distribution of the delay from symptom onset to reporting (in days) for cases
with known onset between `from` and `to`.
"""
function rki_onset_delays(df::DataFrame; from::Date = Date(2020, 2, 1), to::Date = Date(2020, 4, 30))
    return df[.!ismissing.(df.onset), :] |>
        d -> d[from .<= d.onset .<= to, :] |>
        d -> transform(d, [:report, :onset] => ByRow((r, o) -> Dates.value(r - o)) => :delay) |>
        d -> groupby(d, :delay) |>
        d -> combine(d, :cases => sum => :cases) |>
        d -> sort(d, :delay)
end
