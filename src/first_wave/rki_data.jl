###
### FIRST WAVE 2020: DATA AND COMMON HELPERS
###

export load_rki_cases


"""
    RKI_CASES_FILE

Location of the RKI case data (`Aktuell_Deutschland_SarsCov2_Infektionen.csv`).
Downloaded from `RKI_CASES_URL` if it does not exist.
"""
const RKI_CASES_FILE = joinpath("data", "Aktuell_Deutschland_SarsCov2_Infektionen.csv")

"""
    RKI_CASES_URL

Version of the RKI case data used in the paper (14.10.2026): Robert Koch-Institut
(2026). SARS-CoV-2 Infektionen in Deutschland [Data set]. Zenodo.
https://doi.org/10.5281/zenodo.23357370. License: CC BY 4.0.
"""
const RKI_CASES_URL = "https://zenodo.org/records/23357370/files/" *
    "Aktuell_Deutschland_SarsCov2_Infektionen.csv"

"""
    RKI_CASES_SHA256

SHA-256 checksum of the file at `RKI_CASES_URL`.
"""
const RKI_CASES_SHA256 = "8ad25a5f0ac6ef5eaaac600a5b3a01140bca46752edb1bee4342b94bd961a3e4"


"""
    obtain_rki_cases(file = RKI_CASES_FILE)

Downloads the RKI case data from `RKI_CASES_URL` (about 400 MB) if `file` does
not exist, and checks that it is the version used in the paper.
"""
function obtain_rki_cases(file = RKI_CASES_FILE)
    if !isfile(file)
        printinfo("Downloading the RKI case data (about 400 MB) from $RKI_CASES_URL")
        mkpath(dirname(file))
        Downloads.download(RKI_CASES_URL, file)
    end
    checksum = bytes2hex(open(SHA.sha256, file))
    if checksum != RKI_CASES_SHA256
        @warn "$file is not the version of the RKI case data used in the paper " *
            "($RKI_CASES_URL)."
    end
    return file
end

"""
    HEINSBERG

AGS of Kreis Heinsberg (carnival cluster in February 2020).
"""
const HEINSBERG = AGS(5370000)

"""
    RKI_STATES

GEMS population codes of the federal states, in the order of the state numbers in the AGS.
"""
const RKI_STATES = ["SH", "HH", "NI", "HB", "NRW", "HE", "RP", "BW", "BY", "SL",
    "BE", "BB", "MV", "SN", "ST", "TH"]

# name of the pathogen in the GEMS configs
const PATHOGEN = "Covid19"


# Monday and Sunday of ISO week w in 2020, and the ISO week of a date
monday_2020(w::Integer) = Date(2019, 12, 30) + Week(w - 1)
sunday_2020(w::Integer) = monday_2020(w) + Day(6)
week_2020(d::Date) = Dates.value(d - Date(2019, 12, 30)) ÷ 7 + 1

# date as "dd.mm."
ddmm(d::Date) = Dates.format(d, "dd.mm.")

# number of the federal state (1-16) of an AGS
state_number(a::AGS) = a.id ÷ 1_000_000

# table sorted by the column `ags`
sort_by_ags(d::DataFrame) = sort(d, :ags, by = a -> a.id)


"""
    rki_county(id::Integer)

Maps an RKI `IdLandkreis` to a GEMS county id (Berlin districts are merged into 11000).
"""
rki_county(id::Integer) = 11001 <= id <= 11012 ? 11000 : Int(id)


"""
    gems_county(a::AGS)

County of a GEMS AGS as used in the RKI data: Eisenach (16056), still a
separate county in the GEMS population, is merged into the Wartburgkreis (16063).
"""
function gems_county(a::AGS)
    c = county(a)
    return c.id == 16056000 ? AGS(16063000) : c
end


"""
    load_rki_cases(;
        file = RKI_CASES_FILE,
        max_onset_lead = 21,
        max_onset_lag = nothing)

Loads the RKI case data (downloaded if needed, see `obtain_rki_cases`). Columns:

| Name      | Description                                        |
| :-------- | :------------------------------------------------- |
| `county`  | County id (Berlin districts merged)                |
| `ags`     | County `AGS`                                       |
| `report`  | Reporting date (`Meldedatum`)                      |
| `onset`   | Symptom onset, `missing` if unknown                |
| `refdate` | Onset if known, otherwise reporting date           |
| `cases`   | Number of cases (`AnzahlFall`)                     |

Onsets more than `max_onset_lead` days before the reporting date are treated as
unknown (data errors), as are onsets more than `max_onset_lag` days after it
(if given): a case can only be reported before its onset if it tested positive
after infection, so the onset is at most one incubation period after the
report (e.g. `model_incubation(configfile).max_lag`).
"""
function load_rki_cases(;
    file = RKI_CASES_FILE,
    max_onset_lead = 21,
    max_onset_lag = nothing)

    # the default file is downloaded if needed and checked against the version of the paper
    if file == RKI_CASES_FILE
        obtain_rki_cases(file)
    elseif !isfile(file)
        throw(ArgumentError("RKI case data not found at $file."))
    end

    df = CSV.read(file, DataFrame,
        select = [:IdLandkreis, :Meldedatum, :Refdatum, :IstErkrankungsbeginn,
            :NeuerFall, :AnzahlFall],
        dateformat = "yyyy-mm-dd")

    # NeuerFall == -1: removed since the previous day's version
    df = df[df.NeuerFall .!= -1, :]

    # reference date `r` is the onset if it is known (e == 1) and within the bounds around the
    # report `m`
    function onset_of(e, m, r)
        known = e == 1 && Dates.value(m - r) <= max_onset_lead
        plausible = isnothing(max_onset_lag) || Dates.value(r - m) <= max_onset_lag
        return known && plausible ? r : missing
    end
    onset = onset_of.(df.IstErkrankungsbeginn, df.Meldedatum, df.Refdatum)
    county = rki_county.(df.IdLandkreis)

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
    onset_basis(rki)

RKI cases with a known symptom onset, dated by onset (column `report`).
"""
function onset_basis(rki::DataFrame)
    d = rki[.!ismissing.(rki.onset), :]
    d.report = Date.(d.onset)
    return d
end


"""
    county_population(sim::Simulation)

Population per county (as in the RKI data, see `gems_county`), summed up from
the municipalities of the population model.
"""
function county_population(sim::Simulation)
    return sim |> municipalities |>
        m -> DataFrame(ags = gems_county.(ags.(m)), population = size.(m)) |>
        d -> groupby(d, :ags) |>
        d -> combine(d, :population => sum => :population)
end


"""
    rki_onset_delays_by_week(df; from = Date(2020, 2, 1), to = Date(2020, 5, 31), min_cases = 200)

Distributions of the delay from symptom onset to reporting per onset week.
Weeks with fewer than `min_cases` cases are left out; use `delay_distribution`
to look up the distribution for a date.
"""
function rki_onset_delays_by_week(df::DataFrame;
    from::Date = Date(2020, 2, 1),
    to::Date = Date(2020, 5, 31),
    min_cases::Int = 200)

    d = df[.!ismissing.(df.onset), :] |>
        x -> x[from .<= x.onset .<= to, :] |>
        x -> transform(x,
            [:report, :onset] => ByRow((r, o) -> Dates.value(r - o)) => :delay,
            :onset => ByRow(firstdayofweek) => :week) |>
        x -> combine(groupby(x, [:week, :delay]), :cases => sum => :cases)
    dists = Dict(first(g.week) => sort(DataFrames.select(g, :delay, :cases), :delay)
        for g in groupby(d, :week) if sum(g.cases) >= min_cases)
    return (weeks = sort(collect(keys(dists))), dists = dists)
end


"""
    delay_distribution(delays, date::Date)

Delay distribution (columns `delay`, `cases`) of the onset week of `date`,
or of the nearest week available (output of `rki_onset_delays_by_week`).
"""
function delay_distribution(delays, date::Date)
    week = firstdayofweek(date)
    haskey(delays.dists, week) && return delays.dists[week]
    return delays.dists[delays.weeks[argmin(abs.(Dates.value.(delays.weeks .- week)))]]
end

# no reporting delay (comparison by symptom onset date, see `report_offsets` and `gems_reported`)
delay_distribution(::Nothing, ::Date) = DataFrame(delay = [0], cases = [1])


"""
    model_incubation(configfile; q = 0.99)

Time from infection to symptom onset in a GEMS config: 1 + Poisson(`lambda`)
days for symptomatic progressions (`lambda` = exposure to infectiousness +
infectiousness to symptom onset), its rounded mean (`incubation`), its `q`
quantile (`max_lag`), and the mean share of asymptomatic infections over the
age groups (`asymptomatic`).
"""
function model_incubation(configfile; q::Real = 0.99)
    p = TOML.parsefile(configfile)["Pathogens"][PATHOGEN]
    s = p["progressions"]["Symptomatic"]
    lambda = only(s["exposure_to_infectiousness_onset"]["parameters"]) +
        only(s["infectiousness_onset_to_symptom_onset"]["parameters"])
    assignment = p["progression_assignment"]["parameters"]
    k_asymptomatic = findfirst(==("Asymptomatic"), assignment["progression_categories"])
    asymptomatic = mean(row[k_asymptomatic] for row in assignment["stratification_matrix"])

    # q quantile of Poisson(lambda)
    k, pk, cdf = 0, exp(-lambda), exp(-lambda)
    while cdf < q
        k += 1
        pk *= lambda / k
        cdf += pk
    end
    return (
        lambda = lambda,
        incubation = round(Int, 1 + lambda),
        max_lag = 1 + k,
        asymptomatic = asymptomatic)
end


"""
    config_startdate(configfile)

Start date of the simulation in a GEMS config.
"""
config_startdate(configfile) = Date(TOML.parsefile(configfile)["Simulation"]["startdate"])


"""
    config_transmission_rate(configfile)

Transmission rate in a GEMS config.
"""
function config_transmission_rate(configfile)
    p = TOML.parsefile(configfile)["Pathogens"][PATHOGEN]
    return p["transmission_function"]["parameters"]["transmission_rate"]
end
