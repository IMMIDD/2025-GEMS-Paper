###
### FIRST WAVE 2020: OUTCOMES BY AGE AND SETTING, AND THE MIXING COMPARISON
###

# Outcomes that differ between age-based (reference model) and uniform contact
# sampling in municipalities (who meets whom, not how many): age distribution of
# cases and deaths, fatality, setting shares of the transmissions, and the share
# of municipality transmissions within the same age band. `mixing_summary`
# compares the runs of both assumptions (see the `:mixing` part of
# `scenario_first_wave`).

export use_uniform_municipalities!


"""
    use_uniform_municipalities!(sim; contactparameter = nothing)

Replaces the contact sampling of all municipalities of `sim` with uniform
sampling from the members (`ContactparameterSampling`), keeping the mean number
of contacts of the config (or `contactparameter`). Returns the new sampling
method.
"""
function use_uniform_municipalities!(sim; contactparameter = nothing)
    municipalities = settings(sim, Municipality)
    c = something(contactparameter, contact_sampling_method(first(municipalities)).contactparameter)
    csm = ContactparameterSampling(Float64(c))
    foreach(s -> contact_sampling_method!(s, csm), municipalities)
    return csm
end


"""
    RKI_AGE_GROUPS

Age groups of the RKI case data (`Altersgruppe`): lower limit and label.
"""
const RKI_AGE_GROUPS = [0 => "A00-A04", 5 => "A05-A14", 15 => "A15-A34", 35 => "A35-A59",
    60 => "A60-A79", 80 => "A80+"]

# index of the RKI age group of an age in years
rki_age_index(age) = searchsortedlast(first.(RKI_AGE_GROUPS), age)

"""
    SETTING_LAYERS

Setting layers of the transmissions by `setting_type` (see `GEMS.settingchar`).
"""
const SETTING_LAYERS = Dict(
    'h' => :household, 'm' => :municipality,
    'o' => :work, 'd' => :work, 'w' => :work, 'p' => :work,
    'c' => :school, 'y' => :school, 's' => :school, 'x' => :school)

"""
    SAME_AGE_BAND

Width of the age bands (years) for the share of municipality transmissions
between individuals of the same band.
"""
const SAME_AGE_BAND = 10

"""
    MATURE_DAYS

Days before the end of a run after which infections are left out of the
infection fatality rate (their outcome is not yet known), as for the
ascertainment.
"""
const MATURE_DAYS = 21


"""
    rki_age_table(s::RKISettings; file = RKI_CASES_FILE, from = s.startdate, to)

RKI cases and deaths per age group (`RKI_AGE_GROUPS`, unknown ages left out)
reported between `from` and `to` in the counties of `s`.
"""
function rki_age_table(s::RKISettings; file = RKI_CASES_FILE, from::Date = s.startdate, to::Date)
    obtain_rki_cases(file)
    df = CSV.read(file, DataFrame,
        select = [:IdLandkreis, :Altersgruppe, :Meldedatum, :NeuerFall, :NeuerTodesfall,
            :AnzahlFall, :AnzahlTodesfall],
        dateformat = "yyyy-mm-dd")
    df = df[from .<= df.Meldedatum .<= to, :]
    if !isnothing(s.counties)
        df = df[in.(AGS.(rki_county.(df.IdLandkreis) .* 1000), Ref(Set(s.counties))), :]
    end
    # NeuerFall / NeuerTodesfall == -1: removed since the previous day's version
    cases = [sum(df.AnzahlFall[(df.Altersgruppe .== g) .& (df.NeuerFall .!= -1)]; init = 0)
        for g in last.(RKI_AGE_GROUPS)]
    deaths = [sum(df.AnzahlTodesfall[(df.Altersgruppe .== g) .& in.(df.NeuerTodesfall, Ref((0, 1)))]; init = 0)
        for g in last.(RKI_AGE_GROUPS)]
    return DataFrame(age_group = last.(RKI_AGE_GROUPS), cases = cases, deaths = deaths)
end


"""
    infection_age_table(inf, last_tick)

Infections, detected infections, and deaths per RKI age group of a run (`inf`:
`infections` of the post processor), each up to `last_tick`.
"""
function infection_age_table(inf, last_tick::Integer)
    n = length(RKI_AGE_GROUPS)
    infected, detected, deaths = zeros(Int, n), zeros(Int, n), zeros(Int, n)
    for (a, t, d, death) in zip(inf.age_b, inf.tick, inf.first_detected_tick, inf.death)
        (ismissing(a) || t > last_tick) && continue
        g = rki_age_index(Int(a))
        infected[g] += 1
        detected[g] += !ismissing(d) && d <= last_tick
        deaths[g] += 0 <= death <= last_tick
    end
    return DataFrame(age_group = last.(RKI_AGE_GROUPS), infections = infected,
        detected = detected, deaths = deaths)
end


"""
    mixing_outcomes(inf, days)

Outcomes of a run of `days` days that depend on who meets whom: mean age of the
infected, deaths up to the last day (as `infection_age_table`) and the infection
fatality rate of the infections at least `MATURE_DAYS` before the end, share of
the transmissions (seeded infections left out) per setting layer, and the share
of municipality transmissions between individuals of the same `SAME_AGE_BAND`.
"""
function mixing_outcomes(inf, days::Integer)
    last_tick = days - 1
    died = 0 .<= inf.death .<= last_tick
    mature = inf.tick .<= days - MATURE_DAYS
    transmitted = inf.id_a .!= -1 # -1: seeded infection
    layer = [get(SETTING_LAYERS, c, :other) for c in inf.setting_type[transmitted]]
    share(l) = mean(layer .== l)
    municipality = transmitted .& (inf.setting_type .== 'm')
    band(a) = Int(a) ÷ SAME_AGE_BAND
    return (
        mean_age = mean(Int.(skipmissing(inf.age_b))),
        deaths = count(died),
        ifr = count(died[mature]) / max(count(mature), 1),
        share_household = share(:household),
        share_municipality = share(:municipality),
        share_work = share(:work),
        share_school = share(:school),
        municipality_same_age = mean(band.(inf.age_a[municipality]) .== band.(inf.age_b[municipality])))
end


"""
    national_fit(daily, weeks)

Fit of the national reported cases of GEMS to RKI (`daily` from
`rki_daily_summary`) in the ISO `weeks`: ratio GEMS/RKI of the weekly cases
(lowest and highest week), root mean square of the log ratios, and the reported
cases of the whole run.
"""
function national_fit(daily, weeks)
    w = week_2020.(daily.date)
    ratio = [sum(daily.gems[w .== k]) / sum(daily.rki[w .== k]) for k in weeks]
    return (
        ratio_low = minimum(ratio),
        ratio_high = maximum(ratio),
        rms_log_ratio = sqrt(mean(log.(ratio) .^ 2)),
        cases_gems = sum(daily.gems),
        cases_rki = sum(daily.rki))
end


"""
    paired_spearman_difference(tbl_a, tbl_b, window; bootstrap = 2000, seed = 0)

Difference of the county Spearman correlations with RKI (incidence of `window`,
see `rki_cumulative_table`) of two county tables of the same counties, `a` minus
`b`, with its 95% interval and p-value from a bootstrap over states in which both
use the same samples (paired, as `rki_weekly_metrics`).
"""
function paired_spearman_difference(tbl_a::DataFrame, tbl_b::DataFrame, window;
    bootstrap::Int = 2000,
    seed = 0)

    tbl_a.ags == tbl_b.ags || throw(ArgumentError("The tables have different counties."))
    ri = tbl_a[!, county_column("rki", "incidence", window)]
    ga = tbl_a[!, county_column("gems", "incidence", window)]
    gb = tbl_b[!, county_column("gems", "incidence", window)]
    rng = Random.Xoshiro(seed)
    states = sort(unique(tbl_a.state))
    rows = Dict(s => findall(tbl_a.state .== s) for s in states)
    diffs = map(1:bootstrap) do _
        i = reduce(vcat, [rows[s] for s in rand(rng, states, length(states))])
        corspearman(ri[i], ga[i]) - corspearman(ri[i], gb[i])
    end
    diffs = filter(!isnan, diffs)
    return (
        difference = corspearman(ri, ga) - corspearman(ri, gb),
        low = quantile(diffs, 0.025),
        high = quantile(diffs, 0.975),
        p = _bootstrap_p(diffs))
end

"""
    run_incidences(tbls, window)

GEMS incidence per county (`window`, see `rki_cumulative_table`) of each single
run, from the county tables of the runs (`tbls`, see `run_rki_phase_sims`), as a
long table (`run`, `ags`, `incidence`).
"""
function run_incidences(tbls::Vector{DataFrame}, window)
    col = county_column("gems", "incidence", window)
    runs = [DataFrame(run = k, ags = map(a -> a.id, t.ags),
            incidence = rki_cumulative_table(t, RKI_WINDOWS)[!, col])
        for (k, t) in enumerate(tbls)]
    return reduce(vcat, runs; init = DataFrame(run = Int[], ags = Int[], incidence = Float64[]))
end

"""
    mean_run_spearman(a, b)

Mean Spearman correlation of the county incidences of single runs
(`run_incidences`) between the runs in `a` and those in `b`, over all pairs of
runs with different seeds (run `k` of both tables has the same seed). With `a`
and `b` from one model, this is the agreement between its seeds, the reference
for the agreement between two models. `NaN` without such pairs.
"""
function mean_run_spearman(a::DataFrame, b::DataFrame)
    incidences(d) = [sort(d[d.run .== k, :], :ags) for k in sort(unique(d.run))]
    ra, rb = incidences(a), incidences(b)
    pairs = [(i, j) for i in eachindex(ra) for j in eachindex(rb) if i != j]
    isempty(pairs) && return NaN
    all(ra[i].ags == rb[j].ags for (i, j) in pairs) ||
        throw(ArgumentError("The runs have different counties."))
    return mean(corspearman(ra[i].incidence, rb[j].incidence) for (i, j) in pairs)
end


"""
    mixing_age_table(rki, uniform, age_based)

Shares per RKI age group of the reported (RKI) or detected (GEMS) cases, the
deaths, and (GEMS) all infections, for RKI (`rki_age_table`) and the runs with
uniform and age-based municipality contacts (`infection_age_table`, mean of the
runs).
"""
function mixing_age_table(rki, uniform, age_based)
    shares(x) = x ./ sum(x)
    return DataFrame(
        age_group = rki.age_group,
        cases_rki = shares(rki.cases),
        cases_uniform = shares(uniform.detected),
        cases_age_based = shares(age_based.detected),
        deaths_rki = shares(rki.deaths),
        deaths_uniform = shares(uniform.deaths),
        deaths_age_based = shares(age_based.deaths),
        infections_uniform = shares(uniform.infections),
        infections_age_based = shares(age_based.infections))
end


"""
    mixing_summary(; uniform, age_based, ages, run_spearman)

Comparison of uniform and age-based municipality contacts, one row per outcome
(`block`, `outcome`, `uniform`, `age_based`, `rki`, and for the county
correlations the paired difference age-based minus uniform with its 95%
interval). `uniform` and `age_based` are named tuples with the results of the
runs with the RKI introductions (`final`: `daily`, `counties`, `metrics`,
`stats`) and with introductions in random counties (`random`: `counties`,
`metrics`), `ages` comes from `mixing_age_table`, and `run_spearman` holds the
mean county correlations of single runs with the RKI introductions
(`mean_run_spearman`) between the seeds of each model (`uniform`, `age_based`)
and between the models (`between`).
"""
function mixing_summary(; uniform, age_based, ages, run_spearman)
    rows = NamedTuple[]
    add!(block, outcome, u, a; rki = missing, d = nothing) = push!(rows, (
        block = block, outcome = outcome, uniform = u, age_based = a, rki = rki,
        difference = isnothing(d) ? missing : d.difference,
        difference_low = isnothing(d) ? missing : d.low,
        difference_high = isnothing(d) ? missing : d.high,
        difference_p = isnothing(d) ? missing : d.p))

    # national fit in the fitted weeks
    fu, fa = (national_fit(x.final.daily, RKI_FIT_WEEKS) for x in (uniform, age_based))
    add!("national", "weekly cases GEMS/RKI, lowest fitted week", fu.ratio_low, fa.ratio_low)
    add!("national", "weekly cases GEMS/RKI, highest fitted week", fu.ratio_high, fa.ratio_high)
    add!("national", "RMS of the log ratio, fitted weeks", fu.rms_log_ratio, fa.rms_log_ratio)
    add!("national", "reported cases (whole run)", fu.cases_gems, fa.cases_gems, rki = fu.cases_rki)

    # county pattern
    metric(x, part, window, col) = only(x[part].metrics[x[part].metrics.week .== window, col])
    for (part, label) in ((:final, "RKI introductions"), (:random, "introductions in random counties"))
        for window in first.(RKI_WINDOWS)
            d = paired_spearman_difference(age_based[part].counties, uniform[part].counties, window)
            add!("spatial", "county Spearman with RKI, $window, $label",
                metric(uniform, part, window, :spearman), metric(age_based, part, window, :spearman), d = d)
        end
        add!("spatial", "state Spearman with RKI, fit, $label",
            metric(uniform, part, "fit", :spearman_states), metric(age_based, part, "fit", :spearman_states))
    end
    add!("spatial", "gain over the introductions, fit, introductions in random counties",
        metric(uniform, :random, "fit", :gain_seeds), metric(age_based, :random, "fit", :gain_seeds))
    add!("spatial", "mean county Spearman of single runs between seeds of the model, fit",
        run_spearman.uniform, run_spearman.age_based)
    add!("spatial", "mean county Spearman of single runs between both models, fit",
        run_spearman.between, run_spearman.between)

    # age and severity
    stat(x, col) = mean(x.final.stats[!, col])
    add!("age", "mean age of the infected", stat(uniform, :mean_age), stat(age_based, :mean_age))
    older(col) = sum(ages[ages.age_group .∈ Ref(["A60-A79", "A80+"]), col])
    add!("age", "share of cases aged 60+", older(:cases_uniform), older(:cases_age_based),
        rki = older(:cases_rki))
    oldest(col) = only(ages[ages.age_group .== "A80+", col])
    add!("age", "share of deaths aged 80+", oldest(:deaths_uniform), oldest(:deaths_age_based),
        rki = oldest(:deaths_rki))
    add!("age", "deaths (up to the last day)", stat(uniform, :deaths), stat(age_based, :deaths))
    add!("age", "infection fatality rate", stat(uniform, :ifr), stat(age_based, :ifr))

    # settings (check that only who meets whom changed)
    for (col, label) in ((:share_household, "household"), (:share_municipality, "municipality"),
            (:share_work, "work"), (:share_school, "education"))
        add!("settings", "share of transmissions: $label", stat(uniform, col), stat(age_based, col))
    end
    add!("settings", "municipality transmissions within the same $SAME_AGE_BAND-year age band",
        stat(uniform, :municipality_same_age), stat(age_based, :municipality_same_age))
    return DataFrame(rows)
end
