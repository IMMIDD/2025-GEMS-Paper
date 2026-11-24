###
### FIRST WAVE 2020: OUTCOMES BY AGE AND SETTING, AND THE MIXING COMPARISON
###

# Mixing assumptions for the municipality contacts (who meets whom, not how
# many): the POLYMOD matrix of the reference model, and as alternatives uniform
# mixing and contacts within the own age group only. Outcomes that differ
# between them: age distribution of cases and deaths, fatality, setting shares
# of the transmissions, and the share of municipality transmissions within the
# same age band. `mixing_summary` compares the runs of the assumptions (see the
# `:mixing` part of `scenario_first_wave`).

export MUNICIPALITY_MIXING
export use_municipality_mixing!
export municipality_contact_ratio


"""
    MUNICIPALITY_MIXING

Alternative mixing assumptions for the municipality contacts, compared with the
POLYMOD matrix of the reference model (`use_municipality_mixing!`): `:uniform`
(no age preference, the lower bound of age assortativity) and `:same_age`
(contacts within the own age group only, `same_age_matrix`, the upper bound).
"""
const MUNICIPALITY_MIXING = (:uniform, :same_age)

"""
    same_age_matrix()

Age-age matrix of contacts within the own age group only: the identity matrix
over the age groups of the POLYMOD municipality matrix (`MATRIX_GROUPS` groups
of `MATRIX_INTERVAL` years, the last one open).
"""
same_age_matrix() = ContactMatrix{Float64}(
    Matrix{Float64}(LinearAlgebra.I, MATRIX_GROUPS, MATRIX_GROUPS),
    MATRIX_INTERVAL, MATRIX_GROUPS * MATRIX_INTERVAL)

"""
    use_municipality_mixing!(sim, mixing; contactparameter = nothing)

Replaces the contact sampling of all municipalities of `sim` by the mixing
assumption `mixing` (`MUNICIPALITY_MIXING`), keeping the mean number of contacts
of the config (or `contactparameter`): `:uniform` draws the contacts uniformly
from the members (`ContactparameterSampling`), `:same_age` by age with
`same_age_matrix` (`AgeBasedContactSampling`, each municipality with the age
structure of its own members). Returns the new sampling method.
"""
function use_municipality_mixing!(sim, mixing::Symbol; contactparameter = nothing)
    municipalities = settings(sim, Municipality)
    c = Float64(something(contactparameter, contact_sampling_method(first(municipalities)).contactparameter))
    csm = if mixing == :uniform
        ContactparameterSampling(c)
    elseif mixing == :same_age
        AgeBasedContactSampling(c, MATRIX_INTERVAL, same_age_matrix(), Float64[])
    else
        throw(ArgumentError("Unknown municipality mixing :$mixing, use one of $MUNICIPALITY_MIXING."))
    end
    foreach(s -> contact_sampling_method!(s, csm), municipalities)
    return csm
end

"""
    municipality_contact_ratio(sim)

Expected municipality contacts per member relative to the contact parameter,
averaged over all members of the municipalities of `sim`, for their current
contact sampling methods: shows that a mixing assumption changes who meets whom,
not how many.

Uniform sampling gives each member of a municipality with at least two members
the contact parameter c. `AgeBasedContactSampling` draws Poisson(c m_max /
(q ⋅ m_g)) candidates among the other members (q the age structure of the
municipality, m_g the matrix row of the member's age group g, m_max its maximum)
and keeps each with probability m / m_max. With n_j members in group j of a
municipality of N members, a member of group g then has on average
c S / (S + m_gg) N / (N - 1) contacts, S = Σ_j n_j m_gj - m_gg (fewer than c
only if few members of the preferred ages live in the municipality).
"""
function municipality_contact_ratio(sim)
    total, members = 0.0, 0
    for s in settings(sim, Municipality)
        inds = individuals(s)
        N = length(inds)
        members += N
        N < 2 && continue # no one to meet
        csm = contact_sampling_method(s)
        if !(csm isa AgeBasedContactSampling)
            total += N
            continue
        end
        m = csm.contact_matrix.data
        interval, max_age = csm.contact_matrix.interval_steps, csm.contact_matrix.aggregation_bound - 1
        n = zeros(size(m, 1))
        foreach(i -> n[min(Int(age(i)), max_age) ÷ interval + 1] += 1, inds)
        for g in eachindex(n)
            n[g] == 0 && continue
            T = sum(n[j] * m[g, j] for j in eachindex(n))
            T > 0 || continue
            total += n[g] * (T - m[g, g]) / T * N / (N - 1)
        end
    end
    return total / max(members, 1)
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
    mixing_age_table(rki, ages)

Shares per RKI age group of the reported (RKI) or detected (GEMS) cases, the
deaths, and (GEMS) all infections, for RKI (`rki_age_table`) and the runs of
each mixing assumption (`ages`: `model => infection_age_table`, mean of the
runs).
"""
function mixing_age_table(rki, ages)
    shares(x) = x ./ sum(x)
    columns(kind, col) = [Symbol("$(kind)_$m") => shares(a[!, col]) for (m, a) in ages]
    return DataFrame(
        :age_group => rki.age_group,
        :cases_rki => shares(rki.cases),
        columns("cases", :detected)...,
        :deaths_rki => shares(rki.deaths),
        columns("deaths", :deaths)...,
        columns("infections", :infections)...)
end


"""
    mixing_summary(; models, ages, run_spearman)

Comparison of the municipality mixing assumptions, one row per outcome: `block`,
`outcome`, one column per model, `rki`, and for the county correlations the
paired difference of each alternative minus the reference model (the first of
`models`) with its 95% interval and p-value (`difference_<model>`, `_low`,
`_high`, `_p`, see `paired_spearman_difference`).

`models` holds `model => results` in the order of the columns, each with the
runs with the RKI introductions (`final`: `daily`, `counties`, `metrics`,
`stats`) and with introductions in random counties (`random`: `counties`,
`metrics`). `ages` comes from `mixing_age_table`, and `run_spearman` holds the
mean county correlations of single runs with the RKI introductions
(`mean_run_spearman`) between the seeds of each model (`within`) and between
each alternative and the reference model (`between`), as `model => value`.
"""
function mixing_summary(; models, ages, run_spearman)
    names = first.(models)
    ref, alternatives = first(names), names[2:end]
    res = Dict(models)
    diff_cols(a) = [Symbol("difference_$a$suffix") for suffix in ("", "_low", "_high", "_p")]
    cols = vcat([:block, :outcome], names, [:rki], reduce(vcat, diff_cols.(alternatives); init = Symbol[]))
    rows = Dict{Symbol, Any}[]
    function add!(block, outcome, value; rki = missing, difference = nothing)
        row = Dict{Symbol, Any}(c => missing for c in cols)
        row[:block], row[:outcome], row[:rki] = block, outcome, rki
        foreach(m -> row[m] = value(m), names)
        if !isnothing(difference)
            for a in alternatives
                d = difference(a)
                row[Symbol("difference_$a")] = d.difference
                row[Symbol("difference_$(a)_low")] = d.low
                row[Symbol("difference_$(a)_high")] = d.high
                row[Symbol("difference_$(a)_p")] = d.p
            end
        end
        push!(rows, row)
    end

    # national fit in the fitted weeks
    fit = Dict(m => national_fit(res[m].final.daily, RKI_FIT_WEEKS) for m in names)
    add!("national", "weekly cases GEMS/RKI, lowest fitted week", m -> fit[m].ratio_low)
    add!("national", "weekly cases GEMS/RKI, highest fitted week", m -> fit[m].ratio_high)
    add!("national", "RMS of the log ratio, fitted weeks", m -> fit[m].rms_log_ratio)
    add!("national", "reported cases (whole run)", m -> fit[m].cases_gems, rki = fit[ref].cases_rki)

    # county pattern
    metric(m, part, window, col) = only(res[m][part].metrics[res[m][part].metrics.week .== window, col])
    for (part, label) in ((:final, "RKI introductions"), (:random, "introductions in random counties"))
        for window in first.(RKI_WINDOWS)
            add!("spatial", "county Spearman with RKI, $window, $label",
                m -> metric(m, part, window, :spearman),
                difference = a -> paired_spearman_difference(res[a][part].counties, res[ref][part].counties, window))
        end
        add!("spatial", "state Spearman with RKI, fit, $label", m -> metric(m, part, "fit", :spearman_states))
    end
    add!("spatial", "gain over the introductions, fit, introductions in random counties",
        m -> metric(m, :random, "fit", :gain_seeds))
    add!("spatial", "mean county Spearman of single runs between seeds of the model, fit",
        m -> run_spearman.within[m])
    add!("spatial", "mean county Spearman of single runs with the reference model, fit",
        m -> m == ref ? missing : run_spearman.between[m])

    # age and severity
    stat(m, col) = mean(res[m].final.stats[!, col])
    add!("age", "mean age of the infected", m -> stat(m, :mean_age))
    older(col) = sum(ages[ages.age_group .∈ Ref(["A60-A79", "A80+"]), col])
    add!("age", "share of cases aged 60+", m -> older(Symbol("cases_$m")), rki = older(:cases_rki))
    oldest(col) = only(ages[ages.age_group .== "A80+", col])
    add!("age", "share of deaths aged 80+", m -> oldest(Symbol("deaths_$m")), rki = oldest(:deaths_rki))
    add!("age", "deaths (up to the last day)", m -> stat(m, :deaths))
    add!("age", "infection fatality rate", m -> stat(m, :ifr))

    # settings (check that only who meets whom changed)
    add!("settings", "expected municipality contacts per member, relative to the contact parameter",
        m -> stat(m, :municipality_contact_ratio))
    for (col, label) in ((:share_household, "household"), (:share_municipality, "municipality"),
            (:share_work, "work"), (:share_school, "education"))
        add!("settings", "share of transmissions: $label", m -> stat(m, col))
    end
    add!("settings", "municipality transmissions within the same $SAME_AGE_BAND-year age band",
        m -> stat(m, :municipality_same_age))
    return DataFrame([c => [r[c] for r in rows] for c in cols])
end
