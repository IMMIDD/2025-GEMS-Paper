###
### FIRST WAVE 2020: SIMULATIONS (SETTINGS, SEEDING, CONTACT PHASES, VARIANTS)
###

export RKISettings
export run_rki_sim
export run_rki_phase_sims


"""
    RKI_R_WINDOW

First and last infection date of the R per phase: 01.03. to 12.04.2020, one
week before the end of the fitted data (19.04.), as the infections of the last
week before it are not yet completely observed.
"""
const RKI_R_WINDOW = (Date(2020, 3, 1), Date(2020, 4, 12))


"""
    RKISettings

Settings of the simulations of the first wave 2020 (see the keyword constructor).
"""
struct RKISettings
    rki::DataFrame
    days::Int
    configfile::String
    population::String
    counties::Union{Nothing, Vector{AGS}}
    onset_share::Function
    ascertainment::Float64
    care_homes::Int
    allocation::ImportParams
    dates::Vector{Date}
    transitions::Vector{Int}
    school_closure::Union{Nothing, Date}
    r_window::Tuple{Date, Date}
    startdate::Date
end


"""
    RKISettings(;
        rki,
        days,
        configfile = "models/DE_model_first_wave.toml",
        population = "DE",
        counties = population_counties(population),
        onset_share = rki_onset_share(rki, counties = counties),
        ascertainment = ASCERTAINMENT_2020,
        care_homes = CARE_HOME_SIZE,
        allocation = ImportParams(),
        dates = PHASE_DATES,
        transitions = PHASE_TRANSITIONS,
        school_closure = SCHOOL_CLOSURE_2020,
        r_window = RKI_R_WINDOW,
        startdate = config_startdate(configfile))

Settings of the simulations: RKI cases (`load_rki_cases`), simulated days,
model, counties of the population (`nothing`: Germany), share of RKI cases with
a known onset (`rki_onset_share`), ascertainment, care home size, import
allocation, contact phases (start dates, transitions, school closure), the
window of the R per phase, and the start date of the config.
"""
function RKISettings(;
    rki,
    days,
    configfile = "models/DE_model_first_wave.toml",
    population = "DE",
    counties = population_counties(population),
    onset_share = rki_onset_share(rki, counties = counties),
    ascertainment = ASCERTAINMENT_2020,
    care_homes = CARE_HOME_SIZE,
    allocation = ImportParams(),
    dates = PHASE_DATES,
    transitions = PHASE_TRANSITIONS,
    school_closure = SCHOOL_CLOSURE_2020,
    r_window = RKI_R_WINDOW,
    startdate = config_startdate(configfile))

    check_transitions(dates, transitions)
    return RKISettings(rki, days, configfile, population, counties, onset_share, ascertainment,
        care_homes, allocation, dates, transitions, school_closure, r_window, startdate)
end


"""
    with_settings(s::RKISettings; changes...)

Copy of `s` with the fields in `changes` replaced, e.g. `with_settings(s, days = 50)`.
"""
function with_settings(s::RKISettings; changes...)
    fields = fieldnames(RKISettings)
    v = merge(NamedTuple{fields}(getfield.(Ref(s), fields)), values(changes))
    return RKISettings((v[f] for f in fields)...)
end


"""
    population_counties(population)

Counties of a GEMS population: `nothing` for Germany ("DE"), the counties of
the state for a state code (`RKI_STATES`).
"""
function population_counties(population)
    population == "DE" && return nothing
    state = findfirst(==(population), RKI_STATES)
    if isnothing(state)
        throw(ArgumentError("Unknown population '$population', pass its `counties`."))
    end
    geo = county_geography()
    return geo.ags[geo.state .== state]
end


"""
    simulation_days(startdate, last_day)

Number of simulated days from `startdate` to `last_day`.
"""
simulation_days(startdate::Date, last_day::Date) = Dates.value(last_day - startdate) + 1



###
### SINGLE RUNS
###

"""
    run_rki_sim(s::RKISettings;
        seed = 0,
        setup = nothing,
        heinsberg = true,
        imports_from = nothing,
        extra_seeds = nothing,
        transmission_rate = nothing,
        switches = Pair{Date, Float64}[])

Runs a simulation for `s.days` days, seeded with the RKI imports (`import_curve`
divided by `s.ascertainment`, allocated by `tourism_weights`) and the Heinsberg
cluster (its RKI cases, 1:1). `setup` is a function `sim -> x` that is called
before the run (e.g. `phase_setup`).

Variants: `heinsberg = false` drops the Heinsberg cluster, `imports_from` drops
the imports infected before that date, `extra_seeds` is an additional seeding
schedule, `transmission_rate` replaces the rate of the config (the `transmission_scale`
of an `InfectiousnessProfile`), and `switches` (`date => rate`) change it during the run.

Returns `sim`, the post processor `pp`, the seeding `schedule`, the RKI cases
in the counties of the model (`rki`), the population per county (`pop`), the
reporting `delays` and the incubation `model`, `startdate`, `enddate`, and the
return value of `setup`.
"""
function run_rki_sim(s::RKISettings;
    seed = 0,
    setup = nothing,
    heinsberg = true,
    imports_from = nothing,
    extra_seeds = nothing,
    transmission_rate = nothing,
    switches = Pair{Date, Float64}[])

    sim = Simulation(
        configfile = s.configfile,
        population = s.population,
        stop_criterion = TimesUp(limit = s.days),
        seed = seed)
    # keeps the transmission function of the config (Simulation(transmission_rate = ...) would
    # replace it by a ConstantTransmissionRate)
    if !isnothing(transmission_rate)
        pathogen(sim).transmission_function = with_transmission_scale(pathogen(sim).transmission_function, transmission_rate)
    end
    if sim.startdate != s.startdate
        throw(ArgumentError("The config starts on $(sim.startdate), not on $(s.startdate)."))
    end
    enddate = sim.startdate + Day(s.days - 1)
    pop = county_population(sim)
    model = model_incubation(s.configfile)

    # reporting delays per onset week, from all of Germany
    delays = rki_onset_delays_by_week(s.rki)
    imports = _import_seeds(s, pop, delays, model,
        enddate = enddate,
        imports_from = imports_from,
        seed = seed)
    rki = s.rki[in.(s.rki.ags, Ref(Set(pop.ags))), :]

    parts = [imports]
    if heinsberg && HEINSBERG in pop.ags
        cluster = heinsberg_schedule(rki,
            startdate = sim.startdate,
            incubation = model.incubation)
        push!(parts, cluster)
    end
    if !isnothing(extra_seeds)
        push!(parts, extra_seeds)
    end
    schedule = vcat(parts...) |>
        d -> groupby(d, [:tick, :ags]) |>
        d -> combine(d, :count => sum => :count)
    others = sum(schedule.count) - sum(imports.count)
    printinfo("Seeding $(sum(imports.count)) imports and $others other infections")

    add_seeding!(sim, schedule, seed = seed)
    if !isempty(switches)
        _add_rate_switches!(sim, switches)
    end
    setup_info = isnothing(setup) ? nothing : setup(sim)

    run!(sim)
    # we only need the infections; ResultData would compute all default analyses
    pp = PostProcessor(sim)

    return (
        sim = sim,
        pp = pp,
        schedule = schedule,
        rki = rki,
        pop = pop,
        delays = delays,
        model = model,
        startdate = sim.startdate,
        enddate = enddate,
        setup = setup_info)
end


# seeding schedule of the imports until `enddate` (import curve / ascertainment, allocated by the
# tourism model with the county populations of the model, so that the seeds per inhabitant are as
# intended in the model), without the imports infected before `imports_from`
function _import_seeds(s::RKISettings, pop, delays, model; enddate, imports_from, seed)
    curve = import_curve(delays, model, startdate = s.startdate)
    geo = county_geography()
    model_population = Dict(zip(pop.ags, pop.population))
    geo.population = [get(model_population, a, p) for (a, p) in zip(geo.ags, geo.population)]
    imports = import_schedule(curve, tourism_weights_by_date(geo, s.allocation),
        startdate = s.startdate,
        enddate = enddate,
        scale = 1 / s.ascertainment,
        counties = pop.ags,
        seed = seed)
    if !isnothing(imports_from)
        imports = imports[imports.tick .>= Dates.value(imports_from - s.startdate), :]
    end
    return imports
end


# changes the transmission rate to `rate` from `date` on, for each `date => rate` of `switches`
function _add_rate_switches!(sim, switches)
    inner = sim.stepmod
    # stepmod runs at the end of a tick, so a new rate applies from its date
    rate_at = Dict(Dates.value(d - sim.startdate) - 1 => r for (d, r) in switches)
    sim.stepmod = s -> begin
        inner(s)
        r = get(rate_at, Int(tick(s)), nothing)
        if !isnothing(r)
            pathogen(s).transmission_function = with_transmission_scale(pathogen(s).transmission_function, r)
        end
        s
    end
    return sim
end


"""
    phase_setup(s::RKISettings, factors; seed = 0, school_reopening = nothing)

`setup` function for `run_rki_sim` that adds the contact phases of `s` with
the contact `factors` and the testing for the ascertainment of `s`. Returns
the probability that a symptomatic case is tested and the hospitalized share.
"""
function phase_setup(s::RKISettings, factors; seed = 0, school_reopening = nothing)
    return function (sim)
        shares = progression_shares(sim)
        p = symptom_test_probability(s.ascertainment, shares.symptomatic, shares.hospitalized)
        add_contact_phases!(sim, factors,
            dates = s.dates,
            transitions = s.transitions,
            school_closure = s.school_closure,
            school_reopening = school_reopening,
            care_homes = s.care_homes)
        add_pcr_testing!(sim, p_symptomatic = p, seed = seed)
        return (p_symptomatic = p, hospitalized = shares.hospitalized)
    end
end


"""
    gems_reported(pp, delays, model; startdate, reporting_fraction = 1.0, seed = 0)

Reported cases of a simulation: the detected infections, each kept with
probability `reporting_fraction` (a number or a function of the onset date),
reported at symptom onset plus a delay from the RKI delays of the onset week
(`delays = nothing`: at onset). Columns `tick` and `ags`.
"""
function gems_reported(pp, delays, model; startdate::Date, reporting_fraction = 1.0, seed = 0)
    inf = infections(pp)
    rng = Random.Xoshiro(seed)
    n = nrow(inf)
    # infections without symptoms: infection plus the incubation period
    symptomatic = inf.symptom_onset .>= 0
    onset = ifelse.(symptomatic, Int.(inf.symptom_onset), Int.(inf.tick) .+ model.incubation)
    if reporting_fraction isa Real
        p = fill(Float64(reporting_fraction), n)
    else
        p = [reporting_fraction(startdate + Day(o)) for o in onset]
    end
    reported = (rand(rng, n) .< p) .& .!ismissing.(inf.first_detected_tick)

    delay = zeros(Int, n)
    week = firstdayofweek.(startdate .+ Day.(onset))
    for w in unique(week)
        idx = findall(==(w), week)
        dist = delay_distribution(delays, w)
        delay[idx] = sample(rng, dist.delay, Weights(dist.cases), length(idx))
    end
    # household_ags_b: county of the infected individual's household
    return DataFrame(
        tick = (onset .+ delay)[reported],
        ags = gems_county.(inf.household_ags_b[reported]))
end


"""
    phase_shares(pp, startdate; dates, transitions, window)

Share of transmissions outside households per phase, for the infections within
`window` (first and last date).
"""
function phase_shares(pp, startdate::Date; dates, transitions, window)
    inf = infections(pp)
    t = Int.(inf.tick)
    first_tick, last_tick = (Dates.value(d - startdate) for d in window)
    # id_a == -1: seeded infections
    keep = (inf.id_a .!= -1) .& (first_tick .<= t .<= last_tick)
    V = phase_matrix(first_tick:last_tick, [Dates.value(d - startdate) for d in dates], transitions)
    Vt = V[t[keep] .- first_tick .+ 1, :]
    outside = inf.setting_type[keep] .!= 'h' # 'h': household
    return vec(sum(Vt .* outside, dims = 1) ./ sum(Vt, dims = 1))
end



###
### RUNS WITH CONTACT PHASES
###

"""
    rki_daily_summary(dailies; quantiles = nothing)

Mean of the daily GEMS cases of several runs, with their range (`gems_low`,
`gems_high`), or the `quantiles` (low, high) if given.
"""
function rki_daily_summary(dailies::Vector{DataFrame}; quantiles = nothing)
    g = reduce(hcat, [d.gems for d in dailies])
    if isnothing(quantiles)
        low, high = vec(minimum(g, dims = 2)), vec(maximum(g, dims = 2))
    else
        low = [quantile(r, quantiles[1]) for r in eachrow(g)]
        high = [quantile(r, quantiles[2]) for r in eachrow(g)]
    end
    return DataFrame(
        date = dailies[1].date,
        rki = dailies[1].rki,
        gems = vec(mean(g, dims = 2)),
        gems_low = low,
        gems_high = high)
end


# one run with contact phases and testing: per basis (reporting date, and symptom onset with the
# GEMS detected infections thinned to the share of RKI cases with a known onset) the daily cases,
# cases per state, and the county table; the share of transmissions outside households per phase;
# and statistics of the run
function _phase_run(s::RKISettings, factors, seed; map_weeks, stage, school_reopening, sim_options)
    printinfo("$stage: seed $seed, factors $(round.(factors, digits = 4))")
    t0 = time()
    setup = phase_setup(s, factors, seed = seed, school_reopening = school_reopening)
    res = run_rki_sim(s; seed = seed, setup = setup, sim_options...)
    runtime = time() - t0

    by_report = gems_reported(res.pp, res.delays, res.model,
        startdate = res.startdate,
        seed = seed)
    by_onset = gems_reported(res.pp, nothing, res.model,
        startdate = res.startdate,
        reporting_fraction = s.onset_share,
        seed = seed)
    res_onset = merge(res, (rki = onset_basis(res.rki),))
    function basis_results(r, reported)
        tbl = isempty(map_weeks) ? nothing : rki_weekly_table(r, map_weeks, reported = reported)
        return (
            daily = rki_reported_daily(r, reported),
            state = rki_state_daily(r, reported),
            tbl = tbl)
    end
    shares = phase_shares(res.pp, res.startdate,
        dates = s.dates,
        transitions = s.transitions,
        window = s.r_window)

    # ascertainment of the infections at least 21 days before the end (symptoms and hospital
    # admission within the run)
    inf = infections(res.pp)
    mature = inf.tick .<= s.days - 21
    stats = merge(
        (stage = stage,
            seed = seed,
            runtime_s = round(runtime, digits = 1),
            infections = nrow(inf),
            detected = count(!ismissing, inf.first_detected_tick),
            ascertainment = mean(.!ismissing.(inf.first_detected_tick[mature])),
            p_symptomatic = res.setup.p_symptomatic,
            hospitalized = res.setup.hospitalized),
        (; (Symbol("f$k") => x for (k, x) in enumerate(factors))...),
        (; (Symbol("share_outside_$k") => x for (k, x) in enumerate(shares))...))
    printinfo("$stage: seed $seed done after $(round(runtime / 60, digits = 1)) min, " *
        "ascertainment $(round(stats.ascertainment, digits = 3))")
    return (
        report = basis_results(res, by_report),
        onset = basis_results(res_onset, by_onset),
        shares = shares,
        stats = stats,
        model = res.model)
end


"""
    run_rki_phase_sims(s::RKISettings, factors, seeds;
        map_weeks = Int[],
        stage = "",
        school_reopening = nothing,
        sim_options...)

Runs GEMS with the contact phases of `s` and testing, one run per factor
vector and seed. Returns per basis (`report`: by reporting date, `onset`: by
symptom onset) the daily cases as mean and range of the runs (`daily`), the
county tables of the runs for `map_weeks` (`tbls`), and the daily cases per
state (`states`); further the statistics per run (`stats`), the mean share of
transmissions outside households per phase (`shares`), and the incubation
`model`. `sim_options` (seeding variants) are passed to `run_rki_sim`.
"""
function run_rki_phase_sims(s::RKISettings, factors::Vector{<:AbstractVector}, seeds;
    map_weeks = Int[],
    stage = "",
    school_reopening = nothing,
    sim_options...)

    length(factors) == length(seeds) || throw(ArgumentError("Need one factor vector per seed."))
    runs = map(zip(factors, seeds)) do (f, seed)
        r = _phase_run(s, f, seed,
            map_weeks = map_weeks,
            stage = stage,
            school_reopening = school_reopening,
            sim_options = sim_options)
        # free the memory of the simulation before the next one
        GC.gc()
        r
    end
    function summary(b)
        return (
            daily = rki_daily_summary([r[b].daily for r in runs]),
            tbls = isempty(map_weeks) ? DataFrame[] : [r[b].tbl for r in runs],
            states = rki_state_summary([r[b].state for r in runs]))
    end
    return (
        report = summary(:report),
        onset = summary(:onset),
        stats = DataFrame([r.stats for r in runs]),
        shares = vec(mean(reduce(hcat, [r.shares for r in runs]), dims = 2)),
        model = last(runs).model)
end



###
### PROJECTIONS AND COUNTERFACTUALS
###

# label of a projection that returns to the contacts of phase k (of K) after `start`
function _projection_label(k, K, start, dates, schools_open)
    if k == K
        return "contacts kept after $(ddmm(start))"
    elseif k == 1
        return "back to contacts before $(ddmm(dates[1])) after $(ddmm(start)) (no measures)"
    else
        return "back to contacts from $(ddmm(dates[k - 1])) after $(ddmm(start))" *
            (schools_open ? " (schools open)" : "")
    end
end


"""
    run_rki_projections(s::RKISettings, factors, seeds; start, phases)

Runs that return to the contacts of phase k on the day after `start`, for each
k in `phases` (schools reopen if they were open in phase k). Returns the daily
cases by reporting date (`report`) and by onset (`onset`) as `label => daily`
pairs, and the statistics of the runs.
"""
function run_rki_projections(s::RKISettings, factors, seeds; start::Date, phases)
    report, onset, stats = Pair{String, DataFrame}[], Pair{String, DataFrame}[], DataFrame()
    for k in phases
        # phase k starts at dates[k - 1]; schools are open in phases that start before their closure
        schools_open = !isnothing(s.school_closure) && (k == 1 || s.dates[k - 1] < s.school_closure)
        label = _projection_label(k, length(factors), start, s.dates, schools_open)
        # the projection is an additional phase from the day after `start` on, without transition
        projected = with_settings(s,
            dates = vcat(s.dates, start + Day(1)),
            transitions = vcat(s.transitions, 0))
        out = run_rki_phase_sims(projected, fill(vcat(factors, factors[k]), length(seeds)), seeds,
            school_reopening = schools_open ? start + Day(1) : nothing,
            stage = label)
        push!(report, label => out.report.daily)
        push!(onset, label => out.onset.daily)
        stats = vcat(stats, insertcols(out.stats, :phase => k), cols = :union)
    end
    return (report = report, onset = onset, stats = stats)
end


"""
    phase_counterfactuals(factors; dates, school_closure, phases = 1:length(dates))

Counterfactuals in which the contacts of phase k are kept from the start of
phase k + 1 on, for each k in `phases`. Returns `label => (factors,
school_closure)` pairs.
"""
function phase_counterfactuals(factors::AbstractVector{<:Real};
    dates,
    school_closure,
    phases = 1:length(dates))

    return map(phases) do k
        start = k == 1 ? Date(0) : dates[k - 1]
        # schools only close if they closed in phase k
        closure = !isnothing(school_closure) && school_closure <= start ? school_closure : nothing
        if k == 1
            label = "no reduction (phase 1 contacts kept)"
        else
            label = "phase $k contacts kept from $(ddmm(dates[k]))"
        end
        kept = vcat(factors[1:k], fill(factors[k], length(factors) - k))
        label => (factors = kept, school_closure = closure)
    end
end


"""
    rki_seeded_infections(s::RKISettings)

Expected seeded infections per day of `run_rki_sim` in the counties of `s`:
imports from Austria, Italy, and other countries (import curve / ascertainment,
allocated by `tourism_weights`), and the Heinsberg cluster (its RKI cases, 1:1).
"""
function rki_seeded_infections(s::RKISettings)
    model = model_incubation(s.configfile)
    curve = import_curve(rki_onset_delays_by_week(s.rki), model, startdate = s.startdate)
    weights = tourism_weights_by_date(county_geography(), s.allocation)
    imports = expected_imports(curve, weights, counties = s.counties)

    heinsberg = zeros(nrow(imports))
    if isnothing(s.counties) || HEINSBERG in s.counties
        sched = heinsberg_schedule(s.rki, startdate = s.startdate, incubation = model.incubation)
        for r in eachrow(sched)
            if r.tick < length(heinsberg)
                heinsberg[r.tick + 1] += r.count
            end
        end
    end
    return DataFrame(
        date = imports.date,
        austria = imports.austria ./ s.ascertainment,
        italy = imports.italy ./ s.ascertainment,
        other = imports.other ./ s.ascertainment,
        heinsberg = heinsberg)
end
