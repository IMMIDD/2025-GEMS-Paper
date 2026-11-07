###
### FIRST WAVE 2020: SCHEDULED SEEDING (VIA STEPMOD)
###

# Workaround for GEMS v1.0.2: RegionalSeeds only seeds at tick 0 and filters
# counties with in_state instead of in_county, so we seed via stepmod instead.


"""
    add_seeding!(sim::Simulation, schedule::DataFrame; seed = 0)

Attaches a seeding schedule to the simulation (see `seeding_stepmod`), after
the existing `stepmod`. Also adds a tick trigger with an empty strategy that
keeps the simulation awake: GEMS skips ticks without active infections and
does not notice infections added via `stepmod`.
"""
function add_seeding!(sim::Simulation, schedule::DataFrame; seed = 0)
    inner = sim.stepmod
    seeding = seeding_stepmod(schedule, seed = seed)
    sim.stepmod = s -> begin
        inner(s)
        seeding(s)
    end
    add_tick_trigger!(sim, STickTrigger(Municipality, SStrategy("keep awake (seeding)", sim)))
    return sim
end


# draws n susceptible individuals (indices into inds) from the county index by
# rejection sampling (at most 20 tries per individual); falls back to filtering
# the county if that fails
function _sample_susceptible(rng, cand::Vector{Int32}, inds, n::Int, a::AGS, sim)
    ok(k) = !infected(inds[k]) && !dead(inds[k])
    picked = Set{Int32}()
    attempts = 0
    while length(picked) < n && attempts < 20 * n + 100
        attempts += 1
        k = cand[rand(rng, 1:length(cand))]
        if ok(k)
            push!(picked, k)
        end
    end
    length(picked) == n && return collect(picked)

    candidates = filter(ok, cand)
    if length(candidates) < n
        @warn "Seeding: county $(a.id) has only $(length(candidates)) susceptible individuals, " *
            "$n requested at tick $(tick(sim))."
    end
    return gems_sample(rng, candidates, min(n, length(candidates)), replace = false)
end


"""
    seeding_stepmod(schedule::DataFrame; seed = 0)

Returns a `stepmod` function that infects `count` random susceptible
individuals in county `ags` at tick `tick` (columns of `schedule`).
Counties that are not in the population model are skipped.
"""
function seeding_stepmod(schedule::DataFrame; seed = 0)
    for c in ["tick", "ags", "count"]
        c in names(schedule) || throw(ArgumentError("The seeding schedule needs a '$c' column."))
    end
    # county level; city-states (Berlin, Hamburg) count as counties here although is_county is false
    all(a -> county(a) == a, schedule.ags) ||
        throw(ArgumentError("All AGS in the seeding schedule must be counties."))

    # tick -> [(county, count), ...]
    by_tick = Dict{Int, Vector{Tuple{AGS, Int}}}()
    for r in eachrow(schedule)
        if r.count > 0
            push!(get!(by_tick, r.tick, Tuple{AGS, Int}[]), (r.ags, r.count))
        end
    end

    rng = Random.Xoshiro(seed) # infect! requires a Xoshiro RNG
    county_index = Dict{AGS, Vector{Int32}}() # county -> individual indices, built on first call
    skipped = Set{AGS}()

    return function (sim)
        entries = get(by_tick, Int(tick(sim)), nothing)
        isnothing(entries) && return sim

        inds = individuals(sim)
        if isempty(county_index)
            for (k, i) in enumerate(inds)
                push!(get!(county_index, gems_county(ags(household(i, sim))), Int32[]), Int32(k))
            end
        end

        for (a, n) in entries
            if !haskey(county_index, a)
                if !(a in skipped)
                    @warn "Seeding: county $(a.id) is not in the population model and is skipped."
                end
                push!(skipped, a)
                continue
            end

            for k in _sample_susceptible(rng, county_index[a], inds, n, a, sim)
                i = inds[k]
                infect!(i, tick(sim), pathogen(sim), sim = sim, rng = rng)
                # settings of an individual are only simulated while active
                for (type, id) in settings_tuple(i)
                    id != GEMS.DEFAULT_SETTING_ID && activate!(settings(sim, type)[id], sim)
                end
            end
        end
        return sim
    end
end


"""
    merge_schedule(schedule::DataFrame)

Seeding schedule with one row per tick and county, sorted by tick and county.
"""
merge_schedule(schedule::DataFrame) = schedule |>
    d -> groupby(d, [:tick, :ags]) |>
    d -> combine(d, :count => sum => :count) |>
    d -> sort(d, [:tick, order(:ags, by = a -> a.id)])


"""
    MEDIAN_REPORT_DELAY

Median delay (days) from symptom onset to reporting in the RKI data.
"""
const MEDIAN_REPORT_DELAY = 5


"""
    rki_seeding_schedule(df::DataFrame;
        startdate::Date,
        cutoff::Date,
        incubation::Int,
        reported_by = cutoff + Day(7),
        report_delay::Int = MEDIAN_REPORT_DELAY)

Seeding schedule from RKI cases with onset up to `cutoff` that were reported
by `reported_by` (`nothing`: all). Cases are seeded at onset minus `incubation`,
or, if the onset is unknown, at reporting date minus `incubation + report_delay`.
"""
function rki_seeding_schedule(df::DataFrame;
    startdate::Date,
    cutoff::Date,
    incubation::Int,
    reported_by::Union{Nothing, Date} = cutoff + Day(7),
    report_delay::Int = MEDIAN_REPORT_DELAY)

    reported = isnothing(reported_by) ? trues(nrow(df)) : df.report .<= reported_by
    cases = df[(df.refdate .<= cutoff) .& reported, :]

    infection_date = [
        ismissing(o) ? r - Day(incubation + report_delay) : o - Day(incubation)
        for (o, r) in zip(cases.onset, cases.report)]
    sched = DataFrame(
        tick = Dates.value.(infection_date .- startdate),
        ags = cases.ags,
        count = cases.cases) |>
        merge_schedule

    dropped = sum(sched.count[sched.tick .< 0], init = 0)
    if dropped > 0
        printinfo("Seeding: $dropped RKI cases would have been infected before $startdate " *
            "and are not seeded.")
    end
    return sched[(sched.tick .>= 0) .& (sched.count .> 0), :]
end
