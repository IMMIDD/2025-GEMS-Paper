###
### FIRST WAVE 2020: COMPARISON OF GEMS AND RKI CASES (NATIONAL, STATES, COUNTIES)
###

"""
    county_column(source, kind, week)

Column of the county tables: `source` ("rki", "gems", or "seeds"), `kind`
("cases", "incidence", or "" for seeds), and ISO week or window name.
"""
function county_column(source, kind, week)
    isempty(kind) && return "$(source)_w$week"
    return "$(source)_$(kind)_w$week"
end


"""
    rki_reported_daily(res, reported)

Daily cases of RKI (column `report` of `res.rki`) and GEMS (`reported`, from
`gems_reported`) over the simulated period of `res` (from `run_rki_sim`).
"""
function rki_reported_daily(res, reported::DataFrame)
    gems = countmap(reported.tick)
    rki = combine(groupby(res.rki, :report), :cases => sum => :cases)
    rki = Dict(zip(rki.report, rki.cases))
    days = collect(res.startdate:Day(1):res.enddate)
    return DataFrame(
        date = days,
        rki = [get(rki, d, 0) for d in days],
        gems = [get(gems, Dates.value(d - res.startdate), 0) for d in days])
end



###
### COUNTIES
###

# adds the `cases` of a county table as column `name` (0 for counties without cases)
function _add_counts(tbl::DataFrame, counts::DataFrame, name)
    d = leftjoin(tbl, rename(counts, :cases => name), on = :ags) |> sort_by_ags
    d[!, name] = coalesce.(d[!, name], 0)
    return d
end


"""
    rki_weekly_table(res, weeks; reported)

RKI cases (column `report` of `res.rki`) and GEMS cases (`reported`, from
`gems_reported`) per county and ISO week, with incidences per 100,000 (RKI with
the official population, GEMS with the population of the model), and the seeds
up to the end of each week. Also contains the columns used for the baselines of
`rki_weekly_metrics`, e.g. the RKI incidence before the first week
(`rki_incidence_before`).
"""
function rki_weekly_table(res, weeks; reported::DataFrame)
    geo = DataFrames.select(county_geography(),
        :ags, :state, :distance_km, :population => :population_official)
    tbl = disallowmissing!(leftjoin(res.pop, geo, on = :ags)) |> sort_by_ags
    for w in weeks
        t_from = Dates.value(monday_2020(w) - res.startdate)
        rki = rki_cases_by_county(res.rki,
            date = :report,
            from = monday_2020(w),
            to = sunday_2020(w))
        gems = reported[t_from .<= reported.tick .<= t_from + 6, :] |>
            d -> combine(groupby(d, :ags), nrow => :cases)
        seeds = res.schedule[res.schedule.tick .<= t_from + 6, :] |>
            d -> combine(groupby(d, :ags), :count => sum => :cases)

        tbl = _add_counts(tbl, rki, county_column("rki", "cases", w))
        tbl[!, county_column("rki", "incidence", w)] =
            tbl[!, county_column("rki", "cases", w)] .* 100_000 ./ tbl.population_official
        tbl = _add_counts(tbl, gems, county_column("gems", "cases", w))
        tbl[!, county_column("gems", "incidence", w)] =
            tbl[!, county_column("gems", "cases", w)] .* 100_000 ./ tbl.population
        tbl = _add_counts(tbl, seeds, county_column("seeds", "", w))
    end

    before = rki_cases_by_county(res.rki,
        date = :report,
        from = Date(2020, 1, 1),
        to = monday_2020(minimum(weeks)) - Day(1))
    tbl = _add_counts(tbl, before, :rki_incidence_before)
    tbl.rki_incidence_before = tbl.rki_incidence_before .* 100_000 ./ tbl.population_official
    return tbl
end


"""
    rki_average_runs(tbls)

County table with the GEMS columns (and seeds) averaged over the tables of
several runs (`rki_weekly_table`, with the same counties).
"""
function rki_average_runs(tbls::Vector{DataFrame})
    tbl = copy(tbls[1])
    all(t -> t.ags == tbl.ags, tbls) || throw(ArgumentError("The tables have different counties."))
    for c in names(tbl)
        if startswith(c, "gems") || startswith(c, "seeds")
            tbl[!, c] = mean(t[!, c] for t in tbls)
        end
    end
    return tbl
end


"""
    rki_cumulative_table(tbl, windows)

Adds the cases and incidences summed over each window `name => weeks` to a
table from `rki_weekly_table` (seeds: up to the end of the window).
"""
function rki_cumulative_table(tbl::DataFrame, windows)
    out = copy(tbl)
    for (name, weeks) in windows
        for src in ("rki", "gems")
            cases = sum(tbl[!, county_column(src, "cases", w)] for w in weeks)
            out[!, county_column(src, "cases", name)] = cases
        end
        out[!, county_column("rki", "incidence", name)] =
            out[!, county_column("rki", "cases", name)] .* 100_000 ./ tbl.population_official
        out[!, county_column("gems", "incidence", name)] =
            out[!, county_column("gems", "cases", name)] .* 100_000 ./ tbl.population
        out[!, county_column("seeds", "", name)] = tbl[!, county_column("seeds", "", last(weeks))]
    end
    return out
end


# incidence per state from county cases and populations
function state_incidence(cases, pop, state)
    return [sum(cases[state .== s]) / sum(pop[state .== s]) for s in sort(unique(state))]
end


"""
    rki_weekly_metrics(tbl, weeks;
        persistence = tbl.rki_incidence_before,
        bootstrap = 500,
        seed = 0)

Spearman correlation between the RKI and GEMS incidence per county for each
week (or window of `rki_cumulative_table`), with a 95% bootstrap interval over
states, and per state. The baselines are the correlations of the RKI incidence
with the population, the seeds per inhabitant, the closeness to Ischgl, and
`persistence` (by default the RKI incidence before the first week of the table).
"""
function rki_weekly_metrics(tbl::DataFrame, weeks;
    persistence = tbl.rki_incidence_before,
    bootstrap::Int = 500,
    seed = 0)

    rng = Random.Xoshiro(seed)
    res = DataFrame()
    states = sort(unique(tbl.state))
    for w in weeks
        rki_cases = tbl[!, county_column("rki", "cases", w)]
        gems_cases = tbl[!, county_column("gems", "cases", w)]
        ri = tbl[!, county_column("rki", "incidence", w)]
        gi = tbl[!, county_column("gems", "incidence", w)]
        seeds = tbl[!, county_column("seeds", "", w)] ./ tbl.population

        # bootstrap over states
        boot = map(1:bootstrap) do _
            i = reduce(vcat, [findall(tbl.state .== s) for s in rand(rng, states, length(states))])
            corspearman(ri[i], gi[i])
        end
        if length(states) < 3
            by_state = NaN
        else
            by_state = corspearman(
                state_incidence(rki_cases, tbl.population_official, tbl.state),
                state_incidence(gems_cases, tbl.population, tbl.state))
        end

        push!(res, (
            week = w,
            spearman = corspearman(ri, gi),
            spearman_low = quantile(boot, 0.025),
            spearman_high = quantile(boot, 0.975),
            spearman_states = by_state,
            base_population = corspearman(ri, tbl.population_official),
            base_seeds = corspearman(ri, seeds),
            base_distance = corspearman(ri, -tbl.distance_km),
            base_persistence = corspearman(ri, persistence)))
    end
    return res
end


"""
    rki_county_metrics(tbl, windows, weeks)

`rki_weekly_metrics` of the windows (`name => weeks`, see
`rki_cumulative_table`) and of single `weeks`. The persistence baseline of a
window is the RKI incidence in the previous window (first window and weeks:
before the first week of the table).
"""
function rki_county_metrics(tbl::DataFrame, windows, weeks)
    m = rki_weekly_metrics(tbl, [first(windows[1])])
    for k in 2:length(windows)
        previous = tbl[!, county_column("rki", "incidence", first(windows[k - 1]))]
        metrics = rki_weekly_metrics(tbl, [first(windows[k])], persistence = previous)
        m = vcat(m, metrics, cols = :union)
    end
    return vcat(m, rki_weekly_metrics(tbl, weeks), cols = :union)
end



###
### STATES
###

"""
    rki_state_daily(res, reported)

Daily cases per state of RKI (column `report` of `res.rki`) and GEMS
(`reported`, from `gems_reported`) over the simulated period of `res`, for the
states of the model.
"""
function rki_state_daily(res, reported::DataFrame)
    days = collect(res.startdate:Day(1):res.enddate)
    index = Dict(d => i for (i, d) in enumerate(days))
    states = sort(unique(state_number.(res.pop.ags)))
    column = Dict(s => j for (j, s) in enumerate(states))
    rki, gems = zeros(length(days), length(states)), zeros(length(days), length(states))
    for (d, a, c) in zip(res.rki.report, res.rki.ags, res.rki.cases)
        i = get(index, d, 0)
        if i > 0
            rki[i, column[state_number(a)]] += c
        end
    end
    for (t, a) in zip(reported.tick, reported.ags)
        i = get(index, res.startdate + Day(t), 0)
        if i > 0
            gems[i, column[state_number(a)]] += 1
        end
    end
    return DataFrame(
        date = repeat(days, length(states)),
        state = repeat(states, inner = length(days)),
        rki = vec(rki),
        gems = vec(gems))
end


"""
    rki_state_summary(tables)

Mean, minimum (`gems_low`), and maximum (`gems_high`) of the GEMS cases over
runs of `rki_state_daily` tables.
"""
function rki_state_summary(tables::Vector{DataFrame})
    return reduce(vcat, tables) |>
        d -> groupby(d, [:date, :state], sort = true) |>
        d -> combine(d,
            :rki => first => :rki,
            :gems => mean => :gems,
            :gems => minimum => :gems_low,
            :gems => maximum => :gems_high)
end


"""
    rki_state_table(tbl, weeks)

Cases per state and week of a county table (`rki_weekly_table`), with the
official and the model population.
"""
function rki_state_table(tbl::DataFrame, weeks)
    cols = [county_column(s, "cases", w) => sum => county_column(s, "cases", w)
        for s in ["rki", "gems"] for w in weeks]
    return tbl |>
        d -> groupby(d, :state) |>
        d -> combine(d,
            :population => sum => :population,
            :population_official => sum => :population_official,
            cols...) |>
        d -> sort(d, :state)
end


"""
    state_row_incidence(r, source, week)

Incidence per 100,000 of `source` ("rki" or "gems") in a week (or window) of a
row of `rki_state_table`.
"""
function state_row_incidence(r, source, week)
    population = source == "rki" ? r.population_official : r.population
    return 1e5 * r[county_column(source, "cases", week)] / population
end
