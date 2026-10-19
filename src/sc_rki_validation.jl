###
### RKI VALIDATION SCENARIO (FIRST WAVE 2020)
###

export gems_active_by_county
export rki_active_by_county
export run_rki_quick_sim
export rki_comparison_table
export rki_comparison_metrics
export plot_rki_quick_panel
export scenario_rki_quick

export gems_onset_by_county
export run_rki_import_sim
export rki_weekly_table
export rki_weekly_metrics
export plot_rki_weekly_panel
export plot_rki_epicurve
export plot_rki_epicurves
export scenario_rki_imports
export rki_national_weekly
export fit_rki_transmission_rate


"""
    gems_active_by_county(pp, t::Integer)

Number of infected individuals per county at tick `t`.
"""
function gems_active_by_county(pp, t::Integer)
    inf = infections(pp)
    active = (inf.tick .<= t) .& (max.(inf.recovery, inf.death) .> t)
    return DataFrame(ags = county.(inf.household_ags_b[active])) |>
        d -> groupby(d, :ags) |>
        d -> combine(d, nrow => :active)
end


"""
    rki_active_by_county(df::DataFrame, d::Date; days = 14)

Number of RKI cases per county with onset (or reporting date) in the `days`
days up to `d`. The RKI data has no recovery dates.
"""
rki_active_by_county(df::DataFrame, d::Date; days::Int = 14) =
    rki_cases_by_county(df, date = :refdate, from = d - Day(days - 1), to = d) |>
        x -> rename(x, :cases => :active)


"""
    run_rki_quick_sim(;
        configfile = "models/DE_model_RKI_2020.toml",
        population = "DE",
        rki = load_rki_cases(),
        cutoff = Date(2020, 3, 1),
        scale = 1.0,
        transmission_rate = nothing,
        incubation = 7,
        seed = 0)

Runs a simulation seeded with all RKI cases up to `cutoff`. `transmission_rate`
overrides the config, `incubation` has to match the model (see `rki_seeding_schedule`).
"""
function run_rki_quick_sim(;
    configfile = "models/DE_model_RKI_2020.toml",
    population = "DE",
    rki::DataFrame = load_rki_cases(),
    cutoff::Date = Date(2020, 3, 1),
    scale::Real = 1.0,
    transmission_rate = nothing,
    incubation::Int = 7,
    seed = 0)

    sim = Simulation(configfile = configfile, population = population,
        transmission_rate = transmission_rate, seed = seed)
    pop = county_population(sim)
    rki = rki[in.(rki.ags, Ref(Set(pop.ags))), :]

    schedule = rki_seeding_schedule(rki, startdate = sim.startdate, cutoff = cutoff,
        scale = scale, incubation = incubation)
    printinfo("Seeding $(sum(schedule.count)) infections in $(length(unique(schedule.ags))) counties")
    add_seeding!(sim, schedule, seed = seed)

    run!(sim)
    # we only need the infections; ResultData would compute all default analyses
    pp = PostProcessor(sim)

    return (sim = sim, pp = pp, schedule = schedule, rki = rki, pop = pop)
end


# incidence per 100,000 per county (0 for counties without cases)
function _incidence(active::DataFrame, pop::DataFrame)
    return leftjoin(pop, active, on = :ags) |>
        d -> transform(d, :active => ByRow(a -> coalesce(a, 0)) => :active) |>
        d -> transform(d, [:active, :population] => ((a, p) -> a .* 100_000 ./ p) => :incidence) |>
        d -> sort(d, :ags, by = a -> a.id)
end


# averages the GEMS columns (and seeds) of per-run county tables
function _average_runs(tbls::Vector{DataFrame})
    tbl = copy(tbls[1])
    for c in names(tbl)
        (startswith(c, "gems") || c == "seeds") || continue
        tbl[!, c] = mean(t[!, c] for t in tbls)
    end
    return tbl
end


"""
    rki_comparison_table(res, dates; days = 14)

Active and cumulative cases with incidences per county (RKI and GEMS) for all
`dates`. Cumulative: RKI cases and GEMS symptom onsets up to the date.
"""
function rki_comparison_table(res, dates::Vector{Date}; days::Int = 14)
    tbl = sort(res.pop, :ags, by = a -> a.id)
    for d in dates
        t = Dates.value(d - res.sim.startdate)
        tag = Dates.format(d, "mm-dd")
        rki_cum = rki_cases_by_county(res.rki, date = :refdate, from = Date(2020, 1, 1), to = d) |>
            x -> rename(x, :cases => :active)
        gems_cum = gems_onset_by_county(res.pp, 0, t) |> x -> rename(x, :cases => :active)
        for (m, r, g) in [
                ("active", rki_active_by_county(res.rki, d, days = days), gems_active_by_county(res.pp, t)),
                ("cum", rki_cum, gems_cum)]
            r, g = _incidence(r, res.pop), _incidence(g, res.pop)
            tbl[!, "rki_$(m)_$tag"] = r.active
            tbl[!, "rki_$(m)_incidence_$tag"] = r.incidence
            tbl[!, "gems_$(m)_$tag"] = g.active
            tbl[!, "gems_$(m)_incidence_$tag"] = g.incidence
        end
    end
    return tbl
end


"""
    rki_comparison_metrics(tbl, dates; measure = "active")

Spearman correlation of county incidences and R² of log(cases + 1) per date,
with log(population) as baseline. `measure` is "active" or "cum".
"""
function rki_comparison_metrics(tbl::DataFrame, dates::Vector{Date}; measure = "active")
    res = DataFrame(measure = String[], date = Date[], spearman_incidence = Float64[],
        r2_log_cases = Float64[], r2_log_cases_population = Float64[])
    for d in dates
        tag = Dates.format(d, "mm-dd")
        lr = log.(tbl[!, "rki_$(measure)_$tag"] .+ 1)
        push!(res, (
            measure = measure,
            date = d,
            spearman_incidence = corspearman(tbl[!, "rki_$(measure)_incidence_$tag"], tbl[!, "gems_$(measure)_incidence_$tag"]),
            r2_log_cases = cor(lr, log.(tbl[!, "gems_$(measure)_$tag"] .+ 1))^2,
            r2_log_cases_population = cor(lr, log.(tbl.population))^2))
    end
    return res
end


# upper color limit: `cap` quantile of the values (outliers like Heinsberg saturate)
_clim_max(x, cap) = (q = quantile(x, cap); q > 0 ? q : max(maximum(x), 1.0))


"""
    plot_rki_quick_panel(tbl, dates; measure = "active", cap = 0.95)

RKI (left) and GEMS (right) incidence maps, one row per date. Color scales
end at the `cap` quantile of each map.
"""
function plot_rki_quick_panel(tbl::DataFrame, dates::Vector{Date}; measure = "active", cap::Real = 0.95)
    labels = measure == "active" ?
        ("RKI active cases", "GEMS active infections") :
        ("RKI cumulative cases", "GEMS cumulative onsets")
    plts = []
    for d in dates
        tag = Dates.format(d, "mm-dd")
        lbl = Dates.format(d, "d u")
        for (src, title) in [("rki", "$(labels[1]), $lbl"), ("gems", "$(labels[2]), $lbl")]
            df = DataFrame(ags = tbl.ags, incidence = tbl[!, "$(src)_$(measure)_incidence_$tag"])
            push!(plts, agsmap(df,
                title = title,
                clims = (0, _clim_max(df.incidence, cap)),
                titlefontsize = map_titlefontsize(),
                colorbar_titlefontsize = map_legendfontsize(),
                ytickfontsize = map_legendfontsize(),
                linecolor = nothing,
                size = map_dims()))
        end
    end
    return plot(plts...,
        layout = (length(dates), 2),
        size = (2 * map_dims()[1], length(dates) * map_dims()[2]),
        fontfamily = "Times Roman")
end


"""
    scenario_rki_quick(;
        configfile = "models/DE_model_RKI_2020.toml",
        population = "DE",
        cutoff = Date(2020, 3, 1),
        dates = [Date(2020, 3, 1), Date(2020, 3, 15)],
        scale = 1.0,
        transmission_rate = nothing,
        incubation = 7,
        seeds = [0],
        folder = nothing)

Runs the quick RKI comparison (seeds up to `cutoff`) and all analyses.
GEMS results are averaged over one run per seed in `seeds`.
"""
function scenario_rki_quick(;
    configfile = "models/DE_model_RKI_2020.toml",
    population = "DE",
    cutoff::Date = Date(2020, 3, 1),
    dates::Vector{Date} = [Date(2020, 3, 1), Date(2020, 3, 15)],
    scale::Real = 1.0,
    transmission_rate = nothing,
    incubation::Int = 7,
    seeds = [0],
    folder = nothing)

    isnothing(folder) ? folder = result_folder() : nothing

    rki = load_rki_cases()
    tbls = DataFrame[]
    schedule = DataFrame()
    for s in seeds
        res = run_rki_quick_sim(configfile = configfile, population = population, rki = rki,
            cutoff = cutoff, scale = scale, transmission_rate = transmission_rate,
            incubation = incubation, seed = s)
        push!(tbls, rki_comparison_table(res, dates))
        schedule = res.schedule
        res = nothing
        GC.gc()
    end

    tbl = _average_runs(tbls)
    metrics = vcat(rki_comparison_metrics(tbl, dates, measure = "active"),
        rki_comparison_metrics(tbl, dates, measure = "cum"))
    println(metrics)

    CSV.write(joinpath(folder, "sc_rki_quick_counties.csv"),
        transform(tbl, :ags => ByRow(a -> a.id) => :ags))
    CSV.write(joinpath(folder, "sc_rki_quick_metrics.csv"), metrics)
    CSV.write(joinpath(folder, "sc_rki_quick_seeding.csv"),
        transform(schedule, :ags => ByRow(a -> a.id) => :ags))
    png(plot_rki_quick_panel(tbl, dates, measure = "active"), joinpath(folder, "sc_rki_quick_panel.png"))
    png(plot_rki_quick_panel(tbl, dates, measure = "cum"), joinpath(folder, "sc_rki_quick_panel_cumulative.png"))

    return (tbl = tbl, metrics = metrics, folder = folder)
end



###
### IMPORT VERSION
###

"""
    gems_onset_by_county(pp, t_from::Integer, t_to::Integer)

Number of symptom onsets per county between ticks `t_from` and `t_to`.
"""
function gems_onset_by_county(pp, t_from::Integer, t_to::Integer)
    inf = infections(pp)
    sel = (inf.symptom_onset .>= t_from) .& (inf.symptom_onset .<= t_to)
    return DataFrame(ags = county.(inf.household_ags_b[sel])) |>
        d -> groupby(d, :ags) |>
        d -> combine(d, nrow => :cases)
end


"""
    run_rki_import_sim(;
        configfile = "models/DE_model_RKI_2020.toml",
        population = "DE",
        rki = load_rki_cases(),
        decay = 1.0,
        import_scale = 1.0,
        heinsberg = true,
        transmission_rate = nothing,
        incubation = 7,
        seed = 0)

Runs a simulation seeded with the RKI imports (see `import_schedule`) and,
optionally, the Heinsberg cluster. `transmission_rate` overrides the config,
`incubation` has to match the model.
"""
function run_rki_import_sim(;
    configfile = "models/DE_model_RKI_2020.toml",
    population = "DE",
    rki::DataFrame = load_rki_cases(),
    decay::Real = 1.0,
    import_scale::Real = 1.0,
    heinsberg::Bool = true,
    transmission_rate = nothing,
    incubation::Int = 7,
    seed = 0)

    sim = Simulation(configfile = configfile, population = population,
        transmission_rate = transmission_rate, seed = seed)
    pop = county_population(sim)

    # reporting delays from all of Germany
    curve = import_curve(rki_onset_delays(rki), incubation = incubation)
    rki = rki[in.(rki.ags, Ref(Set(pop.ags))), :]

    weights = import_weights(county_geography(), decay = decay)
    imports = import_schedule(curve, weights,
        startdate = sim.startdate, enddate = sim.enddate,
        scale = import_scale, counties = pop.ags, seed = seed)

    schedule = imports
    if heinsberg && AGS(5370000) in pop.ags
        schedule = vcat(imports, heinsberg_schedule(rki, startdate = sim.startdate, incubation = incubation)) |>
            d -> groupby(d, [:tick, :ags]) |>
            d -> combine(d, :count => sum => :count)
    end
    printinfo("Seeding $(sum(imports.count)) imports and $(sum(schedule.count) - sum(imports.count)) Heinsberg infections")
    add_seeding!(sim, schedule, seed = seed)

    run!(sim)
    # we only need the infections; ResultData would compute all default analyses
    pp = PostProcessor(sim)

    # expected imports per county (baseline)
    c = curve[sim.startdate .<= curve.date .<= sim.enddate, :]
    expected = DataFrame(ags = weights.ags,
        expected_imports = import_scale .* (weights.alpine .* sum(c.alpine) .+ weights.other .* sum(c.other)))

    return (sim = sim, pp = pp, schedule = schedule, rki = rki, pop = pop,
        curve = curve, expected_imports = expected)
end


"""
    rki_weekly_table(res, weeks)

RKI cases and GEMS onsets per county and ISO week (2020), with incidences,
seeds, and expected imports.
"""
function rki_weekly_table(res, weeks)
    tbl = sort(res.pop, :ags, by = a -> a.id)
    for w in weeks
        from = _monday_2020(w)
        t_from = Dates.value(from - res.sim.startdate)

        r = rki_cases_by_county(res.rki, date = :refdate, from = from, to = from + Day(6)) |>
            d -> rename(d, :cases => :active)
        g = gems_onset_by_county(res.pp, t_from, t_from + 6) |>
            d -> rename(d, :cases => :active)
        r, g = _incidence(r, res.pop), _incidence(g, res.pop)

        tbl[!, "rki_cases_w$w"] = r.active
        tbl[!, "rki_incidence_w$w"] = r.incidence
        tbl[!, "gems_cases_w$w"] = g.active
        tbl[!, "gems_incidence_w$w"] = g.incidence
    end

    seeds = res.schedule |>
        d -> groupby(d, :ags) |>
        d -> combine(d, :count => sum => :seeds)
    tbl = leftjoin(tbl, seeds, on = :ags)
    tbl.seeds = coalesce.(tbl.seeds, 0)
    if hasproperty(res, :expected_imports)
        tbl = leftjoin(tbl, res.expected_imports, on = :ags)
    end
    return sort(tbl, :ags, by = a -> a.id)
end


"""
    rki_weekly_metrics(tbl, weeks)

Agreement between RKI and GEMS per week: Spearman correlation of incidences,
R² of log(cases + 1), the same for the baselines (population, expected
imports), and the Spearman correlation for counties without seeds.
"""
function rki_weekly_metrics(tbl::DataFrame, weeks)
    res = DataFrame()
    has_imp = hasproperty(tbl, :expected_imports)
    unseeded = tbl.seeds .== 0
    for w in weeks
        rc, gc = tbl[!, "rki_cases_w$w"], tbl[!, "gems_cases_w$w"]
        ri, gi = tbl[!, "rki_incidence_w$w"], tbl[!, "gems_incidence_w$w"]
        lr = log.(rc .+ 1)
        push!(res, (
            week = w,
            spearman_incidence = corspearman(ri, gi),
            r2_log_cases = cor(lr, log.(gc .+ 1))^2,
            r2_log_cases_population = cor(lr, log.(tbl.population))^2,
            r2_log_cases_imports = has_imp ? cor(lr, log.(tbl.expected_imports))^2 : NaN,
            spearman_incidence_imports = has_imp ? corspearman(ri, tbl.expected_imports ./ tbl.population) : NaN,
            n_unseeded = sum(unseeded),
            spearman_incidence_unseeded = sum(unseeded) >= 5 ? corspearman(ri[unseeded], gi[unseeded]) : NaN))
    end
    return res
end


"""
    plot_rki_weekly_panel(tbl, weeks; cap = 0.95)

RKI (top) and GEMS (bottom) incidence maps per week, one color scale per row
that ends at the `cap` quantile of the row.
"""
function plot_rki_weekly_panel(tbl::DataFrame, weeks; cap::Real = 0.95)
    plts = []
    for (src, lbl) in [("rki", "RKI"), ("gems", "GEMS")]
        cmax = _clim_max(reduce(vcat, [tbl[!, "$(src)_incidence_w$w"] for w in weeks]), cap)
        for w in weeks
            df = DataFrame(ags = tbl.ags, incidence = tbl[!, "$(src)_incidence_w$w"])
            push!(plts, agsmap(df,
                title = "$lbl, onsets week $w",
                clims = (0, cmax),
                titlefontsize = map_titlefontsize(),
                colorbar_titlefontsize = map_legendfontsize(),
                ytickfontsize = map_legendfontsize(),
                linecolor = nothing,
                size = map_dims()))
        end
    end
    return plot(plts...,
        layout = (2, length(weeks)),
        size = (length(weeks) * map_dims()[1], 2 * map_dims()[2]),
        fontfamily = "Times Roman")
end


# daily RKI cases and GEMS symptom onsets over the simulated period
function _daily_cases(res)
    onset = infections(res.pp).symptom_onset
    days = collect(res.sim.startdate:Day(1):res.sim.enddate)
    gems = [sum(onset .== Dates.value(d - res.sim.startdate)) for d in days]
    rki = [sum(res.rki.cases[res.rki.refdate .== d]) for d in days]
    return DataFrame(date = days, rki = rki, gems = gems)
end


"""
    plot_rki_epicurve(daily::DataFrame)

Daily RKI cases and GEMS symptom onsets (columns `date`, `rki`, `gems`).
"""
function plot_rki_epicurve(daily::DataFrame)
    return plot(daily.date, [daily.rki daily.gems],
        label = ["RKI (onset or report)" "GEMS (symptom onset)"],
        linewidth = 3,
        ylabel = "daily cases",
        legend = :topleft,
        size = (1600, 800),
        fontfamily = "Times Roman")
end


"""
    plot_rki_epicurves(runs::Vector{<:Pair})

Daily RKI cases and GEMS symptom onsets of several runs (`label => daily`,
see `plot_rki_epicurve`). RKI is taken from the first run.
"""
function plot_rki_epicurves(runs::Vector{<:Pair})
    daily = last(first(runs))
    p = plot(daily.date, daily.rki,
        label = "RKI (onset or report)",
        linewidth = 3,
        ylabel = "daily cases",
        legend = :topleft,
        size = (1600, 800),
        fontfamily = "Times Roman")
    for (lbl, d) in runs
        plot!(p, d.date, d.gems, label = "GEMS ($lbl)", linewidth = 3)
    end
    return p
end


"""
    scenario_rki_imports(;
        configfile = "models/DE_model_RKI_2020.toml",
        population = "DE",
        weeks = 10:12,
        decay = 1.0,
        import_scale = 1.0,
        heinsberg = true,
        transmission_rate = nothing,
        incubation = 7,
        seeds = [0],
        folder = nothing)

Runs the RKI comparison with imports and all analyses. GEMS results are
averaged over one run per seed in `seeds`.
"""
function scenario_rki_imports(;
    configfile = "models/DE_model_RKI_2020.toml",
    population = "DE",
    weeks = 10:12,
    decay::Real = 1.0,
    import_scale::Real = 1.0,
    heinsberg::Bool = true,
    transmission_rate = nothing,
    incubation::Int = 7,
    seeds = [0],
    folder = nothing)

    isnothing(folder) ? folder = result_folder() : nothing

    rki = load_rki_cases()
    tbls = DataFrame[]
    schedules = DataFrame[]
    dailies = DataFrame[]
    curve = DataFrame()
    for s in seeds
        res = run_rki_import_sim(configfile = configfile, population = population, rki = rki,
            decay = decay, import_scale = import_scale, heinsberg = heinsberg,
            transmission_rate = transmission_rate, incubation = incubation, seed = s)
        push!(tbls, rki_weekly_table(res, weeks))
        push!(schedules, insertcols(res.schedule, 1, :seed => s))
        push!(dailies, _daily_cases(res))
        curve = res.curve
        res = nothing
        GC.gc()
    end

    tbl = _average_runs(tbls)
    daily = _average_runs(dailies)
    metrics = rki_weekly_metrics(tbl, weeks)
    println(metrics)

    CSV.write(joinpath(folder, "sc_rki_imports_counties.csv"),
        transform(tbl, :ags => ByRow(a -> a.id) => :ags))
    CSV.write(joinpath(folder, "sc_rki_imports_metrics.csv"), metrics)
    CSV.write(joinpath(folder, "sc_rki_imports_seeding.csv"),
        transform(vcat(schedules...), :ags => ByRow(a -> a.id) => :ags))
    CSV.write(joinpath(folder, "sc_rki_imports_curve.csv"), curve)
    CSV.write(joinpath(folder, "sc_rki_imports_daily.csv"), daily)
    png(plot_rki_weekly_panel(tbl, weeks), joinpath(folder, "sc_rki_imports_panel.png"))
    png(plot_rki_epicurve(daily), joinpath(folder, "sc_rki_imports_epicurve.png"))

    return (tbl = tbl, metrics = metrics, daily = daily, folder = folder)
end


"""
    rki_national_weekly(res, weeks)

Weekly RKI cases and GEMS symptom onsets (summed over all counties).
"""
function rki_national_weekly(res, weeks)
    onset = infections(res.pp).symptom_onset
    rki = [sum(res.rki.cases[_monday_2020(w) .<= res.rki.refdate .<= _monday_2020(w) + Day(6)]) for w in weeks]
    gems = map(weeks) do w
        t0 = Dates.value(_monday_2020(w) - res.sim.startdate)
        sum((onset .>= t0) .& (onset .<= t0 + 6))
    end
    return DataFrame(week = collect(weeks), rki = rki, gems = gems)
end


"""
    fit_rki_transmission_rate(;
        rates,
        weeks = 10:11,
        seeds = [0],
        folder = nothing,
        kwargs...)

Grid search for the transmission rate. Objective: squared difference of
log weekly cases (GEMS vs. RKI, national), averaged over seeds.
Other keyword arguments are passed to `run_rki_import_sim`.
"""
function fit_rki_transmission_rate(;
    rates::AbstractVector{<:Real},
    weeks = 10:11,
    seeds = [0],
    rki::DataFrame = load_rki_cases(),
    folder = nothing,
    kwargs...)

    isnothing(folder) ? folder = result_folder() : nothing
    fits = DataFrame()

    for tr in rates, s in seeds
        res = run_rki_import_sim(; rki = rki, transmission_rate = tr, seed = s, kwargs...)
        nw = rki_national_weekly(res, weeks)
        obj = sum((log.(nw.gems .+ 1) .- log.(nw.rki .+ 1)) .^ 2)
        for r in eachrow(nw)
            push!(fits, (transmission_rate = tr, seed = s, week = r.week, rki = r.rki, gems = r.gems, objective = obj))
        end
        # write after every run to keep partial results
        CSV.write(joinpath(folder, "sc_rki_fit.csv"), fits)

        res = nothing
        GC.gc()
    end

    summary = unique(fits[:, [:transmission_rate, :seed, :objective]]) |>
        d -> groupby(d, :transmission_rate) |>
        d -> combine(d, :objective => mean => :objective)
    best = summary.transmission_rate[argmin(summary.objective)]
    printinfo("Best transmission rate: $best")
    if best == minimum(rates) || best == maximum(rates)
        @warn "Best transmission rate ($best) is at the edge of the grid."
    end

    return (fits = fits, summary = summary, best = best, folder = folder)
end
