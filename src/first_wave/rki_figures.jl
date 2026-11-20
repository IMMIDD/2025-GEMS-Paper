###
### FIRST WAVE 2020: PLOTS
###

export plot_rki_projection
export plot_rki_phase_projections
export plot_rki_counterfactuals
export plot_rki_imports
export plot_rki_seeding
export plot_rki_maps
export plot_rki_scatter
export plot_rki_state_weekly
export plot_rki_state_daily
export plot_rki_state_scatter
export plot_rki_state_maps


# axis and legend labels per basis (`:report` or `:onset`)
function _basis_labels(basis)
    if basis == :onset
        return (
            ylabel = "cases with known symptom onset per day (7-day mean)",
            rki = "RKI (onset, cases with known onset)",
            gems = "GEMS (detected x RKI share with known onset)",
            name = "onset")
    end
    return (
        ylabel = "reported cases per day (7-day mean)",
        rki = "RKI (reporting date)",
        gems = "GEMS",
        name = "reported")
end

# centered 7-day mean (shorter windows at the edges)
rolling7(x) = [mean(x[max(1, i - 3):min(length(x), i + 3)]) for i in eachindex(x)]

# upper color limit: `cap` quantile of the values (outliers like Heinsberg saturate)
function _clim_max(x, cap)
    q = quantile(x, cap)
    return q > 0 ? q : max(maximum(x), 1.0)
end

# x ticks at the given dates, formatted as "dd.mm."
_date_ticks(dates) = (Dates.value.(dates), ddmm.(dates))

# first and last date and label of a week or of a window `name => weeks`
function _week_span(w)
    if w isa Pair
        weeks = last(w)
        return (key = first(w), from = monday_2020(first(weeks)), to = sunday_2020(last(weeks)),
            label = "cumulative, weeks $(first(weeks))-$(last(weeks))")
    end
    return (key = w, from = monday_2020(w), to = sunday_2020(w), label = "week $w")
end

# runs of different length (a projection that ends earlier): all on the dates of the longest run,
# missing days as NaN (not drawn), RKI from the longest run
function _pad_runs(runs)
    ref = last(runs[argmax([nrow(d) for (_, d) in runs])])
    padded = Pair{String, DataFrame}[]
    for (label, d) in runs
        if nrow(d) == nrow(ref)
            push!(padded, label => d)
            continue
        end
        row = Dict(zip(d.date, eachrow(d)))
        p = DataFrame(date = ref.date, rki = ref.rki)
        for c in (:gems, :gems_low, :gems_high)
            p[!, c] = [haskey(row, x) ? Float64(row[x][c]) : NaN for x in ref.date]
        end
        push!(padded, label => p)
    end
    return padded
end



###
### NATIONAL CURVES
###

"""
    plot_rki_projection(runs;
        phase_dates,
        cutoff,
        basis = :report,
        fit_label = "end of fitted data (forecast after)")

Plots the daily cases of RKI and of GEMS runs (`label => daily`, columns
`date`, `rki`, `gems`, `gems_low`, `gems_high`) as 7-day means, with the range
of the runs as a band, by reporting date or symptom onset (`basis`). Vertical
lines at the phase starts and at the end of the fitted data (`cutoff`).
"""
function plot_rki_projection(runs::Vector{<:Pair};
    phase_dates,
    cutoff::Date,
    basis = :report,
    fit_label = "end of fitted data (forecast after)")

    lab = _basis_labels(basis)
    daily = last(first(runs))
    p = scatter(daily.date, daily.rki,
        label = "RKI, daily",
        color = :gray,
        alpha = 0.4,
        markersize = 3,
        markerstrokewidth = 0,
        ylabel = lab.ylabel,
        legend = :topright,
        size = (1600, 800),
        left_margin = 10mm,
        bottom_margin = 10mm,
        fontfamily = "Times Roman")
    plot!(p, daily.date, rolling7(daily.rki),
        color = :black,
        linewidth = 3,
        label = lab.rki)
    for (label, d) in runs
        m = rolling7(d.gems)
        plot!(p, d.date, m,
            ribbon = (m .- rolling7(d.gems_low), rolling7(d.gems_high) .- m),
            fillalpha = 0.25,
            linewidth = 3,
            label = label)
    end
    vline!(p, phase_dates,
        linestyle = :dot,
        color = :black,
        label = "phase start")
    vline!(p, [cutoff],
        linestyle = :dash,
        color = :black,
        label = fit_label)
    # x ticks: first of each month, the phase starts, and the end of the fitted data
    months = [d for d in first(daily.date):Day(1):last(daily.date) if Dates.day(d) == 1]
    plot!(p,
        xticks = _date_ticks(sort(unique(vcat(months, phase_dates, cutoff)))),
        xrotation = 45)
    return p
end


"""
    plot_rki_phase_projections(runs; phase_dates, cutoff, basis = :report, top = 1.2)

Plots projections (`label => daily`, see `plot_rki_projection`) that may end on
different days, on a linear scale up to `top` x the RKI peak.
"""
function plot_rki_phase_projections(runs::Vector{<:Pair};
    phase_dates,
    cutoff::Date,
    basis = :report,
    top::Real = 1.2)

    padded = _pad_runs(runs)
    p = plot_rki_projection(padded,
        phase_dates = phase_dates,
        cutoff = cutoff,
        basis = basis)
    plot!(p, ylims = (0, top * maximum(rolling7(last(first(padded)).rki))))
    return p
end


"""
    plot_rki_counterfactuals(runs; phase_dates, cutoff, basis = :report, top = 3.0)

Plots counterfactual runs (`label => daily`, the fitted run first) against
RKI, on a linear scale up to `top` x the RKI peak and on a log scale.
"""
function plot_rki_counterfactuals(runs::Vector{<:Pair};
    phase_dates,
    cutoff::Date,
    basis = :report,
    top::Real = 3.0)

    lin = plot_rki_projection(runs,
        phase_dates = phase_dates,
        cutoff = cutoff,
        basis = basis,
        fit_label = "end of fitted data")
    plot!(lin,
        ylims = (0, top * maximum(rolling7(last(first(runs)).rki))),
        legend = false)

    # values below 1 are shown as 1
    columns = [:rki, :gems, :gems_low, :gems_high]
    clamp1(d) = transform(d, columns .=> (x -> max.(x, 1.0)), renamecols = false)
    lg = plot_rki_projection([label => clamp1(d) for (label, d) in runs],
        phase_dates = phase_dates,
        cutoff = cutoff,
        basis = basis,
        fit_label = "end of fitted data")
    plot!(lg,
        yscale = :log10,
        ylims = (1, Inf),
        legend = :bottomright)
    return plot(lin, lg,
        layout = (1, 2),
        size = (2400, 800),
        margin = 10mm)
end


"""
    plot_rki_imports(daily, seeded; phase_dates, cutoff, basis = :report)

National curve of GEMS (`daily`, see `plot_rki_projection`) against RKI, with
the seeded infections per day by source below (`seeded`, columns `date`,
`austria`, `italy`, `other`, `heinsberg`, from `rki_seeded_infections`).
"""
function plot_rki_imports(daily::DataFrame, seeded::DataFrame;
    phase_dates,
    cutoff::Date,
    basis = :report)

    top = plot_rki_projection([_basis_labels(basis).gems => daily],
        phase_dates = phase_dates,
        cutoff = cutoff,
        basis = basis)
    xl = Dates.value.(extrema(daily.date))
    plot!(top, xlims = xl)

    stacked = cumsum(hcat(seeded.heinsberg, seeded.other, seeded.italy, seeded.austria), dims = 2)
    labels = ["Heinsberg (RKI cases, 1:1)", "other countries", "Italy", "Austria"]
    colors = [:gray, :steelblue, :darkgreen, :firebrick]
    bottom = plot(
        ylabel = "seeded infections per day",
        legend = :topright,
        xlims = xl,
        fontfamily = "Times Roman")
    for k in reverse(eachindex(labels))
        plot!(bottom, seeded.date, stacked[:, k],
            fillrange = 0,
            fillalpha = 0.85,
            linewidth = 0,
            color = colors[k],
            label = labels[k])
    end
    vline!(bottom, phase_dates,
        linestyle = :dot,
        color = :black,
        label = "")
    plot!(bottom,
        xticks = xticks(top)[1],
        xrotation = 45)
    return plot(top, bottom,
        layout = grid(2, 1, heights = [0.72, 0.28]),
        size = (1600, 1100),
        left_margin = 10mm)
end


"""
    plot_rki_seeding(daily, seeded, infections; phase_dates, cutoff, basis = :report)

Seeding figure: left `plot_rki_imports`, right the infections per day in GEMS
(`infections`, columns `date`, `infections`, `seeded`, from `run_rki_phase_sims`)
as stacked areas of seeded infections and infections by transmission in the model.
"""
function plot_rki_seeding(daily::DataFrame, seeded::DataFrame, infections::DataFrame;
    phase_dates,
    cutoff::Date,
    basis = :report)

    left = plot_rki_imports(daily, seeded,
        phase_dates = phase_dates,
        cutoff = cutoff,
        basis = basis)

    right = plot(infections.date, infections.infections,
        fillrange = 0,
        fillalpha = 0.85,
        linewidth = 0,
        color = :darkorange,
        label = "by transmission in GEMS",
        ylabel = "infections per day",
        legend = :topright,
        size = (1600, 1100),
        left_margin = 10mm,
        bottom_margin = 10mm,
        fontfamily = "Times Roman")
    plot!(right, infections.date, infections.seeded,
        fillrange = 0,
        fillalpha = 0.85,
        linewidth = 0,
        color = :gray,
        label = "seeded (imports and Heinsberg)")
    vline!(right, phase_dates,
        linestyle = :dot,
        color = :black,
        label = "phase start")
    vline!(right, [cutoff],
        linestyle = :dash,
        color = :black,
        label = "end of fitted data")
    share = round(100 * sum(infections.seeded) / sum(infections.infections), digits = 1)
    plot!(right,
        title = "seeded: $share % of all infections",
        xticks = xticks(left[1]),
        xrotation = 45)
    return plot(left, right,
        layout = (1, 2),
        size = (3200, 1100))
end



###
### COUNTIES
###

"""
    plot_rki_maps(tbl, weeks; basis = :report, cap = 0.95)

Maps of the RKI (top) and GEMS (bottom) incidence per county (`tbl`, from
`rki_weekly_table`) for `weeks`, with one color scale up to the `cap` quantile
of all values.
"""
function plot_rki_maps(tbl::DataFrame, weeks; basis = :report, cap::Real = 0.95)
    name = _basis_labels(basis).name
    columns = [county_column(s, "incidence", w) for s in ["rki", "gems"] for w in weeks]
    c = _clim_max(reduce(vcat, [tbl[!, col] for col in columns]), cap)
    plts = []
    for (src, label) in [("rki", "RKI"), ("gems", "GEMS")], w in weeks
        incidence = DataFrame(ags = tbl.ags, incidence = tbl[!, county_column(src, "incidence", w)])
        push!(plts, agsmap(incidence,
            title = "$label, $name week $w ($(ddmm(monday_2020(w))) to $(ddmm(sunday_2020(w))))",
            clims = (0, c),
            titlefontsize = map_titlefontsize(),
            colorbar_titlefontsize = map_legendfontsize(),
            ytickfontsize = map_legendfontsize(),
            linecolor = nothing,
            size = map_dims()))
    end
    return plot(plts...,
        layout = (2, length(weeks)),
        size = (length(weeks) * map_dims()[1], 2 * map_dims()[2]),
        fontfamily = "Times Roman")
end


"""
    plot_rki_scatter(tbl, weeks; basis = :report, highlight = HEINSBERG => "Heinsberg")

Scatter plots of the RKI and GEMS incidence per county (`tbl`) for each entry
of `weeks` (a week, or `name => weeks` for a window of `rki_cumulative_table`),
with the county of `highlight` (`ags => label`) highlighted.
"""
function plot_rki_scatter(tbl::DataFrame, weeks;
    basis = :report,
    highlight = HEINSBERG => "Heinsberg")

    name = _basis_labels(basis).name
    plts = map(weeks) do w
        span = _week_span(w)
        r = tbl[!, county_column("rki", "incidence", span.key)]
        g = tbl[!, county_column("gems", "incidence", span.key)]
        # log scale: values below 0.5 are shown as 0.5
        x, y = max.(r, 0.5), max.(g, 0.5)
        lim = (0.5, 1.2 * max(maximum(x), maximum(y)))
        p = scatter(x, y,
            xscale = :log10,
            yscale = :log10,
            xlims = lim,
            ylims = lim,
            markersize = 3,
            alpha = 0.5,
            markerstrokewidth = 0,
            color = :black,
            label = "",
            title = "$name $(span.label) ($(ddmm(span.from)) to $(ddmm(span.to))), " *
                "Spearman $(round(corspearman(r, g), digits = 2))",
            xlabel = "RKI per 100,000",
            ylabel = "GEMS per 100,000",
            titlefontsize = 10)
        plot!(p, [lim...], [lim...],
            color = :gray,
            linestyle = :dash,
            label = "equal")
        h = tbl.ags .== first(highlight)
        if any(h)
            scatter!(p, x[h], y[h],
                color = 2,
                markersize = 6,
                markerstrokewidth = 0,
                label = last(highlight))
        end
        p
    end
    return plot(plts...,
        layout = (1, length(weeks)),
        size = (520 * length(weeks), 540),
        fontfamily = "Times Roman",
        left_margin = 12mm,
        bottom_margin = 12mm,
        margin = 6mm)
end



###
### STATES
###

"""
    plot_rki_state_weekly(st, weeks; cutoff, basis = :report, states = RKI_STATES)

Weekly incidence per state (`st`, from `rki_state_table`), one panel per state
with the same y scale, and a line at the end of the fitted data (`cutoff`).
"""
function plot_rki_state_weekly(st::DataFrame, weeks;
    cutoff::Date,
    basis = :report,
    states = RKI_STATES)

    sources = ["rki", "gems"]
    incidences = [state_row_incidence(r, s, w) for r in eachrow(st), s in sources, w in weeks]
    ymax = 1.05 * maximum(incidences)
    mondays = monday_2020.(weeks)
    panels = map(enumerate(eachrow(st))) do (k, r)
        p = plot(mondays, [state_row_incidence(r, "rki", w) for w in weeks],
            color = :black,
            linewidth = 2,
            label = "RKI",
            title = states[r.state],
            legend = k == 1 ? :topright : false,
            titlefontsize = 11,
            ylims = (0, ymax))
        plot!(p, mondays, [state_row_incidence(r, "gems", w) for w in weeks],
            color = :darkorange,
            linewidth = 2,
            label = "GEMS")
        vline!(p, [cutoff],
            linestyle = :dash,
            color = :gray,
            label = "")
        plot!(p,
            xticks = _date_ticks(monday_2020.(first(weeks):2:last(weeks))),
            xrotation = 45)
    end
    return plot(panels...,
        layout = length(panels),
        size = (1600, 1300),
        fontfamily = "Times Roman",
        plot_title = "weekly $(_basis_labels(basis).name) cases per 100,000 by state " *
            "(dashed: end of fitted data)")
end


"""
    plot_rki_state_daily(st, daily;
        first_day,
        cutoff,
        last_day,
        basis = :report,
        states = RKI_STATES)

Daily incidence per state (7-day mean, `daily` from `rki_state_summary`) from
`first_day` to `last_day` with the range of the runs, one panel per state of
`st` (`rki_state_table`, for the populations) with the same y scale.
"""
function plot_rki_state_daily(st::DataFrame, daily::DataFrame;
    first_day::Date,
    cutoff::Date,
    last_day::Date,
    basis = :report,
    states = RKI_STATES)

    per_100k(x, population) = 1e5 .* rolling7(x) ./ population
    curves = map(eachrow(st)) do r
        rows = (daily.state .== r.state) .& (first_day .<= daily.date .<= last_day)
        d = sort(daily[rows, :], :date)
        (state = r.state,
            date = d.date,
            rki = per_100k(d.rki, r.population_official),
            gems = per_100k(d.gems, r.population),
            low = per_100k(d.gems_low, r.population),
            high = per_100k(d.gems_high, r.population))
    end
    ymax = 1.05 * maximum(max(maximum(c.rki), maximum(c.high)) for c in curves)
    months = [d for d in first_day:Day(1):last_day if Dates.day(d) == 1]
    ticks = sort(unique(vcat(first_day, months, cutoff, last_day)))
    panels = map(enumerate(curves)) do (k, c)
        p = plot(c.date, c.rki,
            color = :black,
            linewidth = 2,
            label = "RKI",
            title = states[c.state],
            titlefontsize = 11,
            legend = k == 1 ? :topright : false,
            ylims = (0, ymax))
        plot!(p, c.date, c.gems,
            ribbon = (c.gems .- c.low, c.high .- c.gems),
            fillalpha = 0.25,
            color = :darkorange,
            linewidth = 2,
            label = "GEMS")
        vline!(p, [cutoff],
            linestyle = :dash,
            color = :gray,
            label = "")
        plot!(p,
            xticks = _date_ticks(ticks),
            xrotation = 45)
    end
    return plot(panels...,
        layout = length(panels),
        size = (1600, 1300),
        fontfamily = "Times Roman",
        plot_title = "$(_basis_labels(basis).name) cases per day per 100,000 by state " *
            "(7-day mean; dashed: end of fitted data)")
end


"""
    plot_rki_state_scatter(st, windows; basis = :report, states = RKI_STATES)

Cumulative incidence per state (`st`, from `rki_state_table`), GEMS against
RKI, one panel per window (`name => weeks`).
"""
function plot_rki_state_scatter(st::DataFrame, windows; basis = :report, states = RKI_STATES)
    panels = map(windows) do (_, weeks)
        x = [sum(state_row_incidence(r, "rki", w) for w in weeks) for r in eachrow(st)]
        y = [sum(state_row_incidence(r, "gems", w) for w in weeks) for r in eachrow(st)]
        lim = extrema(vcat(x, y)) .* (0.7, 1.4)
        p = scatter(x, y,
            xscale = :log10,
            yscale = :log10,
            xlims = lim,
            ylims = lim,
            label = "",
            color = :darkorange,
            xlabel = "RKI per 100,000",
            ylabel = "GEMS per 100,000",
            markersize = 6,
            title = "$(_basis_labels(basis).name) cumulative, " *
                "weeks $(first(weeks))-$(last(weeks)), " *
                "Spearman $(round(corspearman(x, y), digits = 2))",
            series_annotations = text.(" " .* states[st.state], 9, :left))
        plot!(p, collect(lim), collect(lim),
            color = :gray,
            linestyle = :dash,
            label = "")
    end
    return plot(panels...,
        layout = (1, length(windows)),
        size = (750 * length(windows), 700),
        fontfamily = "Times Roman",
        margin = 8mm)
end


"""
    plot_rki_state_maps(tbl, weeks; basis = :report)

Maps as `plot_rki_maps`, each county with the incidence of its state.
"""
function plot_rki_state_maps(tbl::DataFrame, weeks; basis = :report)
    maps = copy(tbl)
    for w in weeks, s in ["rki", "gems"]
        population = s == "rki" ? :population_official : :population
        x = combine(groupby(maps, :state),
            county_column(s, "cases", w) => sum => :cases,
            population => sum => :population)
        incidence = Dict(zip(x.state, 1e5 .* x.cases ./ x.population))
        maps[!, county_column(s, "incidence", w)] = [incidence[k] for k in maps.state]
    end
    return plot_rki_maps(maps, weeks, basis = basis)
end
