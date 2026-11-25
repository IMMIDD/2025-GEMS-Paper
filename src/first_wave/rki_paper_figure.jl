###
### FIRST WAVE 2020: PAPER FIGURES
###

# The plots are placed in a layout without padding (`remove_padding!`); titles, labels, legends, and
# colorbars are drawn on the canvas of the figure, into the empty cells of the layout. All sizes in inches.

export plot_paper_figure_1
export plot_paper_figure_2

# colors of GEMS and RKI (reference, grey), and of lines that mark dates
const _FIG1_GEMS = RGB(0.24, 0.64, 0.31)
const _FIG1_RKI = RGB(0.45, 0.45, 0.45)
const _FIG1_RKI_DAILY = RGB(0.7, 0.7, 0.7)
const _FIG1_GREY = RGB(0.35, 0.35, 0.35)

# bounds of the county shapes of GEMS (`germanshapes(2)`), to draw the maps without margin
const _FIG1_XLIMS = (3.2804e6, 3.9216e6)
const _FIG1_YLIMS = (5.2375e6, 6.1035e6)
const _FIG1_ASPECT = (_FIG1_XLIMS[2] - _FIG1_XLIMS[1]) / (_FIG1_YLIMS[2] - _FIG1_YLIMS[1])

# line widths of GEMS, RKI, and the lines that mark dates
const _FIG1_LINEWIDTHS = (gems = 1.4, rki = 0.9, marker = 0.7)

# style of the axes of both figures
_fig1_axes() = (
    title = "",
    ylabel = "",
    legend = false,
    fontfamily = FIG_FONT,
    tickfontsize = FIG_SMALL_FONTSIZE,
    grid = true,
    gridalpha = 0.15,
    gridlinewidth = 0.5,
    framestyle = :axes,
    tick_direction = :in,
    widen = false)

# positions on the canvas of data coordinates of the subplot `k`
function _fig1_position(p, k)
    area = plot_corners(p, k)
    scale(v, lims, from, to) = from + (v - lims[1]) / (lims[2] - lims[1]) * (to - from)
    return (x, y) -> (
        scale(x, Plots.xlims(p[k]), area.left, area.right),
        scale(y, Plots.ylims(p[k]), area.bottom, area.top))
end



###
### FIGURE 1: GERMANY
###

# widths of the columns for the row labels (left) and for the colorbar (right), gap between the maps
const _FIG1_LEFT = 0.56
const _FIG1_RIGHT = 0.5
const _FIG1_GAP = 0.06

# heights of the week headers, gap below the maps, week tags and height of the curve, and
# its tick labels
const _FIG1_HEADER = 0.34
const _FIG1_BELOW = FIG_HEADING_ABOVE_IN
const _FIG1_TAGS = 0.17
const _FIG1_CURVE = 1.75
const _FIG1_TICKS = 0.22

# colorbar: width, distance from the maps, and share of the height of the maps
const _FIG1_COLORBAR = (width = 0.1, gap = 0.1, share = 0.7)

# heights of the labels in the curve (shares of its height): first one, and distance between the phases
const _FIG1_LABELS = (top = 0.95, step = 0.085)

function _fig1_map(tbl, column, c, cmap)
    return agsmap(DataFrame(ags = tbl.ags, incidence = tbl[!, column]),
        clims = (0, c),
        fillcolor = cmap,
        colorbar = false,
        linecolor = nothing,
        legend = false,
        framestyle = :none,
        ticks = false,
        xlims = _FIG1_XLIMS,
        ylims = _FIG1_YLIMS,
        widen = false)
end

# series of the curve and their labels in the legend
_fig1_legend(band, rki) = [
    "band" => isnothing(band) ? "GEMS, range of the runs" : "GEMS, 90% calibration band",
    "gems" => "GEMS, 7-day mean",
    "daily" => "$rki, daily",
    "rki" => "$rki, 7-day mean",
    "weeks" => "Weeks of the maps"]

# daily cases of RKI (reference, grey) and GEMS (7-day means; GEMS with the 90% band of
# the calibration runs, `band`, or else the range of the runs), the weeks of the maps
# shaded, the starts of the intervention phases, and the end of the fitted data (its
# labels: `_fig1_curve_labels!`)
function _fig1_curve(daily, band, weeks; phase_dates, cutoff)
    days = Dates.value.(daily.date)
    ribbon_src = isnothing(band) ? daily : band
    low, high = rolling7(ribbon_src.gems_low), rolling7(ribbon_src.gems_high)
    ymax = 1.4 * maximum(vcat(daily.rki, high))
    mondays = [d for d in first(daily.date):Day(1):last(daily.date) if Dates.dayofweek(d) == 1]

    p = plot(; _fig1_axes()...,
        xlims = (first(days), last(days) + 1),
        ylims = (0, ymax),
        xticks = _date_ticks(mondays))
    plot!(p, days, low,
        fillrange = high,
        color = _FIG1_GEMS,
        fillalpha = 0.25,
        linewidth = 0,
        label = "band")
    plot!(p, days, rolling7(daily.gems),
        color = _FIG1_GEMS,
        linewidth = _FIG1_LINEWIDTHS.gems,
        label = "gems")
    scatter!(p, days, daily.rki,
        color = _FIG1_RKI_DAILY,
        alpha = 0.6,
        markersize = 1.8,
        markerstrokewidth = 0,
        label = "daily")
    plot!(p, days, rolling7(daily.rki),
        color = _FIG1_RKI,
        linewidth = _FIG1_LINEWIDTHS.rki,
        label = "rki")

    for (k, w) in enumerate(weeks)
        vspan!(p, Dates.value.([monday_2020(w), sunday_2020(w) + Day(1)]),
            color = RGB(0, 0, 0),
            alpha = 0.07,
            linewidth = 0,
            label = k == 1 ? "weeks" : "")
    end

    # lines that end in a dot, dotted at the starts of the intervention phases (each a bit lower than the
    # one before, as they are only a week apart) and dashed at the end of the fitted data
    for (x, y, style) in _fig1_marks(phase_dates, cutoff)
        _fig1_mark!(p, x, y * ymax, style)
    end
    return p
end

# dates that are marked in the curve, heights of their lines (shares), and line styles; the data are fitted
# up to and including `cutoff`, so its line is at the end of that day
_fig1_marks(phase_dates, cutoff) = [
    [(Dates.value(d), _FIG1_LABELS.top - _FIG1_LABELS.step * (k - 1), :dot) for (k, d) in enumerate(phase_dates)];
    (Dates.value(cutoff) + 1, _FIG1_LABELS.top, :dash)]

# vertical line at `x` up to `y` that ends in a dot
function _fig1_mark!(p, x, y, linestyle)
    plot!(p, [x, x], [0, y],
        linestyle = linestyle,
        linewidth = _FIG1_LINEWIDTHS.marker,
        color = _FIG1_GREY,
        label = "")
    return scatter!(p, [x], [y],
        color = _FIG1_GREY,
        markersize = 1.5,
        markerstrokewidth = 0,
        label = "")
end

# labels of the curve (subplot `k`) on the canvas: weeks of the maps (above it), intervention phases, and
# the end of the fitted data
function _fig1_curve_labels!(p, canvas, k, weeks; phase_dates, cutoff)
    position = _fig1_position(p, k)
    ymax = last(Plots.ylims(p[k]))
    label!(x, y, s; kwargs...) = canvas_text!(p, canvas, x, y, s; size = FIG_SMALL_FONTSIZE, kwargs...)

    for w in weeks
        x, y = position(sum(Dates.value.([monday_2020(w), sunday_2020(w) + Day(1)])) / 2, ymax)
        label!(x, y + _FIG1_TAGS / 2, "W$w", font = FIG_FONT_BOLD)
    end
    marks = [position(x, y * ymax) for (x, y, _) in _fig1_marks(phase_dates, cutoff)]
    for (i, (x, y)) in enumerate(marks[1:end-1])
        label!(x + 0.05, y, "Intervention phase $i", halign = :left)
    end
    x, y = last(marks)
    label!(x - 0.06, y, "← Fitted until $(ddmm(cutoff))", halign = :right)
    label!(x + 0.06, y, "Forecast →", halign = :left)
    return p
end

"""
    plot_paper_figure_1(tbl, daily, weeks;
        band = nothing,
        phase_dates,
        cutoff,
        basis = :report,
        cap = 0.95,
        colormap = :inferno)

Figure of the first wave in Germany.

A: maps of the RKI (top row) and GEMS (bottom row) incidence per county (`tbl`,
from `rki_weekly_table`) for `weeks` (columns), with one color scale up to the
`cap` quantile of all values (higher values saturate) and a colorbar on the right.

B: daily cases of RKI (reference, grey) and GEMS (`daily`, as in
`plot_rki_projection`) as 7-day means, below the maps and aligned with them,
with the 90% band of the calibration runs (`band`, as `daily`; if `nothing`, the
range of the runs of `daily`), the weeks of the maps shaded and tagged, the
starts of the intervention phases (`phase_dates`), and the end of the fitted
data (`cutoff`).
"""
function plot_paper_figure_1(tbl::DataFrame, daily::DataFrame, weeks;
    band = nothing,
    phase_dates,
    cutoff::Date,
    basis = :report,
    cap::Real = 0.95,
    colormap = :inferno)

    n = length(weeks)
    columns = [county_column(s, "incidence", w) for s in ["rki", "gems"] for w in weeks]
    c = _clim_max(reduce(vcat, [tbl[!, col] for col in columns]), cap)
    cmap = cgrad(colormap)
    onset = basis == :onset
    heading_a = "A) " * (onset ? "Weekly cases by symptom onset" : "Weekly reported cases") * " per 100,000 by county"
    heading_b = "B) " * (onset ? "Cases by symptom onset" : "Reported cases") * " in Germany"

    # rows: (heading A and week headers) | maps RKI | gap | maps GEMS | (heading B and tags) | curve | (ticks);
    # columns: (row labels) | maps with gaps between them | (colorbar)
    middle = FIG_FULL_WIDTH_IN - _FIG1_LEFT - _FIG1_RIGHT
    mw = (middle - (n - 1) * _FIG1_GAP) / n
    mh = mw / _FIG1_ASPECT
    heights = [FIG_HEADING_IN + _FIG1_HEADER, mh, _FIG1_GAP, mh, _FIG1_BELOW + FIG_HEADING_IN + _FIG1_TAGS, _FIG1_CURVE,
        _FIG1_TICKS]
    height = sum(heights)
    right = 0.1 # right of the curve
    l = grid(length(heights), 1, heights = heights ./ height)
    for r in (2, 4)
        l[r, 1] = grid(1, 2n + 1, widths = [_FIG1_LEFT; [mw; repeat([_FIG1_GAP, mw], n - 1)]; _FIG1_RIGHT] ./ FIG_FULL_WIDTH_IN)
    end
    l[6, 1] = grid(1, 3, widths = [_FIG1_LEFT, FIG_FULL_WIDTH_IN - _FIG1_LEFT - right, right] ./ FIG_FULL_WIDTH_IN)

    map_row(s) = reduce(vcat, [[empty_panel(), _fig1_map(tbl, county_column(s, "incidence", w), c, cmap)] for w in weeks])
    p = plot(
        empty_panel(),
        map_row("rki")..., empty_panel(),
        empty_panel(),
        map_row("gems")..., empty_panel(),
        empty_panel(),
        empty_panel(), _fig1_curve(daily, band, weeks, phase_dates = phase_dates, cutoff = cutoff), empty_panel(),
        empty_panel(),
        layout = l,
        size = (100FIG_FULL_WIDTH_IN, 100height),
        fontfamily = FIG_FONT)
    remove_padding!(p)
    canvas = add_canvas!(p)

    # subplots of the maps and of the curve
    rki_maps = 1 .+ 2 .* (1:n)
    gems_maps = rki_maps .+ (2n + 2)
    curve = last(gems_maps) + 4
    corners(k) = plot_corners(p, k)
    mid(area) = ((area.left + area.right) / 2, (area.bottom + area.top) / 2)

    canvas_heading!(p, canvas, height - FIG_HEADING_IN / 2, heading_a)
    for (k, w) in zip(rki_maps, weeks)
        x, top = first(mid(corners(k))), corners(k).top
        canvas_text!(p, canvas, x, top + 0.75 * _FIG1_HEADER, "Week $w", font = FIG_FONT_BOLD)
        canvas_text!(p, canvas, x, top + 0.3 * _FIG1_HEADER, "$(ddmm(monday_2020(w)))–$(ddmm(sunday_2020(w)))",
            size = FIG_SMALL_FONTSIZE)
    end
    for (k, s) in ((first(rki_maps), "Robert Koch Institute"), (first(gems_maps), "GEMS"))
        canvas_text!(p, canvas, _FIG1_LEFT / 2, last(mid(corners(k))), s, font = FIG_FONT_BOLD, rotation = 90)
    end
    bottom, top = corners(last(gems_maps)).bottom, corners(last(rki_maps)).top
    bar = _FIG1_COLORBAR.share * (top - bottom)
    county_colorbar!(p, canvas, CountyMap(column = :incidence, title = "", note = ""), (0.0, Float64(c)),
        corners(last(rki_maps)).right + _FIG1_COLORBAR.gap, (top + bottom - bar) / 2, _FIG1_COLORBAR.width, bar,
        colors = cmap)

    canvas_heading!(p, canvas, corners(curve).top + _FIG1_TAGS + FIG_HEADING_IN / 2, heading_b)
    canvas_text!(p, canvas, 0.08, last(mid(corners(curve))), "Cases per day", rotation = 90)
    _fig1_curve_labels!(p, canvas, curve, weeks, phase_dates = phase_dates, cutoff = cutoff)
    canvas_legend!(p, canvas, curve, _fig1_legend(band, onset ? "RKI (known onset)" : "RKI"))
    return p
end



###
### FIGURE 2: STATES
###

# English names of the states (by state number, as `RKI_STATES`)
const _FIG2_STATES = ["Schleswig-Holstein", "Hamburg", "Lower Saxony", "Bremen",
    "North Rhine-Westphalia", "Hesse", "Rhineland-Palatinate", "Baden-Württemberg", "Bavaria",
    "Saarland", "Berlin", "Brandenburg", "Mecklenburg-Vorpommern", "Saxony",
    "Saxony-Anhalt", "Thuringia"]

# columns and rows of the panels, their height, the gaps between them, the widths left (labels) and
# right of them, and the heights of the heading (with the space below it), the tick labels, and the rows of the key
const _FIG2_GRID = 4
const _FIG2_PANEL = 1.3
const _FIG2_HGAP = 0.12
const _FIG2_VGAP = 0.1
const _FIG2_LEFT = 0.45
const _FIG2_RIGHT = 0.1
const _FIG2_TITLE = FIG_HEADING_IN + 0.05
const _FIG2_TICKS = 0.24
const _FIG2_KEY = 0.17

# texts of a panel: distance from its axis, and from its top (name, population, peak height)
const _FIG2_TEXT = (x = 0.05, name = 0.09, population = 0.24, peak = 0.37)

# locator maps: width, distances from the right and top edge of the panel, colors of the
# state of the panel and of the others
const _FIG2_LOCATOR = (width = 0.42, right = 0.02, top = 0.2)
const _FIG2_LOCATOR_ON = RGB(0.35, 0.35, 0.35)
const _FIG2_LOCATOR_OFF = RGB(0.86, 0.86, 0.86)
# city states (Hamburg, Bremen, Berlin): marked with a dot on the locator map
const _FIG2_CITY_STATES = (2, 4, 11)

# curves of a state (`r`, row of `rki_state_table`) per 100,000 as 7-day means from
# `first_day` to `last_day`: RKI (official population), GEMS (model population) of the
# runs (`states`), and the low and high end of the band (`band`, or else of `states`);
# peak height of the 7-day mean of GEMS relative to RKI (ratio); means and peaks on the
# whole data, only then cut to the shown window
function _fig2_curves(r, states, band; first_day, last_day)
    rows(d) = sort(d[d.state .== r.state, :], :date)
    d = rows(states)
    b = isnothing(band) ? d : rows(band)
    keep(x, dates) = x[first_day .<= dates .<= last_day]
    rki = 1e5 .* rolling7(d.rki) ./ r.population_official
    gems = 1e5 .* rolling7(d.gems) ./ r.population
    per_100k(x, population) = keep(1e5 .* rolling7(x) ./ population, b.date)
    return (
        state = r.state,
        date = keep(d.date, d.date),
        rki = keep(rki, d.date),
        gems = keep(gems, d.date),
        low = per_100k(b.gems_low, r.population),
        high = per_100k(b.gems_high, r.population),
        peak_ratio = maximum(gems) / maximum(rki),
        population = r.population_official)
end

# peak height (PH) of GEMS relative to RKI as text ("PH ×0.73"), explained in the key
function _fig2_peak_text(c)
    n = round(Int, 100 * c.peak_ratio)
    return "PH ×$(n ÷ 100).$(lpad(n % 100, 2, '0'))"
end

# one state panel: GEMS (band and 7-day mean), RKI (reference, grey), end of the fitted
# data; tick labels only on the left column and on the bottom row
function _fig2_panel(c, k; ymax, xticks, yticks, cutoff)
    left, bottom = mod1(k, _FIG2_GRID) == 1, k > _FIG2_GRID * (_FIG2_GRID - 1)
    days = Dates.value.(c.date)
    blank(t) = (t, fill("", length(t)))
    ylabels(t) = [isinteger(y) ? string(Int(y)) : string(y) for y in t]

    p = plot(; _fig1_axes()...,
        xlims = (first(days), last(days) + 1),
        ylims = (0, ymax),
        xticks = bottom ? _date_ticks(xticks) : blank(Dates.value.(xticks)),
        yticks = left ? (yticks, ylabels(yticks)) : blank(yticks))
    plot!(p, days, c.low,
        fillrange = c.high,
        color = _FIG1_GEMS,
        fillalpha = 0.25,
        linewidth = 0,
        label = "")
    plot!(p, days, c.gems,
        color = _FIG1_GEMS,
        linewidth = _FIG1_LINEWIDTHS.gems,
        label = "")
    plot!(p, days, c.rki,
        color = _FIG1_RKI,
        linewidth = _FIG1_LINEWIDTHS.rki,
        label = "")
    # end of the fitted data (fitted up to and including `cutoff`: at the end of that day), up to the height of
    # the population
    _fig1_mark!(p, Dates.value(cutoff) + 1, (1 - _FIG2_TEXT.population / _FIG2_PANEL) * ymax, :dash)
    return p
end

# name, population (official, in millions), and peak height of the state of subplot `k` (top left)
function _fig2_panel_text!(p, canvas, k, c; peaks)
    area = plot_corners(p, k)
    x = area.left + _FIG2_TEXT.x
    label!(dy, s; kwargs...) = canvas_text!(p, canvas, x, area.top - dy, s; halign = :left, kwargs...)
    label!(_FIG2_TEXT.name, _FIG2_STATES[c.state], size = FIG_SMALL_FONTSIZE, font = FIG_FONT_BOLD)
    label!(_FIG2_TEXT.population, "pop. $(round(c.population / 1e6, digits = 1)) M", size = FIG_SMALL_FONTSIZE)
    peaks && label!(_FIG2_TEXT.peak, _fig2_peak_text(c), size = FIG_SMALL_FONTSIZE)
    return p
end

# key below the panels, centered at `x`: the series and the end of the fitted data (first row, at `y`) and the
# abbreviations (second row, bold)
function _fig2_key!(p, canvas, x, y, band_label, rki_label; cutoff, peaks)
    sample, gap, sep = 0.22, 0.05, 0.3
    tw(s; bold = false) = 1.05 * text_width_in(s, FIG_SMALL_FONTSIZE, bold = bold)
    label!(x, y, s; kwargs...) = canvas_text!(p, canvas, x, y, s; size = FIG_SMALL_FONTSIZE, halign = :left, kwargs...)
    starts(widths) = x - (sum(widths) + sep * (length(widths) - 1)) / 2 .+ [0; cumsum(widths[1:end-1] .+ sep)]

    series = [(:band, band_label), (:gems, "GEMS"), (:rki, rki_label), (:cutoff, "Fitted until $(ddmm(cutoff))")]
    for ((kind, s), x0) in zip(series, starts([sample + gap + tw(s) for (_, s) in series]))
        line(; kwargs...) = plot!(p, [x0, x0 + sample], [y, y]; subplot = canvas, label = "", kwargs...)
        if kind == :band
            plot!(p, x0 .+ [0, sample, sample, 0], y .+ 0.045 .* [-1, -1, 1, 1], subplot = canvas, label = "",
                seriestype = :shape, fillcolor = _FIG1_GEMS, fillalpha = 0.25, linealpha = 0)
        elseif kind == :gems
            line(color = _FIG1_GEMS, linewidth = _FIG1_LINEWIDTHS.gems)
        elseif kind == :rki
            line(color = _FIG1_RKI, linewidth = _FIG1_LINEWIDTHS.rki)
        else
            line(color = _FIG1_GREY, linewidth = _FIG1_LINEWIDTHS.marker, linestyle = :dash)
        end
        label!(x0 + sample + gap, y, s)
    end

    abbreviations = [("pop.", "= population in millions")]
    peaks && push!(abbreviations, ("PH", "= peak height, GEMS/RKI"))
    widths = [tw(a, bold = true) + 0.03 + tw(s) for (a, s) in abbreviations]
    for ((a, s), x0) in zip(abbreviations, starts(widths))
        label!(x0, y - _FIG2_KEY, a, font = FIG_FONT_BOLD)
        label!(x0 + tw(a, bold = true) + 0.03, y - _FIG2_KEY, s)
    end
    return p
end

"""
    plot_paper_figure_2(st, states;
        band = nothing,
        runs = nothing,
        first_day = nothing,
        cutoff,
        last_day = nothing,
        basis = :report,
        peaks = true,
        locator = true)

Figure of the first wave in the federal states: daily cases per 100,000 (7-day means) per state from
`first_day` to `last_day` (default: the first and last day of `states`), one
panel per state of `st` (`rki_state_table`, for the populations) with a shared
y scale: RKI (reference, grey) and GEMS (mean of the runs of `states`, from
`rki_state_summary`) with the 90% band of the calibration runs (`band`, as
`states`; if `nothing`, the range of the `runs` runs of `states`, with a
warning), and the end of the fitted data (`cutoff`). With `peaks`, each panel
shows the peak height (PH) of the GEMS 7-day mean relative to the RKI one
(ratio), below the name and the population of the state. With `locator`, each
panel has a small map of Germany on the right with its state highlighted. All
values are computed on the whole data, `first_day` and `last_day` only set the
shown window.
"""
function plot_paper_figure_2(st::DataFrame, states::DataFrame;
    band = nothing,
    runs = nothing,
    first_day::Union{Date, Nothing} = nothing,
    cutoff::Date,
    last_day::Union{Date, Nothing} = nothing,
    basis = :report,
    peaks::Bool = true,
    locator::Bool = true)

    first_day = something(first_day, minimum(states.date))
    last_day = something(last_day, maximum(states.date))
    if isnothing(band)
        @warn "No 90% calibration band per state: the figure shows the range of " *
            "$(isnothing(runs) ? "the" : runs) final runs instead."
    end
    band_label = !isnothing(band) ? "GEMS, 90% calibration band" :
        isnothing(runs) ? "GEMS, range of the runs" : "GEMS, range of $runs runs"
    onset = basis == :onset
    rki_label = onset ? "RKI (known onset)" : "RKI"
    heading = (onset ? "Cases by symptom onset" : "Reported cases") * " per day by federal state"
    ylabel = (onset ? "Cases by symptom onset" : "Reported cases") *
        " per day per 100,000 (7-day mean)"

    curves = [_fig2_curves(r, states, band,
        first_day = first_day,
        last_day = last_day) for r in eachrow(st)]
    # shared scale, with room above the curves for the names
    ymax = 1.3 * maximum(max(maximum(c.rki), maximum(c.high)) for c in curves)
    # (ticks only below the names)
    yticks = filter(y -> y < 0.85ymax, first(Plots.PlotUtils.optimize_ticks(0.0, ymax, k_min = 3, k_max = 4)))
    mondays = [d for d in first_day:Day(1):last_day if Dates.dayofweek(d) == 1]
    # every 4th Monday (a subset of the ticks of the figure of Germany), so the labels of neighboring
    # panels do not run into each other
    xticks = mondays[1:4:end]

    # rows: (heading) | panels with gaps between them | (tick labels and key)
    n = _FIG2_GRID
    heights = [_FIG2_TITLE; [_FIG2_PANEL; repeat([_FIG2_VGAP, _FIG2_PANEL], n - 1)]; _FIG2_TICKS + 2 * _FIG2_KEY]
    height = sum(heights)
    pw = (FIG_FULL_WIDTH_IN - _FIG2_LEFT - _FIG2_RIGHT - (n - 1) * _FIG2_HGAP) / n
    widths = [_FIG2_LEFT; [pw; repeat([_FIG2_HGAP, pw], n - 1)]; _FIG2_RIGHT]
    l = grid(length(heights), 1, heights = heights ./ height)
    for r in 1:n
        l[2r, 1] = grid(1, length(widths), widths = widths ./ FIG_FULL_WIDTH_IN)
    end

    items, panels = Any[empty_panel()], Int[]
    for (k, c) in enumerate(curves)
        push!(items, empty_panel())
        push!(items, _fig2_panel(c, k, ymax = ymax, xticks = xticks, yticks = yticks, cutoff = cutoff))
        push!(panels, length(items))
        mod1(k, n) == n && append!(items, [empty_panel(), empty_panel()]) # end of the row, gap below it
    end
    p = plot(items...,
        layout = l,
        size = (100FIG_FULL_WIDTH_IN, 100height),
        fontfamily = FIG_FONT)
    remove_padding!(p)
    if locator
        shapes = _fig2_state_shapes()
        for (k, c) in zip(panels, curves)
            _fig2_locator!(p, k, c.state, shapes)
        end
    end

    canvas = add_canvas!(p)
    top, bottom = plot_corners(p, first(panels)).top, plot_corners(p, last(panels)).bottom
    canvas_heading!(p, canvas, height - FIG_HEADING_IN / 2, heading)
    canvas_text!(p, canvas, 0.08, (top + bottom) / 2, ylabel, rotation = 90)
    for (k, c) in zip(panels, curves)
        _fig2_panel_text!(p, canvas, k, c, peaks = peaks)
    end
    _fig2_key!(p, canvas, (_FIG2_LEFT + FIG_FULL_WIDTH_IN - _FIG2_RIGHT) / 2, bottom - _FIG2_TICKS - _FIG2_KEY / 2,
        band_label, rki_label, cutoff = cutoff, peaks = peaks)
    return p
end

# land areas of the states (shapes of GEMS without water areas), by state number
function _fig2_state_shapes()
    t = DataFrame(germanshapes(1))
    t = t[coalesce.(t.GF, 0) .== 4, :]
    return [parse(Int, a[1:2]) => g for (a, g) in zip(t.AGS_0, t.geometry)]
end

# locator map of Germany as an inset of subplot `parent` (right, below the name), with `state`
# highlighted (the city states, too small at this size, with a dot in addition)
function _fig2_locator!(p, parent, state, shapes)
    inches(x) = 25.4x * mm
    box = bbox(inches(_FIG2_LOCATOR.right), inches(_FIG2_LOCATOR.top), inches(_FIG2_LOCATOR.width),
        inches(_FIG2_LOCATOR.width / _FIG1_ASPECT), :top, :right)
    n = length(p.subplots) + 1
    for (k, (s, g)) in enumerate(shapes)
        args = k == 1 ? (inset = (parent, box), subplot = n) : (subplot = n,)
        plot!(p, g;
            args...,
            fillcolor = s == state ? _FIG2_LOCATOR_ON : _FIG2_LOCATOR_OFF,
            linecolor = RGB(1, 1, 1),
            linewidth = 0.3,
            label = "",
            framestyle = :none,
            ticks = false,
            legend = false,
            aspect_ratio = :equal,
            xlims = _FIG1_XLIMS,
            ylims = _FIG1_YLIMS,
            widen = false,
            margin = -2mm,
            background_color_subplot = RGBA(1, 1, 1, 0),
            background_color_inside = RGBA(1, 1, 1, 0))
    end
    if state in _FIG2_CITY_STATES
        points = reduce(vcat, [g.points for (s, g) in shapes if s == state])
        scatter!(p, [mean(q.x for q in points)], [mean(q.y for q in points)],
            subplot = n,
            color = _FIG2_LOCATOR_ON,
            markersize = 1.5,
            markerstrokewidth = 0,
            label = "")
    end
    return p
end
