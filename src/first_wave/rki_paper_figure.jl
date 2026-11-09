###
### FIRST WAVE 2020: PAPER FIGURE
###

export plot_paper_figure_1


# sizes (px): map (aspect of the bounds below), horizontal gap between the maps (margin
# on each side), width of the row labels, of the colorbar (bar only) and of the empty
# cell right of it for its ticks and label (drawn outside of the bar, so the curve can
# span to the right edge below), the colorbar's share of the height of the maps, height
# of the panel titles, of the week headers, and of the curve
const _FIG1_MAP = (1000, 1350)
const _FIG1_GAP = 10mm
const _FIG1_LABEL_WIDTH = 360
const _FIG1_COLORBAR_WIDTH = 60
const _FIG1_COLORBAR_LABELS = 340
const _FIG1_COLORBAR_SHARE = 0.7
const _FIG1_TITLE_HEIGHT = 190
const _FIG1_HEADER_HEIGHT = 200
const _FIG1_CURVE_HEIGHT = 1300

# fonts as in the other figures of the paper (Times; about 0.009 pt per px of the figure
# width for titles, 0.0075 for ticks): panel letters ("A)", bold, top left), panel
# titles (centered), headers (weeks and rows), and everything else
const _FIG1_FONT = "Times Roman"
const _FIG1_FONT_BOLD = "Times Bold"
const _FIG1_LETTERFONTSIZE = 66
const _FIG1_TITLEFONTSIZE = 52
const _FIG1_HEADERFONTSIZE = 50
const _FIG1_FONTSIZE = 44
const _FIG1_BLACK = RGB(0, 0, 0)
const _FIG1_GREY = RGB(0.35, 0.35, 0.35)
const _FIG1_LIGHTGREY = RGB(0.55, 0.55, 0.55)
const _FIG1_GEMS = RGB(0.24, 0.64, 0.31)
const _FIG1_RKI = RGB(0.45, 0.45, 0.45)
const _FIG1_RKI_DAILY = RGB(0.7, 0.7, 0.7)

# bounds of the county shapes of GEMS (`germanshapes(2)`), to draw the maps without margin
const _FIG1_XLIMS = (3.2804e6, 3.9216e6)
const _FIG1_YLIMS = (5.2375e6, 6.1035e6)

# a text of a text panel: position (x, y from 0 to 1) and style
_fig1_t(s; x = 0.5, y = 0.5, size = _FIG1_FONTSIZE, font = _FIG1_FONT, color = _FIG1_BLACK,
    rotation = 0, halign = :center) = (s, x, y, size, font, color, rotation, halign)

# empty panel with texts (`_fig1_t`)
function _fig1_text(texts...)
    p = plot(framestyle = :none, ticks = false, xlims = (0, 1), ylims = (0, 1), legend = false,
        margin = 0mm)
    for (s, x, y, size, font, color, rotation, halign) in texts
        annotate!(p, x, y, text(s, size, halign, font, color, rotation = rotation))
    end
    return p
end

# panel title as in the other figures of the paper: bold letter ("A)") top left, title
# centered
_fig1_title(letter, title) = _fig1_text(
    _fig1_t("$letter)", x = 0, y = 0.6, size = _FIG1_LETTERFONTSIZE, font = _FIG1_FONT_BOLD,
        halign = :left),
    _fig1_t(title, y = 0.35, size = _FIG1_TITLEFONTSIZE))

_fig1_header(w) = _fig1_text(
    _fig1_t("Week $w", y = 0.68, size = _FIG1_HEADERFONTSIZE, font = _FIG1_FONT_BOLD),
    _fig1_t("$(ddmm(monday_2020(w)))–$(ddmm(sunday_2020(w)))", y = 0.22, color = _FIG1_GREY))

_fig1_rowlabel(s) = _fig1_text(
    _fig1_t(s, size = _FIG1_HEADERFONTSIZE, font = _FIG1_FONT_BOLD, rotation = 90))

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
        widen = false,
        margin = 0mm,
        left_margin = _FIG1_GAP,
        right_margin = _FIG1_GAP)
end

# vertical colorbar (heatmap of the color scale, no ticks)
function _fig1_colorbar(c, cmap)
    v = range(0, c, length = 256)
    return heatmap([0, 1], v, repeat(collect(v), 1, 2),
        clims = (0, c),
        color = cmap,
        colorbar = false,
        ticks = false,
        ylims = (0, c),
        framestyle = :box,
        foreground_color_border = _FIG1_LIGHTGREY,
        margin = 0mm,
        top_margin = 8mm,
        bottom_margin = 8mm,
        left_margin = _FIG1_GAP)
end

# ticks and label of the colorbar, drawn in the empty cell right of it (as their own
# panel: as axis ticks of the colorbar, their padding would widen the whole right
# column of the figure, and the curve below could not span to the right edge)
function _fig1_colorbar_ticks(c, label)
    p = plot(framestyle = :none, ticks = false, xlims = (0, 1), ylims = (0, c), legend = false,
        widen = false, margin = 0mm, top_margin = 8mm, bottom_margin = 8mm)
    ticks = first(Plots.PlotUtils.optimize_ticks(0.0, Float64(c), k_min = 4))
    for t in ticks
        plot!(p, [0, 0.07], [t, t], color = _FIG1_LIGHTGREY, linewidth = 3)
        annotate!(p, 0.11, t, text(isinteger(t) ? string(Int(t)) : string(t), _FIG1_FONTSIZE,
            :left, _FIG1_FONT))
    end
    annotate!(p, 0.82, c / 2, text(label, _FIG1_FONTSIZE, :center, _FIG1_FONT, rotation = 90))
    return p
end

# daily cases of RKI (reference, grey) and GEMS (7-day means; GEMS with the 90% band of
# the calibration runs, `band`, or else the range of the runs), the weeks of the maps
# shaded and tagged below the axis, the starts of the intervention phases, and the end
# of the fitted data; framed legend top right (as in the other figures of the paper)
function _fig1_curve(daily, band, weeks; phase_dates, cutoff, basis)
    days = Dates.value.(daily.date)
    rki7 = rolling7(daily.rki)
    gems7 = rolling7(daily.gems)
    ribbon_src = isnothing(band) ? daily : band
    low, high = rolling7(ribbon_src.gems_low), rolling7(ribbon_src.gems_high)
    ymax = 1.3 * maximum(vcat(daily.rki, high))
    xmin, xmax = first(days), last(days) + 1
    mondays = [d for d in first(daily.date):Day(1):last(daily.date) if Dates.dayofweek(d) == 1]
    rki_label = basis == :onset ? "RKI (known onset)" : "RKI"

    p = plot(
        ylabel = "Cases per day (7-day mean)",
        yguidefontsize = _FIG1_HEADERFONTSIZE,
        tickfontsize = _FIG1_FONTSIZE,
        legendfontsize = _FIG1_FONTSIZE,
        legend = :topright,
        foreground_color_legend = _FIG1_LIGHTGREY,
        background_color_legend = RGB(1, 1, 1),
        grid = true,
        gridalpha = 0.15,
        gridlinewidth = 2,
        framestyle = :axes,
        tick_direction = :in,
        widen = false,
        xlims = (xmin, xmax),
        ylims = (0, ymax),
        xticks = _date_ticks(mondays),
        margin = 0mm,
        top_margin = 26mm,
        bottom_margin = 14mm)
    plot!(p, days, low,
        fillrange = high,
        color = _FIG1_GEMS,
        fillalpha = 0.25,
        linewidth = 0,
        label = isnothing(band) ? "GEMS, range of the runs" : "GEMS, 90% calibration band")
    plot!(p, days, gems7,
        color = _FIG1_GEMS,
        linewidth = 8,
        label = "GEMS, 7-day mean")
    scatter!(p, days, daily.rki,
        label = "$rki_label, daily",
        color = _FIG1_RKI_DAILY,
        alpha = 0.6,
        markersize = 11,
        markerstrokewidth = 0)
    plot!(p, days, rki7,
        color = _FIG1_RKI,
        linewidth = 5,
        label = "$rki_label, 7-day mean")

    # weeks of the maps: shaded, tagged above the plot (like the week headers of the maps)
    for (k, w) in enumerate(weeks)
        span = Dates.value.([monday_2020(w), sunday_2020(w) + Day(1)])
        vspan!(p, span,
            color = _FIG1_BLACK,
            alpha = 0.07,
            linewidth = 0,
            label = k == 1 ? "weeks of the maps" : "")
        annotate!(p, sum(span) / 2, 1.045 * ymax,
            text("W$w", _FIG1_FONTSIZE, :center, _FIG1_FONT_BOLD, _FIG1_GREY))
    end

    # starts of the intervention phases: dotted lines that end in a dot with their label
    # right of it, each label a bit lower than the one before (the lines are only a week
    # apart)
    for (k, d) in enumerate(phase_dates)
        x, y = Dates.value(d), (0.955 - 0.065 * (k - 1)) * ymax
        plot!(p, [x, x], [0, y],
            linestyle = :dot,
            linewidth = 4,
            color = _FIG1_GREY,
            label = "")
        scatter!(p, [x], [y],
            color = _FIG1_GREY,
            markersize = 10,
            markerstrokewidth = 0,
            label = "")
        annotate!(p, x + 0.6, y,
            text("intervention phase $k", _FIG1_FONTSIZE, :left, _FIG1_FONT, _FIG1_GREY))
    end

    # end of the fitted data (fitted up to and including `cutoff`: the line is at the end
    # of that day), ends in a dot between its labels
    c = Dates.value(cutoff) + 1
    plot!(p, [c, c], [0, 0.955 * ymax],
        linestyle = :dash,
        linewidth = 4,
        color = _FIG1_GREY,
        label = "")
    scatter!(p, [c], [0.955 * ymax],
        color = _FIG1_GREY,
        markersize = 10,
        markerstrokewidth = 0,
        label = "")
    annotate!(p, c - 0.8, 0.955 * ymax,
        text("← fitted until $(ddmm(cutoff))", _FIG1_FONTSIZE, :right, _FIG1_FONT, _FIG1_GREY))
    annotate!(p, c + 0.8, 0.955 * ymax,
        text("forecast →", _FIG1_FONTSIZE, :left, _FIG1_FONT, _FIG1_GREY))
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

Figure 1 of the paper.

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
    # titles with the year (and the months of the curve)
    from, to = first(daily.date), last(daily.date)
    year = Dates.year(from)
    title_a = (onset ? "Weekly cases by symptom onset per 100,000 by county" :
        "Weekly reported cases per 100,000 by county") * ", $year"
    title_b = (onset ? "Cases by symptom onset in Germany" : "Reported cases in Germany") *
        ", $(Dates.monthname(from))–$(Dates.monthname(to)) $year"

    # rows: title A | label column, maps (week headers, RKI, GEMS), colorbar (with an
    # empty cell for its ticks and label) | title B | curve (starts below the first map,
    # its tick labels reach into the label column, the cell left of it is almost empty;
    # ends at the right edge)
    mw, mh = _FIG1_MAP
    lw, cw, cl = _FIG1_LABEL_WIDTH, _FIG1_COLORBAR_WIDTH, _FIG1_COLORBAR_LABELS
    th, hh, ch = _FIG1_TITLE_HEIGHT, _FIG1_HEADER_HEIGHT, _FIG1_CURVE_HEIGHT
    width = lw + n * mw + cw + cl
    heights = [th, hh + 2 * mh, th, ch]
    block = [hh, mh, mh] ./ (hh + 2 * mh)
    rest = (1 - _FIG1_COLORBAR_SHARE) * mh
    l = grid(4, 1, heights = heights ./ sum(heights))
    l[2, 1] = grid(1, 3, widths = [lw, n * mw, cw + cl] ./ width)
    l[2, 1][1, 1] = grid(3, 1, heights = block)
    l[2, 1][1, 2] = grid(3, 1, heights = block)
    for r in 1:3
        l[2, 1][1, 2][r, 1] = grid(1, n)
    end
    l[2, 1][1, 3] = grid(3, 1, heights = [hh + rest, 2 * mh - 2 * rest, rest] ./ (hh + 2 * mh))
    l[2, 1][1, 3][2, 1] = grid(1, 2, widths = [cw, cl] ./ (cw + cl))
    l[4, 1] = grid(1, 2, widths = [0.001, 0.999])

    p = plot(
        _fig1_title("A", title_a),
        _fig1_text(), _fig1_rowlabel("Robert Koch Institute"), _fig1_rowlabel("GEMS"),
        _fig1_header.(weeks)...,
        [_fig1_map(tbl, county_column("rki", "incidence", w), c, cmap) for w in weeks]...,
        [_fig1_map(tbl, county_column("gems", "incidence", w), c, cmap) for w in weeks]...,
        _fig1_text(), _fig1_colorbar(c, cmap), _fig1_colorbar_ticks(c, "cases per 100,000"),
        _fig1_text(),
        _fig1_title("B", title_b),
        _fig1_text(), _fig1_curve(daily, band, weeks,
            phase_dates = phase_dates,
            cutoff = cutoff,
            basis = basis),
        layout = l,
        size = (width, sum(heights)),
        fontfamily = _FIG1_FONT)

    sp = p.subplots
    title, label, header = sp[1], sp[3], sp[5]
    first_map = sp[n + 5]
    cb, cblabels, curve = sp[3n + 6], sp[3n + 7], sp[3n + 11]
    w(s) = Measures.width(s.plotarea).value
    h(s) = Measures.height(s.plotarea).value
    mm2px(x) = x * Plots.DPI / 25.4

    # the widths and heights of a layout are shares of the plot areas, the paddings
    # (gaps between the maps, ticks and labels) come on top: enlarge the plot by the
    # paddings so that the maps get their size
    Plots.prepare_output(p)
    used_w = mm2px(w(label) + n * w(first_map) + w(cb) + w(cblabels))
    used_h = mm2px(2 * h(title) + h(header) + 2 * h(first_map) + h(curve))
    plot!(p, size = round.(Int, (2 * width - used_w, 2 * sum(heights) - used_h)))

    # align the left edge of the curve with the first map: move it with its margin (its
    # tick labels reach into the label column; repeat, the move changes the plot areas)
    for _ in 1:3
        Plots.prepare_output(p)
        offset = Plots.left(curve.plotarea) - Plots.left(first_map.plotarea)
        abs(offset.value) < 0.1 && break
        curve[:left_margin] = curve[:left_margin] - offset
    end
    return p
end
