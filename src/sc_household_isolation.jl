export run_household_isolation_sims
export plot_scenario_household_isolation
export plot_r0_by_region
export plot_case_fatality
export plot_combined_readouts

export scenario_household_isolation



"""
    run_household_isolation_sims(;
        configfile = "models/DE_model_first_wave.toml",
        infected_fraction = 0.0001,
        stop_criterion = TimesUp(limit = 180),
        population = "DE")

Runs the household isolation scenario simulations and returns
two ResultData objects: one for the baseline simulation and one for the scenario with household isolation.
And also returns last simulation object for further data extraction.

The baseline uses the default ResultData style, as its readout figure
(`plot_combined_readouts`) needs the raw infections; the scenario uses `"LightRD"`.
"""
function run_household_isolation_sims(;
    configfile = "models/DE_model_first_wave.toml",
    infected_fraction = 0.0001,
    stop_criterion = TimesUp(limit = 180),
    population = "DE")

    # baseline simulation
    sim = Simulation(
        label = "Baseline",
        configfile = configfile,
        population = population,
        seed = 0,
        infected_fraction = infected_fraction,
        stop_criterion = stop_criterion)

    run!(sim)
    rd_b = ResultData(sim)


    # simulation with 14-day household isolation: upon symptom onset, all household members
    # (incl. the symptomatic individual) isolate at home; each new symptomatic member restarts it
    sim = Simulation(
        label = "Scenario",
        configfile = configfile,
        population = population,
        seed = 0,
        infected_fraction = infected_fraction,
        stop_criterion = stop_criterion)

    member_isolation = IStrategy("Isolation", sim)
    add_measure!(member_isolation, SelfIsolation(14))
    member_identification = IStrategy("Find Household Members", sim)
    add_measure!(member_identification,
        FindSettingMembers(Household, member_isolation))
    trigger = SymptomTrigger(member_identification)
    add_symptom_trigger!(sim, trigger)

    run!(sim)
    rd_s = ResultData(sim, style = "LightRD")

    return rd_b, rd_s, sim
end


"""
    plot_scenario_household_isolation(rd_b, rd_s, sim)

Plots the household isolation scenario results.
Requires two ResultData objects as input, one for the baseline and one for the scenario with household isolation,
and the simulation object (for the household sizes and ages per county).
"""
function plot_scenario_household_isolation(rd_b, rd_s, sim)

    # calculate r0 difference
    r0b = r0_per_county(rd_b)
    r0s = r0_per_county(rd_s)
    r0_diff = r0b |>
        df -> leftjoin(df, rename(r0s, :r0 => :r0_new), on = :ags) |>
        df -> transform(df, [:r0, :r0_new] => ByRow(-) => :r0_diff) |>
        df -> DataFrames.select(df, :ags, :r0_diff)

    # mean household sizes per county
    hh_size = households(sim) |>
        vec -> DataFrame(
            ags = county.(ags.(vec)),
            size = size.(vec)) |>
        df -> groupby(df, :ags) |>
        df -> combine(df, :size => mean => :size)

    # mean ages per county
    ind_ages = individuals(sim) |>
        vec -> DataFrame(
            ags = county.(ags.((i -> household(i, sim)).(vec))),
            age = age.(vec)) |>
        df -> groupby(df, :ags) |>
        df -> combine(df, :age => mean => :age)

    # create DataFrame with r0 difference, household size, and mean age
    df = r0_diff |>
        df -> leftjoin(df, hh_size, on = :ags) |>
        df -> leftjoin(df, ind_ages, on = :ags)


    l = @layout [a{0.5w} b{0.5w}]
    l2 = @layout [a{0.75h}
                b{0.25h}]

    return plot(
            plot(
                agsmap(
                    select(df, :ags, :r0_diff),
                    size = 0.5 .* map_dims(),
                    title = "R-Reduction by County",
                    linecolor = nothing,
                ),
                plot(
                    scatter(df.age, df.r0_diff,
                        title = "R-Reduction by Mean Age",
                        xlabel = "Mean Age",
                        ylabel = "R-Reduction",
                        label = "County",
                        markerstrokewidth=0),
                    scatter(df.size, df.r0_diff,
                        title = "R-Reduction by Mean HH-Size",
                        xlabel = "Mean Household Size",
                        ylabel = "R-Reduction",
                        label = "County",
                        markerstrokewidth=0),
                    left_margin = 15mm,
                    bottom_margin = 10mm,
                    layout = (1, 2)
                ),
                layout = l2
            ),
            gemsplot([rd_b, rd_s],
                left_margin = 15mm,
                bottom_margin = 10mm),
            fontfamily = "Times Roman",
            titlefontsize = round(0.7 * map_titlefontsize()),
            tickfontsize = round(0.7 * map_legendfontsize()),
            labelfontsize = round(0.8 * map_legendfontsize()),
            legendfontsize = round(0.6 * map_legendfontsize()),
            layout = l,
            size = 1.5 .* (1800, 1000))
end


###
### GEMS READOUTS OF THE BASELINE (combined figure)
###

# The panels are drawn by GEMS (`gemsplot`, `gemsmap`); only their keyword arguments
# (axes, ticks, fonts, colors) and the layout are set here. Font, panel titles, and
# colorbars in the style of the first wave figures (`_FIG1_*`); width as figure 2.

# sizes (px): figure width, height of a map (its width follows from the aspect of the
# county shapes, `_FIG1_MAP`; two maps and their titles make the figure 7,370 px high,
# a full page of PLOS: 7.5 × 8.75 in), share of the map height that its colorbar
# spans, height of the title rows, width of the y label column of the heatmap (F; its
# tick labels start there; the other rows' y label columns end where their tick labels
# start, see `_ro_combined`), gap between y label and tick labels
const _RO_WIDTH = _FIG2_WIDTH
const _RO_MAP_HEIGHT = 3455
const _RO_COLORBAR_SHARE = 0.7
const _RO_TITLE_HEIGHT = 230
const _RO_YLABEL_WIDTH_HEATMAP = 180
const _RO_YLABEL_GAP = 30

# fonts (pt): about 20 % larger than in figures 1 and 2 (a full page with large
# panels): panel letters, titles, axis labels, ticks and colorbars, legends
const _RO_LETTERFONTSIZE = 80
const _RO_TITLEFONTSIZE = 64
const _RO_GUIDEFONTSIZE = 60
const _RO_FONTSIZE = 54
const _RO_LEGENDFONTSIZE = 46

# legend box width relative to what GR measures for the labels (it measures Times
# wider than it draws it)
const _RO_LEGEND_WIDTH = 0.82

# number as tick label: integers without decimals, else at most two decimals
_ro_num(x) = (r = round(x, digits = 2); isinteger(r) ? string(Int(r)) : string(r))
_ro_millions(x) = _ro_num(x / 1e6)
_ro_percent(x) = _ro_num(100 * x)

# axes of the GEMS plots as in the first wave figures (no titles and y labels: both are
# panels of their own, `_ro_title` and `_ro_ylabel`, so that they line up)
_ro_axes() = (
    title = "",
    plot_title = "",
    ylabel = "",
    fontfamily = _FIG1_FONT,
    tickfontsize = _RO_FONTSIZE,
    guidefontsize = _RO_GUIDEFONTSIZE,
    legendfontsize = _RO_LEGENDFONTSIZE,
    legend_wfactor = _RO_LEGEND_WIDTH,
    extra_kwargs = :subplot,
    foreground_color_legend = _FIG1_LIGHTGREY,
    background_color_legend = RGB(1, 1, 1),
    grid = true,
    gridalpha = 0.15,
    gridlinewidth = 2,
    framestyle = :axes,
    tick_direction = :in,
    margin = 0mm,
    bottom_margin = 18mm)

# map without colorbar, axes, or margin (the colorbar is a panel of its own,
# `_ro_colorbar`), cut to the county shapes as in figure 1
_ro_map() = (
    title = "",
    plot_title = "",
    level = 2,
    fillcolor = :inferno,
    colorbar = false,
    linecolor = nothing,
    legend = false,
    framestyle = :none,
    ticks = false,
    xlims = _FIG1_XLIMS,
    ylims = _FIG1_YLIMS,
    widen = false,
    margin = 0mm,
    left_margin = _FIG1_GAP)

# color range of a plot that GEMS drew (its first subplot)
function _ro_clims(p)
    lo, hi = Plots.update_clims(p[1])
    return isfinite(lo) && isfinite(hi) && lo < hi ? (lo, hi) : (0.0, 1.0)
end

# panel title as `_fig1_title` (bold letter top left, title centered), larger fonts
_ro_title(letter, title) = _fig1_text(
    _fig1_t("$letter)", x = 0, y = 0.6, size = _RO_LETTERFONTSIZE, font = _FIG1_FONT_BOLD,
        halign = :left),
    _fig1_t(title, y = 0.35, size = _RO_TITLEFONTSIZE))

# width (px) of a text in the tick font
_ro_tickwidth(s) = _fig2_text_width(s, _FIG1_FONT, _RO_FONTSIZE)

# y label of a panel of the left column, in a column of its own of `width` px (the
# labels line up, ending `_RO_YLABEL_GAP` before the column's right edge; the row of
# the panel gives it the vertical extent of the panel's plot area)
function _ro_ylabel(s, width)
    half = _RO_GUIDEFONTSIZE / 72 * Plots.DPI / 2
    return _fig1_text(_fig1_t(s, x = 1 - (_RO_YLABEL_GAP + half) / width,
        size = _RO_GUIDEFONTSIZE, rotation = 90))
end

# thin vertical colorbar from `lo` to `hi` (as `_fig1_colorbar`, for any range)
function _ro_colorbar(lo, hi, cmap)
    v = range(lo, hi, length = 256)
    return heatmap([0, 1], v, repeat(collect(v), 1, 2),
        clims = (lo, hi),
        color = cmap,
        colorbar = false,
        ticks = false,
        ylims = (lo, hi),
        framestyle = :box,
        foreground_color_border = _FIG1_LIGHTGREY,
        margin = 0mm,
        left_margin = _FIG1_GAP)
end

# ticks and label of a colorbar in the empty cell right of it (as
# `_fig1_colorbar_ticks`, for any range and tick format). The label is rotated; with
# `sub`, it gets a subscript (e.g. R with "init": GR's Times has no subscript letters),
# centered for a colorbar of `height` px. `labelx`: position of the label (0: at the
# colorbar, 1: right edge of the cell).
function _ro_colorbar_ticks(lo, hi, label; format = _ro_num, sub = nothing,
    height = _RO_COLORBAR_SHARE * _RO_MAP_HEIGHT, labelx = 0.82)

    p = plot(framestyle = :none, ticks = false, xlims = (0, 1), ylims = (lo, hi),
        legend = false, widen = false, margin = 0mm)
    ticks = first(Plots.PlotUtils.optimize_ticks(Float64(lo), Float64(hi), k_min = 4))
    for t in filter(t -> lo <= t <= hi, ticks)
        plot!(p, [0, 0.07], [t, t], color = _FIG1_LIGHTGREY, linewidth = 3)
        annotate!(p, 0.11, t, text(format(t), _RO_FONTSIZE, :left, _FIG1_FONT))
    end
    mid = (lo + hi) / 2
    if isnothing(sub)
        annotate!(p, labelx, mid, text(label, _RO_FONTSIZE, :center, _FIG1_FONT, rotation = 90))
    else
        # label up to `y`, subscript from `y` on (the text runs upwards; the subscript
        # sits lower, i.e. to the right), both centered as one
        subsize = round(Int, 0.7 * _RO_FONTSIZE)
        w = _fig2_text_width(label, _FIG1_FONT, _RO_FONTSIZE)
        ws = _fig2_text_width(sub, _FIG1_FONT, subsize)
        y = mid + (w - ws) / 2 * (hi - lo) / height
        annotate!(p, labelx, y, text(label, _RO_FONTSIZE, :right, _FIG1_FONT, rotation = 90))
        annotate!(p, labelx + 22 / _FIG1_COLORBAR_LABELS, y,
            text(sub, subsize, :left, _FIG1_FONT, rotation = 90))
    end
    return p
end


"""
    plot_case_fatality(rd)

Plots the infection fatality map (deaths per infection, %, by county; GEMS
`CaseFatalityMap`) without colorbar (see `plot_combined_readouts`).
"""
function plot_case_fatality(rd)
    return gemsmap(rd; type = :CaseFatalityMap, _ro_map()...)
end


"""
    plot_r0_by_region(rd)

Plots the map of the initial reproduction number (R_init) by county (GEMS `R0Map`)
without colorbar (see `plot_combined_readouts`).
"""
function plot_r0_by_region(rd)
    return gemsmap(rd; type = :R0Map, _ro_map()...)
end


# effective R plot, set on the drawn plot: shorter legend labels (GEMS's labels span
# most of the plot width; "7-day rolling" goes into the y label instead), and a
# thicker R = 1 line (GEMS draws it 1 px wide)
function _ro_short_labels!(p)
    for s in p[1].series_list
        l = string(s[:label])
        if occursin("R=1", l)
            s[:label] = "R = 1"
            s[:linewidth] = 4
        else
            s[:label] = occursin("In Households", l) ? "In Households" :
                occursin("Outside Households", l) ? "Outside Households" :
                occursin("Effective R", l) ? "Overall" : l
        end
    end
    return p
end


# the six GEMS plots of the left column (time series cut to the first 100 days, where
# the epidemic takes place; counts in millions; legends only where they name series)
function plot_gemsplots(rd)
    axes = _ro_axes()
    return [
        gemsplot(rd; type = :TickCases, axes...,
            series = [:exposed],
            linewidth = 6,
            xlims = (0, 100),
            xlabel = "Days",
            yformatter = _ro_millions,
            legend = false),
        _ro_short_labels!(gemsplot(rd; type = :EffectiveReproduction, axes...,
            linewidth = 5,
            xlims = (0, 100),
            ylims = (0, 4),
            legend = :topright)),
        # attack rate on its full range (0–100 %), a bit below 0 for the dots on 0
        gemsplot(rd; type = :HouseholdAttackRate, axes...,
            ar_only = true,
            markersize = 30,
            xlims = (1.5, 9.5),
            xticks = 2:9,
            ylims = (-0.05, 1.0),
            yticks = 0:0.2:1,
            yformatter = _ro_percent,
            legend = false),
        # 10 % room above the peak (GEMS's data of the plot)
        gemsplot(rd; type = :InfectionDuration, axes...,
            linecolor = :transparent,
            xlims = (0, 41),
            ylims = (0, 1.1 * maximum(aggregated_compartment_periods(rd).total)),
            yformatter = _ro_percent,
            legend = false),
        gemsplot(rd; type = :CumulativeDiseaseProgressions, axes...,
            xlims = (0, 41),
            yformatter = _ro_millions,
            legend = :topright),
        gemsplot(rd; type = :ProgressionCategories, axes...,
            color = :viridis,
            colorbar = false,
            framestyle = :box,
            grid = false)]
end


"""
    plot_combined_readouts(rd, r0map, fatalitymap)

Combines the GEMS readout plots of `rd` (daily infections, effective R, household
attack rate, infection duration, cumulative disease progressions, progression
categories; left column, A–F) with the maps of the initial reproduction number
(R_init) and the infection fatality by county (`r0map`, `fatalitymap` from
`plot_r0_by_region`, `plot_case_fatality`; right column, G–H). `rd` must be a
default-style ResultData object (`ProgressionCategories` and `CaseFatalityMap` need
the raw infections).

All panels are drawn by GEMS; this function only sets their axes and fonts (as in
the first wave figures), adds titles, y labels, and thin colorbars, and aligns the
plot areas of the left column.
"""
plot_combined_readouts(rd, r0map, fatalitymap) =
    _ro_combined(plot_gemsplots(rd), r0map, fatalitymap)

# layout of the combined figure (`panels`: the six plots of the left column, the
# sixth a heatmap)
function _ro_combined(panels, r0map, fatalitymap)

    # titles (the GEMS plot types) and y labels of the left column
    titles = [
        "Daily Infections",
        "Effective Reproduction Number",
        "In-Household Attack Rate",
        "Infection Duration",
        "Cumulative Disease Progressions",
        "Progression Categories"]
    ylabels = [
        "Infections (Millions)",
        "Effective R (7-Day Mean)",
        "Attack Rate (%)",
        "Fraction of Infections (%)",
        "Individuals (Millions)",
        ""]

    # sizes (px): title rows, maps, colorbar with its ticks, y labels; the left column
    # takes the rest of the width, its six panels share the height of the two maps
    W = _RO_WIDTH
    th, mh = _RO_TITLE_HEIGHT, _RO_MAP_HEIGHT
    mw = round(Int, mh * _FIG1_MAP[1] / _FIG1_MAP[2])
    cw, cl = _FIG1_COLORBAR_WIDTH, _FIG1_COLORBAR_LABELS
    gap = round(Int, _FIG1_GAP.value / 25.4 * Plots.DPI)
    rw = mw + gap + cw + cl
    lw = W - rw
    H = 2 * th + 2 * mh
    ph = (H - 6 * th) / 6
    rest = (1 - _RO_COLORBAR_SHARE) / 2

    # y label columns: the heatmap's tick labels (the widest) start after its narrow
    # column, which sets where all plot areas start; the other columns end where the
    # widest of their tick labels ("100" of the attack rate) starts
    heat = panels[6]
    heat_ticks = heat[1][:yaxis][:ticks]
    heat_labels = heat_ticks isa Tuple ? heat_ticks[2] : String[]
    wh = isempty(heat_labels) ? 0.0 : maximum(_ro_tickwidth.(string.(heat_labels)))
    yw_axes = round(Int, _RO_YLABEL_WIDTH_HEATMAP + wh - _ro_tickwidth("100"))
    ywidths = [fill(yw_axes, 5); _RO_YLABEL_WIDTH_HEATMAP]

    # color ranges as GEMS drew them
    heat_lo, heat_hi = _ro_clims(heat)
    r0_lo, r0_hi = _ro_clims(r0map)
    cf_lo, cf_hi = _ro_clims(fatalitymap)
    plot!(r0map, clims = (r0_lo, r0_hi))
    plot!(fatalitymap, clims = (cf_lo, cf_hi))

    # Layout. Plots lines up the plot areas of the cells of a row (top, bottom) and of
    # a column (left, right). Hence: each row of the left column starts with its y
    # label (no padding, so the titles above start at the left edge) and ends with a
    # strip for the heatmap's colorbar (empty for the other panels, so all plot areas
    # end at the same x); the colorbars of the maps sit in a cell of their own with
    # empty space above and below (as in figure 1; margins would shrink the maps); both
    # columns end with an empty row (the heatmap's tick labels would lift the second
    # map otherwise).
    tiny = 1
    l = grid(1, 2, widths = [lw, rw] ./ W)
    l[1, 1] = grid(13, 1, heights = [repeat([th, ph], 6); tiny] ./ (H + tiny))
    for k in 1:6
        yw = ywidths[k]
        l[1, 1][2k, 1] = grid(1, 3, widths = [yw, lw - yw - cw - cl, cw + cl] ./ lw)
    end
    l[1, 1][12, 1][1, 3] = grid(1, 2, widths = [cw, cl] ./ (cw + cl))
    l[1, 2] = grid(5, 1, heights = [th, mh, th, mh, tiny] ./ (H + tiny))
    for r in (2, 4)
        l[1, 2][r, 1] = grid(1, 2, widths = [mw + gap, cw + cl] ./ rw)
        l[1, 2][r, 1][1, 2] = grid(3, 1, heights = [rest, 1 - 2 * rest, rest])
        l[1, 2][r, 1][1, 2][2, 1] = grid(1, 2, widths = [cw, cl] ./ (cw + cl))
    end

    # panels in the order of the layout; remember where the plots of the left column
    # and the heatmap's colorbar are
    items = Any[]
    column = Int[]
    heat_cb = Int[]
    for k in 1:6
        push!(items, _ro_title('A' + (k - 1), titles[k]), _ro_ylabel(ylabels[k], ywidths[k]),
            panels[k])
        push!(column, length(items))
        if k < 6
            push!(items, _fig1_text())
        else
            push!(items,
                _ro_colorbar(heat_lo, heat_hi, cgrad(:viridis)),
                _ro_colorbar_ticks(heat_lo, heat_hi, "Infections (Millions)",
                    format = _ro_millions, labelx = 0.52))
            append!(heat_cb, length(items) - 1:length(items))
        end
    end
    push!(items, _fig1_text())
    for (letter, title, m, lo, hi, label, sub) in (
            ("G", "Initial Reproduction Number by County", r0map, r0_lo, r0_hi, "R", "init"),
            ("H", "Infection Fatality by County", fatalitymap, cf_lo, cf_hi,
                "Deaths per Infection (%)", nothing))
        push!(items,
            _ro_title(letter, title),
            m,
            _fig1_text(),
            _ro_colorbar(lo, hi, cgrad(:inferno)),
            _ro_colorbar_ticks(lo, hi, label, sub = sub),
            _fig1_text())
    end
    push!(items, _fig1_text())

    p = plot(items...,
        layout = l,
        size = (W, H),
        dpi = Plots.DPI,
        fontfamily = _FIG1_FONT)

    # line up the left edges of the plot areas of the left column (tick labels of
    # different widths), and the heatmap's colorbar with the heatmap if the row did not
    # (repeat: each move changes the plot areas)
    sp = p.subplots
    hm = sp[column[6]]
    for _ in 1:4
        Plots.prepare_output(p)
        target = maximum(Plots.left(sp[i].plotarea) for i in column)
        moved = false
        for i in column
            offset = target - Plots.left(sp[i].plotarea)
            if abs(offset.value) > 0.1
                sp[i][:left_margin] = sp[i][:left_margin] + offset
                moved = true
            end
        end
        for i in heat_cb
            top = Plots.top(hm.plotarea) - Plots.top(sp[i].plotarea)
            bottom = Plots.bottom(sp[i].plotarea) - Plots.bottom(hm.plotarea)
            if abs(top.value) > 0.1 || abs(bottom.value) > 0.1
                sp[i][:top_margin] = sp[i][:top_margin] + top
                sp[i][:bottom_margin] = sp[i][:bottom_margin] + bottom
                moved = true
            end
        end
        moved || break
    end
    return p
end


"""
    scenario_household_isolation(;
        configfile = "models/DE_model_first_wave.toml",
        infected_fraction = 0.0001,
        stop_criterion = TimesUp(limit = 180),
        population = "DE",
        folder = nothing,
        folder_tag = nothing)

Runs the household isolation scenario and all analyses.

If `folder` is `nothing`, a new result folder is created, with `folder_tag`
appended to its timestamp-based name (see `result_folder`).
"""
function scenario_household_isolation(;
    configfile = "models/DE_model_first_wave.toml",
    infected_fraction = 0.0001,
    stop_criterion = TimesUp(limit = 180),
    population = "DE",
    folder = nothing,
    folder_tag = nothing)

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder(folder_tag) : nothing

    # run simulations and get result data objects
    rd_b, rd_s, sim = run_household_isolation_sims(
        configfile = configfile,
        infected_fraction = infected_fraction,
        stop_criterion = stop_criterion,
        population = population)

    p_hi = plot_scenario_household_isolation(rd_b, rd_s, sim)
    png(p_hi, joinpath(folder, "sc_household_isolation_r0_reduction.png"))

    # GEMS readouts of the baseline
    p_readouts = plot_combined_readouts(rd_b, plot_r0_by_region(rd_b), plot_case_fatality(rd_b))
    png(p_readouts, joinpath(folder, "sc_household_isolation_baseline_combined.png"))
end