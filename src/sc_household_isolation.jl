export HouseholdIsolationRD
export run_household_isolation_sims
export plot_scenario_household_isolation
export plot_r0_by_region
export plot_case_fatality
export plot_combined_readouts
export save_household_isolation_figure

export scenario_household_isolation
export replot_household_isolation

# label of the reduction of the initial reproduction number (R^init in the paper);
# plain text, as GR's LaTeX rendering is unreliable with the Times font
const _RINIT_REDUCTION_LABEL = "Reduction of initial R"


###
### RESULT DATA STYLE
###

"""
    HouseholdIsolationRD <: ResultDataStyle

Result data of the household isolation scenario: only what the figures need, so that both
runs can be stored (`exportJLD`) and loaded again (`import_resultdata`). Of the raw infections,
it holds their numbers per progression category and age group (`progression_ages`) and the
infection fatality per county (`fatality_per_county`), and of the population the mean
household size and age per county (`county_means`).
"""
mutable struct HouseholdIsolationRD <: ResultDataStyle
    data::Dict{String, Any}
    function HouseholdIsolationRD(pP::PostProcessor)
        sim = simulation(pP)
        return new(Dict(
            "sim_data" => Dict(
                "label" => label(sim),
                "seed" => seed(sim),
                "final_tick" => tick(sim),
                "tick_unit" => tickunit(sim),
                "total_infections" => nrow(infectionsDF(pP)),
                "attack_rate" => attack_rate(pP),
                "r0" => r0(pP)),
            "dataframes" => Dict(
                "tick_cases" => tick_cases(pP),
                "cumulative_cases" => cumulative_cases(pP),
                "effectiveR" => effectiveR(pP),
                "household_attack_rates" => household_attack_rates(pP),
                "aggregated_compartment_periods" => aggregated_compartment_periods(pP),
                "cumulative_disease_progressions" => cumulative_disease_progressions(pP),
                "r0_per_county" => r0_per_county(pP),
                "progression_ages" => progression_ages(infectionsDF(pP)),
                "fatality_per_county" => fatality_per_county(infectionsDF(pP)),
                "county_means" => county_means(sim))))
    end
end

# progression categories by severity and width of the age groups (as in the plot
# `ProgressionCategories` of GEMS)
const HI_CATEGORIES = ["Asymptomatic", "Symptomatic", "Severe", "Hospitalized", "Critical"]
const HI_AGE_GROUP = 5

# number of infections per progression category and age group (named by its youngest age)
function progression_ages(infections)
    counts = DataFrame(
            category = infections.progression_category,
            age_group = infections.age_b .÷ HI_AGE_GROUP .* HI_AGE_GROUP) |>
        df -> combine(groupby(df, [:category, :age_group]), nrow => :count) |>
        df -> transform(df, :category => ByRow(string) => :category)
    issubset(counts.category, HI_CATEGORIES) ||
        throw(ArgumentError("Unknown progression categories: $(setdiff(counts.category, HI_CATEGORIES))."))

    # all combinations of the categories and age groups
    counted = Dict(zip(zip(counts.category, counts.age_group), counts.count))
    return DataFrame([(category = c, age_group = a, count = get(counted, (c, a), 0))
        for c in HI_CATEGORIES if c in counts.category
        for a in 0:HI_AGE_GROUP:maximum(counts.age_group)])
end

# deaths per infection (%) by county of residence (as in the map `CaseFatalityMap` of GEMS)
fatality_per_county(infections) =
    DataFrame(ags = county.(infections.household_ags_b), death = infections.death .>= 0) |>
        df -> combine(groupby(df, :ags), :death => (d -> 100 * mean(d)) => :case_fatality_rate)

# mean household size and mean age per county
function county_means(sim)
    household_sizes = DataFrame(
        ags = county.(ags.(households(sim))),
        size = size.(households(sim)))
    ages = DataFrame(
        ags = county.(ags.((i -> household(i, sim)).(individuals(sim)))),
        age = age.(individuals(sim)))
    return innerjoin(
        combine(groupby(household_sizes, :ags), :size => mean => :size),
        combine(groupby(ages, :ags), :age => mean => :age), on = :ags)
end


###
### SIMULATION
###

"""
    run_household_isolation_sims(;
        configfile = "models/DE_model_first_wave.toml",
        infected_fraction = 0.0001,
        stop_criterion = TimesUp(limit = 180),
        population = "DE")

Runs the household isolation scenario simulations and returns two ResultData objects (style
`HouseholdIsolationRD`): one for the baseline simulation and one for the scenario with household isolation.
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
    rd_b = ResultData(sim, style = "HouseholdIsolationRD")


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
    rd_s = ResultData(sim, style = "HouseholdIsolationRD")

    return rd_b, rd_s
end


###
### FIGURE
###

# height of the plots (inches)
const HI_HEIGHT_IN = 4.17

# share of the width for the time series, and of the height for the map
const HI_TIME_SERIES_SHARE = 0.55
const HI_MAP_SHARE = 0.62

# steps between the ticks of the scatter plots
const HI_AGE_STEP = 2
const HI_SIZE_STEP = 0.1
const HI_R_STEP = 0.1

# size of their points and of the point in their legend, size of the legend and its distance from the axes (inches)
const HI_MARKERSIZE = 1.5
const HI_LEGEND_MARKERSIZE = 2.5
const HI_LEGEND_IN = (0.55, 0.16)
const HI_LEGEND_GAP_IN = 0.05

# width of the plot area of the map relative to the map (room for the colorbar)
const HI_MAP_ROOM = 1.25
const HI_MAP_XLIMS = (_FIG1_XLIMS[1], _FIG1_XLIMS[1] + HI_MAP_ROOM * (_FIG1_XLIMS[2] - _FIG1_XLIMS[1]))

# width of the colorbar and its distance from the map (inches)
const HI_COLORBAR_IN = 0.12
const HI_COLORBAR_GAP_IN = 0.1

# limits of an axis: the ticks in steps of `step` around the values
tick_limits(values, step) = (step * floor(minimum(values) / step), step * ceil(maximum(values) / step))

"""
    plot_scenario_household_isolation(rd_b, rd_s)

Plots the household isolation scenario results.
Requires two ResultData objects (style `HouseholdIsolationRD`) as input, one for the baseline and
one for the scenario with household isolation.
"""
function plot_scenario_household_isolation(rd_b, rd_s)

    # reduction of R0, mean household size, and mean age per county
    df = innerjoin(r0_per_county(rd_b), rename(r0_per_county(rd_s), :r0 => :r0_new), on = :ags) |>
        df -> transform(df, [:r0, :r0_new] => ByRow(-) => :r0_diff) |>
        df -> innerjoin(df, dataframes(rd_b)["county_means"], on = :ags)

    style = (
        titlelocation = :left,
        fontfamily = FIG_FONT,
        titlefontsize = FIG_HEADING_FONTSIZE,
        guidefontsize = FIG_FONTSIZE,
        tickfontsize = FIG_SMALL_FONTSIZE,
        legendfontsize = FIG_SMALL_FONTSIZE)
    axes = (gridalpha = 0.15, gridlinewidth = 0.5)

    county_map = agsmap(select(df, :ags, :r0_diff); style...,
        title = "$_RINIT_REDUCTION_LABEL by county",
        fillcolor = :inferno,
        clims = extrema(df.r0_diff),
        colorbar = false, # drawn on the canvas (thinner)
        legend = false,
        linecolor = nothing,
        framestyle = :none,
        xlims = HI_MAP_XLIMS,
        ylims = _FIG1_YLIMS,
        widen = false)

    # x label as annotation (GR cuts off the lower parts of letters in its own), `y` relative to the plot area
    xlabel!(p, label, y) = annotate!(p, (0.5, y), Plots.text(label, FIG_FONTSIZE, FIG_FONT, :center))
    xlabel_margin = 5mm

    r_limits = tick_limits(df.r0_diff, HI_R_STEP)
    r_ticks = r_limits[1]:HI_R_STEP:r_limits[2]
    r_limits = r_limits .+ (0, HI_R_STEP / 2) # room for a legend
    function by_county(column, step, xlabel; kwargs...)
        limits = tick_limits(df[!, column], step)
        p = scatter(df[!, column], df.r0_diff; style..., axes...,
            xlims = limits, xticks = limits[1]:step:limits[2],
            ylims = r_limits, yticks = r_ticks,
            widen = false,
            label = "",
            markersize = HI_MARKERSIZE,
            markerstrokewidth = 0,
            top_margin = -4mm, # into the empty space below the map
            bottom_margin = xlabel_margin,
            kwargs...)
        return xlabel!(p, xlabel, -0.27)
    end

    # y labels at the same distance from all time series
    function over_time(type, title, ylabel; kwargs...)
        p = gemsplot([rd_b, rd_s]; type = type, style..., axes...,
            plot_title = "",
            title = title,
            xlabel = "",
            ylabel = "",
            left_margin = 5mm,
            linewidth = 1,
            legend = false,
            kwargs...)
        return annotate!(p, (-0.13, 0.5), Plots.text(ylabel, FIG_FONTSIZE, FIG_FONT, :center, rotation = 90))
    end

    l = grid(1, 2, widths = [1 - HI_TIME_SERIES_SHARE, HI_TIME_SERIES_SHARE])
    l[1, 1] = grid(2, 1, heights = [HI_MAP_SHARE, 1 - HI_MAP_SHARE])
    l[1, 1][2, 1] = grid(1, 2, widths = [0.46, 0.54])
    l[1, 2] = grid(3, 1)

    p = plot(
        county_map,
        by_county(:age, HI_AGE_STEP, "Mean age", ylabel = _RINIT_REDUCTION_LABEL,
            xformatter = x -> string(round(Int, x)), left_margin = 1.5mm),
        by_county(:size, HI_SIZE_STEP, "Mean household size",
            xformatter = x -> string(round(x, digits = 1)), yformatter = _ -> "",
            left_margin = -2.5mm, right_margin = 2.5mm), # further left
        over_time(:TickCases, "Daily infections", "Millions",
            yformatter = _ro_millions,
            legend = :topright, foreground_color_legend = FIG_COLORS.lightgray),
        over_time(:CumulativeCases, "Cumulative infections", "Millions",
            yformatter = _ro_millions),
        xlabel!(over_time(:EffectiveReproduction, "Effective reproduction number", "7-day mean",
            bottom_margin = xlabel_margin), "Day", -0.39),
        layout = l,
        size = (100FIG_FULL_WIDTH_IN, 100HI_HEIGHT_IN),
        fontfamily = FIG_FONT)
    return draw_legends!(p, df)
end

# position (left, right, bottom, top; inches from the lower left corner) of the plot area of the subplot `k`
# of a figure that is laid out; `aspect`: width per height, if fixed (centered in the area that Plots reserves)
function plot_corners(p, k; aspect = nothing)
    inches(x) = x.value / 25.4
    area = p[k].plotarea
    w, h = inches(Plots.width(area)), inches(Plots.height(area))
    x, y = inches(Plots.left(area)) + w / 2, p.attr[:size][2] / 100 - inches(Plots.top(area)) - h / 2
    width = isnothing(aspect) ? w : min(w, h * aspect)
    height = isnothing(aspect) ? h : width / aspect
    return (left = x - width / 2, right = x + width / 2, bottom = y - height / 2, top = y + height / 2)
end

# colorbar of the map (subplot 1) and legends of the scatter plots (lower left corner of the first, upper left
# corner of the second, where there are no points; subplots 2 and 3) on a transparent canvas in inches above the figure
function draw_legends!(p, df)
    Plots.prepare_output(p) # lays out the figure
    county_map = plot_corners(p, 1, aspect = (HI_MAP_XLIMS[2] - HI_MAP_XLIMS[1]) / (_FIG1_YLIMS[2] - _FIG1_YLIMS[1]))
    map_right = county_map.left + (county_map.right - county_map.left) / HI_MAP_ROOM
    by_age, by_size = plot_corners(p, 2), plot_corners(p, 3)

    canvas = length(p.subplots) + 1
    plot!(p, inset = bbox(0, 0, 1, 1), subplot = canvas, framestyle = :none, ticks = nothing, legend = false,
        grid = false, widen = false, margin = -2mm, background_color_subplot = :transparent,
        background_color_inside = :transparent, xlims = (0, FIG_FULL_WIDTH_IN), ylims = (0, HI_HEIGHT_IN))

    county_colorbar!(p, canvas, CountyMap(column = :r0_diff, title = "", note = ""), extrema(df.r0_diff),
        map_right + HI_COLORBAR_GAP_IN, county_map.bottom, HI_COLORBAR_IN, county_map.top - county_map.bottom,
        colors = cgrad(:inferno))

    w, h = HI_LEGEND_IN
    for (left, bottom) in (
            (by_age.left + HI_LEGEND_GAP_IN, by_age.bottom + HI_LEGEND_GAP_IN),
            (by_size.left + HI_LEGEND_GAP_IN, by_size.top - HI_LEGEND_GAP_IN - h))
        plot!(p, left .+ [0, w, w, 0, 0], bottom .+ [0, 0, h, h, 0], subplot = canvas, label = "",
            seriestype = :shape, fillcolor = :white, linecolor = FIG_COLORS.lightgray, linewidth = 1)
        scatter!(p, [left + 0.08], [bottom + h / 2], subplot = canvas, label = "",
            markercolor = 1, markersize = HI_LEGEND_MARKERSIZE, markerstrokewidth = 0)
        annotate!(p, [(left + w - 0.04, bottom + h / 2,
            Plots.text("County", FIG_SMALL_FONTSIZE, FIG_FONT, :right, :vcenter))], subplot = canvas)
    end
    return p
end

"""
    save_household_isolation_figure(rd_b, rd_s, folder)

Saves the figure of the paper (`Fig_15_household_isolation.png`: the plots of
`plot_scenario_household_isolation` above the code listing, a screenshot of
`sc_household_isolation_screenshot_file.jl`) and the plots alone (`Fig_15A_household_isolation.png`)
in `folder`.
"""
function save_household_isolation_figure(rd_b, rd_s, folder)
    plots = figure_image(plot_scenario_household_isolation(rd_b, rd_s))
    save_figure(plots, :household_isolation, folder, panel = "A")

    width = size(plots, 2)
    code = load_panel(joinpath(FIG_PANEL_FOLDER, figure_name(:household_isolation, panel = "B", name = "code_listing") * ".png"))
    return save_figure(combine_panels([
        section_strip("A) Scenario comparison in GEMS", width), plots,
        section_strip("B) Scenario code example", width), scale_image(code, width)], direction = :vertical, gap_in = 0.05), :household_isolation, folder)
end


###
### GEMS READOUTS OF THE BASELINE (combined figure)
###

# The panels are drawn by GEMS (`gemsplot`, `gemsmap`); only their keyword arguments
# (axes, ticks, fonts, colors) and the layout are set here. Font, panel titles, and
# colorbars in the style of the first wave figures (`_FIG1_*`); width as figure 2.

# sizes (px): figure width and height (a full page of PLOS: 7.5 × 8.75 in; the maps
# span the plot areas of A–C and D–F, their widths follow from the aspect of the
# county shapes, `_FIG1_MAP`), share of the map height that its colorbar spans, height
# of the title rows, width of the y label column of the heatmap (F; its tick labels
# start there; the other rows' y label columns end where their tick labels start, see
# `_ro_combined`), gap between y label and tick labels
const _RO_WIDTH = _FIG2_WIDTH
const _RO_HEIGHT = 7370

# share of the height of A–C (D–F) by which the maps reach below C (F)
const _RO_MAP_SPILL = 0.05
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
# centered for a colorbar of `height` px (needed only with `sub`). `labelx`: position
# of the label (0: at the colorbar, 1: right edge of the cell).
function _ro_colorbar_ticks(lo, hi, label; format = _ro_num, sub = nothing,
    height = nothing, labelx = 0.82)

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

Plots the infection fatality map (deaths per infection, %, by county; as the GEMS map
`CaseFatalityMap`) without colorbar (see `plot_combined_readouts`).
"""
function plot_case_fatality(rd)
    return agsmap(dataframes(rd)["fatality_per_county"]; _ro_map()...)
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


# heatmap of the progression categories per age group (as the GEMS plot `ProgressionCategories`)
function plot_progression_categories(rd; plotargs...)
    df = dataframes(rd)["progression_ages"]
    return gemsheatmap(df.age_group, Int.(indexin(df.category, HI_CATEGORIES)), df.count;
        xformatter = a -> "$a-$(a + HI_AGE_GROUP - 1)",
        yformatter = k -> HI_CATEGORIES[k],
        xlabel = "Age Group",
        xrotation = 45,
        plotargs...)
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
        # attack rate on its full range (0–100 %), a bit below 0 for the dots on 0; household
        # sizes up to 6 (GEMS samples only a few larger households)
        gemsplot(rd; type = :HouseholdAttackRate, axes...,
            ar_only = true,
            markersize = 30,
            xlims = (1.5, 6.5),
            xticks = 2:6,
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
        plot_progression_categories(rd; axes...,
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
ResultData object of the style `HouseholdIsolationRD`.

All panels are drawn by GEMS; this function only sets their axes and fonts (as in
the first wave figures), adds titles, y labels, and thin colorbars, aligns the plot
areas of the left column, and fits the maps to them (G next to A–C, H next to D–F).
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

    # sizes (px): figure, title rows, panel rows of the left column (six, with their
    # titles, over the full height), colorbar with its ticks. The maps span the plot
    # areas of A–C and D–F, which Plots knows only once it has laid out the figure:
    # the first layout uses an estimate (`mh`), `fit_maps!` below sets the right
    # column to the measured plot areas (rows: title G, map G, space, title H, map H,
    # space; the map titles in the rows of the titles of A and D)
    W, H = _RO_WIDTH, _RO_HEIGHT
    th = _RO_TITLE_HEIGHT
    ph = (H - 6 * th) / 6
    mh = 3 * (th + ph) - th - ph / 4
    aspect = (_FIG1_XLIMS[2] - _FIG1_XLIMS[1]) / (_FIG1_YLIMS[2] - _FIG1_YLIMS[1])
    cw, cl = _FIG1_COLORBAR_WIDTH, _FIG1_COLORBAR_LABELS
    gap = round(Int, _FIG1_GAP.value / 25.4 * Plots.DPI)
    rw = mh * aspect + gap + cw + cl
    lw = W - rw
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
    space = (H + tiny - 2 * th - 2 * mh) / 2
    l[1, 2] = grid(6, 1, heights = [th, mh, space, th, mh, space] ./ (H + tiny))
    for r in (2, 5)
        l[1, 2][r, 1] = grid(1, 2, widths = [rw - cw - cl, cw + cl] ./ rw)
        l[1, 2][r, 1][1, 2] = grid(3, 1, heights = [rest, 1 - 2 * rest, rest])
        l[1, 2][r, 1][1, 2][2, 1] = grid(1, 2, widths = [cw, cl] ./ (cw + cl))
    end

    # panels in the order of the layout; remember where the titles and plots of the left
    # column, the heatmap's colorbar, the maps, and the rows of the right column are
    items = Any[]
    ltitles = Int[]
    column = Int[]
    heat_cb = Int[]
    for k in 1:6
        push!(items, _ro_title('A' + (k - 1), titles[k]))
        push!(ltitles, length(items))
        push!(items, _ro_ylabel(ylabels[k], ywidths[k]), panels[k])
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
    maps = Int[]
    rows = Int[]
    for (letter, title, m, lo, hi, label, sub) in (
            ("G", "Initial Reproduction Number by County", r0map, r0_lo, r0_hi, "R", "init"),
            ("H", "Infection Fatality by County", fatalitymap, cf_lo, cf_hi,
                "Deaths per Infection (%)", nothing))
        push!(items, _ro_title(letter, title))
        push!(rows, length(items))
        push!(items, m)
        push!(maps, length(items))
        push!(rows, length(items))
        push!(items,
            _fig1_text(),
            _ro_colorbar(lo, hi, cgrad(:inferno)),
            _ro_colorbar_ticks(lo, hi, label, sub = sub, height = _RO_COLORBAR_SHARE * mh),
            _fig1_text(),
            _fig1_text())
        push!(rows, length(items))
    end

    p = plot(items...,
        layout = l,
        size = (W, H),
        dpi = Plots.DPI,
        fontfamily = _FIG1_FONT)

    # Right column from the left one as Plots laid it out (px from the top left): the
    # rows of the titles G and H as those of A and D, the plot areas of the maps from
    # the top of A to the bottom of C and from the top of D to the bottom of F (plus
    # `_RO_MAP_SPILL` of that height below), as wide
    # as their height gives (the narrower map moved right, next to its colorbar); the
    # left column takes the rest of the width. Plots splits the height of a column
    # among its rows by their fractions after the rows' padding, hence the padding is
    # taken off. Returns whether anything changed by more than half a pixel.
    sp = p.subplots
    right = p.layout[1, 2]
    topx(m) = m.value / 25.4 * Plots.DPI
    ptop(i; box = false) = topx(Plots.top(box ? sp[i].bbox : sp[i].plotarea))
    pbottom(i; box = false) = topx(Plots.bottom(box ? sp[i].bbox : sp[i].plotarea))
    pwidth(i) = topx(Plots.width(sp[i].plotarea))
    pheight(i; box = false) = pbottom(i, box = box) - ptop(i, box = box)
    function fit_maps!()
        g, h = maps
        span(upper, lower) = (1 + _RO_MAP_SPILL) * (pbottom(lower) - ptop(upper))
        # rows: their borders, padding, and fractions
        ends = [
            pbottom(ltitles[1], box = true),
            ptop(column[1]) + span(column[1], column[3]) + pbottom(g, box = true) - pbottom(g),
            ptop(ltitles[4], box = true),
            pbottom(ltitles[4], box = true),
            ptop(column[4]) + span(column[4], column[6]) + pbottom(h, box = true) - pbottom(h),
            pbottom(rows[end], box = true)]
        sizes = diff([0; ends])
        pads = [pheight(i, box = true) - pheight(i) for i in rows]
        changed = any(abs.(sizes .- [pheight(i, box = true) for i in rows]) .> 0.5)
        fractions = max.(sizes .- pads, 1)
        right.heights = Measure[x / sum(fractions) * Plots.pct for x in fractions]
        # tops of the maps at the tops of A and D (as far as the rows do not set them)
        for (k, start, panel) in ((g, ends[1], column[1]), (h, ends[4], column[4]))
            e = (ptop(panel) - start) - (ptop(k) - ptop(k, box = true))
            if abs(e) > 0.5
                sp[k][:top_margin] = sp[k][:top_margin] + e * Plots.px
                changed = true
            end
        end
        # widths: surplus of each map's plot area over the width its height gives; the
        # column loses the smaller surplus, the other map moves right by the difference
        surplus = [pwidth(k) - span(upper, lower) * aspect
            for (k, upper, lower) in ((g, column[1], column[3]), (h, column[4], column[6]))]
        cut = minimum(surplus)
        if abs(cut) > 0.5
            cwidth = p.layout.widths[2].value * W - cut
            p.layout.widths = Measure[w / W * Plots.pct for w in (W - cwidth, cwidth)]
            for r in (2, 5)
                right[r, 1].widths = Measure[w / cwidth * Plots.pct
                    for w in (cwidth - cw - cl, cw + cl)]
            end
            changed = true
        end
        for (k, s) in zip(maps, surplus)
            if abs(s - cut) > 0.5
                sp[k][:left_margin] = sp[k][:left_margin] + (s - cut) * Plots.px
                changed = true
            end
        end
        return changed
    end

    # fit the maps (see above), line up the left edges of the plot areas of the left
    # column (tick labels of different widths), and the heatmap's colorbar with the
    # heatmap if the row did not (repeat: each move changes the plot areas)
    hm = sp[column[6]]
    for _ in 1:10
        Plots.prepare_output(p)
        moved = fit_maps!()
        moved && Plots.prepare_output(p)
        target = maximum(Plots.left(sp[i].plotarea) for i in column)
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

Runs the household isolation scenario and all analyses. The `ResultData` objects of both runs are
stored in `folder`, from which the figures are drawn (see `replot_household_isolation`).

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

    # run simulations and store result data objects
    rd_b, rd_s = run_household_isolation_sims(
        configfile = configfile,
        infected_fraction = infected_fraction,
        stop_criterion = stop_criterion,
        population = population)
    exportJLD(rd_b, joinpath(folder, HI_BASELINE_FOLDER))
    exportJLD(rd_s, joinpath(folder, HI_SCENARIO_FOLDER))

    replot_household_isolation(folder)
end

# subfolders of the stored result data objects
const HI_BASELINE_FOLDER = joinpath("household_isolation", "baseline")
const HI_SCENARIO_FOLDER = joinpath("household_isolation", "scenario")

"""
    replot_household_isolation(folder)

Draws the figures of the household isolation scenario from the `ResultData` objects that
`scenario_household_isolation()` stored in `folder`, without re-running the simulations (see
`save_household_isolation_figure` and `plot_combined_readouts`).
"""
function replot_household_isolation(folder)
    rd_b, rd_s = [import_resultdata(joinpath(folder, f, "resultdata.jld2"))
        for f in (HI_BASELINE_FOLDER, HI_SCENARIO_FOLDER)]

    save_household_isolation_figure(rd_b, rd_s, folder)

    # GEMS readouts of the baseline
    p_readouts = plot_combined_readouts(rd_b, plot_r0_by_region(rd_b), plot_case_fatality(rd_b))
    png(p_readouts, joinpath(folder, "sc_household_isolation_baseline_combined.png"))
end
