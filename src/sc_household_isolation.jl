export HouseholdIsolationRD
export run_household_isolation_sims
export plot_scenario_household_isolation
export plot_baseline_readouts
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

# layout (inches): left and right edges of the map, the scatter plots, and the time series, heights of the
# titles, the scatter plots, and the space below the plots (with an x label, with tick labels only), and
# distances of the axis labels from the plot areas
const HI_MAP = (0.92, 2.475)
const HI_SCATTER = ((0.5, 1.86), (2.04, 3.4))
const HI_TIME_SERIES = (3.98, 7.42)
const HI_TITLE_IN = 0.26
const HI_SCATTER_IN = 1.28
const HI_BELOW_IN = 0.42
const HI_TICKS_IN = 0.2
const HI_YLABEL_IN = 0.4
const HI_XLABEL_IN = 0.3

# steps between the ticks of the scatter plots
const HI_AGE_STEP = 2
const HI_SIZE_STEP = 0.1
const HI_R_STEP = 0.1

# size of their points and of the point in their legend, size of the legend and its distance from the axes (inches)
const HI_MARKERSIZE = 1.5
const HI_LEGEND_MARKERSIZE = 2.5
const HI_LEGEND_IN = (0.55, 0.16)
const HI_LEGEND_GAP_IN = 0.05

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
        title = "", plot_title = "", xlabel = "", ylabel = "",
        legend = false,
        fontfamily = FIG_FONT,
        tickfontsize = FIG_SMALL_FONTSIZE,
        gridalpha = 0.15, gridlinewidth = 0.5)

    r_limits = tick_limits(df.r0_diff, HI_R_STEP)
    r_ticks = r_limits[1]:HI_R_STEP:r_limits[2]
    r_limits = r_limits .+ (0, HI_R_STEP / 2) # room for a legend
    function by_county(column, step; kwargs...)
        limits = tick_limits(df[!, column], step)
        return scatter(df[!, column], df.r0_diff; style...,
            xlims = limits, xticks = limits[1]:step:limits[2],
            ylims = r_limits, yticks = r_ticks,
            widen = false,
            markersize = HI_MARKERSIZE,
            markerstrokewidth = 0,
            kwargs...)
    end
    over_time(type; kwargs...) = gemsplot([rd_b, rd_s]; type = type, style..., linewidth = 1.2, kwargs...)

    # left column: (title) | map | (space) | scatter plots | (labels); right column: three times (title) | time
    # series | (labels); empty cells around the plots
    width, height = FIG_FULL_WIDTH_IN, HI_HEIGHT_IN
    shares(x) = collect(x) ./ sum(x)
    split = HI_SCATTER[2][2] + 0.1 # between both columns
    map_in = (HI_MAP[2] - HI_MAP[1]) / _FIG1_ASPECT
    space = height - HI_TITLE_IN - map_in - HI_SCATTER_IN - HI_BELOW_IN
    series_in = (height - 3HI_TITLE_IN - 2HI_TICKS_IN - HI_BELOW_IN) / 3
    l = grid(1, 2, widths = [split, width - split] ./ width)
    l[1, 1] = grid(5, 1, heights = shares((HI_TITLE_IN, map_in, space, HI_SCATTER_IN, HI_BELOW_IN)))
    l[1, 1][2, 1] = grid(1, 3, widths = shares((HI_MAP[1], HI_MAP[2] - HI_MAP[1], split - HI_MAP[2])))
    (a, b), (c, d) = HI_SCATTER
    l[1, 1][4, 1] = grid(1, 5, widths = shares((a, b - a, c - b, d - c, split - d)))
    l[1, 2] = grid(7, 1, heights = shares((HI_TITLE_IN, series_in, HI_TICKS_IN + HI_TITLE_IN, series_in,
        HI_TICKS_IN + HI_TITLE_IN, series_in, HI_BELOW_IN)))
    a, b = HI_TIME_SERIES
    for r in (2, 4, 6)
        l[1, 2][r, 1] = grid(1, 3, widths = shares((a - split, b - a, width - b)))
    end

    e = empty_panel
    p = plot(
        e(), e(), plot_county_map(select(df, :ags, :r0_diff)), e(), e(),
        e(), by_county(:age, HI_AGE_STEP, xformatter = x -> string(round(Int, x))), e(),
        by_county(:size, HI_SIZE_STEP, xformatter = x -> string(round(x, digits = 1)), yformatter = _ -> ""), e(), e(),
        e(), e(), over_time(:TickCases, yformatter = _ro_millions), e(), e(),
        e(), over_time(:CumulativeCases, yformatter = _ro_millions), e(), e(),
        e(), over_time(:EffectiveReproduction), e(), e(),
        layout = l,
        size = (100width, 100height),
        fontfamily = FIG_FONT)
    county_map, by_age, by_size, infections, cumulative, reproduction = 3, 7, 9, 14, 18, 22
    remove_padding!(p)
    canvas = add_canvas!(p)

    # titles and axis labels
    corners(k) = plot_corners(p, k)
    for (k, title) in (
            (county_map, "$_RINIT_REDUCTION_LABEL by county"),
            (infections, "Daily infections"),
            (cumulative, "Cumulative infections"),
            (reproduction, "Effective reproduction number"))
        canvas_text!(p, canvas, corners(k).left, corners(k).top + HI_TITLE_IN / 2, title, size = FIG_HEADING_FONTSIZE,
            halign = :left)
    end
    for (k, label) in ((by_age, "Mean age"), (by_size, "Mean household size"), (reproduction, "Day"))
        canvas_text!(p, canvas, (corners(k).left + corners(k).right) / 2, corners(k).bottom - HI_XLABEL_IN, label)
    end
    for (k, label) in ((by_age, _RINIT_REDUCTION_LABEL), (infections, "Millions"), (cumulative, "Millions"),
            (reproduction, "7-day mean"))
        canvas_text!(p, canvas, corners(k).left - HI_YLABEL_IN, (corners(k).bottom + corners(k).top) / 2, label,
            rotation = 90)
    end

    # colorbar, label of R = 1, and legends
    map_colorbar!(p, canvas, county_map, extrema(df.r0_diff))
    r_threshold_label!(p, canvas, reproduction, halign = :right)
    canvas_legend!(p, canvas, infections, ["Baseline" => "Baseline", "Scenario" => "Scenario"])
    county_legends!(p, canvas, by_age, by_size)
    return p
end

# map of the values of the counties in the second column of `df` (GEMS), without its colorbar
plot_county_map(df; plotargs...) = agsmap(df;
    fillcolor = :inferno,
    clims = extrema(df[!, 2]),
    colorbar = false, # drawn on the canvas (thinner)
    legend = false,
    linecolor = nothing,
    framestyle = :none,
    ticks = false,
    xlims = _FIG1_XLIMS,
    ylims = _FIG1_YLIMS,
    widen = false,
    plotargs...)

# colorbar of the map in subplot `k` on the `canvas` of the figure `p` (laid out), right of the map
function map_colorbar!(p, canvas, k, lims)
    area = plot_corners(p, k, aspect = _FIG1_ASPECT)
    x = area.right + HI_COLORBAR_GAP_IN
    return county_colorbar!(p, canvas, CountyMap(column = :value), lims,
        x, area.bottom, HI_COLORBAR_IN, area.top - area.bottom, colors = cgrad(:inferno))
end

# labels the line that GEMS draws at R = 1 in subplot `k` on the `canvas` of the figure `p` (laid out), above
# its left or right end (`halign`); distances from the end and above the line in inches
function r_threshold_label!(p, canvas, k; halign = :left, dx = 0.06, dy = 0.08)
    area = plot_corners(p, k)
    lo, hi = Plots.ylims(p[k])
    x = halign == :left ? area.left + dx : area.right - dx
    y = area.bottom + (area.top - area.bottom) * (1 - lo) / (hi - lo) + dy
    return canvas_text!(p, canvas, x, y, "R = 1", size = FIG_SMALL_FONTSIZE, halign = halign)
end

# legends of the scatter plots (subplots; lower left corner of the first, upper left corner of the second, where
# there are no points) on the `canvas` of the figure `p`
function county_legends!(p, canvas, by_age, by_size)
    by_age, by_size = plot_corners(p, by_age), plot_corners(p, by_size)
    w, h = HI_LEGEND_IN
    for (left, bottom) in (
            (by_age.left + HI_LEGEND_GAP_IN, by_age.bottom + HI_LEGEND_GAP_IN),
            (by_size.left + HI_LEGEND_GAP_IN, by_size.top - HI_LEGEND_GAP_IN - h))
        plot!(p, left .+ [0, w, w, 0, 0], bottom .+ [0, 0, h, h, 0], subplot = canvas, label = "",
            seriestype = :shape, fillcolor = :white, linecolor = FIG_COLORS.lightgray, linewidth = 1)
        scatter!(p, [left + 0.08], [bottom + h / 2], subplot = canvas, label = "",
            markercolor = 1, markersize = HI_LEGEND_MARKERSIZE, markerstrokewidth = 0)
        canvas_text!(p, canvas, left + w - 0.04, bottom + h / 2, "County", size = FIG_SMALL_FONTSIZE, halign = :right)
    end
    return p
end

"""
    save_household_isolation_figure(rd_b, rd_s, folder)

Saves the figure of the paper (`Fig_13_household_isolation.png`: the plots of
`plot_scenario_household_isolation` above the code listing, a screenshot of
`sc_household_isolation_screenshot_file.jl`) and the plots alone (`Fig_13A_household_isolation.png`)
in `folder`.
"""
function save_household_isolation_figure(rd_b, rd_s, folder)
    plots = figure_image(plot_scenario_household_isolation(rd_b, rd_s))
    save_figure(plots, :household_isolation, folder, panel = "A")

    width = size(plots, 2)
    code = load_panel(joinpath(FIG_PANEL_FOLDER, figure_name(:household_isolation, panel = "B", name = "code_listing") * ".png"))
    return save_figure(stack_sections([
        "A) Scenario comparison in GEMS" => plots,
        "B) Scenario code example" => scale_image(code, width)]), :household_isolation, folder)
end


###
### GEMS READOUTS OF THE BASELINE (combined figure)
###

# The panels are drawn by GEMS (`gemsplot`, `gemsheatmap`, `agsmap`); only their keyword arguments (axes, fonts,
# colors) are set here. Titles, axis labels, and colorbars are on the canvas of the figure.

# layout (inches): height (without the additional room for the rotated labels of the heatmap), left and right edges of the plots and of the maps, heights above (title) and below
# (axis labels) the plots, distances of the axis labels from the plots, colorbar of the heatmap (width, distance
# from the heatmap and of its label), and rows and padding of the legends
const RO_HEIGHT_IN = 7.0
const RO_ROTATED_IN = 0.2
const RO_PLOTS = (0.85, 4.22)
const RO_MAPS = (4.62, 6.97)
const RO_ABOVE_IN = 0.2
const RO_BELOW_IN = 0.35
const RO_TITLE_IN = 0.11
const RO_YLABEL_IN = 0.4
const RO_XLABEL_IN = 0.265
const RO_COLORBAR_IN = 0.1
const RO_COLORBAR_GAP_IN = 0.06
const RO_COLORBAR_LABEL_IN = 0.38
const RO_LEGEND = (row_in = 0.118, padding_in = 0.035)

# days shown in the time series (the epidemic takes place in the first 100) and after an exposure
const RO_DAYS = 100
const RO_EXPOSURE_DAYS = 41

# number as tick label: integers without decimals, else at most two decimals
_ro_num(x) = (r = round(x, digits = 2); isinteger(r) ? string(Int(r)) : string(r))
_ro_millions(x) = _ro_num(x / 1e6)
_ro_percent(x) = _ro_num(100 * x)

# heatmap of the progression categories per age group (as the GEMS plot `ProgressionCategories`)
function plot_progression_categories(rd; plotargs...)
    df = dataframes(rd)["progression_ages"]
    return gemsheatmap(df.age_group, Int.(indexin(df.category, HI_CATEGORIES)), df.count;
        xformatter = a -> "$(a)–$(a + HI_AGE_GROUP - 1)",
        xrotation = 45,
        yformatter = k -> HI_CATEGORIES[k],
        plotargs...)
end

# panels of the left column (title, x label, y label, plot)
function readout_panels(rd)
    style = (
        title = "", plot_title = "", xlabel = "", ylabel = "",
        fontfamily = FIG_FONT,
        tickfontsize = FIG_SMALL_FONTSIZE,
        gridalpha = 0.15, gridlinewidth = 0.5)
    return [
        ("Daily infections", "Days", "Millions",
            gemsplot(rd; type = :TickCases, style...,
                series = [:exposed],
                linewidth = 1,
                xlims = (0, RO_DAYS),
                yformatter = _ro_millions,
                legend = false)),
        ("Effective reproduction number", "Days", "7-day mean",
            gemsplot(rd; type = :EffectiveReproduction, style...,
                linewidth = 1,
                xlims = (0, RO_DAYS),
                ylims = (0, 4),
                legend = false)),
        # attack rate on its full range (0–100 %), a bit below 0 for the dots on 0; household
        # sizes up to 6 (GEMS samples only a few larger households)
        ("In-household attack rate", "Household size", "Attack rate (%)",
            gemsplot(rd; type = :HouseholdAttackRate, style...,
                ar_only = true,
                markersize = 3,
                markerstrokewidth = 0,
                xlims = (1.5, 6.5),
                xticks = 2:6,
                ylims = (-0.05, 1.05),
                yticks = 0:0.2:1,
                yformatter = _ro_percent,
                legend = false)),
        # 10 % room above the peak (GEMS's data of the plot)
        ("Infection duration", "Days", "Infections (%)",
            gemsplot(rd; type = :InfectionDuration, style...,
                linecolor = :transparent,
                xlims = (0, RO_EXPOSURE_DAYS),
                ylims = (0, 1.1 * maximum(aggregated_compartment_periods(rd).total)),
                yformatter = _ro_percent,
                legend = false)),
        ("Cumulative disease progressions", "Days after exposure", "Millions",
            gemsplot(rd; type = :CumulativeDiseaseProgressions, style...,
                xlims = (0, RO_EXPOSURE_DAYS),
                yformatter = _ro_millions,
                legend = false)),
        ("Progression categories", "Age group", "",
            plot_progression_categories(rd; style...,
                color = :viridis,
                colorbar = false,
                framestyle = :box,
                grid = false))]
end

# legends of the panels B and E (labels of the series that GEMS draws => labels in the figure)
const RO_R_LEGEND = [
    "(7-day Rolling) Effective R" => "Overall",
    "(7-day Rolling) Effective R (In Households)" => "In households",
    "(7-day Rolling) Effective R (Outside Households)" => "Outside households"]
const RO_PROGRESSION_LEGEND = [
    "Latent" => "Latent",
    "Pre-Symptomatic" => "Pre-symptomatic",
    "Symptomatic" => "Symptomatic",
    "Asymptomatic" => "Asymptomatic"]

"""
    plot_baseline_readouts(rd)

Combines GEMS readout plots of `rd` (daily infections, effective reproduction number, household
attack rate, infection duration, cumulative disease progressions, progression categories; A–F) with
maps of the initial reproduction number and the infection fatality by county (as the GEMS maps
`R0Map` and `CaseFatalityMap`; G, H) in the width of the page. `rd` must be a ResultData object of the
style `HouseholdIsolationRD`.
"""
function plot_baseline_readouts(rd)
    panels = readout_panels(rd)
    maps = [
        ("Initial reproduction number by county", r0_per_county(rd)),
        ("Infection fatality by county (%)", dataframes(rd)["fatality_per_county"])]
    infections = first(panels)[end]
    first(infections.series_list)[:linecolor] = RGBA(parse(RGB, FIG_COLORS.blue)) # GEMS sets a blue of its own

    # left column: (title) | plot | (axis label and title) | plot | ...; right column: (title) | map | (space and
    # title) | map | (space), the maps next to the first and the fourth plot; empty cells around plots and maps
    n = length(panels)
    width, height = FIG_FULL_WIDTH_IN, RO_HEIGHT_IN + RO_ROTATED_IN
    row_in = RO_HEIGHT_IN / n
    plot_in = row_in - RO_ABOVE_IN - RO_BELOW_IN
    map_in = (RO_MAPS[2] - RO_MAPS[1]) / _FIG1_ASPECT
    split = RO_MAPS[1] - 0.05 # between both columns
    shares(x) = collect(x) ./ sum(x)
    l = grid(1, 2, widths = shares((split, width - split)))
    l[1, 1] = grid(2n + 1, 1, heights = shares([RO_ABOVE_IN; repeat([plot_in, RO_BELOW_IN + RO_ABOVE_IN], n - 1);
        plot_in; RO_BELOW_IN + RO_ROTATED_IN]))
    l[1, 2] = grid(5, 1, heights = shares((RO_ABOVE_IN, map_in, n ÷ 2 * row_in - map_in, map_in,
        height - RO_ABOVE_IN - n ÷ 2 * row_in - map_in)))
    for r in 2:2:2n
        l[1, 1][r, 1] = grid(1, 3, widths = shares((RO_PLOTS[1], RO_PLOTS[2] - RO_PLOTS[1], split - RO_PLOTS[2])))
    end
    for r in (2, 4)
        l[1, 2][r, 1] = grid(1, 3, widths = shares((RO_MAPS[1] - split, RO_MAPS[2] - RO_MAPS[1], width - RO_MAPS[2])))
    end
    items, plots = Any[empty_panel()], Int[]
    function add_items!(p)
        push!(items, empty_panel(), p)
        push!(plots, length(items))
        push!(items, empty_panel(), empty_panel()) # right of it and below it
    end
    foreach(add_items!, last.(panels))
    push!(items, empty_panel())
    foreach(m -> add_items!(plot_county_map(m[2])), maps)
    p = plot(items..., layout = l, size = (100width, 100height), fontfamily = FIG_FONT)
    remove_padding!(p)
    canvas = add_canvas!(p)
    corners(k) = plot_corners(p, plots[k])

    for (k, (title, xlabel, ylabel, _)) in enumerate(panels)
        area = corners(k)
        canvas_title!(p, canvas, area.left, area.top + RO_TITLE_IN, 'A' + k - 1, title)
        canvas_text!(p, canvas, (area.left + area.right) / 2, area.bottom - RO_XLABEL_IN - (k == n ? RO_ROTATED_IN : 0), xlabel)
        canvas_text!(p, canvas, area.left - RO_YLABEL_IN, (area.bottom + area.top) / 2, ylabel, rotation = 90)
    end
    canvas_legend!(p, canvas, plots[2], RO_R_LEGEND; RO_LEGEND...)
    r_threshold_label!(p, canvas, plots[2])
    canvas_legend!(p, canvas, plots[5], RO_PROGRESSION_LEGEND; RO_LEGEND...)

    # colorbar of the heatmap
    area = corners(n)
    counts = dataframes(rd)["progression_ages"].count ./ 1e6
    county_colorbar!(p, canvas, CountyMap(column = :value), extrema(counts),
        area.right + RO_COLORBAR_GAP_IN, area.bottom, RO_COLORBAR_IN, area.top - area.bottom, colors = cgrad(:viridis))
    canvas_text!(p, canvas, area.right + RO_COLORBAR_LABEL_IN, (area.bottom + area.top) / 2, "Millions", rotation = 90)

    for (k, (title, df)) in enumerate(maps)
        canvas_title!(p, canvas, corners(n + k).left, corners(n + k).top + RO_TITLE_IN, 'A' + n + k - 1, title)
        map_colorbar!(p, canvas, plots[n + k], extrema(df[!, 2]))
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
`save_household_isolation_figure` and `plot_baseline_readouts`).
"""
function replot_household_isolation(folder)
    rd_b, rd_s = [import_resultdata(joinpath(folder, f, "resultdata.jld2"))
        for f in (HI_BASELINE_FOLDER, HI_SCENARIO_FOLDER)]

    save_household_isolation_figure(rd_b, rd_s, folder)

    save_figure(plot_baseline_readouts(rd_b), :baseline_readouts, folder)
end
