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

    county_map = plot_county_map(select(df, :ags, :r0_diff); style...,
        title = "$_RINIT_REDUCTION_LABEL by county")

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

# map of the values of the counties in the second column of `df` (GEMS), without its colorbar
plot_county_map(df; plotargs...) = agsmap(df;
    fillcolor = :inferno,
    clims = extrema(df[!, 2]),
    colorbar = false, # drawn on the canvas (thinner)
    legend = false,
    linecolor = nothing,
    framestyle = :none,
    xlims = HI_MAP_XLIMS,
    ylims = _FIG1_YLIMS,
    widen = false,
    plotargs...)

# colorbar of the map in subplot `k` on the `canvas` of the figure `p` (laid out), right of the map
function map_colorbar!(p, canvas, k, lims)
    area = plot_corners(p, k, aspect = (HI_MAP_XLIMS[2] - HI_MAP_XLIMS[1]) / (_FIG1_YLIMS[2] - _FIG1_YLIMS[1]))
    x = area.left + (area.right - area.left) / HI_MAP_ROOM + HI_COLORBAR_GAP_IN
    return county_colorbar!(p, canvas, CountyMap(column = :value, title = "", note = ""), lims,
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

# colorbar of the map (subplot 1), label of R = 1 (subplot 6), and legends of the scatter plots (lower left corner of the first, upper left
# corner of the second, where there are no points; subplots 2 and 3) on a transparent canvas in inches above the figure
function draw_legends!(p, df)
    Plots.prepare_output(p) # lays out the figure
    by_age, by_size = plot_corners(p, 2), plot_corners(p, 3)

    canvas = add_canvas!(p)
    map_colorbar!(p, canvas, 1, extrema(df.r0_diff))
    r_threshold_label!(p, canvas, 6, halign = :right)

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

# The panels are drawn by GEMS (`gemsplot`, `gemsheatmap`, `agsmap`); only their keyword arguments (axes, fonts,
# colors) are set here. Titles, axis labels, and colorbars are on the canvas of the figure.

# share of the width for the maps (without labels), distances of the titles and axis labels from the plot areas, and colorbar of the heatmap (width, distance from the
# heatmap and of its label; inches)
const RO_MAP_SHARE = 0.57
const RO_TITLE_IN = 0.12
const RO_YLABEL_IN = 0.4
const RO_XLABEL_IN = 0.28
const RO_COLORBAR_IN = 0.1
const RO_COLORBAR_GAP_IN = 0.07
const RO_COLORBAR_LABEL_IN = 0.38

# days shown in the time series (the epidemic takes place in the first 100) and after an exposure
const RO_DAYS = 100
const RO_EXPOSURE_DAYS = 41

# number as tick label: integers without decimals, else at most two decimals
_ro_num(x) = (r = round(x, digits = 2); isinteger(r) ? string(Int(r)) : string(r))
_ro_millions(x) = _ro_num(x / 1e6)
_ro_percent(x) = _ro_num(100 * x)

# heatmap of the progression categories per age group (as the GEMS plot `ProgressionCategories`), every fourth
# age group labeled
function plot_progression_categories(rd; plotargs...)
    df = dataframes(rd)["progression_ages"]
    return gemsheatmap(df.age_group, Int.(indexin(df.category, HI_CATEGORIES)), df.count;
        xformatter = a -> a % (4 * HI_AGE_GROUP) == 0 ? "$(a)–$(a + HI_AGE_GROUP - 1)" : "",
        yformatter = k -> HI_CATEGORIES[k],
        plotargs...)
end

# panels of the left column (title, x label, y label, plot)
function readout_panels(rd)
    style = (
        title = "", plot_title = "", xlabel = "", ylabel = "",
        fontfamily = FIG_FONT,
        tickfontsize = FIG_SMALL_FONTSIZE,
        gridalpha = 0.15, gridlinewidth = 0.5,
        margin = 0mm, top_margin = 2.5mm, bottom_margin = 2.5mm, right_margin = 10mm)
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
`R0Map` and `CaseFatalityMap`; G, H) on a full page of PLOS. `rd` must be a ResultData object of the
style `HouseholdIsolationRD`.
"""
function plot_baseline_readouts(rd)
    panels = readout_panels(rd)
    maps = [
        ("Initial reproduction number by county", r0_per_county(rd)),
        ("Infection fatality by county (%)", dataframes(rd)["fatality_per_county"])]
    lims(df) = extrema(df[!, 2])

    infections = first(panels)[end]
    plot!(infections, top_margin = 4mm) # room for the first title
    first(infections.series_list)[:linecolor] = RGBA(parse(RGB, FIG_COLORS.blue)) # GEMS sets a blue of its own
    n = length(panels)
    l = grid(1, 2, widths = [1 - RO_MAP_SHARE, RO_MAP_SHARE])
    l[1, 1] = grid(n, 1)
    l[1, 2] = grid(2, 1)
    p = plot(last.(panels)...,
        [plot_county_map(df; ticks = false, margin = 0mm, top_margin = 4mm) for (_, df) in maps]...,
        layout = l,
        size = (100FIG_FULL_WIDTH_IN, 100FIG_MAX_HEIGHT_IN),
        fontfamily = FIG_FONT)
    Plots.prepare_output(p) # lays out the figure
    canvas = add_canvas!(p)

    for (k, (title, xlabel, ylabel, _)) in enumerate(panels)
        area = plot_corners(p, k)
        canvas_title!(p, canvas, area.left, area.top + RO_TITLE_IN, 'A' + k - 1, title)
        canvas_text!(p, canvas, (area.left + area.right) / 2, area.bottom - RO_XLABEL_IN, xlabel)
        canvas_text!(p, canvas, area.left - RO_YLABEL_IN, (area.bottom + area.top) / 2, ylabel, rotation = 90)
    end

    canvas_legend!(p, canvas, 2, RO_R_LEGEND)
    r_threshold_label!(p, canvas, 2)
    canvas_legend!(p, canvas, 5, RO_PROGRESSION_LEGEND)

    # colorbar of the heatmap
    area = plot_corners(p, n)
    counts = dataframes(rd)["progression_ages"].count ./ 1e6
    county_colorbar!(p, canvas, CountyMap(column = :value, title = "", note = ""), extrema(counts),
        area.right + RO_COLORBAR_GAP_IN, area.bottom, RO_COLORBAR_IN, area.top - area.bottom, colors = cgrad(:viridis))
    canvas_text!(p, canvas, area.right + RO_COLORBAR_LABEL_IN, (area.bottom + area.top) / 2, "Millions", rotation = 90)

    # titles of the maps next to those of A and D
    for (k, (title, df)) in enumerate(maps)
        area = plot_corners(p, 1 + (k - 1) * n ÷ 2)
        canvas_title!(p, canvas, plot_corners(p, n + k).left, area.top + RO_TITLE_IN, 'A' + n + k - 1, title)
        map_colorbar!(p, canvas, n + k, lims(df))
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
