export individual_data
export load_reference
export plot_age_per_state
export plot_agemap
export plot_household_sizes
export plot_pop_density
export plot_singles
export plot_elderly
export plot_kids
export plot_combined_demographics
export plot_population_validation

export demographics


# DATA PREPARATION

function individual_data(sim)
    return individuals(sim) |>
        inds -> DataFrame(
            id = id.(inds),
            age = age.(inds),
            sex = sex.(inds),
            ags = (i -> ags(household(i, sim))).(inds)
        )
end


# POPULATION VALIDATION: REFERENCE DATA

const REF_DIR = "reference_data"

"""
    load_reference(key)

Loads the reference table `key` from `reference_data/sources.toml`
and returns the data and its metadata.
"""
function load_reference(key)
    meta = TOML.parsefile(joinpath(REF_DIR, "sources.toml"))[key]
    meta["placeholder"] && @warn "Reference table '$key' contains placeholder data."
    df = CSV.read(joinpath(REF_DIR, meta["file"]), DataFrame)
    return df, meta
end


# POPULATION VALIDATION: MODEL DATA

# Gesyland models care homes as households of 90+ residents (sizes 10-89
# do not occur); official household statistics cover private households only
const PRIVATE_HH_MAX = 9

state_code(a) = id(a) ÷ 10^6

"""
    bin_counts(values, counts, groups)

Sums `counts` per group in `groups` (columns `group_min`, `group_max`),
assigning each entry by its `values`.
"""
function bin_counts(values, counts, groups)
    return map(eachrow(groups)) do g
        hi = ismissing(g.group_max) ? Inf : g.group_max
        sum(counts[g.group_min .<= values .<= hi]; init = 0)
    end
end

household_df(sim) = households(sim) |>
    hhs -> filter(h -> size(h) <= PRIVATE_HH_MAX, hhs) |>
    hhs -> DataFrame(state = state_code.(ags.(hhs)), value = size.(hhs), count = size.(hhs))

age_df(individual_data_df) = individual_data_df |>
    df -> DataFrame(state = state_code.(df.ags), value = df.age, count = 1)

# schools contain school years, which contain classes; a school's size
# is the sum of its classes' sizes (schools include daycare centres and
# higher education, as Gesyland does not label institution types)
school_size_df(sim) = innerjoin(
        DataFrame(year = contained.(schoolclasses(sim)), size = size.(schoolclasses(sim))),
        DataFrame(year = id.(schoolyears(sim)), school = contained.(schoolyears(sim))),
        on = :year) |>
    df -> groupby(df, :school) |>
    df -> combine(df, :size => sum => :value) |>
    df -> filter(r -> r.value > 0, df) |>
    df -> insertcols(df, :count => 1)

# workplaces contain departments, which contain offices;
# a workplace's size is the sum of its offices' sizes
workplace_size_df(sim) = innerjoin(
        DataFrame(dept = contained.(offices(sim)), size = size.(offices(sim))),
        DataFrame(dept = id.(departments(sim)), wp = contained.(departments(sim))),
        on = :dept) |>
    df -> groupby(df, :wp) |>
    df -> combine(df, :size => sum => :value) |>
    df -> insertcols(df, :count => 1)


# POPULATION VALIDATION: COMPARISON

"""
    compare(model, reference)

Joins model and reference counts per group and converts both to shares.
"""
function compare(model, reference)
    groups = select(reference, :group_min, :group_max)
    model_counts = bin_counts(model.value, model.count, groups)
    return DataFrame(
        label = group_label.(groups.group_min, groups.group_max),
        model = sum(model_counts) > 0 ? model_counts ./ sum(model_counts) : fill(NaN, nrow(groups)),
        reference = reference.count ./ sum(reference.count))
end

group_label(lo, hi) = ismissing(hi) ? "$lo+" : lo == hi ? "$lo" : "$(lo)–$(hi)"

# total variation distance: share of persons (or institutions) that would
# have to change groups for the model to match the reference
tvd(df) = sum(abs.(df.model .- df.reference)) / 2


# PLOTTING

function plot_age_per_state(individual_data_df)

    age_group_5(age) = Int(floor(age / 5))
    age_group_10(age) = Int(floor(age / 10))

    function buildplot(ages)

        return ages |>
            vec -> DataFrame(age = ages) |>
            df -> transform(df, :age => (a -> age_group_10.(a)) => :age_group) |>
            df -> groupby(df, :age_group) |>
            df -> combine(df, nrow => :cnt) |>
            df -> transform(df, :cnt => (c -> c ./ sum(c)) => :fraction) |>
            df -> bar(df.age_group, df.fraction,
                xrotation = 45,
                xtick = 0:nrow(df)-1,
                xformatter = (x -> "$(Int(10x))-$(Int(10(x+1)-1))"),
                yformatter = (y -> "$(round(100y, digits = 2)) %"),
                linecolor = nothing,
                fontfamily = "Times Roman",
                ylims = (0, .2),
                tickfontsize = 16,
                titlefontsize = 20,
                legend = nothing)
    end

    state_names = germanshapes(1) |>
        shps -> DataFrame(ags = AGS.(shps.AGS_0), gen = shps.GEN) |>
        # some shapefiles have mulitple entries
        df -> groupby(df, :ags) |>
        df -> combine(df, :gen => first => :gen)

    return individual_data_df |>
        df -> transform(df, :ags => (a -> state.(a)) => :ags) |>
        df -> groupby(df, :ags) |>
        df -> combine(df, :age => buildplot => :age_plot) |>
        # add state names
        df -> leftjoin(df, state_names, on = :ags) |>
        df -> transform(df, [:age_plot, :gen] => ByRow((plt, name) -> plot!(plt, title = name)) => :age_plot) |>
        df -> plot((df.age_plot)...,
            layout = (4, 4),
            size = (2000, 2000))
end


function plot_agemap(sim)
    return gemsmap(sim, type = :AgeMap,
        level = 2,
        size = map_dims(),
        plot_titlefontsize = map_titlefontsize(),
        titlefontsize = map_titlefontsize(),
        ytickfontsize = map_legendfontsize(),
        colorbar_titlefontsize = map_legendfontsize(),
        right_margin=5mm,
        linecolor = nothing,
        fit_lims = true,
        plot_title = "Mean Age Per Region",
        colorbar_title = "Mean Age", 
        fontfamily = "Times Roman")
end

function plot_household_sizes(sim)
    return gemsmap(sim, type = :HouseholdSizeMap,
        level = 2,
        size = map_dims(),
        plot_titlefontsize = map_titlefontsize(),
        titlefontsize = map_titlefontsize(),
        ytickfontsize = map_legendfontsize(),
        colorbar_titlefontsize = map_legendfontsize(),
        right_margin=5mm,
        linecolor = nothing,
        fit_lims = true,
        colorbar_title = "Mean Household Size", 
        fontfamily = "Times Roman")
end

function plot_pop_density(sim)
    return gemsmap(sim, type = :PopDensityMap,
        level = 2,
        size = map_dims(),
        plot_titlefontsize = map_titlefontsize(),
        titlefontsize = map_titlefontsize(),
        ytickfontsize = map_legendfontsize(),
        colorbar_titlefontsize = map_legendfontsize(),
        right_margin=5mm,
        linecolor = nothing,
        fit_lims = true,
        #colorbar_title = "People per km²", 
        fontfamily = "Times Roman")
end

function plot_singles(sim)
    return gemsmap(sim, type = :SinglesMap,
        level = 2,
        size = map_dims(),
        plot_titlefontsize = map_titlefontsize(),
        titlefontsize = map_titlefontsize(),
        ytickfontsize = map_legendfontsize(),
        colorbar_titlefontsize = map_legendfontsize(),
        right_margin=5mm,
        linecolor = nothing,
        fit_lims = true,
        colorbar_title = "Percentage (%)", 
        fontfamily = "Times Roman")
end

function plot_elderly(sim)
    return gemsmap(sim, type = :ElderlyMap,
        level = 2,
        size = map_dims(),
        plot_titlefontsize = map_titlefontsize(),
        titlefontsize = map_titlefontsize(),
        ytickfontsize = map_legendfontsize(),
        colorbar_titlefontsize = map_legendfontsize(),
        right_margin=5mm,
        linecolor = nothing,
        fit_lims = true,
        colorbar_title = "Percentage (%)", 
        fontfamily = "Times Roman")
end

function plot_kids(sim)
    return gemsmap(sim, type = :KidsMap,
        level = 2,
        size = map_dims(),
        plot_titlefontsize = map_titlefontsize(),
        titlefontsize = map_titlefontsize(),
        ytickfontsize = map_legendfontsize(),
        colorbar_titlefontsize = map_legendfontsize(),
        right_margin=5mm,
        linecolor = nothing,
        fit_lims = true,
        colorbar_title = "Percentage (%)", 
        fontfamily = "Times Roman")
end

function plot_combined_demographics(
    p_pop_density, p_mean_hh_size, p_singles,
    p_mean_age, p_kids, p_elderly)

    return plot(
        plot!(p_pop_density,
            plot_title = "",
            title = "",
            xlabel = "Population Density",
            ytickfontsize = map_legendfontsize()*1.5,
            colorbar_titlefontsize = map_legendfontsize()*1.5,
            xlabelfontsize = map_legendfontsize()),
        plot!(p_mean_hh_size,
            plot_title = "",
            title = "",
            xlabel = "Mean Household Size",
            ytickfontsize = map_legendfontsize()*1.5,
            colorbar_titlefontsize = map_legendfontsize()*1.5,
            xlabelfontsize = map_legendfontsize()),
        plot!(p_singles,
            plot_title = "",
            title = "",
            xlabel = "Single-Person Households",
            ytickfontsize = map_legendfontsize()*1.5,
            colorbar_titlefontsize = map_legendfontsize()*1.5,
            xlabelfontsize = map_legendfontsize()),
        plot!(p_mean_age,
            plot_title = "",
            title = "",
            xlabel = "Mean Age",
            ytickfontsize = map_legendfontsize()*1.5,
            colorbar_titlefontsize = map_legendfontsize()*1.5,
            xlabelfontsize = map_legendfontsize()),
        plot!(p_kids,
            plot_title = "",
            title = "",
            xlabel = "Kids (Age 0-14)",
            ytickfontsize = map_legendfontsize()*1.5,
            colorbar_titlefontsize = map_legendfontsize()*1.5,
            xlabelfontsize = map_legendfontsize()),
        plot!(p_elderly, 
            plot_title = "",
            title = "",
            xlabel = "Elderly (Age 65+)",
            ytickfontsize = map_legendfontsize()*1.5,
            colorbar_titlefontsize = map_legendfontsize()*1.5,
            xlabelfontsize = map_legendfontsize()),
        layout = (2, 3),
        top_margin = -50mm,
        right_margin = 10mm,
        bottom_margin = -50mm,
        left_margin = -55mm,
        xlabelfontsize = map_legendfontsize() * 2.5,
        size = (4 * map_dims()[1], 2 * map_dims()[2]))
end


# POPULATION VALIDATION: PLOTTING

const MODEL_COLOR = gemscolors(1)[1]
const REFERENCE_COLOR = :black

const STATE_NAMES = [
    "Schleswig-Holstein", "Hamburg", "Lower Saxony", "Bremen",
    "North Rhine-Westphalia", "Hesse", "Rhineland-Palatinate", "Baden-Württemberg",
    "Bavaria", "Saarland", "Berlin", "Brandenburg",
    "Mecklenburg-Vorpommern", "Saxony", "Saxony-Anhalt", "Thuringia"]

const TICK_FONTSIZE = 42
const LABEL_FONTSIZE = 48
const HEADER_FONTSIZE = 62

"""
    comparison_panel(df; title, ymax, show_x = true, show_y = true,
        legend = false, xrotation = 0, titlefontsize = LABEL_FONTSIZE,
        xlabel = "", xgap = 0.024, ygap = 0.04, xlabel_offset = 0.3, kwargs...)

Bars for the synthetic population, diamonds for the reference, and the
total variation distance in the top-right corner.

Tick labels and the x-axis label are drawn as annotations, as GR places
its own labels at a fixed distance relative to the figure height, which
leaves large gaps in tall figures. `xgap` and `ygap` set the distance
between the axes and their tick labels (as share of the axis span),
`xlabel_offset` the distance of the x-axis label below the x-axis
(as share of `ymax`). The panel's margins need to make room for them.
"""
function comparison_panel(df; title, ymax, show_x = true, show_y = true,
        legend = false, xrotation = 0, titlefontsize = LABEL_FONTSIZE,
        xlabel = "", xgap = 0.024, ygap = 0.04, xlabel_offset = 0.3, kwargs...)

    n = nrow(df)
    x = 1:n
    xmin, xmax = 0.4, n + 0.6
    yticks = 0:(ymax > 0.3 ? 0.1 : 0.05):ymax

    p = bar(x, df.model,
        label = "Gesyland",
        fillcolor = MODEL_COLOR,
        linewidth = 0,
        bar_width = 0.7,
        title = title,
        titlefontsize = titlefontsize,
        xticks = (x, fill("", n)),
        yticks = (yticks, fill("", length(yticks))),
        ylims = (0, ymax),
        xlims = (xmin, xmax),
        grid = :y,
        gridalpha = 0.2,
        legend = legend,
        legendfontsize = LABEL_FONTSIZE,
        foreground_color_legend = nothing,
        # white backdrop hides the grid lines behind the legend
        background_color_legend = :white,
        fontfamily = "Times Roman";
        kwargs...)

    scatter!(p, x, df.reference,
        label = "Reference data",
        markershape = :diamond,
        markersize = 15,
        markercolor = REFERENCE_COLOR,
        markerstrokewidth = 0)

    dev = all(isnan, df.model) ? "no model data" : "TVD $(round(100tvd(df), digits = 1)) %"
    annotate!(p, n + 0.5, 0.86ymax, text(dev, TICK_FONTSIZE, :right, "Times Roman"))

    # tick labels
    if show_y
        for y in yticks
            annotate!(p, xmin - xgap * (xmax - xmin), y,
                text("$(round(Int, 100y)) %", TICK_FONTSIZE, :right, :vcenter, "Times Roman"))
        end
    end
    if show_x
        for (xi, lab) in zip(x, df.label)
            # rotated labels end at their bar, unrotated ones are centered below it
            annotate!(p, xrotation == 0 ? xi : xi + 0.15, -ygap * ymax,
                text(lab, TICK_FONTSIZE, xrotation == 0 ? :center : :right, :top,
                    "Times Roman", rotation = xrotation))
        end
    end
    xlabel != "" && annotate!(p, (xmin + xmax) / 2, -xlabel_offset * ymax,
        text(xlabel, LABEL_FONTSIZE, :center, :top, "Times Roman"))

    return p
end

# figure size in pixels (fonts in the same scale as the other figures)
# and the left margin that holds the y tick labels
const FIG_SIZE = (3000, 5100)
const FIG_WIDTH = FIG_SIZE[1] * Plots.px
const YLABEL_MARGIN = 30mm

"""
    state_grid(model, reference; xrotation = 0, bottom_margin = 20mm)

4×4 grid of comparison panels, one per state, with a shared y-axis.
`bottom_margin` makes room for the x tick labels below the last row.
"""
function state_grid(model, reference; xrotation = 0, bottom_margin = 20mm)
    comps = [compare(model[model.state .== s, :], reference[reference.state .== s, :]) for s in 1:16]
    ymax = 1.15 * maximum(c -> maximum(filter(!isnan, [c.model; c.reference])), comps)

    gap = 2mm
    panels = map(1:16) do s
        comparison_panel(comps[s],
            title = STATE_NAMES[s],
            ymax = ymax,
            show_x = s > 12,
            show_y = s % 4 == 1,
            xrotation = xrotation,
            # only the first column needs room for y tick labels
            left_margin = s % 4 == 1 ? YLABEL_MARGIN : gap,
            right_margin = gap,
            # pull the title closer to the data and the rows closer together
            top_margin = -12mm,
            bottom_margin = s > 12 ? bottom_margin : 0mm)
    end

    # widen the first column so that all panels have the same plot width
    plot_width = (FIG_WIDTH - YLABEL_MARGIN - 7gap) / 4
    w_first = (plot_width + YLABEL_MARGIN + gap) / FIG_WIDTH
    w_other = (1 - w_first) / 3
    return plot(panels...,
        layout = grid(4, 4, widths = [w_first, w_other, w_other, w_other]),
        fontfamily = "Times Roman")
end

"""
    header(txt)

Empty subplot carrying a block title.
"""
header(txt) = plot(framestyle = :none, legend = false, xlims = (0, 1), ylims = (0, 1),
        bottom_margin = -12mm) |>
    p -> annotate!(p, 0.5, 0, text(txt, HEADER_FONTSIZE, :center, :bottom, "Times Roman"))

"""
    plot_population_validation(sim; individual_data_df = individual_data(sim))

Compares the simulation's population with the reference statistics in
`reference_data/`: household sizes and age structure per state, and
workplace and school sizes nationally. Pass `individual_data_df` to reuse
an already computed `individual_data(sim)`.
"""
function plot_population_validation(sim; individual_data_df = individual_data(sim))

    hh_ref, _ = load_reference("household_sizes")
    age_ref, _ = load_reference("age_structure")
    wp_ref, _ = load_reference("workplace_sizes")
    sc_ref, _ = load_reference("school_sizes")

    p_hh = state_grid(household_df(sim), hh_ref, bottom_margin = 16mm)
    p_age = state_grid(age_df(individual_data_df), age_ref, xrotation = 45, bottom_margin = 30mm)

    national = (
        titlefontsize = HEADER_FONTSIZE,
        ymax = 0.5,
        xrotation = 30,
        xgap = 0.012,
        ygap = 0.025,
        xlabel_offset = 0.22,
        left_margin = YLABEL_MARGIN,
        bottom_margin = 45mm)
    p_wp = comparison_panel(compare(workplace_size_df(sim), wp_ref);
        title = "Workplace Sizes (National)",
        xlabel = "Employees per workplace",
        national...)
    p_sc = comparison_panel(compare(school_size_df(sim), sc_ref);
        title = "School Sizes (National)",
        xlabel = "Attendees per school",
        # figure legend, below the TVD in the top-right corner
        legend = (0.72, 0.72),
        national...)

    return plot(
        header("Household Sizes (by State)"), p_hh,
        header("Age Structure (by State)"), p_age,
        p_wp, p_sc,
        layout = @layout([
            a{0.012h}
            b{0.394h}
            c{0.012h}
            d{0.394h}
            [e f]
        ]),
        size = FIG_SIZE,
        fontfamily = "Times Roman")
end


# RUN ALL ANALYSES

"""
    demographics(;
        folder = nothing,
        sim = Simulation(population = "SL"))

Runs all mobilty-related analyses.
"""
function demographics(;
        folder = nothing,
        sim = Simulation(population = "SL", seed = 0))

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder() : nothing

    # calculate individual data
    id = individual_data(sim)

    # plot age structure per state
    p_age_state = plot_age_per_state(id)
    png(p_age_state, joinpath(folder, "age_structure_per_state.png"))

    # plot mean age
    p_mean_age = plot_agemap(sim)   
    png(p_mean_age, joinpath(folder, "mean_age_map.png"))

    # plot mean household size
    p_mean_hh_size = plot_household_sizes(sim)   
    png(p_mean_hh_size, joinpath(folder, "mean_household_size_map.png"))

    # plot population density
    p_pop_density = plot_pop_density(sim)   
    png(p_pop_density, joinpath(folder, "pop_density_map.png"))

    # plot fraction of singles
    p_singles = plot_singles(sim)   
    png(p_singles, joinpath(folder, "fraction_of_single_person_households.png"))

    # plot fraction of elderly
    p_elderly = plot_elderly(sim)
    png(p_elderly, joinpath(folder, "fraction_of_elderly.png"))

    # plot fraction of kids
    p_kids = plot_kids(sim)
    png(p_kids, joinpath(folder, "fraction_of_kids.png"))

    p_combined = plot_combined_demographics(
        p_pop_density, p_mean_hh_size, p_singles,
        p_mean_age, p_kids, p_elderly)
    png(p_combined, joinpath(folder, "demographic_maps.png"))

    # compare population with reference statistics
    p_validation = plot_population_validation(sim, individual_data_df = id)
    png(p_validation, joinpath(folder, "population_validation.png"))
    savefig(p_validation, joinpath(folder, "population_validation.pdf"))
end