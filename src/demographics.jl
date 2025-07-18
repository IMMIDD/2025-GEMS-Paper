export individual_data
export plot_age_per_state
export plot_agemap
export plot_household_sizes
export plot_pop_density
export plot_singles
export plot_elderly
export plot_kids
export plot_combined_demographics

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
            xlabelfontsize = map_legendfontsize()),
        plot!(p_mean_hh_size,
            plot_title = "",
            title = "",
            xlabel = "Mean Household Size",
            xlabelfontsize = map_legendfontsize()),
        plot!(p_singles,
            plot_title = "",
            title = "",
            xlabel = "Single-Person Households",
            xlabelfontsize = map_legendfontsize()),
        plot!(p_mean_age,
            plot_title = "",
            title = "",
            xlabel = "Mean Age",
            xlabelfontsize = map_legendfontsize()),
        plot!(p_kids,
            plot_title = "",
            title = "",
            xlabel = "Kids (Age 0-14)",
            xlabelfontsize = map_legendfontsize()),
        plot!(p_elderly, 
            plot_title = "",
            title = "",
            xlabel = "Elderly (Age 65+)",
            xlabelfontsize = map_legendfontsize()),
        layout = (2, 3),
        top_margin = -50mm,
        right_margin = 10mm,
        bottom_margin = -50mm,
        left_margin = -55mm,
        xlabelfontsize = map_legendfontsize() * 2.5,
        size = (4 * map_dims()[1], 2 * map_dims()[2]))
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
        sim = Simulation(population = "SL"))

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
end