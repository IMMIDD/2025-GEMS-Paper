###
### OUTBREAK SCENARIO
###

export run_outbreak_sim
export plot_weekly_incidence
export plot_case_fatality
export plot_r0_by_region
export plot_combined_outbreak_plots

export scenario_outbreak


# RUN SIMULATION


"""
    run_outbreak_sim(;
        configfile = "models/DE_model_R0_3.28_Heinsberg_Seed.toml",
        population = "DE")

Runs the outbreak scenario simulation and returns the ResultData object. 
"""
function run_outbreak_sim(;
    configfile = "models/DE_model_R0_3.28_Heinsberg_Seed.toml",
    population = "DE")

    sim = Simulation(configfile=configfile, population=population, seed = 0)
    run!(sim)
    rd = ResultData(sim)

   return rd
end


"""
    plot_weekly_incidence(rd; weeks = [10, 18, 20, 24, 26, 30])

Plots weekly incidence for a provided ResultData object.
The `weeks` argument takes a vector of Integers to specify which
weeks should be plotted. In the current plot layout, 
it should have 6 entries.
"""
function plot_weekly_incidence(rd; weeks = [16, 18, 20, 22, 25, 30])
    # wrapper function for gemsmap
    plotweek(rd, week) = gemsmap(rd, type = :WeeklyIncidenceMap, week = week,
        linecolor = nothing,
        titlefontsize = map_titlefontsize(),
        ytickfontsize = map_legendfontsize()*2,
        colorbar_titlefontsize = map_legendfontsize(),
        plot_title = "",
        title = "",
        xlabel = "Week $week",
        xlabelfontsize = map_legendfontsize(),
        size = map_dims(),
    )

    plts = (w -> plotweek(rd, w)).(weeks)

    # build multiplot
    return plot(
        plts[1:6]...,
        layout = (2,3),
        top_margin = -50mm,
        right_margin = 10mm,
        bottom_margin = -50mm,
        left_margin = -80mm,
        xlabelfontsize = map_legendfontsize() * 2.5,
        size = (4 * map_dims()[1], 2 * map_dims()[2])
    )
end


"""
    plot_case_fatality(rd)

Plots case fatality map.
"""
function plot_case_fatality(rd)
    return gemsmap(rd, type = :CaseFatalityMap, level = 2,
        size = map_dims(),
        titlefontsize = map_titlefontsize(),
        ytickfontsize = map_legendfontsize(),
        colorbar_titlefontsize = map_legendfontsize(),
        right_margin=5mm,
        linecolor = nothing,
        tite = "Infection Fatality",
        colorbar_title = "Percentage (%)", 
        fontfamily = "Times Roman")
end


"""
    plot_r0_by_region(rd)

Plots R0 by region map.
"""
function plot_r0_by_region(rd)
    return gemsmap(rd, type = :R0Map, level = 2,
        size = map_dims(),
        titlefontsize = map_titlefontsize(),
        ytickfontsize = map_legendfontsize(),
        colorbar_titlefontsize = map_legendfontsize(),
        right_margin=5mm,
        linecolor = nothing,
        title ="Basic Reproduction Number (R0)",
        plot_title = "",
        colorbar_title = "Basic Reproduction Number (R0)", 
        fontfamily = "Times Roman")
end


function plot_gemsplots(rd)
    return plot(
            gemsplot(rd, type = :TickCases,
                series = [:exposed, :deaths],
                titlefontsize = map_legendfontsize(),
                title = "Daily Infections",
                plot_title = "",
                linewidth = 3,
                top_margin = 10mm),
            gemsplot(rd, type = :EffectiveReproduction,
                linewidth = 3,
                top_margin = -20mm),
            gemsplot(rd, type = :HouseholdAttackRate,
                ar_only = true,
                markersize = 15,
                xlims = (1, 10),
                top_margin = -20mm),
            gemsplot(rd, type = :InfectionDuration,
                linecolor = :transparent,
                xlims = (0, 41),
                top_margin = -20mm),
            gemsplot(rd, type = :CumulativeDiseaseProgressions,
                xlims = (0, 41),
                top_margin = -20mm),
            gemsplot(rd, type = :ProgressionCategories,
                top_margin = -20mm),

            layout = (6, 1),
            tickfontsize = map_legendfontsize(),
            legendfontsize = map_legendfontsize(),
            labelfontsize = map_legendfontsize(),
            titlefontsize = map_titlefontsize(),
            left_margin = 10mm,
            bottom_margin = 20mm
        )
        
        # return gemsplot(rd, type = (
        #     :TickCases,
        #     :EffectiveReproduction,
        #     :HouseholdAttackRate,
        #     :InfectionDuration,
        #     :CumulativeDiseaseProgressions,
        #     :ProgressionCategories,
        #     ),
        #     ar_only = true,
        #     tickfontsize = map_legendfontsize(),
        #     legendfontsize = map_legendfontsize(),
        #     labelfontsize = map_legendfontsize(),
        #     titlefontsize = map_titlefontsize(),
        #     left_margin = 10mm,
        #     top_margin = -20mm,
        #     bottom_margin = 20mm
        # )
end

function plot_combined_outbreak_plots(rd, r0map, fatalitymap)

    l = @layout [a{0.45w} b{0.55w}]

    return plot(
        plot_gemsplots(rd),
        plot(
            plot!(r0map,
                top_margin = -60mm,
                bottom_margin = -30mm),
            plot!(fatalitymap,
                top_margin = -30mm,
                bottom_margin = -50mm),
            left_margin = -40mm,
            layout = (2, 1),
            size = ((map_dims()[1], 2 * map_dims()[2]))
        ),

        layout = l,
        size = 2 .* map_dims(),

        titlefontsize = map_titlefontsize(),
        fontfamily = "Times Roman", 
    )


end


"""
    scenario_outbreak(;
        configfile = "models/DE_model_R0_3.28_Heinsberg_Seed.toml",
        population = "DE",
        folder = nothing)

Runs the outbreak scenario and all analyses.
"""
function scenario_outbreak(;
    configfile = "models/DE_model_R0_3.28_Heinsberg_Seed.toml",
    population = "DE",
    folder = nothing)


    # load, run, and post-process simulation
    rd = run_outbreak_sim(
        configfile = configfile,
        population = population)

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder() : nothing    
    weekly_folder = joinpath(folder, "weekly_incidence")
    mkpath(weekly_folder)

    # weekly incidence
    p_wi = plot_weekly_incidence(rd)
    png(p_wi, joinpath(folder, "sc_outbreak_weekly_incidence.png"))

    #r0
    r0map = plot_r0_by_region(rd)
    png(r0map, joinpath(folder, "sc_outbreak_r0_by_region.png"))

    # case fatality
    p_cf = plot_case_fatality(rd)
    png(p_cf, joinpath(folder, "sc_outbreak_case_fatality.png"))

    # combined plots
    p_combined = plot_combined_outbreak_plots(rd, r0map, p_cf)
    png(p_combined, joinpath(folder, "sc_outbreak_combined.png"))
end