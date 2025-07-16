###
### OUTBREAK SCENARIO
###

export plot_weekly_incidence
export plot_case_fatality

export scenario_outbreak


"""
    plot_weekly_incidence(rd; weeks = [10, 18, 20, 24, 26, 30])

Plots weekly incidence for a provided ResultData object.
The `weeks` argument takes a vector of Integers to specify which
weeks should be plotted. In the current plot layout, 
it should have 6 entries.
"""
function plot_weekly_incidence(rd; weeks = [10, 18, 20, 24, 26, 30])

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
        colorbar_title = "Percentage (%)", 
        fontfamily = "Times Roman")
end


"""
    scenario_outbreak(;
        configfile = "models/DE_model_R0_3.28_Heinsberg_Seed.toml",
        popfile = "data/people.jld2",
        settingsfile = "data/settings.jld2",
        folder = nothing)

Runs the outbreak scenario and all analyses.
"""
function scenario_outbreak(;
    configfile = "models/DE_model_R0_3.28_Heinsberg_Seed.toml",
    popfile = "data/people.jld2",
    settingsfile = "data/settings.jld2",
    folder = nothing)


    # load, run, and post-process simulation
    sim = Simulation(configfile, popfile, settingsfile)
    run!(sim)
    rd = ResultData(sim)

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder() : nothing    

    # weekly incidence
    p_wi = plot_weekly_incidence(rd)
    png(p_wi, joinpath(folder, "sc_outbreak_weekly_incidence.png"))

    p_cf = plot_case_fatality(rd)
    png(p_cf, joinpath(folder, "sc_outbreak_case_fatality.png"))

end