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
        title = "Infection Fatality",
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
end


"""
    plot_combined_readouts(rd, r0map, fatalitymap)

Combines the GEMS readout plots of `rd` (daily infections, effective R, household
attack rate, infection duration, cumulative disease progressions, progression
categories) with the R0 and infection fatality maps. `rd` must be a default-style
ResultData object (`ProgressionCategories` and `CaseFatalityMap` need the raw infections).
"""
function plot_combined_readouts(rd, r0map, fatalitymap)

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