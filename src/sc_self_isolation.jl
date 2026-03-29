export run_self_isolation_sims
export plot_scenario_self_isolation

export scenario_self_isolation



"""
    run_self_isolation_sims(;
        configfile = "models/DE_model_R0_3.28.toml",
        population = "DE")

Runs the self-isolation scenario simulations and returns 
two ResultData objects: one for the baseline simulation and one for the scenario with self-isolation.
And also returns last simulation object for further data extraction.
"""
function run_self_isolation_sims(;
    configfile = "models/DE_model_R0_3.28.toml",
    population = "DE")

    # baseline simultion
    sim = Simulation(configfile=configfile, population = population, seed = 0)
    sim.label = "Baseline"
    run!(sim)
    rd_b = ResultData(sim, style = "LightRD")


    # simulation with 14-day isolation (at home) upon experiencing symptoms
    sim = Simulation(configfile=configfile, population = population, seed = 0)
    sim.label = "Scenario"
    self_isolation = IStrategy("Self Isolation", sim)
    add_measure!(self_isolation, SelfIsolation(14))
    trigger = SymptomTrigger(self_isolation)
    add_symptom_trigger!(sim, trigger)

    run!(sim)
    rd_s = ResultData(sim, style = "LightRD")

    return rd_b, rd_s, sim
end


"""
    plot_scenario_self_isolation(rd_b, rd_s)

Plots the self-isolation scenario results.
Requires two ResultData objects as input, one for the baseline and one for the scenario with self-isolation.
"""
function plot_scenario_self_isolation(rd_b, rd_s, sim)

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
                    # clims = (1.1, 1.2),
                    title = "R0 Reduction by County",
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
    scenario_self_isolation(;
        configfile = "models/DE_model_R0_3.28.toml",
        population = "DE",
        folder = nothing)

Runs the self-isolation scenario and all analyses.
"""
function scenario_self_isolation(;
    configfile = "models/DE_model_R0_3.28.toml",
    population = "DE",
    folder = nothing)

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder() : nothing    

    # run simulations and get result data objects
    rd_b, rd_s, sim = run_self_isolation_sims(
        configfile = configfile,
        population = population)

    p_si_sc = plot_scenario_self_isolation(rd_b, rd_s, sim)
    png(p_si_sc, joinpath(folder, "sc_self_isolation_(R0-comparison).png"))
end