export run_testing_sims
export run_recursive_testing_sims
export plot_scenario_testing

export scenario_testing

"""
    run_testing_sims(;
        configfile = "models/DE_model_R0_3.28_Heinsberg_Seed.toml",
        popfile = "data/people.jld2",
        settingsfile = "data/settings.jld2",
        isolation_duration = 14,
        sensitivity = .8,
        specificity = .99)

Runs three simulation scenarios (isolation-only, isolation-and-test, and isolation-testing-and-tracing).
Returns the three `ResultData` objects and the last `Simulation` object (for any required plotting):
`rd_isolation_only`, `rd_isolation_and_test`, `rd_testing_and_tracing`, `sim`.
"""
function run_testing_sims(;
    configfile = "models/DE_model_R0_3.28_Heinsberg_Seed.toml",
    popfile = "data/people.jld2",
    settingsfile = "data/settings.jld2",
    
    isolation_duration = 14,
    sensitivity = .8,
    specificity = .99)


    ###
    ### ISOLATION ONLY SCENARIO
    ###

    sim = Simulation(configfile, popfile, settingsfile, label = "Isolation Only")

    # STRATEGY & MEASURE DEFINITION
    isolation_only = IStrategy("Self Isolation", sim)
    add_measure!(isolation_only, SelfIsolation(isolation_duration))
    trigger = SymptomTrigger(isolation_only)
    add_symptom_trigger!(sim, trigger)

    run!(sim)
    rd_isolation_only = ResultData(sim, style = "LightRD")



    ###
    ### ISOLATION & TESTING SCENARIO
    ###

    sim = Simulation(configfile, popfile, settingsfile, label = "Isolation & Testing")


    # test type definition
    antigen_test = TestType("Antigen Test", pathogen(sim), sim,
    sensitivity = sensitivity,
    specificity = specificity)

    # STRATEGY & MEASURE DEFINITION
    isolation_and_test = IStrategy("Self Isolation", sim)
    add_measure!(isolation_and_test, SelfIsolation(isolation_duration))
    add_measure!(isolation_and_test,
        Test("symptomatic testing", antigen_test))
    trigger = SymptomTrigger(isolation_and_test)
    add_symptom_trigger!(sim, trigger)

    run!(sim)
    rd_isolation_and_test = ResultData(sim, style = "LightRD")



    ###
    ### TRACING & TESTING SCENARIO
    ###

    sim = Simulation(configfile, popfile, settingsfile, label = "Isolation, Testing & Tracing")

    # test type definition
    antigen_test = TestType("Antigen Test", pathogen(sim), sim,
        sensitivity = sensitivity,
        specificity = specificity)

    # STRATEGY DEFINITIONS
    testing_and_tracing = IStrategy("Isolation & Testing", sim,
        condition = i -> !isdetected(i)) # apply strategy only if individual is not already detected and reported
    contact_identification = IStrategy("Contact Identification", sim)
    isolate_and_test_contacts = IStrategy("Isolate & Test Contacts", sim,
        condition = i -> !isdetected(i)) # apply strategy only if individual is not already detected and reported
    cancel_isolation = IStrategy("Cancel Isolation", sim)


    # MEASURE DEFINITIONS
    # symptomatic isolation and testing strategy
    add_measure!(testing_and_tracing, SelfIsolation(isolation_duration))
    add_measure!(testing_and_tracing,
        Test("symptomatic testing", antigen_test,
            positive_followup = contact_identification,
            negative_followup = cancel_isolation
            ))
    
    # contact identification strategy
    add_measure!(contact_identification,
        FindSettingMembers(Household, isolate_and_test_contacts, nonself = true),
        offset = 1
        )
    add_measure!(contact_identification,
        FindSettingMembers(Office, isolate_and_test_contacts, nonself = true),
        offset = 1,
        condition = is_working)
    add_measure!(contact_identification,
        FindSettingMembers(SchoolClass, isolate_and_test_contacts, nonself = true),
        offset = 1,
        condition = is_student)

    # isolate and test contacts strategy
    add_measure!(isolate_and_test_contacts, SelfIsolation(isolation_duration))
    add_measure!(isolate_and_test_contacts,
        Test("contact testing", antigen_test,
                negative_followup = cancel_isolation))

    # cancel isolation strategy
    add_measure!(cancel_isolation, CancelSelfIsolation())

    # trigger
    trigger = SymptomTrigger(testing_and_tracing)
    add_symptom_trigger!(sim, trigger)

    run!(sim)
    rd_testing_and_tracing = ResultData(sim, style = "LightRD")

    return rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing, sim
end


"""
    plot_scenario_testing(rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing, sim)

Generates a combined plot with multiple analyses of the three scenarios
that are produced in `run_testing_sims()`.

"""
function plot_scenario_testing(rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing, sim)

    # calculate y axis limits for observed reproduction plot
    obs_R_ymax = max(
        observed_R(rd_isolation_only).upper_est_R |>
            vec -> vec[.!ismissing.(vec)] |>
            vec -> isempty(vec) ? 0 : maximum(vec),

        observed_R(rd_isolation_and_test).upper_est_R |>
            vec -> vec[.!ismissing.(vec)] |>
            vec -> isempty(vec) ? 0 : maximum(vec),

        observed_R(rd_testing_and_tracing).upper_est_R |>
            vec -> vec[.!ismissing.(vec)] |>
            vec -> isempty(vec) ? 0 : maximum(vec)
    ) |> ceil |> Int


    # calculate y axis limits for detected cases plot
    det_cases_ymax = max(
        detected_tick_cases(rd_isolation_only) |>
            df -> DataFrames.select(df, :total_reported_cases, :exposed_cnt) |>
            df -> max(maximum(df.total_reported_cases), maximum(df.exposed_cnt)),

        detected_tick_cases(rd_isolation_and_test) |>
            df -> DataFrames.select(df, :total_reported_cases, :exposed_cnt) |>
            df -> max(maximum(df.total_reported_cases), maximum(df.exposed_cnt)),

        detected_tick_cases(rd_testing_and_tracing) |>
            df -> DataFrames.select(df, :total_reported_cases, :exposed_cnt) |>
            df -> max(maximum(df.total_reported_cases), maximum(df.exposed_cnt))
    ) |>
    c -> c * 1.1 |> ceil |> Int

    # calculate y axis limits for dark figure plot 
    dark_fig_ymax = max(
        compartment_fill(rd_isolation_only) |>
            df -> df.exposed_cnt .+ df.infectious_cnt |>
            maximum,

        compartment_fill(rd_isolation_and_test) |>
            df -> df.exposed_cnt .+ df.infectious_cnt |>
            maximum,

        compartment_fill(rd_testing_and_tracing) |>
            df -> df.exposed_cnt .+ df.infectious_cnt |>
            maximum
    ) |>
    c -> c * 1.1 |> ceil |> Int

    return plot(
        gemsplot([rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing],
            type = :TickCases,
            titlefontsize = map_legendfontsize(),
            title = "Daily Infections",
            plot_title = "",
            xlims = (0, 365),
            bottom_margin = 15mm),
        gemsplot([rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing],
            type = :EffectiveReproduction,
            title = "Effective Reproduction Number",
            plot_title = "",
            xlims = (0, 365),
            bottom_margin = 15mm),
        # gemsplot([rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing],
        #     type = :CumulativeIsolations,
        #     title = "Currently Isolated People",
        #     plot_title = "",
        #     xlims = (0, 365),
        #     bottom_margin = 15mm,
        #     left_margin = 30mm,
        #     legend = :topright),
        gemsplot([rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing],
            type = :TotalTests,
            title = "Tests per Day",
            plot_title = "",
            xlims = (0, 365),
            bottom_margin = 15mm,
            left_margin = 30mm,
            legend = :topright),
        gemsplot([rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing],
            type = :ObservedReproduction,
            top_margin = -5mm,
            bottom_margin = 10mm,
            left_margin = 25mm,
            layout = (1, 3),
            xlims = (0, 365),
            ylims = (0, obs_R_ymax)),
        gemsplot([rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing],
            type = :DetectedCases,
            legend = :topright,
            bottom_margin = 10mm,
            left_margin = 30mm,
            layout = (1, 3),
            xlims = (0, 365),
            ylims = (0, det_cases_ymax)),
        plot(
            plot_title = "Dark Figure",
            gemsplot(rd_isolation_only, type = :ActiveDarkFigure,
                title = label(rd_isolation_only),
                plot_title = "",
                legend = :topright,
                xlims = (0, 365),
                ylims = (0, dark_fig_ymax)),
            gemsplot(rd_isolation_and_test, type = :ActiveDarkFigure,
                title = label(rd_isolation_and_test),
                plot_title = "",
                legend = :topright,
                xlims = (0, 365),
                ylims = (0, dark_fig_ymax)),
            gemsplot(rd_testing_and_tracing, type = :ActiveDarkFigure,
                title = label(rd_testing_and_tracing),
                plot_title = "",
                legend = :topright,
                xlims = (0, 365),
                ylims = (0, dark_fig_ymax)),
            bottom_margin = 10mm,
            left_margin = 30mm,
            layout = (1, 3)
        ),
        tickfontsize = map_legendfontsize(),
        legendfontsize = map_legendfontsize(),
        labelfontsize = map_legendfontsize(),
        titlefontsize = map_titlefontsize(),
        layout = (6, 1),
        size = (4500, 4500),
        legend = :topleft
    )
end


"""
    scenario_testing(;
        configfile = "models/DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml",
        popfile = "data/people.jld2",
        settingsfile = "data/settings.jld2",
        folder = nothing)

Runs the testing-scenario and all analyses.
"""
function scenario_testing(;
    configfile = "models/DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml",
    popfile = "data/people.jld2",
    settingsfile = "data/settings.jld2",
    folder = nothing)

    # create result folder
    isnothing(folder) ? folder = result_folder() : nothing

    # run simulations and get result data objects
    rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing, sim = run_testing_sims(
        configfile = configfile,
        popfile = popfile,
        settingsfile = settingsfile)

    # plot results
    p_testing = plot_scenario_testing(rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing, sim)
    png(p_testing, joinpath(folder, "sc_testing_and_observed_progression.png"))
end