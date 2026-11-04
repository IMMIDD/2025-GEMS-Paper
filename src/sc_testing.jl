export calc_outcomes
export outcome_df

export run_testing_sims
export run_recursive_testing_sims
export plot_scenario_testing

export scenario_testing

###
### DATA ANALYSIS
###

"""
    calc_outcomes(rds; kwargs...)

Generates a dataframe from a vector of `ResultData` objects and a series
of keyworded lamda-functions that extract values from each of the `rds`.

This example extracts the label and the number of total infections from each object:

```julia
calc_outcomes([rd1, rd2, rd3],
        label = label,
        total_infections = total_infections)
```
"""
function calc_outcomes(rds; kwargs...)

    res = Dict()

    for (key, val) in kwargs
        # apply function to each rd object, store result under name
        res[key] = val.(rds)
    end

    return DataFrame(res)
end

"""
    outcome_df(rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing)

Builds a dataframe with main simulation outcomes for the three RDs
that come out of `run_testing_sims()`
"""
function outcome_df(rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing)
    return calc_outcomes([rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing],
        label = label,
        total_infections = total_infections,
        total_tests = rd -> try
                total_tests(rd)["Antigen Test"]
            catch
                0
            end,
        total_attack_rate = attack_rate,
        quarantine_days = total_quarantines,
        lost_schooldays = rd -> sum(cumulative_quarantines(rd).students), 
        lost_workdays = rd -> sum(cumulative_quarantines(rd).workers), 
        lost_otherdays = rd -> sum(cumulative_quarantines(rd).other), 
        detection_rate = detection_rate,
        base_reproduction = r0,
        hospitalizations = rd -> sum(tick_hosptitalizations(rd).hospital_admissions),
        icu_admissions = rd -> sum(tick_hosptitalizations(rd).icu_admissions),
        deaths = rd -> sum(tick_deaths(rd).death_cnt),
        time_to_detection = rd -> time_to_detection(rd).mean_time_to_detection |>
            vec -> vec[.!ismissing.(vec)] |>
            vec -> isempty(vec) ? Inf : median(vec),
        max_hospital_load = rd -> maximum(tick_hosptitalizations(rd).current_hospitalized),
        max_icu_load = rd -> maximum(tick_hosptitalizations(rd).current_icu),
        mild_days = rd -> sum(compartment_periods(rd).symptomatic) - sum(compartment_periods(rd).severe),
        severe_days = rd -> sum(compartment_periods(rd).severe),
        hospital_days = rd -> sum(compartment_periods(rd).hospitalized),
        icu_days = rd -> sum(compartment_periods(rd).icu),
        ventilation_days = rd -> sum(compartment_periods(rd).ventilated)
    )
end

"""
    run_testing_sims(;
        configfile = "models/DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml",
        population = "DE",
        start_condition = RegionalSeeds(seeds = Dict(05370008=>10)), # 10 infections in Heinsberg
        stop_criterion = TimesUp(limit = 1095),
        isolation_duration = 14,
        sensitivity = .8,
        specificity = .99)

Runs three simulation scenarios (isolation-only, isolation-and-test, and isolation-testing-and-tracing).
Returns the three `ResultData` objects and the last `Simulation` object (for any required plotting):
`rd_isolation_only`, `rd_isolation_and_test`, `rd_testing_and_tracing`, `sim`.
"""
function run_testing_sims(;
    configfile = "models/DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml",
    population = "DE",

    start_condition = RegionalSeeds(seeds = Dict(05370008=>10)), # 10 infections in Heinsberg
    stop_criterion = TimesUp(limit = 1095),
    
    isolation_duration = 14,
    sensitivity = .8,
    specificity = .99)


    ###
    ### ISOLATION ONLY SCENARIO
    ###

    sim = Simulation(
        label = "Isolation Only",
        configfile = configfile,
        population = population, 
        seed = 0,
        start_condition = start_condition,
        stop_criterion = stop_criterion)

    # STRATEGY & MEASURE DEFINITION
    isolation_only = IStrategy("Self Isolation", sim)
    add_measure!(isolation_only, SelfIsolation(isolation_duration))
    trigger = SymptomTrigger(isolation_only)
    add_symptom_trigger!(sim, trigger)

    run!(sim)
    rd_isolation_only = ResultData(sim)



    ###
    ### ISOLATION & TESTING SCENARIO
    ###

    sim = Simulation(
        label = "Isolation & Testing",
        configfile = configfile,
        population=population, 
        seed = 1,
        start_condition = start_condition,
        stop_criterion = stop_criterion)


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
    rd_isolation_and_test = ResultData(sim)



    ###
    ### TRACING & TESTING SCENARIO
    ###

    sim = Simulation(
        label = "Isolation, Testing & Tracing",
        configfile = configfile,
        population = population, 
        seed = 2,
        start_condition = start_condition,
        stop_criterion = stop_criterion)

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
    rd_testing_and_tracing = ResultData(sim)

    return rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing, sim
end


"""
    plot_scenario_testing(rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing, sim)

Generates a combined plot with multiple analyses of the three scenarios
that are produced in `run_testing_sims()`.

"""
function plot_scenario_testing(rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing, sim)

    max_ticks = final_tick(rd_isolation_only)
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


    # tick conversion function (makes 100_000 to 0.1m)
    tick_to_m(t) = (t / 10^6) |>
        t -> round(t, digits = 1) |>
        t -> t < 10 ? "  $(t)m" : "$(t)m"

    # spacing R plot ticks
    tick_to_R(t) = "        $(round(Int, t))"


    l = @layout [a{0.12h}; b{0.12h}; c{0.12h}; d{0.12h}; e{0.52h}]
    return plot(
        # DAILY INFECTIONS
        gemsplot([rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing],
            type = :TickCases,
            titlefontsize = map_legendfontsize(),
            title = "",
            plot_title = "",
            ylabel = "Infections",
            yformatter = tick_to_m,
            xformatter = x -> "",
            xlims = (0, max_ticks),
            xlabel = "",
            legend = :topright),
            

        # EFFECTIVE REPRODUCTION
        gemsplot([rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing],
            type = :EffectiveReproduction,
            title = "",
            ylabel = "Effective R",
            yformatter = tick_to_R,
            xformatter = x -> "",
            plot_title = "",
            xlims = (0, max_ticks),
            xlabel = "",
            legend = false),

        # TESTS PER DAY
        gemsplot([rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing],
            type = :TotalTests,
            title = "",
            ylabel = "Tests",
            yformatter = tick_to_m,
            xformatter = x -> "",
            plot_title = "",
            xlims = (0, max_ticks),
            xlabel = "",
            #left_margin = 30mm,
            legend = false),

        # ISOLATED PEOPLE
        gemsplot([rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing],
            type = :CumulativeIsolations,
            title = "",
            ylabel = "In Quarantine",
            yformatter = tick_to_m,
            plot_title = "",
            xlims = (0, max_ticks),
            bottom_margin = 15mm,
            #left_margin = 30mm,
            legend = false),

        # GRID
        plot(
            # FIRST COLUMN
            plot(
                # OBSERVED REPRODUCTION
                gemsplot(rd_isolation_only, type = :ObservedReproduction,
                    title = label(rd_isolation_only),
                    plot_title = "",
                    ylabel = "Observed R",
                    yformatter = tick_to_R,
                    ylims = (0, obs_R_ymax),
                    xformatter = x -> "",
                    xlabel = "",
                    legend = :topright),
                # DETECTED CASES
                gemsplot(rd_isolation_only, type = :DetectedCases,
                    title = "",
                    plot_title = "",
                    ylabel = "Detections",
                    yformatter = tick_to_m,
                    ylims = (0, det_cases_ymax),
                    xformatter = x -> "",
                    xlabel = "",
                    legend = :topright),
                # DARK FIGURE
                gemsplot(rd_isolation_only, type = :ActiveDarkFigure,
                    title = "",
                    plot_title = "",
                    ylabel = "Dark Figure",
                    yformatter = tick_to_m,
                    ylims = (0, dark_fig_ymax),
                    xlabel = "Day",
                    legend = :topright),
                layout = (3, 1),
                left_margin = 30mm
            ),

            # SECOND COLUMN
            plot(
                # OBSERVED REPRODUCTION
                gemsplot(rd_isolation_and_test, type = :ObservedReproduction,
                    title = label(rd_isolation_and_test),
                    plot_title = "",
                    ylims = (0, obs_R_ymax),
                    yformatter = x -> "",
                    ylabel = "",
                    xformatter = x -> "",
                    xlabel = "",
                    legend = false),
                # DETECTED CASES
                gemsplot(rd_isolation_and_test, type = :DetectedCases,
                    title = "",
                    plot_title = "",
                    ylims = (0, det_cases_ymax),
                    yformatter = x -> "",
                    ylabel = "",
                    xformatter = x -> "",
                    xlabel = "",
                    legend = false),
                # DARK FIGURE
                gemsplot(rd_isolation_and_test, type = :ActiveDarkFigure,
                    title = "",
                    plot_title = "",
                    ylims = (0, dark_fig_ymax),
                    yformatter = x -> "",
                    ylabel = "",
                    xlabel = "Day",
                    legend = false),
                layout = (3, 1)
            ),

            # THIRD COLUMN
            plot(
                # OBSERVED REPRODUCTION
                gemsplot(rd_testing_and_tracing, type = :ObservedReproduction,
                    title = label(rd_testing_and_tracing),
                    plot_title = "",
                    ylims = (0, obs_R_ymax),
                    yformatter = x -> "",
                    ylabel = "",
                    xformatter = x -> "",
                    xlabel = "",
                    legend = false),         
                # DETECTED CASES
                gemsplot(rd_testing_and_tracing, type = :DetectedCases,
                    title = "",
                    plot_title = "",
                    ylims = (0, det_cases_ymax),
                    yformatter = x -> "",
                    ylabel = "",
                    xformatter = x -> "",
                    xlabel = "",
                    legend = false),
                # DARK FIGURE
                gemsplot(rd_testing_and_tracing, type = :ActiveDarkFigure,
                    title = "",
                    plot_title = "",
                    ylims = (0, dark_fig_ymax),
                    yformatter = x -> "",
                    ylabel = "",
                    xlabel = "Day",
                    legend = false),
                layout = (3, 1),
            ),
            layout = (1, 3),
            bottom_margin = 5mm,
        ),

        tickfontsize = map_legendfontsize(),
        legendfontsize = round(Int, map_legendfontsize() * 0.8),
        xlabelfontsize = map_legendfontsize(),
        ylabelfontsize = map_titlefontsize(),
        titlefontsize = map_legendfontsize(),

        xlims = (0, max_ticks),
        grid = true,
        gridlinewidth = 2.0,
        layout = l,
        size = (2800, 2500),
    )
end


"""
    scenario_testing(;
        configfile = "models/DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml",
        population = "DE",
        folder = nothing)

Runs the testing-scenario and all analyses.
"""
function scenario_testing(;
    configfile = "models/DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml",
    population = "DE",
    folder = nothing)

    # create result folder
    isnothing(folder) ? folder = result_folder() : nothing

    # run simulations and get result data objects
    rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing, sim = run_testing_sims(
        configfile = configfile,
        population = population)

    # plot results
    p_testing = plot_scenario_testing(rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing, sim)
    png(p_testing, joinpath(folder, "sc_testing_and_observed_progression.png"))

    # extract KPIs and store as CSV
    df = outcome_df(rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing)
    CSV.write(joinpath(folder, "sc_testing_outcomes.csv"), df)
end