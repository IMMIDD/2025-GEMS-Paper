export TestingRD
export calc_outcomes
export outcome_df
export comparison_df
export daily_df
export summary_df

export run_testing_sims
export plot_scenario_testing

export scenario_testing
export replot_testing

###
### RESULT DATA STYLE
###

"""
    TestingRD <: ResultDataStyle

Result data of the testing scenarios: only what the figure and the outcome tables
need, so that all runs can be stored (`exportJLD`) and loaded again
(`import_resultdata`). Of the compartment periods, it holds the sums over all
infections. The tests are only held per day, except for the isolations after
false positive tests (see `false_positive_isolations`).
"""
mutable struct TestingRD <: ResultDataStyle
    data::Dict{String, Any}
    function TestingRD(pP::PostProcessor)
        sim = simulation(pP)
        return new(Dict(
            "sim_data" => Dict(
                "label" => label(sim),
                "seed" => seed(sim),
                "final_tick" => tick(sim),
                "tick_unit" => tickunit(sim),
                "total_infections" => nrow(infectionsDF(pP)),
                "attack_rate" => attack_rate(pP),
                "r0" => r0(pP),
                "total_quarantines" => total_quarantines(pP),
                "total_tests" => total_tests(pP),
                "detection_rate" => detection_rate(pP),
                "false_positive_isolations" => false_positive_isolations(pP),
                "members" => count_members(sim)),
            "dataframes" => Dict(
                "tick_cases" => tick_cases(pP),
                "effectiveR" => effectiveR(pP),
                "tick_tests" => tick_tests(pP),
                "cumulative_quarantines" => cumulative_quarantines(pP),
                "observed_R" => observed_R(pP),
                "detected_tick_cases" => detected_tick_cases(pP),
                "compartment_fill" => compartment_fill(pP),
                "tick_hosptitalizations" => hospital_df(pP),
                "tick_deaths" => tick_deaths(pP),
                "time_to_detection" => time_to_detection(pP),
                "compartment_periods" => mapcols(c -> [sum(c)], compartment_periods(pP)))))
    end
end

# number of workers (with an office) and students (with a class of any educational setting)
count_members(sim::Simulation) = Dict(
    "workers" => count(is_working, individuals(sim)),
    "students" => count(is_student, individuals(sim)))

"""
    false_positive_isolations(pP::PostProcessor)

Number of individuals who were isolated after a false positive test, and the days they
spent in isolation because of it. An isolation starts with the test and ends with the next
test of the individual or after the isolation duration. Days on which the individual also
tests negative do not count, since the negative result ends the isolation on the same day.

Assumes that all tested individuals isolate, and that a negative test ends the isolation, as in the
tracing scenario. The isolation duration is taken from the strategies of the simulation.
"""
function false_positive_isolations(pP::PostProcessor)
    tests = testsDF(pP)
    false_positive = tests.test_result .& .!tests.infected
    any(false_positive) || return Dict("individuals" => 0, "days" => 0)

    sim = simulation(pP)
    isolations = [measure(m) for s in strategies(sim) for m in measures(s) if measure(m) isa SelfIsolation]
    isolation_duration = only(unique(GEMS.duration.(isolations)))

    # all tests of the individuals with a false positive test, one row per individual and day
    ids = Set(tests.id[false_positive])
    days = combine(groupby(tests[in.(tests.id, Ref(ids)), [:id, :test_tick, :test_result, :infected]], [:id, :test_tick]),
        :test_result => all => :positive,
        :infected => first => :infected)
    sort!(days, [:id, :test_tick])

    individuals = Set{eltype(days.id)}()
    isolation_days = 0
    for k in eachindex(days.id)
        days.positive[k] && !days.infected[k] || continue
        start = Int(days.test_tick[k])
        stop = min(start + isolation_duration, Int(tick(sim)))
        if k < nrow(days) && days.id[k + 1] == days.id[k]
            stop = min(stop, Int(days.test_tick[k + 1]))
        end
        push!(individuals, days.id[k])
        isolation_days += stop - start
    end
    return Dict("individuals" => length(individuals), "days" => isolation_days)
end

###
### DATA ANALYSIS
###

# first day of the maximum of the daily infections
peak_day(rd) = tick_cases(rd).tick[argmax(tick_cases(rd).exposed_cnt)]

# days from the day on which 5% of all infections have occurred to the day on which 95% have occurred
function wave_duration(rd)
    cases = tick_cases(rd)
    total = cumsum(cases.exposed_cnt)
    day(share) = cases.tick[findfirst(>=(share * total[end]), total)]
    return day(0.95) - day(0.05)
end

# daily antigen tests (`nothing` without tests)
antigen_tests(rd) = get(tick_tests(rd), "Antigen Test", nothing)

# total over all days of a column of the daily counts of individuals, or of the difference of two columns
total_days(rd, col) = sum(compartment_fill(rd)[!, col])
total_days(rd, col, minus) = total_days(rd, col) - total_days(rd, minus)

# `x` per member of the `group` ("workers" or "students"), `missing` if their number is not known
per_member(x, rd, members, group) = x / get(get(sim_data(rd), "members", members), group, missing)

# `key` of the isolations after false positive tests (0 without false positives, `missing` if not stored)
function false_positive_isolation(rd, key)
    sum(detected_tick_cases(rd).false_positives) == 0 && return 0
    haskey(sim_data(rd), "false_positive_isolations") || return missing
    return sim_data(rd)["false_positive_isolations"][key]
end

"""
    calc_outcomes(rds; kwargs...)

Generates a dataframe from a vector of `ResultData` objects and a series
of keyworded lambda functions that extract values from each of the `rds`.

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
    outcome_df(rds; members = Dict())

Builds a dataframe with main simulation outcomes of each of the `rds`
that come out of `run_testing_sims()` (one row per run). Among them:

- `peak_day`, `max_daily_infections`: day and height of the peak of the daily infections
- `wave_duration`: days between 5% and 95% of all infections
- `max_undetected_infections`: maximum of the active infections that are not detected (dark figure)
- `max_quarantined`: maximum of the individuals in isolation or quarantine on one day
- `max_daily_tests`, `test_positivity`: maximum of the daily tests and share of positive results
  (`missing` without tests)
- `false_positives`: number of positive test results of individuals who are not infected
- `false_positive_isolated`, `false_positive_isolation_days`: individuals isolated after a false positive test
  and the days they spent in isolation because of it (see `false_positive_isolations`; `missing` for runs that
  were stored without them)
- `workdays_absent`, `schooldays_absent`: days on which workers (students) are isolated, severely ill, or in hospital.
  They are the sum of `..._isolated_infected` and `..._isolated_not_infected` (isolated, and infected or not infected
  on that day; together `lost_workdays` and `lost_schooldays`) and `..._ill_not_isolated` (severely ill or in
  hospital without isolation). Students are the members of all classes, and every day of the week counts.
- `max_absent_workers`, `max_absent_students`: maximum of the absent workers (students) on one day
- `max_absent_workers_share`, `max_absent_students_share`, `lost_workdays_per_worker`, `lost_schooldays_per_student`:
  the same and the lost days, divided by the number of workers (students) in the population. For runs that were stored
  without these numbers, they are taken from `members` (`missing` if empty).
"""
function outcome_df(rds; members = Dict())
    df = calc_outcomes(rds,
        label = label,
        seed = rd -> sim_data(rd)["seed"],
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
        initial_reproduction = r0,
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
        ventilation_days = rd -> sum(compartment_periods(rd).ventilated),
        peak_day = peak_day,
        max_daily_infections = rd -> maximum(tick_cases(rd).exposed_cnt),
        wave_duration = wave_duration,
        max_undetected_infections = rd -> maximum(compartment_fill(rd).exposed_cnt .+
            compartment_fill(rd).infectious_cnt .- compartment_fill(rd).detected_cnt),
        max_quarantined = rd -> maximum(cumulative_quarantines(rd).quarantined),
        max_daily_tests = rd -> isnothing(antigen_tests(rd)) ? 0 : maximum(antigen_tests(rd).total_tests),
        test_positivity = rd -> isnothing(antigen_tests(rd)) ? missing :
            sum(antigen_tests(rd).positive_tests) / sum(antigen_tests(rd).total_tests),
        false_positives = rd -> sum(detected_tick_cases(rd).false_positives),
        false_positive_isolated = rd -> false_positive_isolation(rd, "individuals"),
        false_positive_isolation_days = rd -> false_positive_isolation(rd, "days"),
        workdays_absent = rd -> total_days(rd, :unable_to_attend_workers),
        workdays_isolated_infected = rd -> total_days(rd, :isolated_workers),
        workdays_isolated_not_infected = rd -> total_days(rd, :quarantined_workers, :isolated_workers),
        workdays_ill_not_isolated = rd -> total_days(rd, :unable_to_attend_workers, :quarantined_workers),
        schooldays_absent = rd -> total_days(rd, :unable_to_attend_students),
        schooldays_isolated_infected = rd -> total_days(rd, :isolated_students),
        schooldays_isolated_not_infected = rd -> total_days(rd, :quarantined_students, :isolated_students),
        schooldays_ill_not_isolated = rd -> total_days(rd, :unable_to_attend_students, :quarantined_students),
        max_absent_workers = rd -> maximum(compartment_fill(rd).unable_to_attend_workers),
        max_absent_students = rd -> maximum(compartment_fill(rd).unable_to_attend_students),
        max_absent_workers_share = rd ->
            per_member(maximum(compartment_fill(rd).unable_to_attend_workers), rd, members, "workers"),
        max_absent_students_share = rd ->
            per_member(maximum(compartment_fill(rd).unable_to_attend_students), rd, members, "students"),
        lost_workdays_per_worker = rd -> per_member(sum(cumulative_quarantines(rd).workers), rd, members, "workers"),
        lost_schooldays_per_student = rd -> per_member(sum(cumulative_quarantines(rd).students), rd, members, "students")
    )
    return DataFrames.select(df, :label, :seed, :) # label and seed first
end

"""
    daily_df(rds)

Builds a dataframe with the daily data of each of the `rds` that come
out of `run_testing_sims()` (one row per run and day).
"""
function daily_df(rds)
    # daily data of one run
    function daily(rd)
        tests = get(tick_tests(rd), "Antigen Test", DataFrame(tick = Int16[], total_tests = Int[]))

        df = reduce((a, b) -> leftjoin(a, b, on = :tick), [
            DataFrames.select(tick_cases(rd), :tick, :exposed_cnt => :infections),
            DataFrames.select(effectiveR(rd), :tick, :rolling_R => :effective_R),
            DataFrames.select(tests, :tick, :total_tests),
            DataFrames.select(cumulative_quarantines(rd), :tick, :quarantined, :students => :quarantined_students,
                :workers => :quarantined_workers, :other => :quarantined_other),
            DataFrames.select(observed_R(rd), :tick, :mean_est_R => :observed_R,
                :lower_est_R => :observed_R_lower, :upper_est_R => :observed_R_upper),
            DataFrames.select(detected_tick_cases(rd), :tick, :new_detections, :double_reports, :false_positives),
            DataFrames.select(compartment_fill(rd), :tick, [:exposed_cnt, :infectious_cnt] => (+) => :active_cases,
                :detected_cnt => :active_detected_cases),
            DataFrames.select(tick_hosptitalizations(rd), :tick, :current_hospitalized, :current_icu),
            DataFrames.select(tick_deaths(rd), :tick, :death_cnt => :deaths)])
        insertcols!(df, 1, :label => label(rd), :seed => sim_data(rd)["seed"])
        df.total_tests = coalesce.(df.total_tests, 0)
        return sort!(df, :tick)
    end
    return reduce(vcat, daily.(rds))
end

"""
    summary_df(outcomes)

Summarizes the `outcomes` of all runs (from `outcome_df()` or `comparison_df()`) per strategy
and metric: number of runs with a value, mean, standard deviation, minimum, maximum, and
95% confidence interval of the mean.
"""
function summary_df(outcomes)
    rows = [(label = lab, metric = metric,
            values = Float64.(collect(skipmissing(outcomes[outcomes.label .== lab, metric]))))
        for metric in names(outcomes, Not(:label, :seed)) for lab in unique(outcomes.label)]
    return DataFrame(map(rows) do r
        stats = isempty(r.values) ? Dict() : GEMS.aggregate_values(r.values)
        stat(key) = get(stats, key, missing)
        (label = r.label, metric = r.metric, runs = length(r.values), mean = stat("mean"), std = stat("std"),
            min = stat("min"), max = stat("max"), lower_95 = stat("lower_95"), upper_95 = stat("upper_95"))
    end)
end

"""
    comparison_df(outcomes; strategy = "Isolation, Testing & Tracing")

Compares the `strategy` with each of the other strategies of the `outcomes` (from `outcome_df()`), one row
per repetition of `run_testing_sims()` (in the order of the random seeds): averted infections, deaths, and
lost workdays and schooldays, additional tests, and additional tests per averted infection, death, lost workday, and
lost schoolday.
"""
function comparison_df(outcomes; strategy = "Isolation, Testing & Tracing")
    runs(lab) = sort(outcomes[outcomes.label .== lab, :], :seed)
    a = runs(strategy)
    return reduce(vcat, map(filter(!=(strategy), unique(outcomes.label))) do lab
        b = runs(lab)
        k = 1:min(nrow(a), nrow(b))
        averted(metric) = b[k, metric] .- a[k, metric]
        additional_tests = -averted(:total_tests)
        DataFrame(
            label = "$strategy vs. $lab",
            seed = a.seed[k],
            averted_infections = averted(:total_infections),
            averted_deaths = averted(:deaths),
            averted_lost_workdays = averted(:lost_workdays),
            averted_lost_schooldays = averted(:lost_schooldays),
            additional_tests = additional_tests,
            tests_per_averted_infection = additional_tests ./ averted(:total_infections),
            tests_per_averted_death = additional_tests ./ averted(:deaths),
            tests_per_averted_lost_workday = additional_tests ./ averted(:lost_workdays),
            tests_per_averted_lost_schoolday = additional_tests ./ averted(:lost_schooldays))
    end)
end

"""
    run_testing_sims(;
        configfile = "models/DE_model_first_wave.toml",
        population = "NRW",
        municipality = 05370008, # Gangelt (Kreis Heinsberg)
        start_condition = RegionalSeeds(seeds = Dict(municipality => 50)),
        stop_criterion = TimesUp(limit = 365),
        municipality_contact_reduction::Float64 = 0.5,
        isolation_duration = 14,
        sensitivity = .8,
        specificity = .99,
        seed = 0,
        after_run = rd -> nothing)

Runs three simulation scenarios (isolation-only, isolation-and-test, and isolation-testing-and-tracing)
on the COVID-19 reference model of the first wave example (`configfile`, without its contact phases)
with the random seeds `seed`, `seed + 1`, and `seed + 2`. `start_condition` and `stop_criterion` replace
those of the config. By default, 50 individuals are infected in the `municipality` (AGS), which must lie within
`population` (e.g., `population = "SL", municipality = 10042111` for Beckingen in the Saarland).
All scenarios assume that a share of the community contacts does not take place
(`municipality_contact_reduction`), which reduces the contact parameter of all municipalities accordingly
(by default, from 4 to 2 in the reference model); contacts in households, workplaces, and schools are not reduced.

In all scenarios, individuals isolate for `isolation_duration` days upon symptom onset.
The isolation-and-test scenario additionally tests them with an antigen test (`sensitivity`, `specificity`)
on the same day, which detects cases without changing the isolation. In the tracing scenario, individuals
who test negative are released, and the household, office, and school class members of a positively tested
individual are identified one day later, isolated, and tested (and released if negative).

Returns a vector with the `ResultData` objects (style `TestingRD`) of the three scenarios. The function `after_run`
is called with each of them as soon as its simulation is done (e.g., to store it).
"""
function run_testing_sims(;
    configfile = "models/DE_model_first_wave.toml",
    population = "NRW",
    municipality = 05370008, # Gangelt (Kreis Heinsberg)

    start_condition = RegionalSeeds(seeds = Dict(municipality => 50)),
    stop_criterion = TimesUp(limit = 365),

    municipality_contact_reduction::Float64 = 0.5,
    isolation_duration = 14,
    sensitivity = .8,
    specificity = .99,

    seed = 0,
    after_run = rd -> nothing)

    0.0 <= municipality_contact_reduction <= 1.0 || throw(ArgumentError("Municipality contact reduction must be between 0 and 1. Got $municipality_contact_reduction."))

    # assume that a share of the community contacts does not take place
    # (scales the contact parameter of the config's municipalities)
    function reduce_municipality_contacts!(sim)
        base_csm = contact_sampling_method(first(settings(sim, Municipality)))
        csm = with_contactparameter(base_csm,
            base_csm.contactparameter * (1.0 - municipality_contact_reduction))
        foreach(s -> contact_sampling_method!(s, csm), settings(sim, Municipality))
    end

    ###
    ### ISOLATION ONLY SCENARIO
    ###

    sim = Simulation(
        label = "Isolation Only",
        configfile = configfile,
        population = population, 
        seed = seed,
        start_condition = start_condition,
        stop_criterion = stop_criterion)

    reduce_municipality_contacts!(sim)

    # STRATEGY & MEASURE DEFINITION
    isolation_only = IStrategy("Self Isolation", sim)
    add_measure!(isolation_only, SelfIsolation(isolation_duration))
    trigger = SymptomTrigger(isolation_only)
    add_symptom_trigger!(sim, trigger)

    run!(sim)
    rd_isolation_only = ResultData(sim, style = "TestingRD")
    after_run(rd_isolation_only)



    ###
    ### ISOLATION & TESTING SCENARIO
    ###

    sim = Simulation(
        label = "Isolation & Testing",
        configfile = configfile,
        population=population,
        seed = seed + 1,
        start_condition = start_condition,
        stop_criterion = stop_criterion)

    reduce_municipality_contacts!(sim)

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
    rd_isolation_and_test = ResultData(sim, style = "TestingRD")
    after_run(rd_isolation_and_test)



    ###
    ### TRACING & TESTING SCENARIO
    ###

    sim = Simulation(
        label = "Isolation, Testing & Tracing",
        configfile = configfile,
        population = population,
        seed = seed + 2,
        start_condition = start_condition,
        stop_criterion = stop_criterion)

    reduce_municipality_contacts!(sim)

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
    rd_testing_and_tracing = ResultData(sim, style = "TestingRD")
    after_run(rd_testing_and_tracing)

    return [rd_isolation_only, rd_isolation_and_test, rd_testing_and_tracing]
end


###
### FIGURE
###

# The panels are drawn by GEMS (`gemsplot`); only their axes and labels are set here.
# Print size, fonts, titles, and section headings as in the population validation
# figure (width of the page of PLOS; the height leaves room for the caption).

# heights (inches) of the figure, its section headings, its column headers, and its title rows
const TS_HEIGHT_IN = 7.4
const TS_HEADER_IN = FIG_HEADING_IN
const TS_COLUMN_IN = 0.28
const TS_TITLE_IN = 0.22

# approximate width of the plot areas (inches), days shown, line width
const TS_PLOT_WIDTH_IN = 6.9
const TS_DAYS = 300
const TS_LINEWIDTH = 0.8

# tick labels (makes 2_500 to 2.5k, 20_000 to 20k, and 1_500_000 to 1.5m; at most one decimal)
ts_number(t) = (r = round(t, digits = 1); isinteger(r) ? "$(Int(r))" : "$r")
ts_tick(t) = abs(t) >= 10^6 ? ts_number(t / 10^6) * "m" : abs(t) >= 10^3 ? ts_number(t / 10^3) * "k" : ts_number(t)

# axes of all panels; tick labels of the x axis (and room for its label) only in the last row of a section
ts_axes(last_row) = (
    title = "", plot_title = "", ylabel = "",
    xlabel = "",
    xformatter = last_row ? (x -> "$(round(Int, x))") : (x -> ""),
    xlims = (0, TS_DAYS),
    legend = false,
    fontfamily = PV_FONT, tickfontsize = PV_SMALL_FONTSIZE, guidefontsize = PV_FONTSIZE,
    grid = true, gridalpha = 0.15, gridlinewidth = 0.5,
    margin = 0mm, right_margin = 1mm, top_margin = -2mm, bottom_margin = last_row ? 1mm : -3mm)

# label of the x axis (close to the tick labels)
ts_xlabel!(p) = annotate!(p, (0.5, -0.39), Plots.text("Day", PV_FONTSIZE, PV_FONT, :center))

# labels of the series of GEMS in the legends of the observed epidemic
const TS_LABELS = Dict(
    "Mean Estimated R" => "Estimated R",
    "95% Confidence Band" => "95% CI of the serial interval",
    "R=1" => "R = 1",
    "New Detections" => "New detections",
    "Double Reports" => "Double reports",
    "False Positives" => "False positives",
    "New True Cases" => "New true cases",
    "Detected Cases" => "Detected cases",
    "Dark Figure" => "Dark figure")

# sets line width (thin outlines of areas) and labels (by GEMS label) of the series that GEMS drew
function ts_style!(p; labels = Dict())
    for s in p[1].series_list
        s[:linewidth] = isnothing(s[:fillrange]) ? TS_LINEWIDTH : 0.3
        s[:label] = get(labels, s[:label], s[:label])
    end
    return p
end

# empty panel for text (x and y from 0 to 1)
ts_strip() = plot(framestyle = :none, ticks = false, legend = false, grid = false,
    xlims = (0, 1), ylims = (0, 1), margin = 0mm, top_margin = -2mm, bottom_margin = -2mm,
    fontfamily = PV_FONT)

# header of a column: name of the strategy above a line in its color
function ts_column(name, color)
    p = plot!(ts_strip(), [0, 1], [0.3, 0.3], linecolor = color, linewidth = 1.5)
    return annotate!(p, 0.5, 0.7, Plots.text(name, PV_FONTSIZE, PV_FONT_BOLD, :center, :vcenter))
end

# title row of a panel: the title (left) and, if `p` is given, the legend of its series
# as GEMS drew them (right), optionally in the `order` of their labels
function ts_title(title, p = nothing; order = nothing)
    strip = ts_strip()
    title_text!(strip, 0, 0.5, title, TS_PLOT_WIDTH_IN, valign = :vcenter)
    isnothing(p) && return strip

    entries = filter(s -> !(string(s[:label]) in ("", "nothing")), p[1].series_list)
    isnothing(order) || sort!(entries, by = s -> something(findfirst(==(s[:label]), order), length(order) + 1))
    x = TS_PLOT_WIDTH_IN # from right to left
    for s in reverse(entries)
        annotate!(strip, x / TS_PLOT_WIDTH_IN, 0.5, Plots.text(s[:label], PV_SMALL_FONTSIZE, PV_FONT, :right, :vcenter))
        x -= drawn_width_in(s[:label])
        key = [x - 0.28, x - 0.06] ./ TS_PLOT_WIDTH_IN
        isnothing(s[:fillrange]) ?
            plot!(strip, key, [0.5, 0.5], linecolor = s[:linecolor], linestyle = s[:linestyle], linewidth = 1.5) :
            plot!(strip, key, [0.7, 0.7], fillrange = 0.3, fillcolor = s[:fillcolor], fillalpha = s[:fillalpha], linewidth = 0)
        x -= 0.42
    end
    return strip
end

# run of the strategy `lab` with the median peak day of infections (lower median)
function typical_run(rds, lab)
    runs = sort(filter(rd -> label(rd) == lab, rds), by = peak_day)
    return runs[cld(length(runs), 2)]
end

"""
    plot_scenario_testing(rds)

Generates a combined plot with multiple analyses of the scenarios in `rds` (any number
of runs of `run_testing_sims()`): all runs in the top panels (one line per run), and for
each strategy the run with the median peak day of infections in the bottom panels.
"""
function plot_scenario_testing(rds)

    strategies = unique(label.(rds))
    typical = [typical_run(rds, lab) for lab in strategies]

    # ALL RUNS
    all_runs(type, formatter; last_row = false) =
        ts_style!(gemsplot(rds; type = type, ts_axes(last_row)..., yformatter = formatter))
    top = [
        all_runs(:TickCases, ts_tick),
        all_runs(:EffectiveReproduction, ts_number),
        all_runs(:TotalTests, ts_tick),
        all_runs(:CumulativeIsolations, ts_tick, last_row = true)]

    # ONE RUN PER STRATEGY (same y axis for all strategies)
    per_strategy(type, formatter, ymax; last_row = false) =
        map(enumerate(typical)) do (k, rd)
            ts_style!(gemsplot(rd; type = type, ts_axes(last_row)...,
                ylims = (0, ymax),
                yformatter = k == 1 ? formatter : (y -> "")); labels = TS_LABELS)
        end
    highest(f) = maximum(rd -> maximum(skipmissing(f(rd)), init = 0.0), typical)
    observed = per_strategy(:ObservedReproduction, ts_number, ceil(highest(rd -> observed_R(rd).upper_est_R)))
    detected = per_strategy(:DetectedCases, ts_tick,
        1.1 * highest(rd -> max.(detected_tick_cases(rd).total_reported_cases, detected_tick_cases(rd).exposed_cnt)))
    darkfigure = per_strategy(:ActiveDarkFigure, ts_tick,
        1.1 * highest(rd -> compartment_fill(rd).exposed_cnt .+ compartment_fill(rd).infectious_cnt),
        last_row = true)

    foreach(ts_xlabel!, [top[end]; darkfigure])

    for (p, rd) in zip(observed, typical)
        sum(detected_tick_cases(rd).new_detections) == 0 && annotate!(p, TS_DAYS / 2, 1.5,
            Plots.text("No detected cases", PV_SMALL_FONTSIZE, PV_FONT, :center))
    end

    # COMBINED FIGURE
    row(panels) = plot(panels..., layout = (1, length(panels)))
    items = [
        section_header("Epidemic dynamics"),
        ts_title("A) Daily infections", top[1], order = strategies), top[1],
        ts_title("B) Effective reproduction number"), top[2],
        ts_title("C) Daily tests"), top[3],
        ts_title("D) Individuals in isolation or quarantine"), top[4],
        section_header("Observed epidemic"),
        row(ts_column.(strategies, gemscolors(length(strategies)))),
        ts_title("E) Observed reproduction number", observed[end]), row(observed),
        ts_title("F) Daily detected cases", detected[end]), row(detected),
        ts_title("G) Active cases", darkfigure[end]), row(darkfigure)]
    panel_in = (TS_HEIGHT_IN - 2 * TS_HEADER_IN - TS_COLUMN_IN - 7 * TS_TITLE_IN) / 7
    heights = [TS_HEADER_IN; repeat([TS_TITLE_IN, panel_in], 4);
        TS_HEADER_IN; TS_COLUMN_IN; repeat([TS_TITLE_IN, panel_in], 3)]
    foreach(p -> keep_margins!(p, text_scaling(PV_WIDTH_IN, TS_HEIGHT_IN)), items)
    p = plot(items...,
        layout = grid(length(items), 1, heights = heights ./ sum(heights)),
        size = (100 * PV_WIDTH_IN, 100 * TS_HEIGHT_IN),
        dpi = PV_DPI,
        fontfamily = PV_FONT)

    # same room for the tick labels of all rows, so that the plot areas of a column (and its
    # header) start and end at the same x (subplots: A-D, then headers and E-G by column)
    layout_figure!(p)
    for ks in ([3, 5, 7, 9, 11, 15, 19, 23], [12, 16, 20, 24], [13, 17, 21, 25]), (i, margin) in ((1, :left_margin), (3, :right_margin))
        room = [p[k].minpad[i] for k in ks]
        foreach((k, r) -> p[k][margin] += (maximum(room) - r) / p.attr[:thickness_scaling], ks, room)
    end

    # line at R = 1 of the effective reproduction number
    layout_figure!(p)
    r_threshold_label!(p, add_canvas!(p), 5, halign = :right)
    return p
end


"""
    scenario_testing(;
        configfile = "models/DE_model_first_wave.toml",
        population = "NRW",
        municipality = 05370008, # Gangelt (Kreis Heinsberg)
        start_condition = RegionalSeeds(seeds = Dict(municipality => 50)),
        stop_criterion = TimesUp(limit = 365),
        municipality_contact_reduction::Float64 = 0.5,
        isolation_duration = 14,
        sensitivity = .8,
        specificity = .99,
        runs = 11,
        folder = nothing,
        folder_tag = nothing)

Runs the testing scenario (see `run_testing_sims` for the parameters) `runs` times with the random seeds
0, 1, 2, ... (three per repetition). The `ResultData` object of every run is stored in `folder`, from which
the figure and the tables are generated (see `replot_testing`). Repetitions that are already stored in `folder`
are skipped, so that an interrupted scenario can be continued (with the same parameters).

If `folder` is `nothing`, a new result folder is created, with `folder_tag`
appended to its timestamp-based name (see `result_folder`).

To test the scenario with the smaller population of the Saarland, starting in Beckingen:

```julia
scenario_testing(population = "SL", municipality = 10042111)
```
"""
function scenario_testing(;
    configfile = "models/DE_model_first_wave.toml",
    population = "NRW",
    municipality = 05370008, # Gangelt (Kreis Heinsberg)

    start_condition = RegionalSeeds(seeds = Dict(municipality => 50)),
    stop_criterion = TimesUp(limit = 365),

    municipality_contact_reduction::Float64 = 0.5,
    isolation_duration = 14,
    sensitivity = .8,
    specificity = .99,

    runs = 11,
    folder = nothing,
    folder_tag = nothing)

    # create result folder
    isnothing(folder) ? folder = result_folder(folder_tag) : nothing
    run_folder(seed) = joinpath(folder, "runs", "seed_$(lpad(seed, 3, '0'))")

    # skip repetitions (three runs each) that are already stored
    is_stored(seed) = isfile(joinpath(run_folder(seed), "resultdata.jld2"))
    seeds = filter(seed -> !all(is_stored, seed:(seed + 2)), range(0, step = 3, length = runs))

    # progress output: estimate the remaining time from the mean duration per run
    total = 3 * length(seeds)
    skipped = 3 * runs - total
    run_count = 0
    start = time()

    # store the result data object of each run
    function store_run(rd)
        exportJLD(rd, run_folder(sim_data(rd)["seed"]))
        run_count += 1
        remaining = round(Int, (time() - start) / run_count * (total - run_count))
        printinfo("Testing run $(skipped + run_count) / $(skipped + total) done, estimated time remaining: $(canonicalize(Second(remaining)))\n")
    end

    # run simulations
    for seed in seeds
        run_testing_sims(
            configfile = configfile,
            population = population,
            start_condition = start_condition,
            stop_criterion = stop_criterion,
            municipality_contact_reduction = municipality_contact_reduction,
            isolation_duration = isolation_duration,
            sensitivity = sensitivity,
            specificity = specificity,
            seed = seed,
            after_run = store_run)
    end

    # plot results and store outcomes
    replot_testing(folder)
end

"""
    replot_testing(folder; population = nothing)

Generates the figure and the tables of the testing scenario from the `ResultData` objects
that `scenario_testing()` stored in `folder`, without re-running the simulations:

- `Fig_18_testing_scenarios.png`: figure (see `plot_scenario_testing`)
- `sc_testing_outcomes.csv`: outcomes of every run (see `outcome_df`)
- `sc_testing_daily.csv`: daily data of every run (see `daily_df`)
- `sc_testing_summary.csv`: statistics of the outcomes over the runs of a strategy (see `summary_df`)
- `sc_testing_comparison.csv`: statistics of the comparison of tracing with the other strategies (see `comparison_df`)

For runs that were stored without the number of workers and students, these are taken from the other runs in `folder`
or counted in the `population` (e.g., `"NRW"`), if provided.
"""
function replot_testing(folder; population = nothing)

    # load the result data objects (sorted by seed)
    rds = [import_resultdata(joinpath(f, "resultdata.jld2"))
        for f in readdir(joinpath(folder, "runs"), join = true)]

    # plot results
    save_figure(plot_scenario_testing(rds), :testing_scenarios, folder)

    # number of workers and students for runs that were stored without them
    stored = [sim_data(rd)["members"] for rd in rds if haskey(sim_data(rd), "members")]
    members = !isempty(stored) ? first(stored) :
        isnothing(population) ? Dict() : count_members(Simulation(population = population))

    # extract KPIs and daily data, store as CSV
    outcomes = outcome_df(rds, members = members)
    CSV.write(joinpath(folder, "sc_testing_outcomes.csv"), outcomes)
    CSV.write(joinpath(folder, "sc_testing_summary.csv"), summary_df(outcomes))
    CSV.write(joinpath(folder, "sc_testing_comparison.csv"), summary_df(comparison_df(outcomes)))
    CSV.write(joinpath(folder, "sc_testing_daily.csv"), daily_df(rds))
end
