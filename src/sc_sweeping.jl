export set_asymptomatic_share!
export run_sweeps
export plot_sweeps

export scenario_sweeping


"""
    set_asymptomatic_share!(sim::Simulation, share::Float64)

Sets the share of asymptomatic progressions in all age groups of the
simulation's pathogen to `share`. All other progression categories are
rescaled proportionally so that each age group still sums up to 1.
Requires the pathogen to use an `AgeBasedProgressionAssignment`.

**Attention**: The initial infections are already assigned their progressions
when the simulation is created. Call `reset!(sim)` afterwards to re-apply the
start condition with the new shares.
"""
function set_asymptomatic_share!(sim::Simulation, share::Float64)

    0.0 <= share <= 1.0 || throw(ArgumentError("Asymptomatic share must be between 0 and 1. Got $share."))

    pa = progression_assignment(pathogen(sim))
    pa isa AgeBasedProgressionAssignment ||
        throw(ArgumentError("The pathogen must use an AgeBasedProgressionAssignment. Got $(typeof(pa))."))

    asymp = findfirst(==(Asymptomatic), pa.progression_categories)
    isnothing(asymp) && throw(ArgumentError("The pathogen has no Asymptomatic progression category."))

    # matrix is stored as categories x age groups
    mtrx = pa.stratification_matrix
    others = setdiff(axes(mtrx, 1), asymp)

    for ag in axes(mtrx, 2)
        rest = 1.0 - mtrx[asymp, ag]
        rest > 0.0 || throw(ArgumentError("Cannot rescale age group $ag as it only contains asymptomatic progressions."))

        # scale all non-asymptomatic categories up (or down) to the remaining share
        mtrx[others, ag] .*= (1.0 - share) / rest
        mtrx[asymp, ag] = share
    end

    return sim
end

"""
    run_sweeps(;
        configfile = "models/DE_model_first_wave.toml",
        population = "SL",
        infected_fraction = 0.001,
        stop_criterion = TimesUp(limit = 365),
        transmission_reduction::Float64 = 0.6,
        rep::Int64 = 3,
        step_size::Float64 = 0.025)

Runs the simulations to sweep parameter space (share of asymptomatic
progressions and self-isolation adherence).
To focus on the regime where self-isolation can be decisive, the scenario assumes
that other measures already prevent a share of all transmissions
(`transmission_reduction`, by default 60 %), such that without isolation R_init
lies slightly above 1. This reduces the `transmission_scale` of the pathogen's
`InfectiousnessProfile` accordingly.
Returns a dataframe with xvals, yvals, and outvals.
"""
function run_sweeps(;
    configfile = "models/DE_model_first_wave.toml",
    population = "SL",
    infected_fraction = 0.001,
    stop_criterion = TimesUp(limit = 365),
    transmission_reduction::Float64 = 0.6,
    rep::Int64 = 3,
    step_size::Float64 = 0.025)

    0.0 <= transmission_reduction <= 1.0 || throw(ArgumentError("Transmission reduction must be between 0 and 1. Got $transmission_reduction."))

    # observables
    xvals = []
    yvals = []
    outvals = []

    # each simulation gets its own (reproducible) seed
    sim_count = 0

    # run all simulations three times
    for _ in 1:rep
        # vary asymptomatic rate
        for asymp in 0.0:step_size:1.0

            # vary adherence
            for adherence in 0.0:step_size:1.0

                sim_count += 1

                sim = Simulation(
                    configfile = configfile,
                    population = population,
                    seed = sim_count,
                    infected_fraction = infected_fraction,
                    stop_criterion = stop_criterion)

                set_asymptomatic_share!(sim, asymp)

                # assume that other measures already prevent a share of all transmissions
                # (scales the transmission_scale of the config's InfectiousnessProfile)
                tf = transmission_function(pathogen(sim))
                pathogen(sim).transmission_function = with_transmission_scale(tf, tf.transmission_scale * (1.0 - transmission_reduction))

                # re-apply the start condition, as the initial infections were
                # already assigned their progressions when the simulation was created
                reset!(sim)

                # draw adherence from the simulation's (seeded, thread-local) RNG
                self_isolation = IStrategy("Self Isolation", sim, condition = ind -> rand(rng(sim)) <= adherence)
                add_measure!(self_isolation, SelfIsolation(14))
                trigger = SymptomTrigger(self_isolation)
                add_symptom_trigger!(sim, trigger)

                run!(sim)

                # extract data for heatmap (R0 directly from the post processor,
                # as building the full ResultData object is not needed)
                push!(xvals, asymp)
                push!(yvals, adherence)
                push!(outvals, r0(PostProcessor(sim)))
            end
        end
    end

    return DataFrame(
        xvals = xvals,
        yvals = yvals,
        outvals = outvals)
end

"""
    plot_sweeps(df)

Generates a GEMS heatmap from a dataframe with xvals, yvals, and outvals.
Requires the input of `run_sweeps()`.
"""
function plot_sweeps(df)

    formattick(x) = Int(round(100x)) |>
        x -> x % 5 == 0 ? "$x %" : ""

    return gemsheatmap(df.xvals, df.yvals, df.outvals,
        xlabel = "Fraction of Asymptomatic Cases",
        ylabel = "Probability of Adherence",
        colorbar_title = "\nInitial Effective Reproduction Number",
        xrev = false,
        yrev = true,
        xformatter = formattick,
        yformatter = formattick,
        color = :r0,
        size = (2000, 1200),
        tickfontsize = 20,
        labelfontsize = 24,
        colorbar_fontsize = 24,
        colorbar_titlefontsize = 24,
        xrotation = 45,
        margin = 15mm,
        right_margin = 18mm)

end

"""
    scenario_sweeping(;
        configfile = "models/DE_model_first_wave.toml",
        population = "SL",
        infected_fraction = 0.001,
        stop_criterion = TimesUp(limit = 365),
        transmission_reduction::Float64 = 0.6,
        rep = 3,
        step_size::Float64 = 0.025,
        folder = nothing)

Runs the parameter sweep scenario and all analyses.
"""
function scenario_sweeping(;
    configfile = "models/DE_model_first_wave.toml",
    population = "SL",
    infected_fraction = 0.001,
    stop_criterion = TimesUp(limit = 365),
    transmission_reduction::Float64 = 0.6, # share of transmissions prevented by other measures
    rep = 3, # repetitions per combination
    step_size::Float64 = 0.025,
    folder = nothing)

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder() : nothing

    # run simulations
    sweeps = run_sweeps(
        configfile = configfile,
        population = population,
        infected_fraction = infected_fraction,
        stop_criterion = stop_criterion,
        transmission_reduction = transmission_reduction,
        rep = rep,
        step_size = step_size)

    # store df as CSV
    CSV.write(joinpath(folder, "sc_sweep_data.csv"), sweeps)

    # generate plot
    p_sw = plot_sweeps(sweeps)
    png(p_sw, joinpath(folder, "sc_sweep_heatmap.png"))
end
