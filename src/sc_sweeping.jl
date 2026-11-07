export set_asymptomatic_share!
export run_sweeps
export plot_sweeps
export sweep_variability
export write_sweep_variability

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
    sweep_variability(df)

Summarizes the stochastic variability of the outcome between the replicates
of each parameter combination (grid cell). Requires the input of `run_sweeps()`
with at least two replicates per combination.

Per cell, the sample standard deviation (n - 1) across replicates is computed.
These are summarized across the grid as the pooled SD (root of the mean
per-cell variance), the mean, 95th percentile, and maximum SD. The standard
error of the plotted cell means is the pooled SD divided by √n.
Also counts the cells whose replicates lie on both sides of R = 1.

Returns a NamedTuple with the results.
"""
function sweep_variability(df)

    cells = combine(groupby(df, [:xvals, :yvals]),
        :outvals => length => :n,
        :outvals => std => :sd,
        :outvals => (r -> minimum(r) < 1 < maximum(r)) => :straddles)

    n = minimum(cells.n)
    n >= 2 || throw(ArgumentError("At least two replicates per combination are required to compute variability. Got $n."))
    n == maximum(cells.n) || @warn "Number of replicates differs between combinations ($n to $(maximum(cells.n)))."

    pooled_sd = sqrt(mean(cells.sd .^ 2))

    return (
        cells = nrow(cells),
        replicates = n,
        pooled_sd = pooled_sd,
        mean_sd = mean(cells.sd),
        p95_sd = quantile(cells.sd, 0.95),
        max_sd = maximum(cells.sd),
        se_mean = pooled_sd / sqrt(n),
        straddling_cells = count(cells.straddles))
end

"""
    write_sweep_variability(df, path)

Computes the variability summary of `sweep_variability(df)` and writes it as
a plain text file to `path`.
"""
function write_sweep_variability(df, path)

    v = sweep_variability(df)
    r(x) = round(x, digits = 3)

    open(path, "w") do io
        println(io, "Parameter sweep: stochastic variability of R_init between replicates")
        println(io, "")
        println(io, "Grid cells (parameter combinations):     $(v.cells)")
        println(io, "Replicates per cell (n):                 $(v.replicates)")
        println(io, "")
        println(io, "Per-cell sample SD (n - 1) across replicates, summarized over all cells:")
        println(io, "  Pooled SD (sqrt of mean variance):     $(r(v.pooled_sd))")
        println(io, "  Mean SD:                               $(r(v.mean_sd))")
        println(io, "  95th percentile SD:                    $(r(v.p95_sd))")
        println(io, "  Maximum SD:                            $(r(v.max_sd))")
        println(io, "")
        println(io, "Standard error of cell means (pooled SD / sqrt(n)): $(r(v.se_mean))")
        println(io, "Cells with replicates on both sides of R = 1:       $(v.straddling_cells) ($(round(100 * v.straddling_cells / v.cells, digits = 1)) %)")
    end

    return v
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

    # summarize variability between replicates (requires rep >= 2)
    rep >= 2 ?
        write_sweep_variability(sweeps, joinpath(folder, "sc_sweep_variability.txt")) :
        @warn "Skipping variability summary, as it requires at least two repetitions per combination."

    # generate plot
    p_sw = plot_sweeps(sweeps)
    png(p_sw, joinpath(folder, "sc_sweep_heatmap.png"))
end
