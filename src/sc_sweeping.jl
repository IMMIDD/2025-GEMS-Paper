export run_sweeps
export plot_sweeps

export scenario_sweeping


"""
    run_sweeps(;
        rep::Int64 = 3,
        step_size::Float64 = 0.1)

Runs the simulations to sweep parameter space.
Returns a dataframe with xvals, yvals, and outvals.
"""
function run_sweeps(;
    rep::Int64 = 3,
    step_size::Float64 = 0.1)

    # observables
    xvals = []
    yvals = []
    outvals = []

    # run all simulations three times
    for i in 1:rep
        # vary symptomatic rate
        for symp in 0.0:step_size:1.0
            # vary adherence
            for adherence in 0.0:step_size:1.0
                # set up scenario
                sim = Simulation(
                    population = "SL",
                    transmission_rate = 0.04,
                    infectious_offset = 0,
                    progression_categories = [1.0-symp, symp, 0, 0])

                self_isolation = IStrategy("Self Isolation", sim, condition = i -> rand() <= adherence)
                add_measure!(self_isolation, SelfIsolation(14))
                trigger = SymptomTrigger(self_isolation)
                add_symptom_trigger!(sim, trigger)

                run!(sim)
                rd = ResultData(sim, style = "LightRD")

                # extract data for heatmap
                push!(xvals, symp)
                push!(yvals, adherence)
                push!(outvals, r0(rd))
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

    return gemsheatmap(df.xvals, df.yvals, df.outvals,
        xlabel = "Rate of Symptomatic Cases",
        ylabel = "Adherence Rate",
        colorbar_title = "\nBasic Reproduction Number",
        xrev = true,
        yrev = true,
        xformatter = (x -> "$(Int(round(100x)))%"),
        yformatter = (x -> "$(Int(round(100x)))%"),
        color = :r0,
        size = (1600, 1000),
        tickfontsize = 24,
        labelfontsize = 24,
        colorbar_fontsize = 24,
        colorbar_titlefontsize = 24,
        margin = 10mm,
        right_margin = 18mm)

end

"""
    scenario_sweeping(;
        rep = 3,
        folder = nothing)

Runs the parameter sweep scenario and all analyses.
"""
function scenario_sweeping(;
    rep = 3, # repetitions per combination
    step_size::Float64 = 0.1,
    folder = nothing)

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder() : nothing    

    # run simulations
    sweeps = run_sweeps(rep = rep, step_size = step_size)
        
    # store df as CSV
    CSV.write(joinpath(folder, "sc_sweep_data.csv"), sweeps)

    # generate plot
    p_sw = plot_sweeps(sweeps)
    png(p_sw, joinpath(folder, "sc_sweep_heatmap.png"))
end