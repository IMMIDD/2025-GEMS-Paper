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
    step_size::Float64 = 0.025)

    # observables
    xvals = []
    yvals = []
    outvals = []

    dp_a = Asymptomatic(
        exposure_to_infectiousness_onset = 0,
        infectiousness_onset_to_recovery = 7  
    )

    dp_s = Symptomatic(
        exposure_to_infectiousness_onset = 0,
        infectiousness_onset_to_symptom_onset = 2, 
        symptom_onset_to_recovery = 7             
    )

    # run all simulations three times
    for i in 1:rep
        # vary symptomatic rate
        for symp in 0.0:step_size:1.0
            
            pa = AgeBasedProgressionAssignment(
                age_groups = ["0-120"], 
                progression_categories = ["Asymptomatic", "Symptomatic"],
                stratification_matrix = [[1.0 - symp, symp]] 
            )

            my_pathogen = Pathogen(
                name = "My Pathogen",
                progressions = [dp_a, dp_s],
                progression_assignment = pa,
                transmission_function = ConstantTransmissionRate(transmission_rate = 0.04)
            )

            # vary adherence
            for adherence in 0.0:step_size:1.0
                
                sim = Simulation(
                    population = "SL",
                    pathogen = my_pathogen,
                    seed = i
                )

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

    formattick(x) = Int(round(100x)) |>
        x -> x % 5 == 0 ? "$x %" : ""

    return gemsheatmap(df.xvals, df.yvals, df.outvals,
        xlabel = "Rate of Symptomatic Cases",
        ylabel = "Adherence Rate",
        colorbar_title = "\nBasic Reproduction Number",
        xrev = true,
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
        rep = 3,
        folder = nothing)

Runs the parameter sweep scenario and all analyses.
"""
function scenario_sweeping(;
    rep = 3, # repetitions per combination
    step_size::Float64 = 0.025,
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