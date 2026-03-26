export run_sims
export plot_timings
export plot_benchmarks
export benchmark

# DATA GENERATION

"""
    mean_infection_duration(rd)

Calculates mean duration of infection (exposure to removal).
"""
mean_infection_duration(rd) = rd |> infections |> df -> df.removed_tick .- df.tick |> mean

"""
    run_sims()

Runs benchmarking simulations for each German state and the complete German model
as well as synthetic populations of various sizes.

Returns a DataFrame with timings.
"""
function run_sims(num_runs::Int = 11)
    # run everything once to pre-compile simulation
    sim = Simulation()
    run!(sim)
    rd = ResultData(sim)

    res = DataFrame(
        model = String[],
        type = String[],
        pop_size = Float64[],       # Changed to Float64 for averages
        infections = Float64[],     # Changed to Float64 for averages
        mean_infection_duration = Float64[],
        initialization_time = Float64[],
        simulation_time = Float64[],
        post_processing_time = Float64[])

    pop_codes = [
        "BW", # Baden-Württemberg
        "BY", # Bayern (Bavaria)
        "BE", # Berlin
        "BB", # Brandenburg
        "HB", # Bremen
        "HH", # Hamburg
        "HE", # Hessen
        "MV", # Mecklenburg-Vorpommern
        "NI", # Niedersachsen (Lower Saxony)
        "NRW", # Nordrhein-Westfalen (North Rhine–Westphalia)
        "RP", # Rheinland-Pfalz (Rhineland-Palatinate)
        "SL", # Saarland
        "SN", # Sachsen (Saxony)
        "ST", # Sachsen-Anhalt (Saxony-Anhalt)
        "SH", # Schleswig-Holstein
        "TH", # Thüringen (Thuringia)
        "DE"  # Germany
    ]

    synth_sizes = [
        100_000,
        1_000_000,
        5_000_000,
        10_000_000,
        20_000_000,
        50_000_000,
        #85_000_000,
    ]

    # run German models
    for p in pop_codes
        # Initialize accumulators for this model
        pop_sizes = Float64[]
        infs = Float64[]
        durs = Float64[]
        init_ts = Float64[]
        sim_ts = Float64[]
        proc_ts = Float64[]

        for i in 1:num_runs
            init_t = @timed Simulation(population = p, seed = i)
            sim = init_t.value # get sim object
            sim_t = @timed run!(sim)
            proc_t = @timed ResultData(sim)

            push!(pop_sizes, length(individuals(sim)))
            push!(infs, total_infections(proc_t.value))
            push!(durs, mean_infection_duration(proc_t.value))

            push!(init_ts, init_t.time)
            push!(sim_ts, sim_t.time)
            push!(proc_ts, proc_t.time)
        end

        push!(res, (
            model = p,
            type = "Gesyland",
            pop_size = median(pop_sizes),
            infections = median(infs),
            mean_infection_duration = median(durs),
            initialization_time = median(init_ts),
            simulation_time = median(sim_ts),
            post_processing_time = median(proc_ts)))
    end

    # run synthetic models
    for s in synth_sizes
        # Initialize accumulators for this model
        pop_sizes = Float64[]
        infs = Float64[]
        durs = Float64[]
        init_ts = Float64[]
        sim_ts = Float64[]
        proc_ts = Float64[]

        for i in 1:num_runs
            init_t = @timed Simulation(pop_size = s, seed = i)
            sim = init_t.value # get sim object
            sim_t = @timed run!(sim)
            proc_t = @timed ResultData(sim)

            push!(pop_sizes, length(individuals(sim)))
            push!(infs, total_infections(proc_t.value))
            push!(durs, mean_infection_duration(proc_t.value))

            push!(init_ts, init_t.time)
            push!(sim_ts, sim_t.time)
            push!(proc_ts, proc_t.time)
        end

        push!(res, (
            model = "Size_$s",
            type = "Synthetic",
            pop_size = median(pop_sizes),
            infections = median(infs),
            mean_infection_duration = median(durs),
            initialization_time = median(init_ts),
            simulation_time = median(sim_ts),
            post_processing_time = median(proc_ts)))
    end

    return res
end

# PLOTTING

"""
    plot_timings(x, y, title, xlab, ylab)

generates a single benchmarking plot.
"""
function plot_timings(x, y, types, title, xlab, ylab)
    return scatter(
        log10.(x), log10.(y),
        title = title,
        color = (t -> (t == "Gesyland" ? :tomato : :deepskyblue)).(types),
        markerstrokewidth=0,
        xticks = ([5, 6, 7, 8], ["100k", "1m", "10m", "100m"]),
        yticks = ([-1, 0, 1, 2, 3], ["0.1s", "1s", "10s", "100s", "1000s"]),
        xlims = (4.9, 8.1),
        ylims = (-1, 3),
        xlabel = xlab,
        ylabel = ylab,
        xtickfontsize = 14,
        ytickfontsize = 14,
        xlabelfontsize = 14,
        ylabelfontsize = 14,
        titlefontsize = 20,
        markersize=7,
        left_margin = 10mm,
        right_margin = 5mm,
        bottom_margin = 10mm,
        fontfamily = "Times Roman",
        legend = nothing)
end

"""
    plot_benchmarks(df)

Returns a multiplot for initialization, simulation, and post-processing benchmarks.
Requires a timings dataframe as input (which is the output of the `run_sims()` function).
"""
function plot_benchmarks(df)
    init_times = plot_timings(
        df.pop_size, df.initialization_time, df.type,
        "Initialization Time", "Population Size", "Runtime")
    sim_times = plot_timings(
        df.pop_size, df.simulation_time, df.type,
        "Simulation Time", "Population Size", "Runtime")
    proc_times = plot_timings(
        df.pop_size, df.post_processing_time, df.type,
        "Post-Processing Time", "Population Size", "Runtime")
        
    return plot(init_times, sim_times, proc_times, layout = (1, 3), size = (1800, 600))
end


# RUN ALL ANALYSES

"""
    benchmark(; folder = nothing)

Runs all benchmark-related analyses.
"""
function benchmark(; folder = nothing)
    # create result folder if not passed
    isnothing(folder) ? folder = result_folder() : nothing
    
    # run benchmark simulations
    df = run_sims()

    # export dataframe
    CSV.write(joinpath(folder, "benchmarks.csv"), df)

    # plot benchmarks
    bench_p = plot_benchmarks(df)
    png(bench_p, joinpath(folder, "benchmarks.png"))
end