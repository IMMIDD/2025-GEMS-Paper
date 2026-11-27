export run_sims
export run_thread_speedup
export run_peak_memory
export summarize_benchmarks
export plot_benchmarks
export plot_thread_speedup
export plot_runtime_per_agent
export replot_benchmarks
export benchmark

# DATA GENERATION

"""
    BENCHMARK_POPULATIONS

GEMS population models of the 16 German states and Germany.
"""
const BENCHMARK_POPULATIONS = ["BW", "BY", "BE", "BB", "HB", "HH", "HE", "MV", "NI", "NRW", "RP", "SL", "SN", "ST", "SH", "TH", "DE"]

"""
    BENCHMARK_SYNTHETIC_SIZES

Sizes of the on-the-fly populations.
"""
const BENCHMARK_SYNTHETIC_SIZES = [100_000, 1_000_000, 5_000_000, 10_000_000, 20_000_000, 50_000_000]

"""
    BENCHMARK_CONFIG

Config of the benchmark: the reference model, run for 365 days with 0.1 % initial infections.
"""
const BENCHMARK_CONFIG = "models/DE_model_first_wave.toml"

"""
    machine_info()

Returns a one-row DataFrame with the hardware, Julia and GEMS versions, and Julia threads of this process.
"""
function machine_info()
    return DataFrame(
        cpu = strip(Sys.cpu_info()[1].model),
        physical_cores = try GEMS.CpuId.cpucores() catch; missing end,
        cpu_threads = Sys.CPU_THREADS,
        memory_gb = round(Sys.total_memory() / 2^30, digits = 1),
        os = string(Sys.KERNEL),
        julia_version = string(VERSION),
        gems_version = string(pkgversion(GEMS)),
        julia_threads = Threads.nthreads())
end

# keyword arguments, name, and type of a population: a population model code or an on-the-fly size
_population_kwargs(p::AbstractString) = (population = p,)
_population_kwargs(p::Integer) = (pop_size = p,)
_population_name(p::AbstractString) = p
_population_name(p::Integer) = "Size_$p"
_population_type(::AbstractString) = "Gesyland"
_population_type(::Integer) = "Synthetic"

# reference model with 0.1 % initial infections, run for 365 days
_benchmark_simulation(; kwargs...) = Simulation(; configfile = BENCHMARK_CONFIG, infected_fraction = 0.001,
    stop_criterion = TimesUp(limit = 365), kwargs...)

# compiles both population paths (on-the-fly and population model) before anything is timed
function _warm_up()
    for kwargs in ((;), (population = "SL",))
        sim = _benchmark_simulation(; kwargs...)
        run!(sim)
        ResultData(sim)
    end
    GC.gc()
end

# one timed run; initialization covers reading the (already downloaded) population and settings
# files up to the ready simulation; live memory is measured after a full GC outside the timings
function _timed_run(p, seed)
    GC.gc()
    init_t = @timed _benchmark_simulation(; seed = seed, _population_kwargs(p)...)
    sim = init_t.value
    sim_t = @timed run!(sim)
    GC.gc()
    sim_live = Base.gc_live_bytes()
    proc_t = @timed ResultData(sim)
    GC.gc()
    post_live = Base.gc_live_bytes()
    # post-processing with the LightRD style (fewer analyses), timed separately
    light_t = @timed ResultData(sim; style = "LightRD")
    return (
        model = _population_name(p),
        type = _population_type(p),
        run = seed,
        threads = Threads.nthreads(),
        pop_size = length(individuals(sim)),
        infections = total_infections(proc_t.value),
        initialization_time = init_t.time,
        simulation_time = sim_t.time,
        post_processing_time = proc_t.time,
        post_processing_light_time = light_t.time,
        compile_time = init_t.compile_time + sim_t.compile_time + proc_t.compile_time,
        gc_time = init_t.gctime + sim_t.gctime + proc_t.gctime,
        allocated_gb = (init_t.bytes + sim_t.bytes + proc_t.bytes) / 2^30,
        live_gb_simulation = sim_live / 2^30,
        live_gb_post_processing = post_live / 2^30)
end

"""
    run_sims(; populations = BENCHMARK_POPULATIONS, synthetic_sizes = BENCHMARK_SYNTHETIC_SIZES, num_runs = 11)

Runs the benchmark `num_runs` times (random seeds 1 to `num_runs`) for each population model and
on-the-fly population size, after a warm-up. Returns a DataFrame with one row per run.
"""
function run_sims(; populations = BENCHMARK_POPULATIONS, synthetic_sizes = BENCHMARK_SYNTHETIC_SIZES, num_runs::Int = 11)
    # download population files before timing anything
    foreach(GEMS.obtain_remote_files, populations)
    _warm_up()
    rows = [_timed_run(p, seed) for p in vcat(populations, synthetic_sizes) for seed in 1:num_runs]
    return DataFrame(rows)
end

# timed runs of one Julia process for `run_thread_speedup`, written to `out`
function _speedup_worker(populations, num_runs, out)
    foreach(GEMS.obtain_remote_files, filter(p -> p isa AbstractString, populations))
    _warm_up()
    CSV.write(out, DataFrame([_timed_run(p, seed) for p in populations for seed in 1:num_runs]))
end

"""
    run_thread_speedup(; populations = ["RP"], threads = [1, 2, 4, 8, 16, 32, 64, 128], num_runs = 5, folder = mktempdir())

Runs the benchmark for each population with each number of threads, in one Julia process per
thread count (the thread count is fixed at startup). Returns a DataFrame with one row per run.
"""
function run_thread_speedup(; populations = ["RP"], threads = [1, 2, 4, 8, 16, 32, 64, 128], num_runs::Int = 5, folder = mktempdir())
    project = dirname(Base.active_project())
    files = String[]
    for n in threads
        out = joinpath(folder, "speedup_$(n)_threads.csv")
        code = "using GEMSPaper; GEMSPaper._speedup_worker($(repr(populations)), $num_runs, $(repr(out)))"
        run(`$(Base.julia_cmd()) --project=$project --threads=$n -e $code`)
        push!(files, out)
    end
    return reduce(vcat, [CSV.read(f, DataFrame) for f in files])
end

# peak memory of one run in a fresh process, written to `out`
function _peak_memory_worker(population, out)
    sim = _benchmark_simulation(; _population_kwargs(population)...)
    run!(sim)
    rss_sim = Sys.maxrss()
    ResultData(sim)
    CSV.write(out, DataFrame(
        model = _population_name(population),
        pop_size = length(individuals(sim)),
        threads = Threads.nthreads(),
        peak_rss_gb_simulation = rss_sim / 2^30,
        peak_rss_gb_post_processing = Sys.maxrss() / 2^30))
end

"""
    run_peak_memory(population = "DE"; threads = Threads.nthreads(), folder = mktempdir())

Peak resident memory of one run (simulation and post-processing) of `population` in a fresh
Julia process, as the RAM a machine needs. Returns a one-row DataFrame.
"""
function run_peak_memory(population = "DE"; threads::Int = Threads.nthreads(), folder = mktempdir())
    population isa AbstractString && GEMS.obtain_remote_files(population)
    out = joinpath(folder, "peak_memory.csv")
    code = "using GEMSPaper; GEMSPaper._peak_memory_worker($(repr(population)), $(repr(out)))"
    run(`$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) --threads=$threads -e $code`)
    return CSV.read(out, DataFrame)
end

"""
    summarize_benchmarks(runs)

Median and interquartile range of each stage per population and thread count.
"""
function summarize_benchmarks(runs::DataFrame)
    q(x, p) = quantile(x, p)
    return combine(groupby(runs, [:model, :type, :threads]),
        :pop_size => median => :pop_size,
        :infections => median => :infections,
        nrow => :runs,
        [s => f => Symbol("$(s)_$(n)") for s in [:initialization_time, :simulation_time, :post_processing_time, :post_processing_light_time]
            for (n, f) in [("median", median), ("q25", x -> q(x, 0.25)), ("q75", x -> q(x, 0.75))]]...,
        :live_gb_simulation => median => :live_gb_simulation,
        :live_gb_post_processing => median => :live_gb_post_processing)
end

# PLOTTING

# stages in the GEMS colors, as in the other figures
const _STAGES = [
    (:initialization_time, "Initialization", FIG_COLORS.blue),
    (:simulation_time, "Simulation", FIG_COLORS.green),
    (:post_processing_time, "Post-processing", FIG_COLORS.red)]

# fonts and resolution of the shared figure style (`figure_style.jl`)
_bench_style() = (
    fontfamily = FIG_FONT,
    tickfontsize = FIG_SMALL_FONTSIZE,
    guidefontsize = FIG_FONTSIZE,
    titlefontsize = FIG_HEADING_FONTSIZE,
    titlelocation = :left,
    legendfontsize = FIG_SMALL_FONTSIZE,
    framestyle = :axes,
    grid = true,
    gridalpha = 0.15,
    gridlinewidth = 0.5,
    tick_direction = :in,
    dpi = FIG_DPI)
const _BENCH_MARKER = (markersize = 3.5, markerstrokewidth = 1)

# log ticks within [lo, hi] with readable labels; adds 2 and 5 if fewer than three decades fit
function _decade_ticks(lo, hi, label)
    ks = floor(Int, log10(lo)):ceil(Int, log10(hi))
    ticks = filter(t -> lo <= t <= hi, [10.0^k for k in ks])
    if length(ticks) < 3
        ticks = filter(t -> lo <= t <= hi, [m * 10.0^k for k in ks for m in (1, 2, 5)])
    end
    return (ticks, label.(ticks))
end
_size_label(x) = x >= 1e6 ? "$(round(Int, x / 1e6))M" : "$(round(Int, x / 1e3))k"
_seconds_label(x) = x >= 1 ? "$(round(Int, x))" : "$(round(x, sigdigits = 1))"

# least-squares slope of log(y) against log(x)
function _loglog_slope(x, y)
    lx, ly = log.(x), log.(y)
    return sum((lx .- mean(lx)) .* (ly .- mean(ly))) / sum((lx .- mean(lx)) .^ 2)
end

# drops the on-the-fly populations if `on_the_fly` is false
_select_populations(summary, on_the_fly) = sort(on_the_fly ? summary : summary[summary.type .!= "Synthetic", :], :pop_size)

# legend entries for the marker fill (population models filled, on-the-fly open)
function _fill_legend!(p, summary)
    any(summary.type .== "Gesyland") && scatter!(p, [NaN], [NaN]; _BENCH_MARKER..., markercolor = :gray, markerstrokecolor = :gray, label = "German states")
    any(summary.type .== "Synthetic") && scatter!(p, [NaN], [NaN]; _BENCH_MARKER..., markercolor = :white, markerstrokecolor = :gray, label = "On-the-fly")
end

# runtime of one stage against the population size, with shared axis limits
function _stage_panel(d, stage, title, color, xlims, ylims; legend = false, xlabel = "Population size", ylabel = "Runtime (s)")
    p = plot(;
        xscale = :log10, yscale = :log10, xlims = xlims, ylims = ylims,
        xticks = _decade_ticks(xlims..., _size_label), yticks = _decade_ticks(ylims..., _seconds_label),
        xlabel = xlabel, ylabel = ylabel, title = title, legend = legend ? :topleft : false,
        _bench_style()...)
    col = "$(stage)_median"
    for (type, name, filled) in (("Gesyland", "German states", true), ("Synthetic", "On-the-fly", false))
        g = d[d.type .== type, :]
        nrow(g) == 0 && continue
        y = g[!, col]
        scatter!(p, g.pop_size, y; _BENCH_MARKER...,
            yerror = (y .- g[!, "$(stage)_q25"], g[!, "$(stage)_q75"] .- y),
            markercolor = filled ? color : :white, markerstrokecolor = color, linecolor = color, label = name)
    end

    # slope of the population models (of all populations if there are none)
    fit = any(d.type .== "Gesyland") ? d[d.type .== "Gesyland", :] : d
    nrow(fit) > 1 && annotate!(p, xlims[2] / 1.6, ylims[1] * 1.6,
        text("slope $(round(_loglog_slope(fit.pop_size, fit[!, col]), digits = 2))", FIG_SMALL_FONTSIZE, FIG_FONT, :right))

    de = d[d.model .== "DE", :]
    nrow(de) > 0 && annotate!(p, de.pop_size[1] / 1.1, de[1, col] * 1.5, text("Germany", FIG_SMALL_FONTSIZE, FIG_FONT, :right))
    return p
end

"""
    plot_benchmarks(summary; speedup = nothing, physical_cores = missing, on_the_fly = true)

Panels A–C: runtime of each stage against the population size (log-log, shared axes; median,
interquartile range, and the log-log slope of the population models). Panel D: `plot_thread_speedup`
if `speedup` is given. Requires the output of `summarize_benchmarks`.
"""
function plot_benchmarks(summary::DataFrame; speedup = nothing, physical_cores = missing, on_the_fly::Bool = true)
    d = _select_populations(summary, on_the_fly)
    ys = vcat([vcat(d[!, "$(s)_q25"], d[!, "$(s)_q75"]) for (s, _, _) in _STAGES]...)
    xlims = (minimum(d.pop_size) / 1.4, maximum(d.pop_size) * 1.4)
    ylims = (minimum(ys) / 1.6, maximum(ys) * 1.8)
    grid = !isnothing(speedup)

    # axis titles only where the neighbouring panel does not repeat them
    ps = [_stage_panel(d, s, "$('A' + i - 1)) $title", color, xlims, ylims;
            legend = i == 1,
            xlabel = grid && i <= 2 ? "" : "Population size",
            ylabel = (grid ? i == 2 : i > 1) ? "" : "Runtime (s)")
        for (i, (s, title, color)) in enumerate(_STAGES)]

    # 2x2 at text width (about the size of four agreement plots of the population validation figure);
    # without the speedup panel, a row at page width
    grid && push!(ps, plot_thread_speedup(speedup; physical_cores = physical_cores, on_the_fly = on_the_fly,
        title = "D) Speedup by threads"))

    # panel letters in bold, as in the other figures (plot areas about 2 in square)
    foreach(p -> (panel_title!(p, 2.0, 2.0); plot!(p, top_margin = 5mm)), ps)
    grid || return plot(ps..., layout = (1, 3), size = (100FIG_FULL_WIDTH_IN, 260), dpi = FIG_DPI, right_margin = 2mm)
    return plot(ps..., layout = (2, 2), size = (100FIG_TEXT_WIDTH_IN, 500), dpi = FIG_DPI, right_margin = 2mm)
end

# TEMPORARY: single-panel alternative to `plot_benchmarks` for comparison; delete once the figure is decided
function plot_benchmarks_single_panel(summary::DataFrame; on_the_fly::Bool = true)
    d = _select_populations(summary, on_the_fly)
    stage_cols = ["$(s)_median" for (s, _, _) in _STAGES]
    ymin = minimum(minimum(d[!, c]) for c in stage_cols)
    ymax = maximum(maximum(d[!, c]) for c in stage_cols)
    p = plot(;
        xscale = :log10, yscale = :log10,
        xticks = _decade_ticks(minimum(d.pop_size), maximum(d.pop_size), _size_label),
        yticks = _decade_ticks(ymin, ymax, _seconds_label), ylims = (ymin / 1.5, ymax * 2.5),
        xlabel = "Population size", ylabel = "Runtime (s)",
        title = "Runtime by population size",
        legend = :topleft, size = (375, 280),
        _bench_style()...)

    # the three stages of each population joined by a grey line
    for r in eachrow(d)
        y = [r[c] for c in stage_cols]
        plot!(p, [r.pop_size, r.pop_size], [minimum(y), maximum(y)], color = :gray70, linewidth = 1, label = "")
    end

    for (s, label, color) in _STAGES
        for (type, filled) in (("Gesyland", true), ("Synthetic", false))
            g = d[d.type .== type, :]
            nrow(g) == 0 && continue
            y = g[!, "$(s)_median"]
            scatter!(p, g.pop_size, y; _BENCH_MARKER...,
                yerror = (y .- g[!, "$(s)_q25"], g[!, "$(s)_q75"] .- y),
                markercolor = filled ? color : :white,
                markerstrokecolor = color, linecolor = color,
                label = (filled || !any(d.type .== "Gesyland")) ? label : "")
        end
    end
    _fill_legend!(p, d)

    # label Germany
    de = d[d.model .== "DE", :]
    if nrow(de) > 0
        annotate!(p, de.pop_size[1] / 1.1, de.simulation_time_median[1] * 1.6, text("Germany", FIG_SMALL_FONTSIZE, FIG_FONT, :right))
    end
    return p
end

"""
    plot_thread_speedup(summary; physical_cores = missing, on_the_fly = true, title = "Speedup by number of threads")

Speedup of each stage (colour) against the number of threads for each population (population
models filled and solid, on-the-fly populations open and dashed), with the linear speedup and, if
given, a vertical line at the number of physical cores. Requires the output of
`summarize_benchmarks` for runs with several thread counts.
"""
function plot_thread_speedup(summary::DataFrame; physical_cores = missing, on_the_fly::Bool = true,
    title = "Speedup by number of threads")
    summary = _select_populations(summary, on_the_fly)
    threads = sort(unique(summary.threads))
    p = plot(;
        xscale = :log2, yscale = :log2,
        xticks = (threads, string.(threads)), yticks = (threads, string.(threads)),
        xlabel = "Threads", ylabel = "Speedup factor", title = title,
        legend = :topleft, size = (375, 280),
        _bench_style()...)
    # reference lines labelled directly, in the free space below the linear speedup
    n = maximum(threads)
    plot!(p, threads, threads, linestyle = :dot, color = :gray, label = "")
    annotate!(p, n / 1.05, n / 2.1, text("linear", FIG_SMALL_FONTSIZE, FIG_FONT, :gray, :right))
    if !ismissing(physical_cores) && physical_cores < n
        vline!(p, [physical_cores], color = :black, linestyle = :dashdot, linewidth = 1, label = "")
        annotate!(p, physical_cores * 1.08, n / 3.2, text("$physical_cores cores", FIG_SMALL_FONTSIZE, FIG_FONT, :left))
    end
    for (_, label, color) in _STAGES
        plot!(p, [NaN], [NaN], color = color, linewidth = 1.5, label = label)
    end

    # one marker shape per population, filled for population models and open for on-the-fly
    shapes = [:circle, :rect, :diamond, :utriangle, :dtriangle]
    groups = groupby(sort(summary, :threads), :model)
    for (i, g) in enumerate(groups)
        filled = g.type[1] == "Gesyland"
        style = (; _BENCH_MARKER..., linestyle = filled ? :solid : :dash, marker = shapes[mod1(i, length(shapes))])
        base = g[g.threads .== minimum(g.threads), :]
        for (s, _, color) in _STAGES
            col = "$(s)_median"
            plot!(p, g.threads, base[1, col] ./ g[!, col]; style...,
                color = color, linewidth = 1, markercolor = filled ? color : :white, markerstrokecolor = color, label = "")
        end
        # legend entry per population only if there are several
        name = filled ? g.model[1] : "$(_size_label(g.pop_size[1])) on-the-fly"
        length(groups) > 1 && plot!(p, [NaN], [NaN]; style..., color = :gray, markercolor = filled ? :gray : :white,
            markerstrokecolor = :gray, label = name)
    end
    return p
end

"""
    plot_runtime_per_agent(summary; on_the_fly = true)

Runtime per million agents of each stage against the population size (log-log); linear scaling is
a flat line. Requires the output of `summarize_benchmarks`.
"""
function plot_runtime_per_agent(summary::DataFrame; on_the_fly::Bool = true)
    d = _select_populations(summary, on_the_fly)
    per_agent(s) = d[!, "$(s)_median"] ./ (d.pop_size ./ 1e6)
    ys = reduce(vcat, [per_agent(s) for (s, _, _) in _STAGES])
    p = plot(;
        xscale = :log10, yscale = :log10,
        xticks = _decade_ticks(minimum(d.pop_size), maximum(d.pop_size), _size_label),
        yticks = _decade_ticks(minimum(ys), maximum(ys), _seconds_label),
        xlabel = "Population size", ylabel = "Runtime per million agents (s)",
        title = "Runtime per agent",
        legend = :topright, size = (375, 280),
        _bench_style()...)
    for (s, label, color) in _STAGES
        for (type, filled) in (("Gesyland", true), ("Synthetic", false))
            g = d[d.type .== type, :]
            nrow(g) == 0 && continue
            scatter!(p, g.pop_size, g[!, "$(s)_median"] ./ (g.pop_size ./ 1e6); _BENCH_MARKER...,
                markercolor = filled ? color : :white, markerstrokecolor = color,
                label = (filled || !any(d.type .== "Gesyland")) ? label : "")
        end
    end
    _fill_legend!(p, d)
    return p
end

"""
    replot_benchmarks(folder; on_the_fly = true)

Draws the benchmark figures from the CSV files that `benchmark()` wrote to `folder`, without
running anything: `Fig_19_benchmarks.png` (`save_figure`; with the speedup panel if its summary
exists), `benchmarks_per_agent.png`, and, temporarily, `benchmarks_single_panel.png`.
"""
function replot_benchmarks(folder; on_the_fly::Bool = true)
    summary = CSV.read(joinpath(folder, "benchmark_summary.csv"), DataFrame)
    machine_file = joinpath(folder, "benchmark_machine.csv")
    machine = isfile(machine_file) ? CSV.read(machine_file, DataFrame) : DataFrame()
    cores = hasproperty(machine, :physical_cores) ? machine.physical_cores[1] : missing
    speedup_file = joinpath(folder, "benchmark_speedup_summary.csv")
    speedup = isfile(speedup_file) ? CSV.read(speedup_file, DataFrame) : nothing
    save_figure(plot_benchmarks(summary; speedup = speedup, physical_cores = cores, on_the_fly = on_the_fly),
        :benchmarks, folder)
    png(set_text_size!(plot_runtime_per_agent(summary; on_the_fly = on_the_fly)), joinpath(folder, "benchmarks_per_agent.png"))

    # TEMPORARY: single-panel alternative, see `plot_benchmarks_single_panel`
    p_single = plot_benchmarks_single_panel(summary; on_the_fly = on_the_fly)
    isnothing(speedup) || (p_single = plot(p_single, plot_thread_speedup(speedup; physical_cores = cores, on_the_fly = on_the_fly),
        layout = (1, 2), size = (100FIG_FULL_WIDTH_IN, 280), dpi = FIG_DPI, left_margin = 3mm, bottom_margin = 3mm))
    png(set_text_size!(p_single), joinpath(folder, "benchmarks_single_panel.png"))
end


# RUN ALL ANALYSES

"""
    benchmark(;
        folder = nothing,
        folder_tag = nothing,
        populations = BENCHMARK_POPULATIONS,
        synthetic_sizes = BENCHMARK_SYNTHETIC_SIZES,
        num_runs = 11,
        speedup_populations = ["RP"],
        speedup_threads = [1, 2, 4, 8, 16, 32, 64, 128],
        speedup_runs = 5,
        peak_memory_population = "DE")

Runs all benchmark-related analyses: runtime by population size (with the threads of this
process), speedup by number of threads (one process per thread count; skipped if `speedup_populations` or
`speedup_threads` is empty), and the peak memory of one run of `peak_memory_population` in a fresh
process (skipped if `nothing`). Writes all runs, their summaries, the machine information, and the
figures.

If `folder` is `nothing`, a new result folder is created, with `folder_tag`
appended to its timestamp-based name (see `result_folder`).
"""
function benchmark(;
    folder = nothing,
    folder_tag = nothing,
    populations = BENCHMARK_POPULATIONS,
    synthetic_sizes = BENCHMARK_SYNTHETIC_SIZES,
    num_runs::Int = 11,
    speedup_populations = ["RP"],
    speedup_threads = [1, 2, 4, 8, 16, 32, 64, 128],
    speedup_runs::Int = 5,
    peak_memory_population = "DE")

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder(folder_tag) : nothing

    CSV.write(joinpath(folder, "benchmark_machine.csv"), machine_info())

    # runtime by population size
    runs = run_sims(populations = populations, synthetic_sizes = synthetic_sizes, num_runs = num_runs)
    CSV.write(joinpath(folder, "benchmark_runs.csv"), runs)
    CSV.write(joinpath(folder, "benchmark_summary.csv"), summarize_benchmarks(runs))

    # speedup by number of threads
    if !isempty(speedup_populations) && !isempty(speedup_threads)
        speedup = run_thread_speedup(populations = speedup_populations, threads = speedup_threads,
            num_runs = speedup_runs, folder = mkpath(joinpath(folder, "thread_speedup")))
        CSV.write(joinpath(folder, "benchmark_speedup_runs.csv"), speedup)
        CSV.write(joinpath(folder, "benchmark_speedup_summary.csv"), summarize_benchmarks(speedup))
    end

    # peak memory
    if !isnothing(peak_memory_population)
        CSV.write(joinpath(folder, "benchmark_peak_memory.csv"),
            run_peak_memory(peak_memory_population, folder = mkpath(joinpath(folder, "peak_memory"))))
    end

    replot_benchmarks(folder)
end
