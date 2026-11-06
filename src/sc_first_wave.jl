###
### FIRST WAVE 2020 SCENARIO (COMPARISON WITH RKI CASES)
###

export FIRST_WAVE_PARTS
export scenario_first_wave
export calibrate_first_wave


"""
    RKI_CUTOFF

Last onset date of the fitted data (19.04.2020; reopenings from 20.04. on).
"""
const RKI_CUTOFF = Date(2020, 4, 19)

"""
    RKI_LAST_DAY

Last simulated day (six weeks of forecast after `RKI_CUTOFF`).
"""
const RKI_LAST_DAY = RKI_CUTOFF + Week(6)

"""
    RKI_SHORT_LAST_DAY

Last simulated day of the calibration runs and of the projection without measures.
"""
const RKI_SHORT_LAST_DAY = Date(2020, 5, 17)

"""
    RKI_FIT_WEEKS

ISO weeks of the fitted data, from week 10 (02.-08.03.) to the week of `RKI_CUTOFF`.
"""
const RKI_FIT_WEEKS = 10:week_2020(RKI_CUTOFF)

"""
    RKI_FORECAST_WEEKS

ISO weeks of the forecast, after `RKI_CUTOFF` up to `RKI_LAST_DAY`.
"""
const RKI_FORECAST_WEEKS = (week_2020(RKI_CUTOFF) + 1):week_2020(RKI_LAST_DAY)

"""
    RKI_WEEKS

All ISO weeks of the county comparison.
"""
const RKI_WEEKS = first(RKI_FIT_WEEKS):last(RKI_FORECAST_WEEKS)

"""
    RKI_WINDOWS

Windows of the cumulative county comparison (`rki_cumulative_table`).
"""
const RKI_WINDOWS = ["fit" => RKI_FIT_WEEKS, "forecast" => RKI_FORECAST_WEEKS]

"""
    RKI_MAP_WEEKS

ISO weeks shown in the maps and scatter plots of single weeks.
"""
const RKI_MAP_WEEKS = [10, 11, 13, 15, 17]

"""
    RKI_STATE_PLOT_START

First day of the daily plots per state.
"""
const RKI_STATE_PLOT_START = Date(2020, 3, 1)

"""
    RKI_BAND_SIZE

Number of runs of the band (factor samples of the calibration).
"""
const RKI_BAND_SIZE = 20

"""
    BASES

The two ways of comparing cases: by reporting date (`:report`) and by symptom onset (`:onset`).
"""
const BASES = (:report, :onset)

"""
    FIRST_WAVE_PARTS

Parts of `scenario_first_wave`, in the order they are run. All parts use the
calibrated contact factors of the config (`first_wave_input`, see
`calibrate_first_wave`):

- `:final`: runs with the calibrated factors until `RKI_LAST_DAY`
- `:projections`: return to the contacts of earlier phases after `RKI_CUTOFF`
- `:counterfactuals`: contacts of a phase kept from the next phase on
- `:band`: runs with the factor samples of the calibration
- `:seeds_robustness`: all RKI cases infected before the first phase change as seeds
- `:figures`: figures of all parts whose results exist (counterfactuals: with `:final`)
"""
const FIRST_WAVE_PARTS = (
    :final,
    :projections,
    :counterfactuals,
    :band,
    :seeds_robustness,
    :figures)


"""
    first_wave_file(part, kind = "", basis = nothing; ext = "csv")

Name of a result file of `scenario_first_wave`, e.g.
`first_wave_file("final", "daily", :onset)` is `sc_first_wave_final_daily_onset.csv`.
"""
function first_wave_file(part, kind = "", basis = nothing; ext = "csv")
    name = isnothing(basis) ? "" : _basis_labels(basis).name
    return "sc_first_wave_" * join(filter(!isempty, [part, kind, name]), "_") * ".$ext"
end

"""
    first_wave_input(kind, configfile)

File in `inputs/` with the calibrated contact factors (`kind = "factors"`) or
their samples for the band (`kind = "band_factors"`) of a config, written by
`calibrate_first_wave`.
"""
function first_wave_input(kind, configfile)
    name = splitext(basename(configfile))[1]
    return joinpath("inputs", "first_wave_$(kind)_$name.csv")
end

# AGS as integers (for CSV files)
_ags_ids(d::DataFrame) = transform(d, :ags => ByRow(a -> a.id) => :ags)

# calibrated contact factors of the config
function _read_factors(s::RKISettings)
    return CSV.read(first_wave_input("factors", s.configfile), DataFrame).factor
end

# runs (`label => daily`) in one file with the column `scenario`, and back
function _write_labelled(file, runs)
    labelled = [insertcols(copy(d), :scenario => label) for (label, d) in runs]
    CSV.write(file, reduce(vcat, labelled))
end
function _read_labelled(d::DataFrame)
    return [s => select(d[d.scenario .== s, :], Not(:scenario)) for s in unique(d.scenario)]
end

# county table with AGS
_read_counties(file) = transform(CSV.read(file, DataFrame), :ags => ByRow(AGS) => :ags)

# statistics, daily curves, and county and state tables of `run_rki_phase_sims` (both bases)
function _write_runs(folder, part, out)
    CSV.write(joinpath(folder, first_wave_file(part, "stats")), out.stats)
    for b in BASES
        tbl = rki_cumulative_table(rki_average_runs(out[b].tbls), RKI_WINDOWS)
        CSV.write(joinpath(folder, first_wave_file(part, "daily", b)), out[b].daily)
        CSV.write(joinpath(folder, first_wave_file(part, "counties", b)), _ags_ids(tbl))
        CSV.write(joinpath(folder, first_wave_file(part, "states", b)), out[b].states)
    end
end



###
### PARTS
###

# forward generation time of the model: Saarland run with a constant transmission rate and 20 extra
# seeds per day in Regionalverband Saarbruecken (ticks 7-27) for a large early cohort; the rate is
# 0.05 per contact and day, i.e. a scale of 0.35 (0.05 over 7 days) for an infectiousness profile
function _rki_generation_time(s::RKISettings)
    seeds = DataFrame(tick = 7:27, ags = AGS(10041000), count = 20)
    rate = config_transmission_type(s.configfile) == "InfectiousnessProfile" ? 0.35 : 0.05
    saarland = RKISettings(
        rki = s.rki,
        days = 70,
        configfile = s.configfile,
        population = "SL",
        ascertainment = 1.0)
    res = run_rki_sim(saarland,
        seed = 2,
        heinsberg = false,
        extra_seeds = seeds,
        transmission_rate = rate)
    w = gems_generation_time(res.pp)
    printinfo("Forward generation time: mean $(round(sum((1:length(w)) .* w), digits = 2)) days")
    CSV.write(generation_time_file(s.configfile), DataFrame(s = 1:length(w), w = w))
end


# calibration on the RKI cases by onset up to `RKI_CUTOFF` (runs until `RKI_SHORT_LAST_DAY`):
# writes the diagnostics and the factors to `folder` and returns the factors and `band_size`
# samples of them
function _rki_calibration(folder, s::RKISettings; seeds, band_size)
    generation_time = CSV.read(generation_time_file(s.configfile), DataFrame).w
    short = with_settings(s, days = simulation_days(s.startdate, RKI_SHORT_LAST_DAY))
    cal = calibrate_rki_phases(short,
        generation_time = generation_time,
        cutoff = RKI_CUTOFF,
        weeks = RKI_FIT_WEEKS,
        seeds = seeds)

    phases = 1:length(cal.factors)
    factors = DataFrame(phase = phases, factor = cal.factors, converged = cal.converged)

    # the band runs get seeds different from those of the other runs
    S = curve_factor_samples(cal, n = band_size)
    band_factors = hcat(
        DataFrame(seed = 1000 .+ (1:size(S, 1))),
        DataFrame(S, ["f$k" for k in 1:size(S, 2)]))

    file(kind) = joinpath(folder, first_wave_file("calibration", kind))
    CSV.write(file(""), factors)
    CSV.write(file("rt"), cal.rt)
    CSV.write(file("targets"), cal.targets)
    CSV.write(file("surrogate"), cal.surrogate)
    CSV.write(file("history"), cal.history)
    CSV.write(file("runs"), cal.runs)
    CSV.write(file("jacobian"), cal.jacobian)
    CSV.write(joinpath(folder, first_wave_file("band", "factors")), band_factors)
    return (factors = factors, band_factors = band_factors)
end


# runs with the calibrated factors, and the expected seeded infections
function _rki_final(folder, s::RKISettings; seeds)
    f = _read_factors(s)
    out = run_rki_phase_sims(s, fill(f, length(seeds)), seeds,
        map_weeks = RKI_WEEKS,
        stage = "final")
    _write_runs(folder, "final", out)
    CSV.write(joinpath(folder, first_wave_file("seeded_infections")), rki_seeded_infections(s))
end


# returns to the contacts of the `phases` after `RKI_CUTOFF`, and to the contacts before the
# measures (phase 1, until `RKI_SHORT_LAST_DAY` only)
function _rki_projections(folder, s::RKISettings; seeds, phases)
    f = _read_factors(s)
    short = with_settings(s, days = simulation_days(s.startdate, RKI_SHORT_LAST_DAY))
    none = run_rki_projections(short, f, seeds,
        start = RKI_CUTOFF,
        phases = [1])
    kept = run_rki_projections(s, f, seeds,
        start = RKI_CUTOFF,
        phases = phases)
    stats = vcat(none.stats, kept.stats, cols = :union)
    CSV.write(joinpath(folder, first_wave_file("projections", "stats")), stats)
    for b in BASES
        file = joinpath(folder, first_wave_file("projections", "daily", b))
        _write_labelled(file, vcat(none[b], kept[b]))
    end
end


# contacts of phase k kept from the start of phase k + 1 on, for each k in `phases`
function _rki_counterfactuals(folder, s::RKISettings; seeds, phases)
    variants = phase_counterfactuals(_read_factors(s),
        dates = s.dates,
        school_closure = s.school_closure,
        phases = phases)
    outs = map(variants) do (label, v)
        out = run_rki_phase_sims(with_settings(s, school_closure = v.school_closure),
            fill(v.factors, length(seeds)), seeds,
            stage = label)
        label => out
    end
    CSV.write(joinpath(folder, first_wave_file("counterfactuals", "stats")),
        reduce(vcat, [out.stats for (_, out) in outs], cols = :union))
    for b in BASES
        _write_labelled(joinpath(folder, first_wave_file("counterfactuals", "daily", b)),
            [label => out[b].daily for (label, out) in outs])
    end
end


# runs with the factor samples (the first `runs`, all if `nothing`) and their 5-95% band; each run
# is written to its own folder, so a restart skips the finished runs
function _rki_band(folder, s::RKISettings; runs = nothing)
    samples = CSV.read(first_wave_input("band_factors", s.configfile), DataFrame)
    if !isnothing(runs)
        samples = samples[1:runs, :]
    end
    factor_columns = ["f$k" for k in 1:(length(s.dates) + 1)]
    dir = mkpath(joinpath(folder, "sc_first_wave_band_runs"))
    run_file(seed, b) = joinpath(dir, "run_$(seed)_$(_basis_labels(b).name).csv")
    stats_file(seed) = joinpath(dir, "stats_$seed.csv")

    for r in eachrow(samples)
        all(isfile(run_file(r.seed, b)) for b in BASES) && continue
        out = run_rki_phase_sims(s, [collect(r[factor_columns])], [r.seed], stage = "band")
        CSV.write(stats_file(r.seed), out.stats)
        for b in BASES
            CSV.write(run_file(r.seed, b), select(out[b].daily, :date, :rki, :gems))
        end
    end

    stats = [CSV.read(stats_file(seed), DataFrame) for seed in samples.seed]
    CSV.write(joinpath(folder, first_wave_file("band", "stats")),
        reduce(vcat, stats, cols = :union))
    for b in BASES
        dailies = [CSV.read(run_file(seed, b), DataFrame) for seed in samples.seed]
        band = rki_daily_summary(dailies, quantiles = (0.05, 0.95))
        CSV.write(joinpath(folder, first_wave_file("band", "daily", b)), band)
    end
end


# robustness: all RKI cases infected before the first phase change as seeds (divided by the
# ascertainment, Heinsberg 1:1), no transmission before it
function _rki_seeds_robustness(folder, s::RKISettings; seeds)
    switch = first(s.dates)
    model = model_incubation(s.configfile)
    # cases without onset are infected `MEDIAN_REPORT_DELAY` days plus the incubation period before
    # their report
    sched = rki_seeding_schedule(s.rki,
        startdate = s.startdate,
        cutoff = switch + Day(model.incubation + MEDIAN_REPORT_DELAY),
        incubation = model.incubation,
        reported_by = nothing)
    sched = sched[sched.tick .< Dates.value(switch - s.startdate), :]
    if !isnothing(s.counties)
        sched = sched[in.(sched.ags, Ref(Set(s.counties))), :]
    end
    scale = [a == HEINSBERG ? 1.0 : 1 / s.ascertainment for a in sched.ags]
    sched.count = round.(Int, sched.count .* scale)
    CSV.write(joinpath(folder, first_wave_file("seeds_robustness", "seeds")), _ags_ids(sched))

    out = run_rki_phase_sims(s, fill(_read_factors(s), length(seeds)), seeds,
        map_weeks = RKI_WEEKS,
        stage = "seeds until $(ddmm(switch))",
        heinsberg = false,
        imports_from = switch,
        extra_seeds = sched,
        transmission_rate = 0.0,
        switches = [switch => config_transmission_rate(s.configfile)])
    _write_runs(folder, "seeds_robustness", out)
end



###
### FIGURES
###

# national curves, county and state comparison of the final runs
function _figures_final(folder, s::RKISettings, b)
    file(kind) = joinpath(folder, first_wave_file("final", kind, b))
    isfile(file("daily")) || return
    png_file(part, kind) = joinpath(folder, first_wave_file(part, kind, b, ext = "png"))
    daily = CSV.read(file("daily"), DataFrame)
    seeded = CSV.read(joinpath(folder, first_wave_file("seeded_infections")), DataFrame)

    runs = [_basis_labels(b).gems => daily]
    png(plot_rki_projection(runs, phase_dates = s.dates, cutoff = RKI_CUTOFF, basis = b),
        png_file("final", "projection"))
    png(plot_rki_imports(daily, seeded, phase_dates = s.dates, cutoff = RKI_CUTOFF, basis = b),
        png_file("final", "projection_imports"))

    tbl = _read_counties(file("counties"))
    _figures_counties(folder, "final", tbl, b)

    st = rki_state_table(tbl, RKI_WEEKS)
    CSV.write(joinpath(folder, first_wave_file("states", "weekly", b)), st)
    png(plot_rki_state_weekly(st, RKI_WEEKS, cutoff = RKI_CUTOFF, basis = b),
        png_file("states", "weekly"))
    png(plot_rki_state_daily(st, CSV.read(file("states"), DataFrame),
            first_day = RKI_STATE_PLOT_START,
            cutoff = RKI_CUTOFF,
            last_day = RKI_LAST_DAY,
            basis = b),
        png_file("states", "daily"))
    png(plot_rki_state_scatter(st, RKI_WINDOWS, basis = b), png_file("states", "scatter"))
    png(plot_rki_state_maps(tbl, RKI_MAP_WEEKS, basis = b), png_file("states", "maps"))
end

# county metrics, scatter plots of the windows and of single weeks, and maps
function _figures_counties(folder, part, tbl, b)
    metrics = rki_county_metrics(tbl, RKI_WINDOWS, RKI_MAP_WEEKS)
    CSV.write(joinpath(folder, first_wave_file(part, "metrics", b)), metrics)
    png_file(kind) = joinpath(folder, first_wave_file(part, kind, b, ext = "png"))
    png(plot_rki_scatter(tbl, RKI_WINDOWS, basis = b), png_file("cumulative_scatter"))
    png(plot_rki_scatter(tbl, RKI_MAP_WEEKS, basis = b), png_file("scatter"))
    png(plot_rki_maps(tbl, RKI_MAP_WEEKS, basis = b), png_file("maps"))
end

function _figures_band(folder, s::RKISettings, b)
    file = joinpath(folder, first_wave_file("band", "daily", b))
    isfile(file) || return
    runs = ["GEMS, 90% band (calibration uncertainty)" => CSV.read(file, DataFrame)]
    png(plot_rki_projection(runs, phase_dates = s.dates, cutoff = RKI_CUTOFF, basis = b),
        joinpath(folder, first_wave_file("band", "", b, ext = "png")))
end

function _figures_projections(folder, s::RKISettings, b)
    file = joinpath(folder, first_wave_file("projections", "daily", b))
    isfile(file) || return
    runs = _read_labelled(CSV.read(file, DataFrame))
    png(plot_rki_phase_projections(runs, phase_dates = s.dates, cutoff = RKI_CUTOFF, basis = b),
        joinpath(folder, first_wave_file("projections", "", b, ext = "png")))
end

function _figures_counterfactuals(folder, s::RKISettings, b)
    file = joinpath(folder, first_wave_file("counterfactuals", "daily", b))
    final = joinpath(folder, first_wave_file("final", "daily", b))
    (isfile(file) && isfile(final)) || return
    runs = vcat(["fitted" => CSV.read(final, DataFrame)], _read_labelled(CSV.read(file, DataFrame)))
    png(plot_rki_counterfactuals(runs, phase_dates = s.dates, cutoff = RKI_CUTOFF, basis = b),
        joinpath(folder, first_wave_file("counterfactuals", "", b, ext = "png")))
end

function _figures_seeds_robustness(folder, s::RKISettings, b)
    file(part, kind) = joinpath(folder, first_wave_file(part, kind, b))
    isfile(file("seeds_robustness", "daily")) || return
    runs = ["GEMS, all cases infected before $(ddmm(first(s.dates))) as seeds" =>
        CSV.read(file("seeds_robustness", "daily"), DataFrame)]
    if isfile(file("final", "daily"))
        final = CSV.read(file("final", "daily"), DataFrame)
        push!(runs, "GEMS, main design (imports + Heinsberg)" => final)
    end
    png(plot_rki_projection(runs, phase_dates = s.dates, cutoff = RKI_CUTOFF, basis = b),
        joinpath(folder, first_wave_file("seeds_robustness", "projection", b, ext = "png")))
    tbl = _read_counties(file("seeds_robustness", "counties"))
    _figures_counties(folder, "seeds_robustness", tbl, b)
end

# all figures (and the county metrics) of the parts whose results exist in `folder`
function _rki_figures(folder, s::RKISettings)
    for b in BASES
        _figures_final(folder, s, b)
        _figures_band(folder, s, b)
        _figures_projections(folder, s, b)
        _figures_counterfactuals(folder, s, b)
        _figures_seeds_robustness(folder, s, b)
    end
end



###
### SCENARIO
###

# seeds of n runs (test: one run)
_seeds(n, test) = test ? [0] : collect(0:(n - 1))

# settings of the scenario; onsets more than the 99% quantile of the model's incubation period
# after the report are treated as unknown
function _first_wave_settings(configfile, population)
    rki = load_rki_cases(max_onset_lag = model_incubation(configfile).max_lag)
    return RKISettings(
        rki = rki,
        days = simulation_days(config_startdate(configfile), RKI_LAST_DAY),
        configfile = configfile,
        population = population)
end


"""
    calibrate_first_wave(;
        configfile = "models/DE_model_first_wave.toml",
        population = "DE",
        folder = nothing,
        generation_time = true,
        test = false)

Calibrates the contact factors of `scenario_first_wave` on the RKI cases by
symptom onset up to `RKI_CUTOFF` (after measuring the generation time of the
model, if `generation_time`). Writes the results to `folder`, and the factors for
the scenario to `inputs/` (`first_wave_input`) if calibrated on Germany and not
as a `test` (one run per evaluation).
"""
function calibrate_first_wave(;
    configfile = "models/DE_model_first_wave.toml",
    population = "DE",
    folder = nothing,
    generation_time = true,
    test = false)

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder() : mkpath(folder)

    s = _first_wave_settings(configfile, population)
    if generation_time
        _rki_generation_time(s)
    end
    cal = _rki_calibration(folder, s,
        seeds = _seeds(2, test),
        band_size = test ? 1 : RKI_BAND_SIZE)

    # only a full calibration on Germany replaces the factors of the scenario
    if population == "DE" && !test
        CSV.write(first_wave_input("factors", configfile), cal.factors)
        CSV.write(first_wave_input("band_factors", configfile), cal.band_factors)
        printinfo("Contact factors written to $(first_wave_input("factors", configfile))")
    end
    return folder
end


"""
    scenario_first_wave(;
        configfile = "models/DE_model_first_wave.toml",
        population = "DE",
        folder = nothing,
        parts = FIRST_WAVE_PARTS,
        test = false)

First wave 2020 in Germany compared with the RKI cases: runs with the contact
factors calibrated on the cases by symptom onset up to `RKI_CUTOFF` (see
`calibrate_first_wave`) until `RKI_LAST_DAY`, phase projections,
counterfactuals, a band of the calibration uncertainty, a robustness run, and
all figures (see `FIRST_WAVE_PARTS`).

Single `parts` can be run in a folder of an earlier run. For a state population
(e.g. "SL"), all data are restricted to its counties. With `test`, every part
runs only once (one seed, one variant). The RKI case data are downloaded on the
first run (`obtain_rki_cases`).
"""
function scenario_first_wave(;
    configfile = "models/DE_model_first_wave.toml",
    population = "DE",
    folder = nothing,
    parts = FIRST_WAVE_PARTS,
    test = false)

    unknown = setdiff(parts, FIRST_WAVE_PARTS)
    isempty(unknown) || throw(ArgumentError("Unknown parts $unknown, see FIRST_WAVE_PARTS."))

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder() : mkpath(folder)

    s = _first_wave_settings(configfile, population)

    K = length(s.dates) + 1

    if :final in parts
        _rki_final(folder, s, seeds = _seeds(3, test))
    end
    if :projections in parts
        _rki_projections(folder, s, seeds = _seeds(2, test), phases = test ? [K] : 2:K)
    end
    if :counterfactuals in parts
        _rki_counterfactuals(folder, s, seeds = [0], phases = test ? [K - 1] : 1:(K - 1))
    end
    if :band in parts
        _rki_band(folder, s, runs = test ? 1 : nothing)
    end
    if :seeds_robustness in parts
        _rki_seeds_robustness(folder, s, seeds = _seeds(3, test))
    end
    if :figures in parts
        _rki_figures(folder, s)
    end
    return folder
end
