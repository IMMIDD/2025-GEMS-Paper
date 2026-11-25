###
### FIRST WAVE 2020 SCENARIO (COMPARISON WITH RKI CASES)
###

export FIRST_WAVE_PARTS
export FIRST_WAVE_PAPER_PARTS
export scenario_first_wave
export calibrate_first_wave
export replot_first_wave


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
- `:random_seeds`: the same number of seeds per day in counties drawn in proportion to their
  population
- `:mixing`: uniform and same-age instead of POLYMOD contacts in the municipalities
  (`MUNICIPALITY_MIXING`, same number of contacts, no recalibration), with the RKI
  introductions and with random seed counties, compared with the runs of `:final` and
  `:random_seeds` (run if their outputs are not in the folder): national fit, county
  correlations, age groups, fatality, and settings (`mixing_summary`)
- `:figures`: the paper figures and the outputs of the paper text; with
  `diagnostics`, also the figures of all parts whose results exist (counterfactuals:
  with `:final`)

`scenario_first_wave` runs only the `FIRST_WAVE_PAPER_PARTS` by default.
"""
const FIRST_WAVE_PARTS = (
    :final,
    :projections,
    :counterfactuals,
    :band,
    :seeds_robustness,
    :random_seeds,
    :mixing,
    :figures)

"""
    FIRST_WAVE_PAPER_PARTS

Parts of `scenario_first_wave` whose results go into the paper (Sec. 3.1.1,
Figs. 10-11): the final runs, the band, the runs with random seed counties, the
comparison of the municipality mixing, and the figures.
"""
const FIRST_WAVE_PAPER_PARTS = (:final, :band, :random_seeds, :mixing, :figures)


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

# statistics, infections by age group, daily curves, county and state tables of
# `run_rki_phase_sims` (of the `bases`), and the county incidences of the single runs in the
# fit window (`run_incidences`)
function _write_runs(folder, part, out; bases = BASES)
    CSV.write(joinpath(folder, first_wave_file(part, "stats")), out.stats)
    CSV.write(joinpath(folder, first_wave_file(part, "ages")), out.ages)
    for b in bases
        tbl = rki_cumulative_table(rki_average_runs(out[b].tbls), RKI_WINDOWS)
        CSV.write(joinpath(folder, first_wave_file(part, "daily", b)), out[b].daily)
        CSV.write(joinpath(folder, first_wave_file(part, "counties", b)), _ags_ids(tbl))
        CSV.write(joinpath(folder, first_wave_file(part, "states", b)), out[b].states)
        CSV.write(joinpath(folder, first_wave_file(part, "run_incidences", b)),
            run_incidences(out[b].tbls, "fit"))
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


# runs with the calibrated factors, the expected seeded infections, and the infections per day
# of the runs
function _rki_final(folder, s::RKISettings; seeds, bases = BASES)
    f = _read_factors(s)
    out = run_rki_phase_sims(s, fill(f, length(seeds)), seeds,
        map_weeks = RKI_WEEKS,
        stage = "final")
    _write_runs(folder, "final", out, bases = bases)
    CSV.write(joinpath(folder, first_wave_file("seeded_infections")), rki_seeded_infections(s))
    CSV.write(joinpath(folder, first_wave_file("final", "infections")), out.infections)
end


# runs with the same number of seeds per day, in counties drawn in proportion to their population
function _rki_random_seeds(folder, s::RKISettings; seeds, bases = BASES)
    out = run_rki_phase_sims(s, fill(_read_factors(s), length(seeds)), seeds,
        map_weeks = RKI_WEEKS,
        stage = "random seeds",
        random_seeding = true)
    _write_runs(folder, "random_seeds", out, bases = bases)
end


# outputs of a part of `_write_runs` needed by the mixing comparison
_mixing_files(folder, part) = [joinpath(folder, first_wave_file(part, k, b))
    for (k, b) in (("daily", :report), ("counties", :report), ("run_incidences", :report),
        ("stats", nothing), ("ages", nothing))]

# municipality contacts of the reference model (POLYMOD matrix) and of the alternative mixing
# assumptions (`MUNICIPALITY_MIXING`), with the same number of contacts and the calibrated factors
# of the config (no recalibration), each with the RKI introductions and with introductions in
# random counties; the reference runs are those of `:final` and `:random_seeds` (reused if their
# outputs exist in `folder`, run and written as these parts otherwise); writes the runs of each
# alternative (parts `mixing_<mixing>` and `mixing_<mixing>_random_seeds`), the county metrics of
# all runs, the shares by age group (`mixing_ages`), and the comparison (`mixing_summary`, see
# `mixing_summary`)
function _rki_mixing(folder, s::RKISettings; seeds, bases = BASES)
    municipality = TOML.parsefile(s.configfile)["Settings"]["Municipality"]["contact_sampling_method"]
    if municipality["type"] != "AgeBasedContactSampling"
        printinfo("Mixing comparison skipped: the municipalities of $(s.configfile) do not use " *
            "age-based sampling (AgeBasedContactSampling)")
        return nothing
    end
    factors = fill(_read_factors(s), length(seeds))
    file(part, kind, b = nothing) = joinpath(folder, first_wave_file(part, kind, b))

    models = map((:polymod, MUNICIPALITY_MIXING...)) do model
        mixing = model == :polymod ? nothing : model
        parts = map(((:final, false), (:random, true))) do (key, random)
            part = isnothing(mixing) ? (random ? "random_seeds" : "final") :
                "mixing_$mixing" * (random ? "_random_seeds" : "")
            if !isnothing(mixing) || !all(isfile, _mixing_files(folder, part))
                out = run_rki_phase_sims(s, factors, seeds,
                    map_weeks = RKI_WEEKS,
                    stage = replace(part, "_" => " "),
                    random_seeding = random,
                    municipality_mixing = mixing)
                _write_runs(folder, part, out, bases = bases)
            end
            counties = _read_counties(file(part, "counties", :report))
            metrics = rki_county_metrics(counties, RKI_WINDOWS, RKI_MAP_WEEKS)
            CSV.write(file(part, "metrics", :report), metrics)
            key => (
                daily = CSV.read(file(part, "daily", :report), DataFrame),
                counties = counties,
                metrics = metrics,
                runs = CSV.read(file(part, "run_incidences", :report), DataFrame),
                stats = CSV.read(file(part, "stats"), DataFrame),
                ages = CSV.read(file(part, "ages"), DataFrame))
        end
        model => (; parts...)
    end
    runs(m) = Dict(models)[m].final.runs
    reference = first(first(models))

    rki = rki_age_table(s, to = s.startdate + Day(s.days - 1))
    ages = mixing_age_table(rki, [m => r.final.ages for (m, r) in models])
    CSV.write(file("mixing", "ages"), ages)
    summary = mixing_summary(
        models = collect(models),
        ages = ages,
        run_spearman = (
            within = Dict(m => mean_run_spearman(runs(m), runs(m)) for (m, _) in models),
            between = Dict(m => mean_run_spearman(runs(m), runs(reference)) for (m, _) in models)))
    CSV.write(file("mixing", "summary"), summary)
    printinfo("Mixing comparison written to $(file("mixing", "summary"))")
    return summary
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


# runs with the factor samples (the first `runs`, all if `nothing`) and their 5-95% band, national
# and per state (of the `bases`); each run is written to its own folder (both bases), so a restart
# skips the finished runs (runs without the state file, from before it was written, are repeated)
function _rki_band(folder, s::RKISettings; runs = nothing, bases = BASES)
    samples = CSV.read(first_wave_input("band_factors", s.configfile), DataFrame)
    if !isnothing(runs)
        samples = samples[1:runs, :]
    end
    factor_columns = ["f$k" for k in 1:(length(s.dates) + 1)]
    dir = mkpath(joinpath(folder, "sc_first_wave_band_runs"))
    run_file(seed, b) = joinpath(dir, "run_$(seed)_$(_basis_labels(b).name).csv")
    state_file(seed, b) = joinpath(dir, "run_$(seed)_$(_basis_labels(b).name)_states.csv")
    stats_file(seed) = joinpath(dir, "stats_$seed.csv")

    for r in eachrow(samples)
        all(isfile(run_file(r.seed, b)) && isfile(state_file(r.seed, b)) for b in BASES) && continue
        out = run_rki_phase_sims(s, [collect(r[factor_columns])], [r.seed], stage = "band")
        CSV.write(stats_file(r.seed), out.stats)
        for b in BASES
            CSV.write(run_file(r.seed, b), select(out[b].daily, :date, :rki, :gems))
            CSV.write(state_file(r.seed, b), select(out[b].states, :date, :state, :rki, :gems))
        end
    end

    stats = [CSV.read(stats_file(seed), DataFrame) for seed in samples.seed]
    CSV.write(joinpath(folder, first_wave_file("band", "stats")),
        reduce(vcat, stats, cols = :union))
    for b in bases
        dailies = [CSV.read(run_file(seed, b), DataFrame) for seed in samples.seed]
        band = rki_daily_summary(dailies, quantiles = (0.05, 0.95))
        CSV.write(joinpath(folder, first_wave_file("band", "daily", b)), band)
        states = [CSV.read(state_file(seed, b), DataFrame) for seed in samples.seed]
        CSV.write(joinpath(folder, first_wave_file("band", "states", b)),
            rki_state_summary(states, quantiles = (0.05, 0.95)))
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
    infections_file = joinpath(folder, first_wave_file("final", "infections"))
    if isfile(infections_file)
        infections = CSV.read(infections_file, DataFrame)
        p = plot_rki_seeding(daily, seeded, infections,
            phase_dates = s.dates,
            cutoff = RKI_CUTOFF,
            basis = b)
        png(p, png_file("final", "seeding"))
    end

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

"""
    replot_first_wave(folder; phase_dates = PHASE_DATES)

Draws the paper figures of the final runs (reported cases only) that `scenario_first_wave()` stored in
`folder`, without re-running the simulations: `Fig_11_first_wave.png` (see `plot_paper_figure_1`) and
`Fig_12_first_wave_states.png` (see `plot_paper_figure_2`), with the band of the calibration runs if it
exists (else with the range of the final runs).
"""
function replot_first_wave(folder; phase_dates = PHASE_DATES)
    b = :report
    file(part, kind) = joinpath(folder, first_wave_file(part, kind, b))
    isfile(file("final", "daily")) || return
    read_if_exists(f) = isfile(f) ? CSV.read(f, DataFrame) : nothing
    tbl = _read_counties(file("final", "counties"))

    p = plot_paper_figure_1(tbl, CSV.read(file("final", "daily"), DataFrame), RKI_MAP_WEEKS,
        band = read_if_exists(file("band", "daily")),
        phase_dates = phase_dates,
        cutoff = RKI_CUTOFF,
        basis = b)
    save_figure(p, :first_wave, folder)

    stats = read_if_exists(joinpath(folder, first_wave_file("final", "stats")))
    p = plot_paper_figure_2(rki_state_table(tbl, RKI_WEEKS), CSV.read(file("final", "states"), DataFrame),
        band = read_if_exists(file("band", "states")),
        runs = isnothing(stats) ? nothing : nrow(stats),
        cutoff = RKI_CUTOFF,
        last_day = RKI_LAST_DAY,
        basis = b)
    return save_figure(p, :first_wave_states, folder)
end

# county metrics (Spearman correlations of the windows and single weeks, see `rki_county_metrics`)
function _county_metrics(folder, part, tbl, b)
    metrics = rki_county_metrics(tbl, RKI_WINDOWS, RKI_MAP_WEEKS)
    CSV.write(joinpath(folder, first_wave_file(part, "metrics", b)), metrics)
end

# county metrics of the paper (final runs and random seed counties, by reporting date)
function _paper_metrics(folder)
    for part in ("final", "random_seeds")
        file = joinpath(folder, first_wave_file(part, "counties", :report))
        isfile(file) && _county_metrics(folder, part, _read_counties(file), :report)
    end
end

# county metrics, scatter plots of the windows and of single weeks, and maps
function _figures_counties(folder, part, tbl, b)
    _county_metrics(folder, part, tbl, b)
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

function _figures_random_seeds(folder, s::RKISettings, b)
    file(part, kind) = joinpath(folder, first_wave_file(part, kind, b))
    isfile(file("random_seeds", "daily")) || return
    runs = ["GEMS, seeds in random counties (by population)" =>
        CSV.read(file("random_seeds", "daily"), DataFrame)]
    if isfile(file("final", "daily"))
        final = CSV.read(file("final", "daily"), DataFrame)
        push!(runs, "GEMS, main design (tourism model + Heinsberg)" => final)
    end
    png(plot_rki_projection(runs, phase_dates = s.dates, cutoff = RKI_CUTOFF, basis = b),
        joinpath(folder, first_wave_file("random_seeds", "projection", b, ext = "png")))
    tbl = _read_counties(file("random_seeds", "counties"))
    _figures_counties(folder, "random_seeds", tbl, b)
end

# seeded (introduced) infections of the final runs, for the paper text: infections and seeded
# infections (mean of the runs) for the whole run, the fitted period, and each phase (from its
# start date, no transition), and the planned seeds per source
function _seeding_summary(folder, s::RKISettings)
    file = joinpath(folder, first_wave_file("final", "infections"))
    isfile(file) || return
    inf = CSV.read(file, DataFrame)
    starts = vcat(first(inf.date), s.dates)
    ends = vcat(s.dates .- Day(1), last(inf.date))
    periods = vcat(
        ["whole run" => (first(inf.date), last(inf.date)),
         "fitted period" => (first(inf.date), RKI_CUTOFF)],
        ["phase $k" => (starts[k], ends[k]) for k in eachindex(starts)])
    rows = map(periods) do (name, (from, to))
        i = from .<= inf.date .<= to
        (period = name, from = from, to = to, infections = sum(inf.infections[i]),
            seeded = sum(inf.seeded[i]), seeded_share = sum(inf.seeded[i]) / sum(inf.infections[i]))
    end
    CSV.write(joinpath(folder, first_wave_file("final", "seeding_summary")), DataFrame(rows))

    seeds = CSV.read(joinpath(folder, first_wave_file("seeded_infections")), DataFrame)
    sources = setdiff(names(seeds), ["date"])
    totals = [sum(seeds[!, c]) for c in sources]
    CSV.write(joinpath(folder, first_wave_file("final", "seeding_sources")),
        DataFrame(source = sources, seeded = totals, share = totals ./ sum(totals)))
end

# paper figures, seeding summary, and county metrics of the paper; with `diagnostics`, all figures
# (and the county metrics) of the parts whose results exist in `folder`
function _rki_figures(folder, s::RKISettings; diagnostics = false)
    _seeding_summary(folder, s)
    # before all other figures: GR maps each new color to an entry of a color table that lasts for the whole
    # Julia session, and once the county maps of the other figures have filled it, new colors are drawn wrong
    replot_first_wave(folder, phase_dates = s.dates)
    if !diagnostics
        _paper_metrics(folder)
        return
    end
    for b in BASES
        _figures_final(folder, s, b)
        _figures_band(folder, s, b)
        _figures_projections(folder, s, b)
        _figures_counterfactuals(folder, s, b)
        _figures_seeds_robustness(folder, s, b)
        _figures_random_seeds(folder, s, b)
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
        folder_tag = nothing,
        generation_time = true,
        test = false)

Calibrates the contact factors of `scenario_first_wave` on the RKI cases by
symptom onset up to `RKI_CUTOFF` (after measuring the generation time of the
model, if `generation_time`). Writes the results to `folder`, and the factors for
the scenario to `inputs/` (`first_wave_input`) if calibrated on Germany and not
as a `test` (one run per evaluation).

If `folder` is `nothing`, a new result folder is created, with `folder_tag`
appended to its timestamp-based name (see `result_folder`).
"""
function calibrate_first_wave(;
    configfile = "models/DE_model_first_wave.toml",
    population = "DE",
    folder = nothing,
    folder_tag = nothing,
    generation_time = true,
    test = false)

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder(folder_tag) : mkpath(folder)

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
        folder_tag = nothing,
        parts = FIRST_WAVE_PAPER_PARTS,
        diagnostics = false,
        test = false)

First wave 2020 in Germany compared with the RKI cases: runs with the contact
factors calibrated on the cases by symptom onset up to `RKI_CUTOFF` (see
`calibrate_first_wave`) until `RKI_LAST_DAY`, a band of the calibration
uncertainty, runs with random seed counties, and the paper figures (see
`FIRST_WAVE_PAPER_PARTS`). Further parts (phase projections, counterfactuals, a
robustness run of the seeding) can be selected from `FIRST_WAVE_PARTS`.

By default, only the outputs of the paper are written: cases by reporting date,
the paper figures, the seeding summary, and the county metrics. With
`diagnostics`, the comparisons by symptom onset and the figures of all parts are
written as well.

Single `parts` can be run in a folder of an earlier run. For a state population
(e.g. "SL"), all data are restricted to its counties. With `test`, every part
runs only once (one seed, one variant). The RKI case data are downloaded on the
first run (`obtain_rki_cases`).

If `folder` is `nothing`, a new result folder is created, with `folder_tag`
appended to its timestamp-based name (see `result_folder`).
"""
function scenario_first_wave(;
    configfile = "models/DE_model_first_wave.toml",
    population = "DE",
    folder = nothing,
    folder_tag = nothing,
    parts = FIRST_WAVE_PAPER_PARTS,
    diagnostics = false,
    test = false)

    unknown = setdiff(parts, FIRST_WAVE_PARTS)
    isempty(unknown) || throw(ArgumentError("Unknown parts $unknown, see FIRST_WAVE_PARTS."))

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder(folder_tag) : mkpath(folder)

    s = _first_wave_settings(configfile, population)

    K = length(s.dates) + 1
    # the paper compares the cases by reporting date
    bases = diagnostics ? BASES : (:report,)

    if :final in parts
        _rki_final(folder, s, seeds = _seeds(3, test), bases = bases)
    end
    if :projections in parts
        _rki_projections(folder, s, seeds = _seeds(2, test), phases = test ? [K] : 2:K)
    end
    if :counterfactuals in parts
        _rki_counterfactuals(folder, s, seeds = [0], phases = test ? [K - 1] : 1:(K - 1))
    end
    if :band in parts
        _rki_band(folder, s, runs = test ? 1 : nothing, bases = bases)
    end
    if :seeds_robustness in parts
        _rki_seeds_robustness(folder, s, seeds = _seeds(3, test))
    end
    if :random_seeds in parts
        _rki_random_seeds(folder, s, seeds = _seeds(3, test), bases = bases)
    end
    if :mixing in parts
        _rki_mixing(folder, s, seeds = _seeds(3, test), bases = bases)
    end
    if :figures in parts
        _rki_figures(folder, s, diagnostics = diagnostics)
    end
    return folder
end
