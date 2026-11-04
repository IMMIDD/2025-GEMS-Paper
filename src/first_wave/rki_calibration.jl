###
### FIRST WAVE 2020: CALIBRATION OF THE CONTACT FACTORS
###

export calibrate_rki_phases


"""
    MIN_FACTOR

Lower bound of the contact factors in the calibration and its samples.
"""
const MIN_FACTOR = 0.01


# sums of daily values y on `days` per ISO week of 2020
weekly_sums(days, y, weeks) = [sum(y[monday_2020(w) .<= days .<= sunday_2020(w)]) for w in weeks]

# 2.5% and 97.5% quantiles of each column
function interval95(samples)
    low = [quantile(c, 0.025) for c in eachcol(samples)]
    high = [quantile(c, 0.975) for c in eachcol(samples)]
    return (low, high)
end

# phase settings of `estimate_phases` and related functions
phase_keywords(s::RKISettings) = (dates = s.dates, transitions = s.transitions, window = s.r_window)

# root mean square of the residuals per calibration week
rms(r, problem) = sqrt(sum(abs2, r) / length(problem.obs))


# Heinsberg seeds of `run_rki_sim` (if in the counties of `s`) in observed cases per infection date
function _heinsberg_reported(s::RKISettings, model)
    seeds = Dict{Date, Float64}()
    if !isnothing(s.counties) && !(HEINSBERG in s.counties)
        return seeds
    end
    sched = heinsberg_schedule(s.rki, startdate = s.startdate, incubation = model.incubation)
    for r in eachrow(sched)
        d = s.startdate + Day(r.tick)
        seeds[d] = get(seeds, d, 0.0) + r.count * s.ascertainment
    end
    return seeds
end


"""
    rki_phase_targets(s::RKISettings; generation_time, cutoff, bootstrap = 1000, block = 7)

RKI R per phase (`estimate_phases`) from the cases with a known onset up to
`cutoff`, dated by onset, in the counties of `s`. Returns the R(t) table
(`rt`), the R per phase (`R`), its bootstrap samples (`bootstrap_phases`), and
a table of the targets with their 95% interval (`targets`).
"""
function rki_phase_targets(s::RKISettings;
    generation_time,
    cutoff::Date,
    bootstrap::Int = 1000,
    block::Int = 7)

    model = model_incubation(s.configfile)
    data = s.rki[s.rki.refdate .<= cutoff, :]
    imports = imports_by_infection_date(rki_onset_delays_by_week(data), model,
        allocation = s.allocation,
        counties = s.counties)
    rt = rki_reproduction(onset_basis(data),
        generation_time = generation_time,
        delays = nothing,
        model = model,
        counties = s.counties,
        imports = onset_scaled(imports, s.onset_share, model.incubation),
        to = cutoff)
    R = estimate_phases(rt; phase_keywords(s)...)
    samples = bootstrap_phases(rt; n = bootstrap, block = block, phase_keywords(s)...)
    low, high = interval95(samples)
    targets = DataFrame(
        phase = 1:length(R),
        start = vcat(missing, s.dates),
        transition_days = vcat(missing, s.transitions),
        R = R,
        low = low,
        high = high)
    return (rt = rt, R = R, samples = samples, targets = targets)
end



###
### SURROGATE AND LEAST SQUARES
###

"""
    curve_surrogate(rt;
        generation_time,
        delays,
        model,
        imports,
        dates,
        transitions,
        startdate,
        seeds = Dict{Date, Float64}())

Returns a function `R -> y`: the expected daily cases on the days of `rt` for
the R per phase `R`, from the renewal equation with the imports and `seeds`
from `startdate` on. Used as a cheap model of GEMS in the calibration.
"""
function curve_surrogate(rt::DataFrame;
    generation_time,
    delays,
    model,
    imports,
    dates,
    transitions,
    startdate::Date,
    seeds = Dict{Date, Float64}())

    w = generation_time
    days, n = rt.date, nrow(rt)
    m = [d >= startdate ? get(imports, d, 0.0) : 0.0 for d in days]
    h = [get(seeds, d, 0.0) for d in days]
    V = phase_matrix(days, dates, transitions)
    offsets = [report_offsets(d, delays, model) for d in days]
    return function (R)
        x = zeros(n)
        Rt = V * R
        for i in 1:n
            if days[i] >= startdate
                past = sum(w[j] * x[i - j] for j in 1:min(length(w), i - 1); init = 0.0)
                x[i] = Rt[i] * past + m[i] + h[i]
            end
        end
        y = zeros(n)
        for i in 1:n, (o, p) in offsets[i]
            if 1 <= i + o <= n
                y[i + o] += x[i] * p
            end
        end
        return y
    end
end


"""
    curve_residuals(v; obs, low, high, weights = ones(length(obs)), penalty = 10.0)

Residuals of the calibration for the observables `v` (log weekly cases, then
the log ratio of GEMS and RKI R per phase): the weighted log ratio to the RKI
weekly cases `obs`, and a penalty for R ratios outside [`low`, `high`]. Also
returns the derivatives `d` of the residuals.
"""
function curve_residuals(v::AbstractVector;
    obs,
    low,
    high,
    weights = ones(length(obs)),
    penalty::Real = 10.0)

    W = length(obs)
    ℓ = v[(W + 1):end]
    up, down = ℓ .> log.(high), ℓ .< log.(low)
    r = vcat(sqrt.(weights) .* (v[1:W] .- log.(obs)),
        penalty .* (up .* (ℓ .- log.(high)) .+ down .* (ℓ .- log.(low))))
    return (r = r, d = vcat(sqrt.(weights), penalty .* (up .| down)))
end


"""
    curve_step(v, J; obs, max_step = 0.3, damping = 0.0, kwargs...)

Levenberg-Marquardt step of the parameters for the observables `v` with the
Jacobian `J` (steps at most `max_step` per parameter). `kwargs` are passed to
`curve_residuals`. Returns the step `δ` and the predicted residuals `r`.
"""
function curve_step(v::AbstractVector, J::AbstractMatrix;
    obs,
    max_step::Real = 0.3,
    damping::Real = 0.0,
    kwargs...)

    P = size(J, 2)
    lo, hi = fill(-max_step, P), fill(max_step, P)
    D = LinearAlgebra.Diagonal(vec(sum(abs2, J[1:length(obs), :], dims = 1)))
    δ = zeros(P)
    for _ in 1:50
        (; r, d) = curve_residuals(v .+ J * δ; obs = obs, kwargs...)
        A = d .* J
        M = A' * A + damping * D
        b = A' * A * δ - A' * r # normal equations for the new δ, active penalties fixed
        # solve with the parameters outside [lo, hi] fixed at the bound
        fixed = falses(P)
        new = zeros(P)
        for _ in 1:P
            free = .!fixed
            new[free] = M[free, free] \ (b[free] - M[free, fixed] * new[fixed])
            out = free .& ((new .< lo) .| (new .> hi))
            any(out) || break
            new[out] = clamp.(new[out], lo[out], hi[out])
            fixed .|= out
        end
        done = maximum(abs.(new .- δ)) < 1e-10
        δ = new
        done && break
    end
    return (δ = δ, r = curve_residuals(v .+ J * δ; obs = obs, kwargs...).r)
end


"""
    fit_curve_surrogate(observables, p; obs, iterations = 100, h = 1e-3, max_step = 0.5, kwargs...)

Levenberg-Marquardt fit of the log parameters `p` of a cheap model
(`observables(p)` returns the observables, `h` is the step of the numerical
derivatives). Returns `p`, the observables `v`, and their Jacobian `J`.
"""
function fit_curve_surrogate(observables::Function, p::AbstractVector;
    obs,
    iterations::Int = 100,
    h::Real = 1e-3,
    max_step::Real = 0.5,
    kwargs...)

    function jacobian(p, v)
        columns = [(observables(p .+ h .* (eachindex(p) .== j)) .- v) ./ h for j in eachindex(p)]
        return reduce(hcat, columns)
    end
    loss(v) = sum(abs2, curve_residuals(v; obs = obs, kwargs...).r)
    p = Float64.(p)
    # damping: start at 0.01, divided by 3 after a successful step, multiplied by 5 otherwise
    v, damping = observables(p), 1e-2
    for _ in 1:iterations
        J = jacobian(p, v)
        δ = curve_step(v, J; obs = obs, damping = damping, max_step = max_step, kwargs...).δ
        v2 = observables(p .+ δ)
        if loss(v2) < loss(v)
            p, v, damping = p .+ δ, v2, damping / 3
            maximum(abs, δ) < 1e-5 && break
        else
            damping *= 5
            damping > 1e8 && break
        end
    end
    return (p = p, v = v, J = jacobian(p, v))
end



###
### CALIBRATION
###

# fit of the true R per phase of `curve_surrogate` (two passes, the second restarts the damping);
# returns the true R, the R of the surrogate fit as computed from the data (`R_pipeline`), its
# observables and their Jacobian, and a table of the fit
function _fit_surrogate(s::RKISettings, rt, problem, R_rki;
    generation_time,
    imports,
    model,
    seeds,
    weeks)

    K, W = length(R_rki), length(weeks)
    sur = curve_surrogate(rt,
        generation_time = generation_time,
        delays = nothing,
        model = model,
        imports = imports,
        dates = s.dates,
        transitions = s.transitions,
        startdate = s.startdate,
        seeds = seeds)
    function observe(y)
        R = gems_phase_R(DataFrame(date = rt.date, gems = y);
            generation_time = generation_time,
            delays = nothing,
            model = model,
            imports = imports,
            phase_keywords(s)...)
        return vcat(log.(weekly_sums(rt.date, y, weeks)), log.(R) .- log.(R_rki))
    end
    observe_p(p) = observe(sur(exp.(p)))
    # start: R falling from 3 to 0.75 over the phases
    start = collect(range(log(3.0), log(0.75), length = K))
    fit = fit_curve_surrogate(observe_p, start; problem...)
    fit = fit_curve_surrogate(observe_p, fit.p; problem...)
    R_true, R_pipeline = exp.(fit.p), exp.(fit.v[(W + 1):end]) .* R_rki
    ratios = exp.(fit.v[1:W]) ./ problem.obs
    printinfo("Surrogate: true R $(round.(R_true, digits = 3)), " *
        "R $(round.(R_pipeline, digits = 3)) (RKI $(round.(R_rki, digits = 3))), " *
        "weekly ratios $(round.(ratios, digits = 2))")
    table = DataFrame(
        phase = 1:K,
        R_true = R_true,
        R_pipeline = R_pipeline,
        R_rki = R_rki,
        low = problem.low .* R_rki,
        high = problem.high .* R_rki)
    return (R_true = R_true, R_pipeline = R_pipeline, v = fit.v, J = fit.J, table = table)
end


# GEMS runs with the factors f (mean over `seeds`): observables (log weekly cases by onset, log
# ratio of the GEMS and RKI R per phase), residuals, R, weekly ratios, shares of transmissions
# outside households, and run statistics
function _evaluate_factors(s::RKISettings, f, it;
    seeds,
    generation_time,
    imports,
    problem,
    R_rki,
    weeks)

    out = run_rki_phase_sims(s, fill(f, length(seeds)), seeds, stage = "calibration $it")
    daily = out.onset.daily
    R = gems_phase_R(daily;
        generation_time = generation_time,
        delays = nothing,
        model = out.model,
        imports = imports,
        phase_keywords(s)...)
    weekly_rki = weekly_sums(daily.date, daily.rki, weeks)
    weekly_gems = weekly_sums(daily.date, daily.gems, weeks)
    if !isapprox(weekly_rki, problem.obs, rtol = 1e-6)
        @warn "RKI weekly cases of the runs differ from the calibration targets."
    end
    vg = vcat(log.(max.(weekly_gems, 1.0)), log.(R))
    v = vg .- vcat(zeros(length(weeks)), log.(R_rki))
    return (
        θ = log.(f),
        v = v,
        vg = vg,
        r = curve_residuals(v; problem...).r,
        R = R,
        ratio = weekly_gems ./ weekly_rki,
        shares = out.shares,
        stats = out.stats)
end


# one Newton step of the factors f from the run `e0` towards the R of the surrogate fit `sur`; d log
# R_k / d log f_k is the share of transmissions outside households in phase k (households are not
# scaled); the factors change at most by `max_change`
function _warm_start(f, e0, sur, R_rki; max_change)
    K = length(f)
    W = length(sur.v) - K
    B = sur.J[W .+ (1:K), :] .* sur.R_pipeline ./ sur.R_true'
    Jf = B .* (e0.shares .* e0.R ./ f)'
    step = f .+ Jf \ (exp.(sur.v[(W + 1):end]) .* R_rki .- e0.R)
    return max.(clamp.(step, f ./ max_change, f .* max_change), MIN_FACTOR)
end


# row of the calibration history for the run `e`
function _history_row(it, kind, e;
    R_rki,
    problem,
    weeks,
    accepted = missing,
    damping = NaN,
    rms_next = NaN)

    K, W = length(R_rki), length(weeks)
    inside = [problem.low[k] <= exp(e.v[W + k]) <= problem.high[k] for k in 1:K]
    return merge(
        (iteration = it,
            kind = kind,
            accepted = accepted,
            rms = rms(e.r, problem),
            rms_predicted = rms_next,
            damping = damping),
        (; (Symbol("f$k") => x for (k, x) in enumerate(exp.(e.θ)))...),
        (; (Symbol("R$k") => x for (k, x) in enumerate(e.R))...),
        (; (Symbol("R_rki$k") => x for (k, x) in enumerate(R_rki))...),
        (; (Symbol("inside_$k") => x for (k, x) in enumerate(inside))...),
        (; (Symbol("ratio_w$wk") => x for (wk, x) in zip(weeks, e.ratio))...))
end


# damped Gauss-Newton steps of the log factors θ from the best run so far, with Broyden updates of
# the Jacobian (starts from the surrogate Jacobian `J_sur` x the shares of transmissions outside
# households of the first run); stops when the predicted improvement of the RMS is below `tol` or
# after `iterations` runs; the damping starts at 0.01 and is divided by 3 after an accepted run,
# multiplied by 4 otherwise
function _gauss_newton(θ, J_sur, evaluate, record; iterations, tol, max_step, problem)
    best, J, damping, converged = nothing, nothing, 0.01, false
    for it in 1:iterations
        e = evaluate(exp.(θ), it)
        accepted = missing
        if isnothing(best)
            J, best = J_sur .* e.shares', e
        else
            Δθ = e.θ .- best.θ
            J .+= (e.vg .- best.vg .- J * Δθ) * Δθ' ./ sum(abs2, Δθ)
            accepted = sum(abs2, e.r) < sum(abs2, best.r)
            if accepted
                best = e
                damping /= 3
            else
                damping *= 4
            end
        end
        step = curve_step(best.v, J; max_step = max_step, damping = damping, problem...)
        kind = it == 1 ? "warm start" : "Gauss-Newton"
        record(it, kind, e; accepted = accepted, damping = damping, rms_next = rms(step.r, problem))
        if rms(best.r, problem) - rms(step.r, problem) < tol
            converged = true
            break
        end
        if it < iterations
            θ = best.θ .+ step.δ
        end
    end
    return (best = best, J = J, converged = converged)
end


"""
    calibrate_rki_phases(s::RKISettings;
        generation_time,
        cutoff,
        weeks,
        factors = nothing,
        seeds = [0],
        weights = [k == 1 ? 0.5 : 1.0 for k in eachindex(weeks)],
        penalty = 10.0,
        iterations = 6,
        tol = 0.015,
        max_step = 0.3,
        max_change = 3.0,
        bootstrap = 1000,
        block = 7)

Calibrates the contact factors of the phases on the RKI cases by symptom onset
in the counties of `s`: the weekly cases in `weeks` (ending by `cutoff`), with
the R per phase of GEMS kept inside the 95% interval of the RKI R per phase
(`rki_phase_targets`, data up to `cutoff`).

The factors start from `factors` or, if `nothing`, from a fit of
`curve_surrogate`, are moved by one Newton step towards the surrogate fit, and
are then updated by damped Gauss-Newton steps, with `iterations + 1` evaluations
at most, each with one GEMS run per seed of `seeds`. The calibration stops when
the predicted improvement of the residuals is below `tol`.

Returns the `factors`, whether the calibration `converged`, the R(t) table
(`rt`), the `targets`, the surrogate fit (`surrogate`), the `history` of the
evaluations, the statistics of their runs (`runs`), the `jacobian` table, and
the quantities needed by `curve_factor_samples`.
"""
function calibrate_rki_phases(s::RKISettings;
    generation_time,
    cutoff::Date,
    weeks,
    factors = nothing,
    seeds = [0],
    weights = [k == 1 ? 0.5 : 1.0 for k in eachindex(weeks)],
    penalty::Real = 10.0,
    iterations::Int = 6,
    tol::Real = 0.015,
    max_step::Real = 0.3,
    max_change::Real = 3.0,
    bootstrap::Int = 1000,
    block::Int = 7)

    if !all(w -> sunday_2020(w) <= cutoff, weeks)
        throw(ArgumentError("All weeks must end by $cutoff."))
    end
    model = model_incubation(s.configfile)
    K = length(s.dates) + 1

    # 1. RKI targets: weekly cases and the 95% interval of the R per phase, relative to the RKI R
    data = rki_phase_targets(s,
        generation_time = generation_time,
        cutoff = cutoff,
        bootstrap = bootstrap,
        block = block)
    R_rki = data.R
    low, high = interval95(data.samples)
    problem = (
        obs = Float64.(weekly_sums(data.rt.date, data.rt.reported, weeks)),
        weights = Float64.(weights),
        penalty = penalty,
        low = low ./ R_rki,
        high = high ./ R_rki)

    # imports for the R of GEMS: delays of all data, as for the simulated cases
    imports = imports_by_infection_date(rki_onset_delays_by_week(s.rki), model,
        allocation = s.allocation,
        counties = s.counties)
    imports = onset_scaled(imports, s.onset_share, model.incubation)

    # 2. surrogate fit (by onset: no reporting delay, Heinsberg seeds in cases with a known onset)
    heinsberg = onset_scaled(_heinsberg_reported(s, model), s.onset_share, model.incubation)
    sur = _fit_surrogate(s, data.rt, problem, R_rki,
        generation_time = generation_time,
        imports = imports,
        model = model,
        seeds = heinsberg,
        weeks = weeks)

    history, runs = DataFrame(), DataFrame()
    function evaluate(f, it)
        return _evaluate_factors(s, f, it,
            seeds = seeds,
            generation_time = generation_time,
            imports = imports,
            problem = problem,
            R_rki = R_rki,
            weeks = weeks)
    end
    function record(it, kind, e; kwargs...)
        row = _history_row(it, kind, e; R_rki = R_rki, problem = problem, weeks = weeks, kwargs...)
        push!(history, row, cols = :union)
        append!(runs, insertcols(e.stats, :iteration => it), cols = :union)
        printinfo("Calibration $it ($kind): RMS $(round(rms(e.r, problem), digits = 4)), " *
            "R $(round.(e.R, digits = 3)) (RKI $(round.(R_rki, digits = 3))), " *
            "weekly ratios $(round.(e.ratio, digits = 2))")
    end

    # 3. start run and one Newton step towards the surrogate fit
    f = Float64.(something(factors, sur.R_pipeline ./ sur.R_pipeline[1]))
    e0 = evaluate(f, 0)
    record(0, "start", e0)
    f = _warm_start(f, e0, sur, R_rki, max_change = max_change)
    printinfo("Warm start: factors $(round.(f, digits = 4))")

    # 4. damped Gauss-Newton steps
    (; best, J, converged) = _gauss_newton(log.(f), sur.J, evaluate, record,
        iterations = iterations,
        tol = tol,
        max_step = max_step,
        problem = problem)
    printinfo("Calibration: factors $(round.(exp.(best.θ), digits = 4)), " *
        "RMS $(round(rms(best.r, problem), digits = 4)), " *
        "R $(round.(best.R, digits = 3)) (RKI $(round.(R_rki, digits = 3))), converged $converged")

    jacobian = hcat(
        DataFrame(observable = vcat(["log_week$x" for x in weeks], ["log_R$k" for k in 1:K])),
        DataFrame(J, ["dlog_f$k" for k in 1:K]))
    return (
        factors = exp.(best.θ),
        converged = converged,
        rt = data.rt,
        targets = data.targets,
        surrogate = sur.table,
        history = history,
        runs = runs,
        jacobian = jacobian,
        v = best.v,
        r = best.r,
        J = J,
        problem = problem)
end


"""
    curve_factor_samples(cal; n, seed = 0)

`n` samples of the contact factors (rows) from the approximate covariance of
the calibration `cal` (`calibrate_rki_phases`).
"""
function curve_factor_samples(cal; n::Int, seed = 0)
    W = length(cal.problem.obs)
    A = curve_residuals(cal.v; cal.problem...).d .* cal.J
    σ = sqrt(sum(abs2, cal.r[1:W]) / W)
    C = σ^2 .* inv(A' * A)
    L = LinearAlgebra.cholesky(LinearAlgebra.Symmetric((C .+ C') ./ 2)).L
    θ = log.(cal.factors)
    rng = Random.Xoshiro(seed)
    S = reduce(vcat, [permutedims(exp.(θ .+ L * randn(rng, length(θ)))) for _ in 1:n])
    return max.(S, MIN_FACTOR)
end
