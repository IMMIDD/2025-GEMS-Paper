###
### FIRST WAVE 2020: REPRODUCTION NUMBER AND R PER PHASE
###

"""
    BACK_CALCULATION_START

First infection date of the back-calculations of R(t).
"""
const BACK_CALCULATION_START = Date(2020, 2, 1)

"""
    generation_time_file(configfile)

File of the forward generation time of the model of a config
(`gems_generation_time`, columns `s` and `w`), written by the part
`:generation_time` of `scenario_first_wave`.
"""
function generation_time_file(configfile)
    name = splitext(basename(configfile))[1]
    return joinpath("inputs", "generation_time_$name.csv")
end


"""
    gems_generation_time(pp; cohort = 10:22, max_days = 30)

Forward generation time distribution (`w[s]`: share of `s` days) of the
infections in the ticks `cohort`. The run must last at least
`last(cohort) + max_days` ticks.
"""
function gems_generation_time(pp; cohort = 10:22, max_days::Int = 30)
    inf = infections(pp)
    tick_of = Dict(zip(inf.infection_id, inf.tick))
    ids = Set(inf.infection_id[in.(inf.tick, Ref(cohort))])
    # transmissions (id_a == -1: seeded) from the cohort
    tx = inf[(inf.id_a .!= -1) .& in.(inf.source_infection_id, Ref(ids)), :]
    gt = Int.(tx.tick) .- [tick_of[s] for s in tx.source_infection_id]
    return [count(==(s), gt) for s in 1:max_days] ./ count(<=(max_days), gt)
end


"""
    rki_onset_share(rki; counties = nothing, weeks = 10:20, lag = 5, min_cases = 50)

Share of RKI cases with a known symptom onset as a function of the onset date:
the share among the cases reported in the week of onset + `lag` days. Weeks
with fewer than `min_cases` reports get the share of the nearest week.
"""
function rki_onset_share(rki::DataFrame;
    counties = nothing,
    weeks = 10:20,
    lag::Int = 5,
    min_cases::Int = 50)

    d = isnothing(counties) ? rki : rki[in.(rki.ags, Ref(Set(counties))), :]
    counts = map(weeks) do w
        x = d[monday_2020(w) .<= d.report .<= sunday_2020(w), :]
        (sum(x.cases[.!ismissing.(x.onset)]; init = 0), sum(x.cases; init = 0))
    end
    valid = [w for (w, (_, n)) in zip(weeks, counts) if n >= min_cases]
    isempty(valid) && throw(ArgumentError("No week with at least $min_cases reported cases."))
    share = Dict(w => k / n for (w, (k, n)) in zip(weeks, counts) if n >= min_cases)
    nearest(w) = valid[argmin(abs.(valid .- w))]
    return date -> share[nearest(week_2020(date + Day(lag)))]
end


"""
    onset_scaled(x::Dict, onset_share, incubation)

Infections per infection date (`Dict` date => value) in cases with a known
onset: each value times the onset share at the onset date (infection plus
`incubation` days).
"""
function onset_scaled(x::Dict, onset_share, incubation::Int)
    return Dict(d => v * onset_share(d + Day(incubation)) for (d, v) in x)
end


"""
    rki_reproduction(rki;
        generation_time,
        delays,
        model,
        imports = Dict{Date, Float64}(),
        counties = nothing,
        from = BACK_CALCULATION_START,
        to = Date(2020, 5, 10),
        iterations = 200)

Reproduction number by infection date from the cases in `rki` (column `report`,
in `counties` or all of Germany): the infections are back-calculated with the
delays of the model (`report_offsets`), and R(t) is the local infections
(without `imports`) divided by the past infections weighted with the generation
time (renewal equation, column `lambda`).
"""
function rki_reproduction(rki::DataFrame;
    generation_time::Vector{Float64},
    delays,
    model,
    imports = Dict{Date, Float64}(),
    counties = nothing,
    from::Date = BACK_CALCULATION_START,
    to::Date = Date(2020, 5, 10),
    iterations::Int = 200)

    d = rki[from .<= rki.report .<= to, :]
    if !isnothing(counties)
        d = d[in.(d.ags, Ref(Set(counties))), :]
    end
    days = collect(from:Day(1):to)
    y = zeros(length(days))
    for r in eachrow(d)
        y[Dates.value(r.report - from) + 1] += r.cases
    end

    K = zeros(length(days), length(days))
    for (i, t) in enumerate(days), (offset, p) in report_offsets(t, delays, model)
        if 1 <= i + offset <= length(days)
            K[i, i + offset] += p
        end
    end
    x = back_calculate(K, y, iterations)

    w = generation_time
    m = [get(imports, t, 0.0) for t in days]
    lambda = zeros(length(x))
    for i in eachindex(x), s in 1:min(length(w), i - 1)
        lambda[i] += w[s] * x[i - s]
    end
    return DataFrame(
        date = days,
        reported = y,
        infections = x,
        imports = m,
        lambda = lambda,
        R = [l > 0 ? max(xi - mi, 0.0) / l : missing for (xi, mi, l) in zip(x, m, lambda)])
end



###
### R PER PHASE
###

# R per phase by Poisson maximum likelihood of the local infections, loc ~ Poisson((V R) lambda)
# with the phase weights V (days x phases); with steps R_k = sum(loc) / sum(lambda) over the days of
# phase k, otherwise iteratively reweighted least squares (200 iterations at most)
function phase_R(loc, lambda, V)
    X = V .* lambda
    R = [sum(loc .* V[:, k]) / sum(X[:, k]) for k in axes(V, 2)]
    if all(x -> x == 0 || x == 1, V) || !all(isfinite, R)
        return R
    end
    for _ in 1:200
        mu = max.(X * R, 1e-9)
        new = max.((X' * (X ./ mu)) \ (X' * (loc ./ mu)), 1e-6)
        done = maximum(abs.(new ./ R .- 1)) < 1e-10
        R = new
        done && break
    end
    return R
end

# local infections, lambda, and phase weights of the days of an R(t) table within `window`
function phase_data(rt::DataFrame; dates, transitions, window)
    check_transitions(dates, transitions)
    from, to = window
    d = rt[from .<= rt.date .<= to, :]
    V = phase_matrix(d.date, dates, transitions)
    if !all(>(0), sum(V, dims = 1))
        throw(ArgumentError("Every phase needs weight on a day between $from and $to."))
    end
    return (loc = max.(d.infections .- d.imports, 0.0), lambda = Float64.(d.lambda), V = V)
end


"""
    estimate_phases(rt; dates, transitions, window)

R per phase from an R(t) table (`rki_reproduction`) on the days of `window`
(first and last date), by Poisson maximum likelihood of the local infections.
"""
function estimate_phases(rt::DataFrame; dates, transitions, window)
    (; loc, lambda, V) = phase_data(rt, dates = dates, transitions = transitions, window = window)
    return phase_R(loc, lambda, V)
end


"""
    bootstrap_phases(rt;
        dates,
        transitions,
        window,
        n = 1000,
        block = 7,
        seed = 0)

`n` samples (rows) of `estimate_phases` from a moving-block bootstrap of the
daily residuals (blocks of `block` days, as the back-calculation smooths over
several days and the reports have a weekly cycle).
"""
function bootstrap_phases(rt::DataFrame;
    dates,
    transitions,
    window,
    n::Int = 1000,
    block::Int = 7,
    seed = 0)

    (; loc, lambda, V) = phase_data(rt, dates = dates, transitions = transitions, window = window)
    mu = (V * phase_R(loc, lambda, V)) .* lambda
    e = [m > 0 ? l / m : 1.0 for (l, m) in zip(loc, mu)]
    N = length(e)
    rng = Random.Xoshiro(seed)
    samples = zeros(n, size(V, 2))
    for b in 1:n
        idx = Int[]
        while length(idx) < N
            s = rand(rng, 1:N)
            append!(idx, mod1.(s:(s + block - 1), N))
        end
        samples[b, :] = phase_R(mu .* e[idx[1:N]], lambda, V)
    end
    return samples
end


"""
    gems_phase_R(daily;
        generation_time,
        delays,
        model,
        imports,
        dates,
        transitions,
        window)

R per phase of daily cases (columns `date` and `gems`), computed in the same way
as for the RKI data (back-calculation from `BACK_CALCULATION_START`).
"""
function gems_phase_R(daily::DataFrame;
    generation_time,
    delays,
    model,
    imports,
    dates,
    transitions,
    window)

    g = DataFrame(report = daily.date, cases = Float64.(daily.gems))
    rt = rki_reproduction(g,
        generation_time = generation_time,
        delays = delays,
        model = model,
        imports = imports,
        to = last(daily.date))
    return estimate_phases(rt, dates = dates, transitions = transitions, window = window)
end
