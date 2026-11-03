###
### FIRST WAVE 2020: IMPORT MODEL
###

# Only national import figures are used (there is no public county-level data); imports are
# distributed across counties by tourism and travel data (`tourism_weights`).

export ImportParams


"""
    RKI_IMPORTS

Cases with infection country abroad per ISO reporting week 2020.

- `abroad`: weeks 10-17 from the RKI situation report of 27.04.2020 (Fig. 3).
  Week 9 is extrapolated from the reports of 04.-08.03. (64, 130, 223 cases).
- `alpine_share`: share of Austria and Italy, from the country tables of the
  reports up to 08.03. (0.90), 08.-22.03. (0.90), and 22.-26.03. (0.84).
- `italy_share`: share of Italy in the Alpine imports, from the country tables
  of the reports of 08.03. (all Alpine imports up to week 10), 22.03. (995 of
  4,849 in weeks 11-12), and 27.04. (474 of 6,443 in weeks 13-17).

The RKI counts mentions; a case can name several countries.
"""
const RKI_IMPORTS = DataFrame(
    week = 9:17,
    abroad = [35, 262, 2781, 6454, 3814, 1095, 267, 113, 31],
    alpine_share = [0.90, 0.90, 0.90, 0.90, 0.84, 0.84, 0.84, 0.84, 0.84],
    italy_share = [1.0, 1.0, 0.21, 0.21, 0.07, 0.07, 0.07, 0.07, 0.07])


"""
    ISCHGL_LONLAT

Coordinates of Ischgl (reference point for ski imports from Austria).
"""
const ISCHGL_LONLAT = (10.2919, 47.0128)

"""
    BOZEN_LONLAT

Coordinates of Bozen (reference point for ski imports from Italy, mostly South Tyrol).
"""
const BOZEN_LONLAT = (11.354, 46.498)

"""
    COUNTY_INCOME_FILE

Disposable income of private households per inhabitant 2019 (EUR) per county:
VGR der Länder, Reihe 2 Band 3, Tab. 2.4 (Berechnungsstand Februar 2025,
statistikportal.de), boundaries as of 2023. License: Datenlizenz Deutschland -
Namensnennung - Version 2.0 (see the header of the file).
"""
const COUNTY_INCOME_FILE = joinpath("inputs", "county_income_2019.csv")


# area-weighted centroid of a shapefile polygon (all parts)
function _centroid(poly)
    pts = poly.points
    bounds = vcat(poly.parts .+ 1, length(pts) + 1)
    a = cx = cy = 0.0
    for k in 1:(length(bounds) - 1)
        for i in bounds[k]:(bounds[k + 1] - 2)
            p, q = pts[i], pts[i + 1]
            c = p.x * q.y - q.x * p.y
            a += c
            cx += (p.x + q.x) * c
            cy += (p.y + q.y) * c
        end
    end
    return (cx / (3a), cy / (3a))
end


# income per county (AGS id); city-states have short keys in the file; Eisenach
# (16056, merged into the Wartburgkreis) gets the value of the Wartburgkreis in
# case older shapes still contain it
function _county_income()
    d = CSV.read(COUNTY_INCOME_FILE, DataFrame, comment = "#")
    key(a) = (a < 100 ? a * 1000 : a) * 1000
    income = Dict(key(a) => Float64(i) for (a, i) in zip(d.ags, d.income_2019))
    income[16056000] = income[16063000]
    return income
end


"""
    county_geography()

Population (`EWZ` of the shapefile), state, distance to Ischgl and Bozen (km,
from the centroid), and disposable income per inhabitant (`COUNTY_INCOME_FILE`)
for all German counties. Does not need a population model.
"""
function county_geography()
    # GF 4: land areas (GF 2: water areas)
    shp = DataFrame(germanshapes(2)) |>
        d -> d[coalesce.(d.GF .== 4, false), :]

    to_gk3 = Proj.Transformation("EPSG:4326", "EPSG:31467", always_xy = true)
    ix, iy = to_gk3(ISCHGL_LONLAT)
    bx, by = to_gk3(BOZEN_LONLAT)
    centroids = _centroid.(shp.geometry)
    ags = AGS.(String.(shp.AGS_0))
    income = _county_income()

    return DataFrame(
        ags = ags,
        state = state_number.(ags),
        population = Int.(shp.EWZ),
        distance_km = [hypot(x - ix, y - iy) / 1000 for (x, y) in centroids],
        distance_bozen_km = [hypot(x - bx, y - by) / 1000 for (x, y) in centroids],
        income = [income[a.id] for a in ags]) |>
        sort_by_ags
end


"""
    report_offsets(infection::Date, delays, model; max_j = 30)

Probabilities of the delay (in days) from an infection on `infection` to its
reporting, as in `gems_reported`: onset after 1 + Poisson(`lambda`) days
(asymptomatic: `model.incubation`), plus the reporting delay of the onset week.
"""
function report_offsets(infection::Date, delays, model; max_j::Int = 30)
    probs = Dict{Int, Float64}()
    function add_onset!(k, w)
        dist = delay_distribution(delays, infection + Day(k))
        total = sum(dist.cases)
        for r in eachrow(dist)
            probs[k + r.delay] = get(probs, k + r.delay, 0.0) + w * r.cases / total
        end
    end
    pj = exp(-model.lambda) # Poisson probability of j
    for j in 0:max_j
        add_onset!(1 + j, (1 - model.asymptomatic) * pj)
        pj *= model.lambda / (j + 1)
    end
    add_onset!(model.incubation, model.asymptomatic)
    return probs
end


"""
    back_calculate(K, y, iterations)

Infections per day from the observations `y` by EM back-calculation with a
smoothing step (Becker et al. 1991); `K[i, j]` is the probability that an
infection on day i is observed in period j.
"""
function back_calculate(K::Matrix{Float64}, y::Vector{Float64}, iterations::Int)
    n = size(K, 1)
    s = vec(sum(K, dims = 2))
    x = fill(sum(y) / n, n)
    for _ in 1:iterations
        mu = K' * x
        new = zeros(n)
        for i in 1:n
            s[i] > 0 || continue
            t = 0.0
            for j in eachindex(y)
                if mu[j] > 0
                    t += K[i, j] * y[j] / mu[j]
                end
            end
            new[i] = x[i] / s[i] * t
        end
        x = copy(new)
        for i in 2:(n - 1)
            x[i] = (new[i - 1] + 2new[i] + new[i + 1]) / 4
        end
    end
    return x
end


"""
    import_curve(delays, model; startdate, enddate = Date(2020, 4, 26), iterations = 200)

Expected imported infections per day between `startdate` and `enddate`,
estimated such that, after the model's time to onset and the reporting delays
(`rki_onset_delays_by_week`, see `report_offsets`), they reproduce the weekly
reported imports (`RKI_IMPORTS`). Columns `austria`, `italy`, and `other`.
"""
function import_curve(delays, model;
    startdate::Date,
    enddate::Date = Date(2020, 4, 26),
    iterations::Int = 200)

    imports = RKI_IMPORTS
    days = collect(startdate:Day(1):enddate)
    col = Dict(w => j for (j, w) in enumerate(imports.week))
    K = zeros(length(days), nrow(imports))
    for (i, d) in enumerate(days), (offset, p) in report_offsets(d, delays, model)
        j = get(col, week_2020(d + Day(offset)), 0)
        if j > 0
            K[i, j] += p
        end
    end

    alpine = Float64.(imports.abroad .* imports.alpine_share)
    other = Float64.(imports.abroad .* (1 .- imports.alpine_share))
    return DataFrame(
        date = days,
        austria = back_calculate(K, alpine .* (1 .- imports.italy_share), iterations),
        italy = back_calculate(K, alpine .* imports.italy_share, iterations),
        other = back_calculate(K, other, iterations))
end


"""
    TIROL_NIGHTS

Overnight stays of guests from Germany in Tyrol, winter season 2018/19, by
German origin region (Amt der Tiroler Landesregierung, "Der Tourismus im Winter
18/19", Tab. 6). Regions as defined by Statistik Austria (`TIROL_REGION`).
"""
const TIROL_NIGHTS = Dict(:BY => 2_876_436, :NW => 2_655_420, :MD => 2_442_360, :BW => 1_994_070,
    :Ost => 1_878_668, :Nord => 1_484_185, :BE => 535_800)

"""
    TIROL_REGION

German origin region (`TIROL_NIGHTS`) of each state (state number of the AGS):
MD = Hessen, Rhineland-Palatinate, Saarland; Nord = Schleswig-Holstein, Hamburg,
Lower Saxony, Bremen; Ost = Brandenburg, Mecklenburg-Western Pomerania, Saxony,
Saxony-Anhalt, Thuringia.
"""
const TIROL_REGION = Dict(1 => :Nord, 2 => :Nord, 3 => :Nord, 4 => :Nord, 5 => :NW, 6 => :MD,
    7 => :MD, 8 => :BW, 9 => :BY, 10 => :MD, 11 => :BE, 12 => :Ost, 13 => :Ost, 14 => :Ost,
    15 => :Ost, 16 => :Ost)

"""
    SCHOOL_HOLIDAYS_2020

School holidays from mid-February to March 2020 by state number (incl. adjacent
weekends). Holidays that end before 15.02. are not listed. In Baden-Württemberg,
each school sets its carnival holidays (share: `ImportParams.bw_coverage`).
"""
const SCHOOL_HOLIDAYS_2020 = Dict(
    9 => (Date(2020, 2, 22), Date(2020, 3, 1)),     # Bavaria
    8 => (Date(2020, 2, 22), Date(2020, 3, 1)),     # Baden-Württemberg
    2 => (Date(2020, 2, 29), Date(2020, 3, 15)),    # Hamburg
    7 => (Date(2020, 2, 15), Date(2020, 2, 23)),    # Rhineland-Palatinate
    10 => (Date(2020, 2, 15), Date(2020, 2, 25)),   # Saarland
    14 => (Date(2020, 2, 8), Date(2020, 2, 23)),    # Saxony
    13 => (Date(2020, 2, 8), Date(2020, 2, 23)))    # Mecklenburg-Western Pomerania


"""
    ImportParams(;
        distance_decay = 0.4,
        income_elasticity = 1.0,
        holiday_factor = 1.0,
        bw_coverage = 0.5,
        italy_ski_share = 0.7)

Parameters of the tourism-based import allocation (`tourism_weights`):

- `distance_decay`: distance decay within a Tyrol origin region (slope of the
  Tyrol overnight stays across the regions)
- `income_elasticity`: income elasticity of ski trips (T-MONA guest survey)
- `holiday_factor`: additional ski trips during school holidays (assumption)
- `bw_coverage`: share of schools with carnival holidays in Baden-Württemberg (assumption)
- `italy_ski_share`: share of ski trips in the imports from Italy (RKI situation reports)
"""
Base.@kwdef struct ImportParams
    distance_decay::Float64 = 0.4
    income_elasticity::Float64 = 1.0
    holiday_factor::Float64 = 1.0
    bw_coverage::Float64 = 0.5
    italy_ski_share::Float64 = 0.7
end


# share of schools on holiday in a state on a date
function _holiday(state::Int, date::Date, bw_coverage::Real)
    haskey(SCHOOL_HOLIDAYS_2020, state) || return 0.0
    from, to = SCHOOL_HOLIDAYS_2020[state]
    from <= date <= to || return 0.0
    return state == 8 ? bw_coverage : 1.0
end

# income per county relative to the national mean
function _relative_income(geo::DataFrame)
    mean_income = sum(geo.income .* geo.population) / sum(geo.population)
    return geo.income ./ mean_income
end

# ski weights: Tyrol overnight stays per origin region, distributed within the
# region by population x income^elasticity x distance^(-decay)
function _ski_weights(geo::DataFrame, distance::Vector{Float64}, p::ImportParams)
    income = _relative_income(geo) .^ p.income_elasticity
    base = geo.population .* income .* (distance ./ 100) .^ (-p.distance_decay)
    region = [TIROL_REGION[s] for s in geo.state]
    total = Dict(r => sum(base[region .== r]) for r in keys(TIROL_NIGHTS))
    return [TIROL_NIGHTS[r] * b / total[r] for (r, b) in zip(region, base)]
end


"""
    tourism_weights(geo, date; params = ImportParams())

Probability that an import infected on `date` lands in a county, separately for
imports from Austria, Italy, and other countries (columns `austria`, `italy`,
`other`), using no case data:

- Austria: ski weights x (1 + holiday factor x school holidays of the state);
  ski weights = Tyrol overnight stays of the origin region (`TIROL_NIGHTS`),
  distributed within the region by population, income, and distance to Ischgl
- Italy: a share like Austria with Bozen as reference point, the rest by
  population and income
- other countries: by population and income
"""
function tourism_weights(geo::DataFrame, date::Date; params::ImportParams = ImportParams())
    p = params
    holiday = 1 .+ p.holiday_factor .* [_holiday(s, date, p.bw_coverage) for s in geo.state]
    normalize(x) = x ./ sum(x)
    general = normalize(geo.population .* _relative_income(geo) .^ p.income_elasticity)
    austria = normalize(_ski_weights(geo, geo.distance_km, p) .* holiday)
    bozen = normalize(_ski_weights(geo, geo.distance_bozen_km, p) .* holiday)
    return DataFrame(
        ags = geo.ags,
        austria = austria,
        italy = p.italy_ski_share .* bozen .+ (1 - p.italy_ski_share) .* general,
        other = general)
end

# `tourism_weights` as a function of the date, computed once per date
function tourism_weights_by_date(geo::DataFrame, params::ImportParams)
    cache = Dict{Date, DataFrame}()
    return d -> get!(() -> tourism_weights(geo, d, params = params), cache, d)
end


"""
    expected_imports(curve, weights; counties = nothing)

Expected imported infections per day in `counties` (all if `nothing`) by source
(columns `date`, `austria`, `italy`, `other`), from the import curve and the
county weights per date (e.g. `tourism_weights_by_date`).
"""
function expected_imports(curve::DataFrame, weights::Function; counties = nothing)
    keep = isnothing(counties) ? nothing : Set(counties)
    shares = map(eachrow(curve)) do r
        w = weights(r.date)
        k = isnothing(keep) ? trues(nrow(w)) : in.(w.ags, Ref(keep))
        (sum(w.austria[k]) * r.austria, sum(w.italy[k]) * r.italy, sum(w.other[k]) * r.other)
    end
    return DataFrame(
        date = curve.date,
        austria = getindex.(shares, 1),
        italy = getindex.(shares, 2),
        other = getindex.(shares, 3))
end


"""
    imports_by_infection_date(delays, model;
        allocation = ImportParams(),
        counties = nothing,
        startdate = BACK_CALCULATION_START)

Reported imports per infection date in `counties` (all if `nothing`), as
`Dict` date => imports for `rki_reproduction`.
"""
function imports_by_infection_date(delays, model;
    allocation::ImportParams = ImportParams(),
    counties = nothing,
    startdate::Date = BACK_CALCULATION_START)

    curve = import_curve(delays, model, startdate = startdate)
    weights = tourism_weights_by_date(county_geography(), allocation)
    e = expected_imports(curve, weights, counties = counties)
    return Dict(zip(e.date, e.austria .+ e.italy .+ e.other))
end


# floor plus Bernoulli draw for the fractional part
_stochastic_round(rng, x) = floor(Int, x) + (rand(rng) < x - floor(x) ? 1 : 0)


"""
    import_schedule(curve, weights;
        startdate::Date,
        enddate::Date,
        scale = 1.0,
        counties = nothing,
        seed = 0)

Draws a seeding schedule of imported infections (`scale` x the import curve)
with the county weights per date (e.g. `tourism_weights_by_date`). If
`counties` is given, only imports into these counties are kept (e.g. when
simulating a single state).
"""
function import_schedule(curve::DataFrame, weights::Function;
    startdate::Date,
    enddate::Date,
    scale::Real = 1.0,
    counties = nothing,
    seed = 0)

    rng = Random.Xoshiro(seed)
    keep = isnothing(counties) ? nothing : Set(counties)
    rows = DataFrame(tick = Int[], ags = AGS[], count = Int[])

    for r in eachrow(curve[startdate .<= curve.date .<= enddate, :])
        t = Dates.value(r.date - startdate)
        w = weights(r.date)
        for (n, p) in [(r.austria, w.austria), (r.italy, w.italy), (r.other, w.other)]
            k = _stochastic_round(rng, n * scale)
            k == 0 && continue
            for a in sample(rng, w.ags, Weights(p), k)
                if isnothing(keep) || a in keep
                    push!(rows, (t, a, 1))
                end
            end
        end
    end
    return merge_schedule(rows)
end


"""
    heinsberg_schedule(rki::DataFrame; startdate::Date, cutoff = Date(2020, 2, 29), kwargs...)

Seeding schedule for the Heinsberg cluster: RKI cases of the county with
onset up to `cutoff` (see `rki_seeding_schedule`).
"""
function heinsberg_schedule(rki::DataFrame;
    startdate::Date,
    cutoff::Date = Date(2020, 2, 29),
    kwargs...)

    cases = rki[rki.ags .== HEINSBERG, :]
    return rki_seeding_schedule(cases; startdate = startdate, cutoff = cutoff, kwargs...)
end
