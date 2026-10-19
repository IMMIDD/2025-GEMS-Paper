###
### IMPORT MODEL (FIRST WAVE 2020)
###

# Only national import figures are used (there is no public county-level data);
# imports are distributed across counties by a distance kernel.

export RKI_IMPORTS
export ISCHGL_LONLAT
export county_geography
export import_curve
export import_weights
export import_schedule
export heinsberg_schedule


"""
    RKI_IMPORTS

Cases with infection country abroad per ISO reporting week 2020.

- `abroad`: weeks 10-17 from the RKI situation report of 27.04.2020 (Fig. 3).
  Week 9 is extrapolated from the reports of 04.-08.03. (64, 130, 223 cases).
- `alpine_share`: share of Austria and Italy, from the country tables of the
  reports up to 08.03. (0.90), 08.-22.03. (0.90), and 22.-26.03. (0.84).

The RKI counts mentions; a case can name several countries.
"""
const RKI_IMPORTS = DataFrame(
    week = 9:17,
    abroad = [35, 262, 2781, 6454, 3814, 1095, 267, 113, 31],
    alpine_share = [0.90, 0.90, 0.90, 0.90, 0.84, 0.84, 0.84, 0.84, 0.84])


"""
    ISCHGL_LONLAT

Coordinates of Ischgl (reference point of the import kernel).
"""
const ISCHGL_LONLAT = (10.2919, 47.0128)


# Monday of ISO week w in 2020
_monday_2020(w::Integer) = Date(2019, 12, 30) + Week(w - 1)


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


"""
    county_geography()

Population (`EWZ` of the shapefile) and distance to Ischgl (km, from the
centroid) for all German counties. Does not need a population model.
"""
function county_geography()
    shp = DataFrame(germanshapes(2)) |>
        d -> d[coalesce.(d.GF .== 4, false), :]   # GF 2 = water areas

    to_gk3 = Proj.Transformation("EPSG:4326", "EPSG:31467", always_xy = true)
    ix, iy = to_gk3(ISCHGL_LONLAT)
    centroids = _centroid.(shp.geometry)

    return DataFrame(
        ags = AGS.(String.(shp.AGS_0)),
        population = Int.(shp.EWZ),
        distance_km = [hypot(x - ix, y - iy) / 1000 for (x, y) in centroids]) |>
        d -> sort(d, :ags, by = a -> a.id)
end


"""
    import_curve(delays::DataFrame; incubation = 7, imports = RKI_IMPORTS)

Expected imported infections per day (`alpine` and `other`). Weekly imports
are spread evenly over the reporting week and shifted back by the reporting
delay (`rki_onset_delays`) and the incubation period.
"""
function import_curve(delays::DataFrame; incubation::Int = 7, imports::DataFrame = RKI_IMPORTS)
    p = delays.cases ./ sum(delays.cases)
    curve = Dict{Date, Tuple{Float64, Float64}}()
    for r in eachrow(imports), day in 0:6
        report = _monday_2020(r.week) + Day(day)
        for (delay, pd) in zip(delays.delay, p)
            d = report - Day(delay + incubation)
            n = r.abroad / 7 * pd
            a, o = get(curve, d, (0.0, 0.0))
            curve[d] = (a + n * r.alpine_share, o + n * (1 - r.alpine_share))
        end
    end
    return DataFrame(date = collect(keys(curve)),
            alpine = first.(values(curve)),
            other = last.(values(curve))) |>
        d -> sort(d, :date)
end


"""
    import_weights(geo::DataFrame; decay = 1.0)

Probability that an import lands in a county: `alpine` ∝ population ×
(distance to Ischgl / 100 km)^(-decay), `other` ∝ population.
"""
function import_weights(geo::DataFrame; decay::Real = 1.0)
    alpine = geo.population .* (geo.distance_km ./ 100) .^ (-decay)
    return DataFrame(ags = geo.ags,
        alpine = alpine ./ sum(alpine),
        other = geo.population ./ sum(geo.population))
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

Draws a seeding schedule of imported infections from the import curve and
the county weights. If `counties` is given, only imports into these counties
are kept (e.g. when simulating a single state).
"""
function import_schedule(curve::DataFrame, weights::DataFrame;
    startdate::Date,
    enddate::Date,
    scale::Real = 1.0,
    counties = nothing,
    seed = 0)

    rng = GEMS.Xoshiro(seed)
    keep = isnothing(counties) ? nothing : Set(counties)
    rows = DataFrame(tick = Int[], ags = AGS[], count = Int[])

    for r in eachrow(curve[startdate .<= curve.date .<= enddate, :])
        t = Dates.value(r.date - startdate)
        for (n, w) in [(r.alpine, weights.alpine), (r.other, weights.other)]
            k = _stochastic_round(rng, n * scale)
            k == 0 && continue
            for a in sample(rng, weights.ags, Weights(w), k)
                (isnothing(keep) || a in keep) && push!(rows, (t, a, 1))
            end
        end
    end

    return rows |>
        d -> groupby(d, [:tick, :ags]) |>
        d -> combine(d, :count => sum => :count) |>
        d -> sort(d, [:tick, order(:ags, by = a -> a.id)])
end


"""
    heinsberg_schedule(rki::DataFrame; startdate::Date, cutoff = Date(2020, 2, 29), kwargs...)

Seeding schedule for the Heinsberg cluster: RKI cases of the county with
onset up to `cutoff` (see `rki_seeding_schedule`).
"""
heinsberg_schedule(rki::DataFrame; startdate::Date, cutoff::Date = Date(2020, 2, 29), kwargs...) =
    rki_seeding_schedule(rki[rki.county .== 5370, :]; startdate = startdate, cutoff = cutoff, kwargs...)
