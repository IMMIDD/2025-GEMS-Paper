###
### COMMUTER VALIDATION: GEMS COMMUTES VS. BA COMMUTER STATISTICS
###
# Commutes between counties (Kreise) in GEMS: the commuting balance per
# county against its population density (map), and the comparison with the
# commuter statistics of the Bundesagentur für Arbeit (BA) on 30.06.2019:
# where the commuters of a county go (flows between pairs of counties) and
# how far people commute. Part of the mobility analysis (`mobility`), as
# the population validation is part of the demographics.
# The code is split into four parts that can be reused separately:
#   1. reading the BA tables in `inputs/commuters/` and the county centres
#      of the GEMS shapefile
#   2. aggregating the simulation's commutes per pair of counties and its
#      population density per county
#   3. drawing the panels and the map
#   4. composing and saving the figure
# Style, panels and layout helpers come from the population validation
# (`population_validation.jl`); sizes are given in print units as there.

export commuter_validation
export commuter_validation_data, model_county_density, save_commuter_data, load_commuter_data
export commuter_figure, balance_map
export ba_county_commuters, ba_county_flows, ba_commuters
export county_centres, county_names, DISTANCE_BANDS, distance_band


# 1. REFERENCE DATA AND GEOGRAPHY

const COMMUTER_DIR = joinpath("inputs", "commuters")
const BA_AREA_FILE = joinpath(COMMUTER_DIR, "ba_commuters_by_area_2019.csv")
const BA_FLOW_FILE = joinpath(COMMUTER_DIR, "ba_commuter_flows_by_county_2019.csv")

"""
    BA_MIN_FLOW

Smallest flow between two counties that the BA lists by name (10
commuters). Smaller flows are only contained in aggregates ("Übrige
Kreise", "Übrige Bundesländer"). The comparison of flows and distances
applies the same threshold to GEMS.
"""
const BA_MIN_FLOW = 10

"""
    COUNTY_MERGERS

Counties of 2019 (the boundaries of the Gesyland population and the BA data)
merged since, as county key (5 digits) => key of the county they belong to
now (the boundaries of the GEMS shapefile): Eisenach (16056) into the
Wartburgkreis (16063) on 01.07.2021.
"""
const COUNTY_MERGERS = Dict(16056 => 16063)

merged_county(k::Integer) = get(COUNTY_MERGERS, k, k)

# county key (5 digits) of an AGS
county_key(a::AGS) = Int(id(a) ÷ 1000)

"""
    ba_count(s)

Value of a cell of the BA extracts: empty or "-" = none (0), "*" =
anonymized (1 or 2, `missing`). No value at county level and no listed flow
between counties is anonymized.
"""
ba_count(s) = ismissing(s) || s == "-" ? 0 : s == "*" ? missing : parse(Int, s)

read_ba_csv(file) = CSV.read(file, DataFrame, comment = "#", types = String)

"""
    ba_county_commuters(file = BA_AREA_FILE)

Employees subject to social insurance per county of 2019 from the BA
"Gemeindedaten": by place of residence (`residence`) and place of work
(`workplace`), and the commuters across the county boundary (`in_commuters`,
`out_commuters`). Landkreise are rows with 5 digits, urban districts
(kreisfreie Städte) appear only as municipality (8 digits).

BA notes: commutes abroad are not included; the districts 16063, 16068 and
16073 (Thüringen) are incomplete by place of residence in the first half of
2019.
"""
function ba_county_commuters(file = BA_AREA_FILE)
    d = read_ba_csv(file)
    districts = Set(d.ags[length.(d.ags) .== 5])
    is_county(a) = length(a) == 5 || (length(a) == 8 && first(a, 5) ∉ districts)
    c = d[is_county.(d.ags), :]
    return DataFrame(
        county = parse.(Int, first.(c.ags, 5)),
        residence = ba_count.(c.residence),
        workplace = ba_count.(c.workplace),
        in_commuters = ba_count.(c.in_commuters),
        out_commuters = ba_count.(c.out_commuters))
end

"""
    ba_county_flows(file = BA_FLOW_FILE)

Out-commuters per county of residence (`home`) and county of work
(`destination`) of 2019, for all flows the BA lists by name (at least
`BA_MIN_FLOW` commuters).
"""
function ba_county_flows(file = BA_FLOW_FILE)
    d = read_ba_csv(file)
    d = d[length.(d.destination) .== 5, :]
    return DataFrame(
        home = parse.(Int, d.home),
        destination = parse.(Int, d.destination),
        commuters = ba_count.(d.total))
end

"""
    merge_counties(counties, flows; county = merged_county)

The BA data in the boundaries of the GEMS shapefile (`COUNTY_MERGERS`):
counts of merged counties are added, and their mutual flows become commutes
within the county, so they leave the in- and out-commuters and the flows.
"""
function merge_counties(counties, flows; county = merged_county)
    f = transform(flows, :home => ByRow(county) => :home, :destination => ByRow(county) => :destination)
    internal = combine(groupby(f[f.home .== f.destination, :], :home), :commuters => sum => :internal)
    c = transform(counties, :county => ByRow(county) => :county)
    c = combine(groupby(c, :county), [:residence, :workplace, :in_commuters, :out_commuters] .=> sum, renamecols = false)
    c = leftjoin(c, rename(internal, :home => :county), on = :county)
    c.internal = coalesce.(c.internal, 0)
    c.in_commuters .-= c.internal
    c.out_commuters .-= c.internal
    f = combine(groupby(f[f.home .!= f.destination, :], [:home, :destination]), :commuters => sum => :commuters)
    return sort!(select!(c, Not(:internal)), :county), f
end

"""
    ba_commuters(; area_file = BA_AREA_FILE, flow_file = BA_FLOW_FILE)

County totals (`ba_county_commuters`) and flows between counties
(`ba_county_flows`) of the BA, in the boundaries of the GEMS shapefile
(`merge_counties`).
"""
ba_commuters(; area_file = BA_AREA_FILE, flow_file = BA_FLOW_FILE) =
    merge_counties(ba_county_commuters(area_file), ba_county_flows(flow_file))

"""
    county_names()

Names of the counties (`GEN` of the GEMS shapefile). Names that two counties
share (e.g. the city and the district of München) are completed by "(city)"
or "(district)".
"""
function county_names()
    shp = DataFrame(germanshapes(2))
    shp = shp[coalesce.(shp.GF .== 4, false), :]
    names = String.(shp.GEN)
    city = in.(String.(shp.BEZ), Ref(("Kreisfreie Stadt", "Stadtkreis", "Stadt")))
    shared = [count(==(n), names) > 1 for n in names]
    label = [s ? "$n ($(c ? "city" : "district"))" : n for (n, c, s) in zip(names, city, shared)]
    return DataFrame(county = parse.(Int, first.(String.(shp.AGS), 5)), name = label)
end

"""
    county_centres()

Population-weighted centre of each county (coordinates of the GEMS
shapefile, Gauss-Krüger zone 3, metres): the centroids of its
municipalities (`_centroid`), weighted by their population (`EWZ`), so that
the centre lies where most people live. Counties in the boundaries of the
shapefile (`COUNTY_MERGERS`).
"""
function county_centres()
    shp = DataFrame(germanshapes(3))
    shp = shp[coalesce.(shp.GF .== 4, false), :]
    xy = _centroid.(shp.geometry)
    d = DataFrame(
        county = merged_county.(parse.(Int, first.(String.(shp.AGS), 5))),
        x = first.(xy), y = last.(xy),
        population = Int.(coalesce.(shp.EWZ, 0)))
    weighted(v, w) = sum(v .* w) / sum(w)
    return combine(groupby(d, :county),
        [:x, :population] => weighted => :x,
        [:y, :population] => weighted => :y,
        :population => sum => :population)
end

"""
    DISTANCE_BANDS

Bands of the distance between the centres of two counties
(`county_centres`), each up to (not including) its `upper` limit in km.
"""
const DISTANCE_BANDS = [
    (label = "<25 km", upper = 25.0),
    (label = "25–50 km", upper = 50.0),
    (label = "50–100 km", upper = 100.0),
    (label = "100–200 km", upper = 200.0),
    (label = "200+ km", upper = Inf),
]

distance_band(km) = findfirst(b -> km < b.upper, DISTANCE_BANDS)


# 2. MODEL DATA

"""
    model_county_flows(county_commutes_df)

Workers per county of residence (`home`) and county of work (`destination`)
in GEMS, from the commutes between counties (`county_commutes`), in the
boundaries of the GEMS shapefile (`COUNTY_MERGERS`). Includes the workers
who work in their home county (`home == destination`).
"""
function model_county_flows(county_commutes_df)
    d = DataFrame(
        home = merged_county.(county_key.(county_commutes_df.hh_ags)),
        destination = merged_county.(county_key.(county_commutes_df.off_ags)),
        commuters = county_commutes_df.commutes)
    return combine(groupby(d, [:home, :destination]), :commuters => sum => :commuters)
end

"""
    model_county_totals(flows)

Workers per county by place of residence and place of work, and the
commuters across the county boundary, from `model_county_flows`.
"""
function model_county_totals(flows)
    crossing = flows[flows.home .!= flows.destination, :]
    by(d, key, name) = rename(combine(groupby(d, key), :commuters => sum => name), key => :county)
    t = outerjoin(by(flows, :home, :residence), by(flows, :destination, :workplace),
        by(crossing, :destination, :in_commuters), by(crossing, :home, :out_commuters), on = :county)
    return sort!(coalesce.(t, 0), :county)
end

"""
    model_county_density(sim)

Population density per county in GEMS (individuals per km²), as GEMS
computes it for its map of the population density (`PopDensityMap`): the
individuals per municipality over the cadastral area of the municipality in
the GEMS shapefile (`region_info`), summed per county, in the boundaries of
the shapefile (`COUNTY_MERGERS`). As there, municipalities without an area
in the shapefile are left out: those merged since 2019, such as Eisenach,
about 0.2 % of the population.
"""
function model_county_density(sim)
    r = region_info(sim)
    r = r[.!ismissing.(r.area) .& coalesce.(r.area .!= 0, false), :]
    d = DataFrame(county = merged_county.(county_key.(r.ags)), individuals = r.pop_size, area = r.area)
    d = combine(groupby(d, :county), :individuals => sum => :individuals, :area => sum => :area)
    return DataFrame(county = d.county, density = d.individuals ./ d.area)
end

"""
    commuter_validation_data(county_commutes_df, density)

Data of the commuter validation, GEMS (from `county_commutes`) and the BA
data (`ba_commuters`) in the counties of the GEMS shapefile:
- `counties`: per county its `name`, the population density in GEMS
  (`model_density`, from `density`, the output of `model_county_density`)
  and, for `model` and `reference`, the workers by place of residence and
  place of work, and the in- and out-commuters (`model_residence`, ...,
  `reference_out`)
- `pairs`: commuters per pair of different counties (`home`,
  `destination`) in GEMS (`model`) and in the BA data (`reference`, 0 if not
  listed, i.e. below `BA_MIN_FLOW`), and the distance between the county
  centres (`distance_km`)
"""
function commuter_validation_data(county_commutes_df, density)
    model_flows = model_county_flows(county_commutes_df)
    reference_counties, reference_flows = ba_commuters()

    cols = [:residence, :workplace, :in_commuters, :out_commuters]
    short = Dict(:residence => "residence", :workplace => "workplace", :in_commuters => "in", :out_commuters => "out")
    prefixed(d, p) = rename(d[:, [:county; cols]], (c => Symbol(p, "_", short[c]) for c in cols)...)
    counties = innerjoin(county_names(), rename(density[:, [:county, :density]], :density => :model_density),
        prefixed(model_county_totals(model_flows), "model"), prefixed(reference_counties, "reference"), on = :county)
    nrow(counties) == length(unique(model_flows.home)) == nrow(reference_counties) ||
        @warn "The counties of GEMS, the BA data and the shapefile differ."

    pairs = outerjoin(
        rename(model_flows[model_flows.home .!= model_flows.destination, :], :commuters => :model),
        rename(reference_flows, :commuters => :reference), on = [:home, :destination])
    pairs.model = coalesce.(pairs.model, 0)
    pairs.reference = coalesce.(pairs.reference, 0)
    centres = Dict(r.county => (r.x, r.y) for r in eachrow(county_centres()))
    distance(a, b) = hypot((centres[a] .- centres[b])...) / 1000
    pairs.distance_km = distance.(pairs.home, pairs.destination)
    return (counties = sort!(counties, :county), pairs = sort!(pairs, [:home, :destination]))
end

const COMMUTER_TABLES = (:counties, :pairs)
commuter_file(folder, table) = joinpath(folder, "commuter_validation_$(table).csv")

"""
    save_commuter_data(data, folder)

Saves the tables of `commuter_validation_data` as CSV files to `folder`.
"""
save_commuter_data(data, folder) =
    foreach(t -> CSV.write(commuter_file(folder, t), getfield(data, t)), COMMUTER_TABLES)

"""
    load_commuter_data(folder)

Reads the tables saved by `save_commuter_data`, to redraw the figure without
loading the population.
"""
load_commuter_data(folder) =
    NamedTuple{COMMUTER_TABLES}(Tuple(CSV.read(commuter_file(folder, t), DataFrame) for t in COMMUTER_TABLES))


# 3. PANELS

# name of the reference data in the figure: all of them come from the BA (the tables are named in the caption)
const CV_REFERENCE_NAME = "Federal Employment Agency (BA)"

# percentage with one significant digit below 1 % (axis ticks spanning several orders of magnitude)
share_label(x) = (100x >= 1 ? string(round(Int, 100x)) : string(round(100x, sigdigits = 1))) * "%"

"""
    compared_pairs(data; min_flow = BA_MIN_FLOW)

The pairs of counties (`data.pairs`) with at least `min_flow` commuters in
GEMS and in the BA data, as the BA lists smaller flows only in aggregates.
"""
compared_pairs(data; min_flow = BA_MIN_FLOW) =
    data.pairs[(data.pairs.model .>= min_flow) .& (data.pairs.reference .>= min_flow), :]

"""
    pair_shares(data; min_flow = BA_MIN_FLOW)

Flows between pairs of counties (`compared_pairs`) as shares of all
out-commuters of the home county, for GEMS (`model`) and the BA data
(`reference`): where the commuters of a county go.
"""
function pair_shares(data; min_flow = BA_MIN_FLOW)
    p, c = compared_pairs(data, min_flow = min_flow), data.counties
    out(source) = getindex.(Ref(Dict(c.county .=> c[!, "$(source)_out"])), p.home)
    return (model = p.model ./ out("model"), reference = p.reference ./ out("reference"),
        distance_km = p.distance_km)
end

"""
    common_part_of_commuters(model, reference)

Common Part of Commuters (CPC) of two sets of flows between the same pairs
of counties (Lenormand et al. 2012, PLoS ONE 7(10): e45985; Lenormand et
al. 2016, J. Transp. Geogr. 51: 158–169): the share of commuters on which
both agree, from 0 (no overlap) to 1 (identical flows). Both sets are
scaled to the same total first, since GEMS counts more workers than the BA
(civil servants and self-employed persons included), so the CPC measures
the distribution of the commuters over the pairs.
"""
common_part_of_commuters(model, reference) = sum(min.(model ./ sum(model), reference ./ sum(reference)))

"""
    distance_shares(data; min_flow = BA_MIN_FLOW)

Shares of workers who work in their home county and who commute over the
distances of `DISTANCE_BANDS` (between the county centres), for GEMS and
the BA data. For both, the shares refer to the workers in their home county
and the flows of at least `min_flow` commuters, the flows the BA lists by
name.
"""
function distance_shares(data; min_flow = BA_MIN_FLOW)
    p, c = data.pairs, data.counties
    band = distance_band.(p.distance_km)
    function shares(source)
        flows = p[!, source]
        same = sum(c[!, "$(source)_residence"]) - sum(c[!, "$(source)_out"])
        by_band = [sum(flows[(flows .>= min_flow) .& (band .== b)]) for b in eachindex(DISTANCE_BANDS)]
        v = [same; by_band]
        return v ./ sum(v)
    end
    return (model = shares("model"), reference = shares("reference"),
        labels = ["Same county"; getfield.(DISTANCE_BANDS, :label)])
end

"""
    commuter_panels(data)

Validation panels of the commuter figure:
- flows between counties as shares of the out-commuters of the home county
  (`pair_shares`), GEMS against the BA data (`plot_agreement`), logarithmic
  axes, colored by distance band (legend in the upper left corner, which
  the points leave free), with the Common Part of Commuters of the flows
  (`common_part_of_commuters`)
- commuting distance (`distance_shares`), as pairs of bars
"""
function commuter_panels(data)
    s = pair_shares(data)
    p = compared_pairs(data)
    cpc = common_part_of_commuters(p.model, p.reference)
    p_pairs = plot_agreement(s.reference, s.model,
        nlabels = 0,
        groups = distance_band.(s.distance_km),
        log_axes = true,
        ticks = 10.0 .^ (-5:0),
        tick_label = share_label,
        summary = "CPC = $(round(cpc, digits = 2))",
        reference_label = CV_REFERENCE_NAME,
        markersize = 1.2,
        markeralpha = 0.35,
        # the many small flows over long distances below the fewer, larger ones nearby
        reverse_groups = true,
        title = "B) Flows between counties")

    d = distance_shares(data)
    p_distance = plot_paired_bars(d.reference, d.model, d.labels,
        title = "C) Commuting distance",
        # short, as the axis of the flow panel names the source in full
        reference_label = "BA",
        legend = :topright,
        xrotation = PV_TICK_ROTATION)
    return p_pairs, p_distance
end

"""
    inset_legend!(p, labels, colors, width_in, height_in; x = 0.04, y = 0.96)

Color legend inside the plot area of `p` (`width_in` × `height_in` inches),
one entry per line from the relative position (`x`, `y`) downward: a marker
in each of `colors` followed by its label. Works on linear and logarithmic
axes.
"""
function inset_legend!(p, labels, colors, width_in, height_in; x = 0.04, y = 0.96)
    axis_value(lims, scale, f) = scale == :log10 ? lims[1] * (lims[2] / lims[1])^f : lims[1] + f * (lims[2] - lims[1])
    sp = p[1]
    data_x(f) = axis_value(xlims(p), sp[:xaxis][:scale], f)
    data_y(f) = axis_value(ylims(p), sp[:yaxis][:scale], f)
    line = line_height_in(PV_SMALL_FONTSIZE) / height_in
    for (k, (l, c)) in enumerate(zip(labels, colors))
        yk = y - (k - 0.5) * line
        # a marker, not a text symbol: GR draws text in the nearest color of its color table
        scatter!(p, [data_x(x + 0.03 / width_in)], [data_y(yk)], markercolor = c, markerstrokewidth = 0,
            markersize = 3, label = "")
        annotate!(p, (x + PV_MARKER_GAP_IN / width_in, yk), Plots.text(l, PV_SMALL_FONTSIZE, PV_FONT, :left, :vcenter))
    end
    return p
end


# MAP OF THE COMMUTING BALANCE

# colors of the map: the population density in grey, from light to a medium grey, so that the peaks
# stand out on the densest counties; net in-commuting in a darker shade of the GEMS blue, net
# out-commuting in the GEMS red
const CV_DENSITY_COLORS = cgrad(["#f4f4f4", "#555555"])
const CV_IN_COLOR = "#00456F"
const CV_OUT_COLOR = gemscolors(3)[3]

# range of the density (individuals per km², logarithmic color scale; values outside are shown at
# its ends) and the values labeled in the legend
const CV_DENSITY_LIMS = (30, 5000)
const CV_DENSITY_TICKS = (30, 300, 3000)

"""
    CV_PEAK

Peaks of the net commuters on the map, in inches, independent of the size
of the map: `height` of the highest peak, `base` height of every peak (so
that small ones stay visible; the height grows linearly with the net
commuters above it), `width` at the foot. The opacity grows from 50 % to
solid at `solid` net commuters, so that the many small peaks recede.
"""
const CV_PEAK = (height = 0.42, base = 0.035, width = 0.042, solid = 100_000)

# net commuters per county in GEMS: workers at the place of work minus workers at the place of
# residence (positive: more in- than out-commuters)
net_commuters(counties) = counties.model_workplace .- counties.model_residence

peak_height(n, max_net) = CV_PEAK.base + (CV_PEAK.height - CV_PEAK.base) * abs(n) / max_net
peak_opacity(n) = 0.5 + 0.5 * min(1.0, abs(n) / CV_PEAK.solid)
peak_color(n) = n > 0 ? CV_IN_COLOR : CV_OUT_COLOR

# profile of a peak: a raised cosine (bell), smooth at the foot and rounded at the top
const PEAK_PROFILE = range(-1, 1, length = 25)

"""
    peak_shape(x, y, n, width, height)

Peak of `width` and `height` with its foot centered at (`x`, `y`), upward
for net in-commuting (`n > 0`) and downward for net out-commuting.
"""
peak_shape(x, y, n, width, height) =
    Shape(x .+ (width / 2) .* PEAK_PROFILE, y .+ sign(n) * height .* (1 .+ cospi.(PEAK_PROFILE)) ./ 2)

"""
    land_shapes()

Counties of the GEMS shapefile with their land areas only (GF 4), keyed by
county (`county`). `agsmap` also draws the water areas of a county, such as
the coastal waters of Hamburg, in its color.
"""
function land_shapes()
    shp = DataFrame(germanshapes(2))
    shp = shp[coalesce.(shp.GF .== 4, false), :]
    return DataFrame(county = parse.(Int, first.(String.(shp.AGS), 5)), geometry = shp.geometry)
end

"""
    balance_map(data, width_in; extent = map_extent(), centres = county_centres())

Map of the commuting balance in GEMS, for a plot area `width_in` inches
wide (with the aspect ratio of `extent`): the counties shaded by their
population density (`model_density`, logarithmic, `CV_DENSITY_COLORS`), and
on the centre of each county (`county_centres`) a peak of its net commuters
(`net_commuters`, `CV_PEAK`), upward in blue for net in-commuting, downward
in red for net out-commuting. Peaks are drawn from north to south, so that
the southern ones lie on top of those they reach.
"""
function balance_map(data, width_in; extent = map_extent(), centres = county_centres())
    c = innerjoin(data.counties, centres[:, [:county, :x, :y]], on = :county)
    net = net_commuters(c)
    max_net = maximum(abs, net)
    lims = log10.(CV_DENSITY_LIMS)
    p = plot(clims = lims, colorbar = false, aspect_ratio = :auto, xlims = extent[1:2], ylims = extent[3:4],
        framestyle = :none, ticks = nothing, legend = false, grid = false, fontfamily = PV_FONT)
    shade = Dict(zip(c.county, clamp.(log10.(c.model_density), lims...)))
    for r in eachrow(land_shapes())
        haskey(shade, r.county) &&
            plot!(p, r.geometry, fill_z = shade[r.county], fillcolor = CV_DENSITY_COLORS, linecolor = nothing, label = "")
    end
    m = (extent[2] - extent[1]) / width_in       # map units (metres) per inch
    for i in sortperm(c.y, rev = true)
        n = net[i]
        n == 0 && continue
        plot!(p, peak_shape(c.x[i], c.y[i], n, m * CV_PEAK.width, m * peak_height(n, max_net)),
            fillcolor = peak_color(n), fillalpha = peak_opacity(n), linecolor = nothing, label = "")
    end
    return p
end

"""
    balance_legend!(figure, k, width_in, height_in, max_net)

Legend of the map in the (inset) subplot `k` of `figure`, in inches
(`width_in` × `height_in`), in one row whose items all end on the line of
the lowest labels: the two directions of the peaks, three peak heights
(10k, 100k, 200k net commuters, at the scale of the map with the highest
peak of `max_net`) and the density scale. Series and annotations are added
to the figure with the subplot index: added through the subplot itself,
they land in another subplot of the combined figure.
"""
function balance_legend!(figure, k, width_in, height_in, max_net)
    plot!(figure, subplot = k, xlims = (0, width_in), ylims = (0, height_in), framestyle = :none, ticks = nothing,
        legend = false, grid = false)
    label!(x, y, s; halign = :left, valign = :vcenter) =
        annotate!(figure, [(x, y, Plots.text(s, PV_SMALL_FONTSIZE, PV_FONT, halign, valign))], subplot = k)
    peak!(x, y, n) = plot!(figure, peak_shape(x, y, n, CV_PEAK.width, peak_height(n, max_net)), subplot = k,
        fillcolor = peak_color(n), fillalpha = peak_opacity(n), linecolor = nothing, label = "")

    # baseline of the peak heights and of the density scale; below it the line of the lowest labels
    base = 0.17
    row = base - PV_TICK_GAP_IN - line_height_in(PV_SMALL_FONTSIZE) / 2
    for (y, n, s) in ((row + 0.21, 50_000, "net in-commuting"), (row, -50_000, "net out-commuting"))
        peak!(0.05, y - sign(n) * peak_height(n, max_net) / 2, n)
        label!(0.13, y, s)
    end
    for (x, n) in zip((1.42, 1.78, 2.14), (10_000, 100_000, 200_000))
        peak!(x, base, n)
        label!(x, base - PV_TICK_GAP_IN, "$(n ÷ 1000)k", halign = :center, valign = :top)
    end

    x0 = 2.48
    label!(x0, base + 0.24, "Individuals per km²")
    lims = log10.(CV_DENSITY_LIMS)
    bar_width, n = width_in - x0 - 0.02, 100
    edges = range(lims..., length = n + 1)
    for i in 1:n
        col = get(CV_DENSITY_COLORS, ((edges[i] + edges[i + 1]) / 2 - lims[1]) / (lims[2] - lims[1]))
        xa = x0 + bar_width * (i - 1) / n
        plot!(figure, Shape([xa, xa + bar_width / n, xa + bar_width / n, xa], [base, base, base + 0.1, base + 0.1]),
            subplot = k, fillcolor = col, linecolor = col, linewidth = 0.5, label = "")
    end
    for v in CV_DENSITY_TICKS
        label!(x0 + bar_width * (log10(v) - lims[1]) / (lims[2] - lims[1]), base - PV_TICK_GAP_IN, string(v),
            halign = :center, valign = :top)
    end
    return figure
end

"""
    map_extent(; pad = 0.01)

Extent of the land area of Germany in the coordinates of the GEMS shapefile
(x0, x1, y0, y1), with a margin of `pad` of its size on each side.
"""
function map_extent(; pad = 0.01)
    # GF 4: land areas (GF 2: water areas, including coastal waters)
    shp = DataFrame(germanshapes(1))
    shp = shp[coalesce.(shp.GF .== 4, false), :]
    xs = [pt.x for g in shp.geometry for pt in g.points]
    ys = [pt.y for g in shp.geometry for pt in g.points]
    dx, dy = maximum(xs) - minimum(xs), maximum(ys) - minimum(ys)
    return (minimum(xs) - pad * dx, maximum(xs) + pad * dx, minimum(ys) - pad * dy, maximum(ys) + pad * dy)
end


# 4. FIGURE

"""
    COMMUTER_PANEL_NOTES

What the values of each panel (A–C) mean, printed in one line below its
title (`panel_header!`), so that the axes only need to name the source.
Details (e.g. the distance between the county centres) belong in the
caption.
"""
const COMMUTER_PANEL_NOTES = [
    "Net commuters and population density per county",
    "Share of home county's out-commuters",
    "Share of workers by distance",
]
const COMMUTER_MAP_TITLE = "A) Commuting balance (GEMS)"

# width of the commuter figure (inches), between the text width (5.2 in) and the print width of PLOS
const CV_WIDTH_IN = 6.0

# height of a panel header: title and one line of its note (inches)
const CV_HEADER_IN = line_height_in(PV_HEADING_FONTSIZE) + PV_TITLE_GAP_IN / 2 +
    line_height_in(PV_SMALL_FONTSIZE) + PV_TITLE_GAP_IN

# legend below the map: its height and its distance from the map (inches)
const CV_LEGEND_IN = 0.56
const CV_LEGEND_GAP_IN = 0.04

"""
    panel_header!(p, title, note, width_in, height_in; offset_in = 0.0)

Draws the header of a panel above its plot area (`width_in` × `height_in`
inches): the `title` (`title_text!`) at the top of the room `CV_HEADER_IN`,
and below it the `note` in small font, both starting `offset_in` inches left
of the plot area.
"""
function panel_header!(p, title, note, width_in, height_in; offset_in = 0.0)
    x = -offset_in / width_in
    title_text!(p, x, 1 + CV_HEADER_IN / height_in, title, width_in, valign = :top)
    y = CV_HEADER_IN - line_height_in(PV_HEADING_FONTSIZE) - PV_TITLE_GAP_IN / 2 - line_height_in(PV_SMALL_FONTSIZE) / 2
    return annotate!(p, (x, 1 + y / height_in), Plots.text(note, PV_SMALL_FONTSIZE, PV_FONT, :left, :vcenter))
end

"""
    CommuterLayout

Sizes of the commuter figure (`commuter_figure`), in inches, measured from
the plot areas that Plots computes (fitted by `fit_commuter_layout`):
- `width_in`: width of the figure (default `CV_WIDTH_IN`)
- `map_padding_in`, `panel_padding_in`: width of the left and right column
  besides the plot area of the map and of the panels
- `row_padding_in`: height of the rows of B and C besides their plot areas
- `map_top_in`: room above the plot area of the map (its header)
- `below_in`: room below the plot area of C (its axis text), which the map
  leaves to its legend
- `map_bottom_mm`: bottom margin of the map, so that its plot area ends
  where the one of C ends
- `header_offset_in`: distance of the plot areas of B and C from the left
  edge of their column, where their headers start
"""
Base.@kwdef struct CommuterLayout
    width_in::Float64 = CV_WIDTH_IN
    map_padding_in::Float64 = 0.05
    panel_padding_in::Float64 = 0.8
    row_padding_in::Vector{Float64} = [0.95, 1.1]
    map_top_in::Float64 = 0.45
    below_in::Float64 = 0.65
    map_bottom_mm::Float64 = 16.5
    header_offset_in::Float64 = 0.6
end

with_layout(l::CommuterLayout; changes...) = CommuterLayout(;
    (f => get(changes, f, getfield(l, f)) for f in fieldnames(CommuterLayout))...)

"""
    commuter_sizes(l::CommuterLayout, ratio)

Sizes that follow from the layout `l` and the aspect ratio of the map
(`ratio`, width / height): the side of the square panels B and C (`side`),
the width of the map (`map_width`) and the height of the figure (`height`),
such that the figure is `l.width_in` wide and the map, from the top of its
plot area to the bottom of the plot area of C, is as high as its width
allows.
"""
function commuter_sizes(l::CommuterLayout, ratio)
    w = l.width_in - l.map_padding_in - l.panel_padding_in
    side = (w / ratio + l.map_top_in + l.below_in - sum(l.row_padding_in)) / (2 + 1 / ratio)
    return (side = side, map_width = w - side, height = 2side + sum(l.row_padding_in))
end

"""
    build_commuter_figure(data, layout::CommuterLayout; extent = map_extent(), dpi = PV_DPI)

Composes the commuter figure with the sizes of `layout` (`commuter_sizes`),
in two columns:
- left: the map of the commuting balance (`balance_map`), its plot area
  ending where the one of C ends, with its legend below it (an inset,
  `balance_legend!`) in the room that the axis text of C takes
- right: the flows between counties (B) above the commuting distance (C)
  (`commuter_panels`), with square plot areas

Above each plot area its header (`panel_header!`), for B and C starting at
the left edge of their column; titles, notes and axis text are placed at
fixed distances (`axis_text!`), as in the population validation. The panels
are marked A–C (`mark_panel!`).
"""
function build_commuter_figure(data, layout::CommuterLayout; extent = map_extent(), dpi = PV_DPI)
    l = layout
    ratio = (extent[2] - extent[1]) / (extent[4] - extent[3])
    z = commuter_sizes(l, ratio)
    s, map_w, height_in = z.side, z.map_width, z.height
    map_h = map_w / ratio

    p_map = balance_map(data, map_w, extent = extent)
    p_pairs, p_distance = commuter_panels(data)
    titles = [COMMUTER_MAP_TITLE; [p[1][:title] for p in (p_pairs, p_distance)]]
    foreach(p -> plot!(p, title = ""), (p_pairs, p_distance))

    # axis text of the panels, legend of the distance bands
    # the y axis label of the flows closer to its tick labels: their width ("0.01%") is overestimated
    # by about 0.05 in (`text_width_in`)
    rooms = [axis_text!(p_pairs, s, s, guide_gap_in = PV_GUIDE_GAP_IN - 0.05), axis_text!(p_distance, s, s)]
    inset_legend!(p_pairs, getfield.(DISTANCE_BANDS, :label), group_colors(length(DISTANCE_BANDS)), s, s)

    for (k, (p, w, h, offset)) in enumerate(((p_map, map_w, map_h, 0.0), (p_pairs, s, s, l.header_offset_in),
            (p_distance, s, s, l.header_offset_in)))
        panel_header!(p, titles[k], COMMUTER_PANEL_NOTES[k], w, h, offset_in = offset)
        mark_panel!(p, string('A' + k - 1))
    end

    # margins: Plots adds 2 mm above a plot area without title; the rotated labels of C get some extra
    # room below (their length is estimated), the right edge some room for the axis label of B
    top = (25.4CV_HEADER_IN - 2) * mm
    plot!(p_map, top_margin = top, bottom_margin = l.map_bottom_mm * mm, left_margin = 0mm, right_margin = 0mm)
    for (k, p) in enumerate((p_pairs, p_distance))
        plot!(p, top_margin = top, right_margin = 4mm)
        reserve!(p, rooms[k], l.width_in, height_in, cushion_mm = 3, bottom_cushion_mm = 4.5)
    end

    right = plot(p_pairs, p_distance, layout = grid(2, 1, heights = (s .+ l.row_padding_in) ./ height_in))
    figure = plot(p_map, right,
        layout = grid(1, 2, widths = [map_w + l.map_padding_in, s + l.panel_padding_in] ./ l.width_in),
        size = (100l.width_in, 100height_in),
        dpi = dpi,
        fontfamily = PV_FONT)

    # legend below the map
    a = findfirst(sp -> any(sr -> sr[:label] == panel_mark("A"), sp.series_list), figure.subplots)
    k = length(figure.subplots) + 1
    plot!(figure, inset = (a, bbox(0.0, 1 + CV_LEGEND_GAP_IN / map_h, 1.0, CV_LEGEND_IN / map_h)), subplot = k)
    balance_legend!(figure, k, map_w, CV_LEGEND_IN, maximum(abs, net_commuters(data.counties)))
    return figure
end

"""
    fit_commuter_layout(data; layout = CommuterLayout(), extent = map_extent(), iterations = 12, tolerance_in = 0.005)

Adjusts `layout` to the plot areas that Plots computes for the figure
(`plot_areas`), until all deviations are below `tolerance_in`: the plot
areas have the sizes of `commuter_sizes`, B and C start at the same
position, and the map ends where C ends. Warns if a note is (by the
estimate of `text_width_in`) wider than the room up to the edge of its
column, or if the legend of the map does not fit below it.
"""
function fit_commuter_layout(data; layout = CommuterLayout(), extent = map_extent(), iterations = 12, tolerance_in = 0.005)
    l = layout
    ratio = (extent[2] - extent[1]) / (extent[4] - extent[3])
    for _ in 1:iterations
        z = commuter_sizes(l, ratio)
        a = plot_areas(build_commuter_figure(data, l, extent = extent))
        A, B, C = (a[k] for k in ("A", "B", "C"))
        left_width = z.map_width + l.map_padding_in
        map_bottom, c_bottom = A.top + A.height, C.top + C.height
        deviation = maximum(abs, [A.width - z.map_width, A.height - z.map_width / ratio,
            B.width - z.side, B.height - z.side, C.width - z.side, C.height - z.side, B.left - C.left,
            map_bottom - c_bottom, B.left - left_width - l.header_offset_in])
        if deviation < tolerance_in
            rooms = [z.map_width, l.width_in - left_width, l.width_in - left_width]
            for (k, n) in enumerate(COMMUTER_PANEL_NOTES)
                text_width_in(n, PV_SMALL_FONTSIZE) <= rooms[k] - 0.1 ||
                    @warn "The note of panel $('A' + k - 1) may not fit into one line: \"$n\""
            end
            CV_LEGEND_GAP_IN + CV_LEGEND_IN <= l.below_in ||
                @warn "The legend of the map does not fit into the room below it."
            return l
        end
        l = with_layout(l,
            map_padding_in = left_width - A.width,
            panel_padding_in = (z.side + l.panel_padding_in) - max(B.width, C.width),
            row_padding_in = (z.side .+ l.row_padding_in) .- [B.height, C.height],
            map_top_in = A.top,
            below_in = z.height - c_bottom,
            map_bottom_mm = l.map_bottom_mm + 25.4 * (map_bottom - c_bottom),
            header_offset_in = B.left - left_width)
    end
    @warn "The layout of the commuter figure did not converge in $iterations iterations."
    return l
end

"""
    commuter_figure(data; layout = fit_commuter_layout(data), dpi = PV_DPI)

Figure of the commuting in GEMS and its validation (`CV_WIDTH_IN` wide), in
two columns: the map of the commuting balance and the population density per
county in GEMS (A, `balance_map`), and beside it the flows between counties
(B) and the commuting distance (C), both against the BA data. The BA data
are named `CV_REFERENCE_NAME`; the source belongs in the caption.
"""
commuter_figure(data; layout = fit_commuter_layout(data), dpi = PV_DPI) =
    build_commuter_figure(data, layout, dpi = dpi)

"""
    commuter_validation(;
        folder = nothing,
        folder_tag = nothing,
        sim = nothing,
        data = nothing,
        county_commutes_df = nothing)

Compares the commutes between counties in GEMS with the BA commuter
statistics of 30.06.2019 and saves the figure (`commuter_figure`, PNG) and
its data (CSV) to `folder`. Pass `data` (from `commuter_validation_data` or
`load_commuter_data`) to redraw the figure without loading the population.
Otherwise the data come from `sim` (default: the German population), with
its population density (`model_county_density`); pass `county_commutes_df`
to reuse already computed `county_commutes` (as `mobility` does).

If `folder` is `nothing`, a new result folder is created, with
`folder_tag` appended to its timestamp-based name (see `result_folder`).
"""
function commuter_validation(;
        folder = nothing,
        folder_tag = nothing,
        sim = nothing,
        data = nothing,
        county_commutes_df = nothing)

    isnothing(folder) ? folder = result_folder(folder_tag) : mkpath(folder)

    if isnothing(data)
        sim = isnothing(sim) ? Simulation(population = "DE", seed = 0) : sim
        isnothing(county_commutes_df) && (county_commutes_df = county_commutes(mun_commutes(sim)))
        data = commuter_validation_data(county_commutes_df, model_county_density(sim))
    end
    save_commuter_data(data, folder)

    p = commuter_figure(data)
    png(p, joinpath(folder, "commuter_validation.png"))
    return p
end
