###
### COMMUTER VALIDATION: GEMS COMMUTES VS. BA COMMUTER STATISTICS
###
# Commutes between counties (Kreise) in GEMS, compared with the commuter
# statistics of the Bundesagentur für Arbeit (BA) on 30.06.2019: where the
# commuters of a county go (flows between pairs of counties) and how far
# people commute. Part of the mobility analysis (`mobility`), as
# the population validation is part of the demographics.
# The code is split into four parts that can be reused separately:
#   1. reading the BA tables in `inputs/commuters/` and the county centres
#      of the GEMS shapefile
#   2. aggregating the simulation's commutes per pair of counties
#   3. drawing the panels
#   4. composing and saving the figure
# Style, panels and layout helpers come from the population validation
# (`population_validation.jl`); sizes are given in print units as there.

export commuter_validation
export commuter_validation_data, save_commuter_data, load_commuter_data
export commuter_figure
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
    commuter_validation_data(county_commutes_df)

Data of the commuter validation, GEMS (from `county_commutes`) and the BA
data (`ba_commuters`) in the counties of the GEMS shapefile:
- `counties`: per county its `name` and, for `model` and `reference`, the
  workers by place of residence and place of work, and the in- and
  out-commuters (`model_residence`, ..., `reference_out`)
- `pairs`: commuters per pair of different counties (`home`,
  `destination`) in GEMS (`model`) and in the BA data (`reference`, 0 if not
  listed, i.e. below `BA_MIN_FLOW`), and the distance between the county
  centres (`distance_km`)
"""
function commuter_validation_data(county_commutes_df)
    model_flows = model_county_flows(county_commutes_df)
    reference_counties, reference_flows = ba_commuters()

    cols = [:residence, :workplace, :in_commuters, :out_commuters]
    short = Dict(:residence => "residence", :workplace => "workplace", :in_commuters => "in", :out_commuters => "out")
    prefixed(d, p) = rename(d[:, [:county; cols]], (c => Symbol(p, "_", short[c]) for c in cols)...)
    counties = innerjoin(county_names(), prefixed(model_county_totals(model_flows), "model"),
        prefixed(reference_counties, "reference"), on = :county)
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
        title = "C) Flows between counties")

    d = distance_shares(data)
    p_distance = plot_paired_bars(d.reference, d.model, d.labels,
        title = "D) Commuting distance",
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


# MAPS OF THE COUNTIES

# colors of the maps: inferno (black to yellow), the default color scale of Plots, which the other
# maps of the paper use (`agsmap` without `fillcolor`)
const CV_MAP_COLORS = cgrad(:inferno)

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

"""
    county_map(counties, values, title; clims, extent = map_extent())

Map of the counties (`agsmap` of GEMS), colored by `values` (one per county
key in `counties`) within `clims`, without color bar and axes. The map fills
its plot area (limits `extent`); the layout gives the plot area the aspect
ratio of the map.
"""
county_map(counties, values, title; clims, extent = map_extent()) =
    agsmap(DataFrame(ags = AGS.(1000 .* counties), value = values),
        fillcolor = CV_MAP_COLORS,
        clims = clims,
        colorbar = false,
        linecolor = nothing,
        aspect_ratio = :auto,
        xlims = extent[1:2],
        ylims = extent[3:4],
        framestyle = :none,
        # no room for the (hidden) tick labels of the coordinates
        ticks = nothing,
        legend = false,
        title = title,
        fontfamily = PV_FONT)

"""
    commuter_maps(data)

Maps of the in- and out-commuter shares per county in GEMS (`county_map`):
the share of a county's workers who live in another county, and the share
of its employed residents who work in another county, on a common color
scale (`clims`, from 0 to the largest share, rounded up to 10 %). Drawn
once, as they take long; the figure uses copies.
"""
function commuter_maps(data)
    c = data.counties
    in_share = c.model_in ./ c.model_workplace
    out_share = c.model_out ./ c.model_residence
    clims = (0.0, ceil(10 * maximum([in_share; out_share])) / 10)
    extent = map_extent()
    return (
        in = county_map(c.county, in_share, "A) In-commuters (GEMS)", clims = clims, extent = extent),
        out = county_map(c.county, out_share, "B) Out-commuters (GEMS)", clims = clims, extent = extent),
        clims = clims,
        ratio = (extent[2] - extent[1]) / (extent[4] - extent[3]))
end

"""
    colorbar!(p, clims, width_in, height_in; colors = CV_MAP_COLORS, n = 100, step = 0.2)

Horizontal color bar of the maps in the (inset) subplot `p` (`width_in` ×
`height_in` inches): `n` rectangles from `clims[1]` (left) to `clims[2]`
(right), with tick labels (percent, every `step`) below it.
"""
function colorbar!(p, clims, width_in, height_in; colors = CV_MAP_COLORS, n = 100, step = 0.2)
    lo, hi = clims
    edges = range(lo, hi, length = n + 1)
    for i in 1:n
        c = get(colors, (edges[i] + edges[i + 1]) / 2 - lo, (0, hi - lo))
        plot!(p, rectangles([edges[i]], [edges[i + 1]], [0.0], [1.0]), seriestype = :shape,
            fillcolor = c, linecolor = c, linewidth = 0.5, label = "")
    end
    for v in lo:step:(hi + 1e-9)
        annotate!(p, ((v - lo) / (hi - lo), -PV_TICK_GAP_IN / height_in),
            Plots.text(percent_label(v), PV_SMALL_FONTSIZE, PV_FONT, :center, :top))
    end
    return p
end

# room below the maps for the color bar: gap, bar, gap, tick labels, gap (inches)
colorbar_room_in(height_in) = 0.08 + height_in + PV_TICK_GAP_IN + line_height_in(PV_SMALL_FONTSIZE) + 0.1


# 4. FIGURE

"""
    COMMUTER_PANEL_NOTES

What the values of each panel (A–D) mean, printed in one line below its
title (`panel_header!`), so that the axes only need to name the source.
Details (e.g. the distance between the county centres) belong in the
caption.
"""
const COMMUTER_PANEL_NOTES = [
    "Share of workers who live in another county",
    "Share of residents who work in another county",
    "Per county pair: share of home county's out-commuters",
    "Share of workers by distance to work county",
]

# width of the text of PLOS (inches): the width of the commuter figure
const CV_WIDTH_IN = 5.2

# height of a panel header: title and one line of its note (inches)
const CV_HEADER_IN = line_height_in(PV_HEADING_FONTSIZE) + PV_TITLE_GAP_IN / 2 +
    line_height_in(PV_SMALL_FONTSIZE) + PV_TITLE_GAP_IN

"""
    panel_header!(p, title, note, width_in, height_in)

Draws the header of a panel above its plot area (`width_in` × `height_in`
inches): the `title` (`title_text!`) at the top of the room `CV_HEADER_IN`,
and below it the `note` in small font.
"""
function panel_header!(p, title, note, width_in, height_in)
    title_text!(p, 0.0, 1 + CV_HEADER_IN / height_in, title, width_in, valign = :top)
    y = CV_HEADER_IN - line_height_in(PV_HEADING_FONTSIZE) - PV_TITLE_GAP_IN / 2 - line_height_in(PV_SMALL_FONTSIZE) / 2
    return annotate!(p, (0.0, 1 + y / height_in), Plots.text(note, PV_SMALL_FONTSIZE, PV_FONT, :left, :vcenter))
end

"""
    CommuterLayout

Sizes of the commuter figure (`commuter_figure`), in inches:
- `width_in`: width of the figure (default: the text width, `CV_WIDTH_IN`)
- `map_width_in`: width of the maps; the panels below them are squares of
  the same width (fitted to `width_in`)
- `widths_in`: widths of the two columns; `heights_in`: heights of the two
  rows below their headers, the upper one including the color bar (fitted)
- `colorbar_width_in`, `colorbar_height_in`: size of the color bar, centered
  below the two maps
"""
Base.@kwdef struct CommuterLayout
    width_in::Float64 = CV_WIDTH_IN
    map_width_in::Float64 = 2.0
    widths_in::Vector{Float64} = [2.7, 2.5]
    heights_in::Vector{Float64} = [3.3, 2.7]
    colorbar_width_in::Float64 = 2.2
    colorbar_height_in::Float64 = 0.09
end

with_layout(l::CommuterLayout; changes...) = CommuterLayout(;
    (f => get(changes, f, getfield(l, f)) for f in fieldnames(CommuterLayout))...)

"""
    build_commuter_figure(data, maps, layout::CommuterLayout; dpi = PV_DPI)

Composes the commuter figure with the sizes of `layout`, as a grid of two
columns, so that the plot areas of a column start and end at the same
positions:
- upper row: the maps of the in- and out-commuters in GEMS (copies of
  `maps`, from `commuter_maps`), with their common color bar centered below
  them (an inset, `colorbar!`)
- lower row: the flows between counties and the commuting distance
  (`commuter_panels`), with square plot areas as wide as the maps

Above each plot area its header (`panel_header!`); titles, notes and axis
text are placed at fixed distances (`axis_text!`), as in the population
validation. The panels are marked A–D (`mark_panel!`).
"""
function build_commuter_figure(data, maps, layout::CommuterLayout; dpi = PV_DPI)
    l = layout
    p_pairs, p_distance = commuter_panels(data)
    p_in, p_out = deepcopy(maps.in), deepcopy(maps.out)
    panels = [p_in, p_out, p_pairs, p_distance]
    titles = [p[1][:title] for p in panels]
    foreach(p -> plot!(p, title = ""), panels)
    s = l.map_width_in
    map_h = s / maps.ratio
    width_in = sum(l.widths_in)
    heights = l.heights_in .+ CV_HEADER_IN
    height_in = sum(heights)

    # axis text of the lower panels, legend of the distance bands
    # the y axis label of the flows closer to its tick labels: their width ("0.01%") is overestimated
    # by about 0.05 in (`text_width_in`)
    rooms = [axis_text!(p_pairs, s, s, guide_gap_in = PV_GUIDE_GAP_IN - 0.05), axis_text!(p_distance, s, s)]
    bottom_in = maximum(r -> r.bottom, rooms)
    inset_legend!(p_pairs, getfield.(DISTANCE_BANDS, :label), group_colors(length(DISTANCE_BANDS)), s, s)

    sizes = [(s, map_h), (s, map_h), (s, s), (s, s)]
    for (k, p) in enumerate(panels)
        panel_header!(p, titles[k], COMMUTER_PANEL_NOTES[k], sizes[k]...)
        mark_panel!(p, string('A' + k - 1))
    end

    # margins: Plots adds 2 mm above a plot area without title; the maps leave room for the color
    # bar below them, the lower panels for their axis text (the grid aligns the columns)
    for p in (p_in, p_out)
        plot!(p, top_margin = (25.4CV_HEADER_IN - 2) * mm, bottom_margin = 25.4colorbar_room_in(l.colorbar_height_in) * mm,
            left_margin = 0mm, right_margin = 1mm)
    end
    for (k, p) in enumerate((p_pairs, p_distance))
        plot!(p, top_margin = (25.4CV_HEADER_IN - 2) * mm, right_margin = 1mm)
        # the rotated labels of the distance panel get some extra room below (their length is estimated)
        reserve!(p, (left = rooms[k].left, bottom = bottom_in), width_in, height_in, bottom_cushion_mm = 3)
    end

    figure = plot(panels...,
        layout = grid(2, 2, widths = l.widths_in ./ width_in, heights = heights ./ height_in),
        size = (100width_in, 100height_in),
        dpi = dpi,
        fontfamily = PV_FONT)

    # color bar: an inset below map A, centered below both maps (whose distance is the column width
    # minus the plot width, as the columns are aligned)
    a = findfirst(sp -> any(sr -> sr[:label] == panel_mark("A"), sp.series_list), figure.subplots)
    k = length(figure.subplots) + 1
    x0 = (l.widths_in[1] + s - l.colorbar_width_in) / 2
    plot!(figure, inset = (a, bbox(x0 / s, 1 + 0.08 / map_h, l.colorbar_width_in / s, l.colorbar_height_in / map_h)),
        subplot = k, framestyle = :none, ticks = nothing, legend = false, grid = false,
        xlims = maps.clims, ylims = (0, 1))
    colorbar!(figure[k], maps.clims, l.colorbar_width_in, l.colorbar_height_in)
    return figure
end

"""
    fit_commuter_layout(data, maps; layout = CommuterLayout(), iterations = 12, tolerance_in = 0.005)

Adjusts `layout` to the plot areas that Plots computes for the figure
(`plot_areas`), until all deviations are below `tolerance_in`: the figure
is `width_in` wide, the plot areas of both columns are `map_width_in` wide,
the maps have the aspect ratio of the map and the lower panels are square.
Warns if a note is (by the estimate of `text_width_in`) wider than the room
up to the next header or the edge of the figure.
"""
function fit_commuter_layout(data, maps; layout = CommuterLayout(), iterations = 12, tolerance_in = 0.005)
    l = layout
    for i in 1:iterations
        a = plot_areas(build_commuter_figure(data, maps, l))
        A, B, C, D = (a[k] for k in ("A", "B", "C", "D"))
        s = l.map_width_in
        map_h = s / maps.ratio
        width_in = sum(l.widths_in)
        deviation = maximum(abs, [width_in - l.width_in; A.width - s; B.width - s; C.width - s; D.width - s;
            A.height - map_h; C.height - s; D.height - s; A.left - C.left; B.left - D.left])
        if deviation < tolerance_in
            rooms = [B.left - A.left, width_in - B.left, D.left - C.left, width_in - D.left]
            for (k, n) in enumerate(COMMUTER_PANEL_NOTES)
                text_width_in(n, PV_SMALL_FONTSIZE) <= rooms[k] - 0.1 ||
                    @warn "The note of panel $('A' + k - 1) may not fit into one line: \"$n\""
            end
            return l
        end
        # columns: padding plus plot area; the plot areas share what the padding leaves of the width
        padding = l.widths_in .- [C.width, D.width]
        s_new = (l.width_in - sum(padding)) / 2
        l = with_layout(l,
            map_width_in = s_new,
            widths_in = padding .+ s_new,
            heights_in = l.heights_in .+ [s_new / maps.ratio - A.height, s_new - C.height])
    end
    @warn "The layout of the commuter figure did not converge in $iterations iterations."
    return l
end

"""
    commuter_figure(data; maps = commuter_maps(data), layout = fit_commuter_layout(data, maps), dpi = PV_DPI)

Figure of the commuting in GEMS and its validation, at the text width of
PLOS (`CV_WIDTH_IN`), in two aligned columns: maps of the in- and
out-commuters per county in GEMS (A, B) with their common color bar, and
below them the flows between counties (C) and the commuting distance (D),
both against the BA data. The BA data are named `CV_REFERENCE_NAME`; the
source belongs in the caption.
"""
commuter_figure(data; maps = commuter_maps(data), layout = fit_commuter_layout(data, maps), dpi = PV_DPI) =
    build_commuter_figure(data, maps, layout, dpi = dpi)

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
`load_commuter_data`) to redraw the figure without loading the population,
or `county_commutes_df` to reuse already computed `county_commutes` (as
`mobility` does); otherwise they are computed from `sim` (default: the
German population).

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
        if isnothing(county_commutes_df)
            sim = isnothing(sim) ? Simulation(population = "DE", seed = 0) : sim
            county_commutes_df = county_commutes(mun_commutes(sim))
        end
        data = commuter_validation_data(county_commutes_df)
    end
    save_commuter_data(data, folder)

    p = commuter_figure(data)
    png(p, joinpath(folder, "commuter_validation.png"))
    return p
end
