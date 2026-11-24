###
### DEMOGRAPHIC MAPS: POPULATION OF GEMS PER COUNTY
###
# Maps of the population model per county (population, households, age). Part
# of the demographics (`demographics`).
# The code is split into three parts that can be reused separately:
#   1. aggregating the population per county
#   2. drawing a map and its colorbar
#   3. composing and saving the figure
# Style and text helpers come from the population validation
# (`population_validation.jl`), the county shapes and the header of a panel from
# the commuter validation (`commuter_validation.jl`). Sizes are given in print
# units as there.

export demographic_maps
export county_demographics, load_county_demographics
export demographic_maps_figure
export CountyMap, DEMOGRAPHIC_MAPS


# 1. DATA

# households with more members (e.g. care homes) are left out of the mean household
# size, as in the map of GEMS (`HouseholdSizeMap`)
const DM_MAX_HOUSEHOLD_SIZE = 10

"""
    county_demographics(sim; individual_data_df = individual_data(sim))

Population of `sim` per county, with the definitions of the maps of GEMS
(`PopDensityMap`, `HouseholdSizeMap`, `SinglesMap`, `AgeMap`, `KidsMap`,
`ElderlyMap`) in the boundaries of the GEMS shapefile (`COUNTY_MERGERS`):
- `density`: individuals per km² (`model_county_density`)
- `household_size`: mean size of the households with up to
  `DM_MAX_HOUSEHOLD_SIZE` members
- `single_households`: share of single-person households
- `mean_age`: mean age (years)
- `children`, `elderly`: shares of the individuals aged 0–14 and 65 or older
"""
function county_demographics(sim; individual_data_df = individual_data(sim))
    hh = DataFrame(county = merged_county.(county_key.(ags.(households(sim)))), size = Int.(size.(households(sim))))
    hh = combine(groupby(hh, :county),
        :size => (s -> mean(filter(<=(DM_MAX_HOUSEHOLD_SIZE), s))) => :household_size,
        :size => (s -> count(==(1), s) / length(s)) => :single_households)

    ind = DataFrame(county = merged_county.(county_key.(individual_data_df.ags)), age = Int.(individual_data_df.age))
    ind = combine(groupby(ind, :county),
        :age => mean => :mean_age,
        :age => (a -> count(<=(14), a) / length(a)) => :children,
        :age => (a -> count(>=(65), a) / length(a)) => :elderly)

    return sort!(innerjoin(model_county_density(sim), hh, ind, on = :county), :county)
end

const DM_DATA_FILE = "demographic_maps.csv"

"""
    load_county_demographics(folder)

Reads the table of `county_demographics` saved by `demographic_maps`, to redraw
the figure without loading the population.
"""
load_county_demographics(folder) = CSV.read(joinpath(folder, DM_DATA_FILE), DataFrame)


# 2. MAPS

"""
    CountyMap

A map of the figure: the `column` of `county_demographics`, the panel `title`
and its `note` (what the values mean), whether the values are shares, labeled in
percent (`share`), and whether the color scale is logarithmic (`log_scale`).
"""
Base.@kwdef struct CountyMap
    column::Symbol
    title::String
    note::String
    share::Bool = false
    log_scale::Bool = false
end

const DEMOGRAPHIC_MAPS = [
    CountyMap(column = :density, title = "Population density", note = "Individuals per km² (log scale)",
        log_scale = true),
    CountyMap(column = :household_size, title = "Household size", note = "Mean persons per household"),
    CountyMap(column = :single_households, title = "Single-person households", note = "Share of households",
        share = true),
    CountyMap(column = :mean_age, title = "Age", note = "Mean age in years"),
    CountyMap(column = :children, title = "Children", note = "Share of population aged 0–14", share = true),
    CountyMap(column = :elderly, title = "Older people", note = "Share of population aged 65+", share = true),
]

# colors of the population model: grey, as the population density in the commuter figure (inferno is left
# to the outcomes of simulations), up to a darker grey than there, where it is the background of the peaks
const DM_COLORS = cgrad([first(color_list(CV_DENSITY_COLORS)), "#1c1c1c"])

# steps between the ticks of a linear color scale (1, 2 or 5 times a power of ten) and their preference
const DM_TICK_STEPS = [(1.0, 1.0), (2.0, 0.9), (5.0, 0.9)]

color_value(m::CountyMap, v) = m.log_scale ? log10(v) : Float64(v)
color_limits(m::CountyMap, data) = extrema(color_value.(Ref(m), data[!, m.column]))

# values labeled at the colorbar of a map with the color limits `lims`; on a logarithmic scale 1, 2 and 5
# times a power of ten or, if fewer than three of them lie within the limits, all multiples
function color_ticks(m::CountyMap, lims)
    inside(ticks) = [t for t in ticks if lims[1] <= color_value(m, t) <= lims[2]]
    m.log_scale || return inside(first(Plots.PlotUtils.optimize_ticks(lims..., k_min = 3, k_max = 7, Q = DM_TICK_STEPS)))
    ticks(steps) = inside([k * 10.0^e for e in 0:4 for k in steps])
    return length(ticks((1, 2, 5))) >= 3 ? ticks((1, 2, 5)) : ticks(1:9)
end

# labels of the `ticks` of a map, with as many decimals as they need
function color_labels(m::CountyMap, ticks)
    v = m.share ? 100 .* ticks : ticks
    digits = something(findfirst(d -> all(isapprox.(v, round.(v, digits = d), rtol = 1e-6)), 0:3), 4) - 1
    number(x) = digits == 0 ? string(round(Int, x)) : string(round(x, digits = digits))
    return number.(v) .* (m.share ? "%" : "")
end

# extent of shapes in map coordinates (x0, x1, y0, y1)
function shape_extent(geometries)
    xs = [pt.x for g in geometries for pt in g.points]
    ys = [pt.y for g in geometries for pt in g.points]
    return (minimum(xs), maximum(xs), minimum(ys), maximum(ys))
end

"""
    county_map!(figure, k, shapes, values, lims, extent; colors = DM_COLORS)

Draws the counties of `shapes` (`land_shapes`) into subplot `k` of `figure`,
each in the color of its value in `values` (by county), with the color limits
`lims`.
"""
function county_map!(figure, k, shapes, values, lims, extent; colors = DM_COLORS)
    plot!(figure, subplot = k, clims = lims, colorbar = false, aspect_ratio = :auto, xlims = extent[1:2],
        ylims = extent[3:4], framestyle = :none, ticks = nothing, legend = false, grid = false)
    for r in eachrow(shapes)
        plot!(figure, r.geometry, subplot = k, fill_z = values[r.county], fillcolor = colors, linecolor = nothing,
            label = "")
    end
    return figure
end

"""
    county_colorbar!(figure, k, m::CountyMap, lims, x, y, width, height; colors = DM_COLORS)

Draws a vertical colorbar with labeled ticks for the map `m` with the color
limits `lims` into subplot `k` of `figure` (whose data coordinates are inches), at
the position (`x`, `y`, lower left corner).
"""
function county_colorbar!(figure, k, m::CountyMap, lims, x, y, width, height; colors = DM_COLORS)
    # a heatmap: shapes of a single color each would fill up the color table of GR over several figures
    edges = range(0, 1, length = 257)
    heatmap!(figure, [x, x + width], y .+ height .* edges, reshape((edges[1:end-1] .+ edges[2:end]) ./ 2, :, 1),
        subplot = k, color = colors, clims = (0, 1), colorbar = false)
    ticks = color_ticks(m, lims)
    for (v, label) in zip(ticks, color_labels(m, ticks))
        t = y + height * (color_value(m, v) - lims[1]) / (lims[2] - lims[1])
        plot!(figure, [x + width, x + width + DM_TICK_IN], [t, t], subplot = k, linecolor = :black, linewidth = 0.5,
            label = "")
        annotate!(figure, [(x + width + DM_TICK_IN + PV_TICK_GAP_IN / 2, t,
            Plots.text(label, PV_SMALL_FONTSIZE, PV_FONT, :left, :vcenter))], subplot = k)
    end
    return figure
end


# 3. FIGURE

# layout (inches): maps per row, gaps between the columns and rows, distance of the colorbar from the
# map, its width, its ticks and the room for their labels, and its share of the height of the map
const DM_COLUMNS = 3
const DM_COLUMN_GAP_IN = 0.12
const DM_ROW_GAP_IN = 0.1
const DM_COLORBAR_GAP_IN = 0.03
const DM_COLORBAR_IN = 0.09
const DM_TICK_IN = 0.03
const DM_LABEL_IN = 0.33
const DM_COLORBAR_SHARE = 0.75

# padding that Plots adds around the plot area of a subplot without axes (removed, so that the plot areas
# have the given sizes)
const DM_PADDING_MM = 2

"""
    demographic_maps_figure(data; maps = DEMOGRAPHIC_MAPS, width_in = PV_WIDTH_IN, dpi = PV_DPI)

Figure of the `maps` of `data` (`county_demographics`), `DM_COLUMNS` per row
and `width_in` inches wide (fonts as in the population validation figure).
Above each map its letter, title and note (`panel_header!`), and to its right
its colorbar, which spans the values of all counties. Shows the counties in
`data` with their land areas (`land_shapes`).
"""
function demographic_maps_figure(data; maps = DEMOGRAPHIC_MAPS, width_in = PV_WIDTH_IN, dpi = PV_DPI)
    shapes = land_shapes()
    shapes = shapes[in.(shapes.county, Ref(Set(data.county))), :]
    extent = shape_extent(shapes.geometry)

    column = (width_in - (DM_COLUMNS - 1) * DM_COLUMN_GAP_IN) / DM_COLUMNS
    map_w = column - DM_COLORBAR_GAP_IN - DM_COLORBAR_IN - DM_TICK_IN - DM_LABEL_IN
    map_h = map_w * (extent[4] - extent[3]) / (extent[2] - extent[1])
    rows = cld(length(maps), DM_COLUMNS)
    height_in = rows * (CV_HEADER_IN + map_h) + (rows - 1) * DM_ROW_GAP_IN

    # the maps are insets of an empty subplot that fills the figure; the headers and colorbars are drawn
    # last, into a canvas in inches above them (GR draws text in the nearest color of its color table,
    # which the maps change)
    blank = (framestyle = :none, ticks = nothing, margin = -DM_PADDING_MM * mm)
    area = (; blank..., legend = false, grid = false, colorbar = false, widen = false,
        xlims = (0, width_in), ylims = (0, height_in), background_color_subplot = :transparent,
        background_color_inside = :transparent)
    figure = plot(; area..., size = (100width_in, 100height_in), dpi = dpi, fontfamily = PV_FONT)
    # the position of a map on the canvas: lower left corner
    corner(i) = ((column + DM_COLUMN_GAP_IN) * ((i - 1) % DM_COLUMNS),
        height_in - (CV_HEADER_IN + map_h) - ((i - 1) ÷ DM_COLUMNS) * (CV_HEADER_IN + map_h + DM_ROW_GAP_IN))
    for i in eachindex(maps)
        x, y = corner(i)
        plot!(figure; inset = (1, bbox(x / width_in, 1 - (y + map_h) / height_in, map_w / width_in,
            map_h / height_in)), subplot = i + 1, blank...)
    end
    canvas = length(maps) + 2
    plot!(figure; inset = (1, bbox(0, 0, 1, 1)), subplot = canvas, area...)

    for (i, m) in enumerate(maps)
        x, y = corner(i)
        lims = color_limits(m, data)
        values = Dict(zip(data.county, color_value.(Ref(m), data[!, m.column])))
        county_map!(figure, i + 1, shapes, values, lims, extent)
        county_colorbar!(figure, canvas, m, lims, x + map_w + DM_COLORBAR_GAP_IN, y + map_h * (1 - DM_COLORBAR_SHARE) / 2,
            DM_COLORBAR_IN, map_h * DM_COLORBAR_SHARE)
        county_map_header!(figure, canvas, m, 'A' + i - 1, x, y + map_h)
    end
    return figure
end

# header of a map (as `panel_header!`) in subplot `k` of `figure` (in inches), above its upper left
# corner (x, y)
function county_map_header!(figure, k, m::CountyMap, letter, x, y)
    top = y + CV_HEADER_IN
    title_y = top - line_height_in(PV_HEADING_FONTSIZE) / 2
    note_y = top - line_height_in(PV_HEADING_FONTSIZE) - PV_TITLE_GAP_IN / 2 - line_height_in(PV_SMALL_FONTSIZE) / 2
    dx = text_width_in("$letter)", PV_HEADING_FONTSIZE, bold = true) + 0.3PV_HEADING_FONTSIZE / 72
    text!(x, y, s, size, font) = annotate!(figure, [(x, y, Plots.text(s, size, font, :left, :vcenter))], subplot = k)
    text!(x, title_y, "$letter)", PV_HEADING_FONTSIZE, PV_FONT_BOLD)
    text!(x + dx, title_y, m.title, PV_HEADING_FONTSIZE, PV_FONT)
    text!(x, note_y, m.note, PV_SMALL_FONTSIZE, PV_FONT)
    return figure
end

"""
    demographic_maps(;
        folder = nothing,
        folder_tag = nothing,
        sim = nothing,
        data = nothing,
        individual_data_df = nothing)

Draws the demographic maps of the population (`demographic_maps_figure`) and
saves the figure (`demographic_maps.png`) and its data (`demographic_maps.csv`)
to `folder`. Pass `data` (from `county_demographics` or
`load_county_demographics`) to redraw the figure without loading the
population. Otherwise the data come from `sim` (default: the German
population). Pass `individual_data_df` to reuse an already computed
`individual_data(sim)`.

If `folder` is `nothing`, a new result folder is created, with `folder_tag`
appended to its timestamp-based name (see `result_folder`).
"""
function demographic_maps(;
        folder = nothing,
        folder_tag = nothing,
        sim = nothing,
        data = nothing,
        individual_data_df = nothing)

    isnothing(folder) ? folder = result_folder(folder_tag) : mkpath(folder)

    if isnothing(data)
        sim = isnothing(sim) ? Simulation(population = "DE", seed = 0) : sim
        isnothing(individual_data_df) && (individual_data_df = individual_data(sim))
        data = county_demographics(sim, individual_data_df = individual_data_df)
    end
    CSV.write(joinpath(folder, DM_DATA_FILE), data)

    p = demographic_maps_figure(data)
    png(p, joinpath(folder, "demographic_maps.png"))
    return p
end
