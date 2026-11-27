###
### SETTINGS MAP: HOUSEHOLDS, SCHOOLS, AND WORKPLACES OF GEMS
###
# Map of all settings of the population model. Part of the demographics
# (`demographics`).
# Instead of drawing millions of points, the settings are counted per pixel of
# the figure, and the pixels are colored.
# The code is split into three parts that can be reused separately:
#   1. counting the settings per pixel
#   2. coloring the pixels
#   3. composing and saving the figure
# Style and images come from `figure_style.jl`, the projection of the locations
# from `settings.jl`, the county shapes from the commuter validation
# (`commuter_validation.jl`) and the demographic maps (`demographic_maps.jl`).

export settings_map
export settings_per_pixel, load_settings_per_pixel
export settings_map_figure
export household_gradient, household_opacity


# 1. DATA

# width of the figure (pixels), margin around the land area (share of its size)
const SM_WIDTH = fig_px(FIG_TEXT_WIDTH_IN)
const SM_PADDING = 0.01

# locations projected at a time
const SM_CHUNK = 1_000_000

const SM_COUNTS = (:households, :schoolclasses, :offices)

# extent of the land area of the `counties` with a margin (x0, x1, y0, y1)
function settings_extent(counties)
    shapes = land_shapes()
    x0, x1, y0, y1 = shape_extent(shapes.geometry[in.(shapes.county, Ref(counties))])
    dx, dy = SM_PADDING * (x1 - x0), SM_PADDING * (y1 - y0)
    return (x0 - dx, x1 + dx, y0 - dy, y1 + dy)
end

# number of `settings` per pixel (row, column) of a map of `extent` with the given size; settings
# outside of `counties` are left out
function count_per_pixel(settings, counties, extent, height, width)
    counts = zeros(Int, height, width)
    outside = 0
    for chunk in Iterators.partition(settings, SM_CHUNK)
        s = filter(s -> merged_county(county_key(ags(s))) in counties && !isnan(lon(s)) && !isnan(lat(s)), chunk)
        isempty(s) && continue
        points = convert_points(DataFrame(lon = lon.(s), lat = lat.(s)))
        for (x, y) in zip(points.lon, points.lat)
            column = floor(Int, (x - extent[1]) / (extent[2] - extent[1]) * width) + 1
            row = floor(Int, (extent[4] - y) / (extent[4] - extent[3]) * height) + 1
            checkbounds(Bool, counts, row, column) ? (counts[row, column] += 1) : (outside += 1)
        end
    end
    outside > 0 && @warn "$outside settings lie outside of the map."
    return counts
end

"""
    settings_per_pixel(sim)

Number of households, school classes, and offices of `sim` per pixel of the map (matrices, rows from
north to south), `SM_WIDTH` pixels wide and as high as the land area of the counties with households
requires. Uses all settings in these counties (smaller populations also contain schools and offices
outside of them).
"""
function settings_per_pixel(sim)
    counties = Set(merged_county(county_key(a)) for a in unique(ags.(households(sim))))
    extent = settings_extent(counties)
    height = round(Int, SM_WIDTH * (extent[4] - extent[3]) / (extent[2] - extent[1]))
    count_settings(settings) = count_per_pixel(settings, counties, extent, height, SM_WIDTH)
    return (households = count_settings(households(sim)), schoolclasses = count_settings(schoolclasses(sim)),
        offices = count_settings(offices(sim)))
end

const SM_DATA_FILE = "settings_map.csv.gz"

# saves the counts of all pixels with settings (and of the last pixel, which holds the size of the map)
function save_settings_per_pixel(data, folder)
    last_pixel = lastindex(data.households)
    pixels = [i for i in eachindex(data.households) if i == last_pixel || any(c -> data[c][i] > 0, SM_COUNTS)]
    index = CartesianIndices(data.households)[pixels]
    table = DataFrame(row = getindex.(index, 1), column = getindex.(index, 2))
    for c in SM_COUNTS
        table[!, c] = data[c][pixels]
    end
    return CSV.write(joinpath(folder, SM_DATA_FILE), table, compress = true)
end

"""
    load_settings_per_pixel(folder)

Reads the counts of `settings_per_pixel` saved by `settings_map`, to redraw the figure without
loading the population.
"""
function load_settings_per_pixel(folder)
    table = CSV.read(joinpath(folder, SM_DATA_FILE), DataFrame)
    function count_matrix(c)
        counts = zeros(Int, maximum(table.row), maximum(table.column))
        counts[CartesianIndex.(table.row, table.column)] .= table[!, c]
        return counts
    end
    return NamedTuple{SM_COUNTS}(count_matrix.(SM_COUNTS))
end


# 2. COLORS

# colors of schools and workplaces (red, green, blue)
function sm_rgb(hex)
    c = parse(Plots.Colors.RGB{Float64}, hex)
    return (c.r, c.g, c.b)
end
const SM_SCHOOL_COLOR = sm_rgb(FIG_COLORS.blue)
const SM_WORKPLACE_COLOR = sm_rgb(FIG_COLORS.green)
sm_setting_color(schoolclasses) = schoolclasses > 0 ? SM_SCHOOL_COLOR : SM_WORKPLACE_COLOR

# households per pixel from which on they cover schools and workplaces
const SM_MANY_HOUSEHOLDS = 150

# gray (0: black) of a pixel with a single household, number of households from which on pixels are
# black, and exponent of the gradient (larger: sparse pixels differ less)
const SM_GRAY = 0.3
const SM_BLACK_HOUSEHOLDS = 1500
const SM_EXPONENT = 2.5
household_gray(n) = SM_GRAY * (1 - min(1.0, log(n) / log(SM_BLACK_HOUSEHOLDS))^SM_EXPONENT)

"""
    household_gradient(households, schoolclasses, offices)

Color of a pixel with the given numbers of settings: dark gray to black for households (darker for
more households, with the largest differences in dense areas), or blue for schools and green for
workplaces, if the pixel contains any and fewer than `SM_MANY_HOUSEHOLDS` households.
"""
function household_gradient(households, schoolclasses, offices)
    (schoolclasses > 0 || offices > 0) && households < SM_MANY_HOUSEHOLDS &&
        return FIG_PIXEL(sm_setting_color(schoolclasses)...)
    households == 0 && return FIG_WHITE
    gray = household_gray(households)
    return FIG_PIXEL(gray, gray, gray)
end

# opacity of the households above schools and workplaces, and number of households from which on it is
# reached
const SM_MAX_OPACITY = 0.7
const SM_OPAQUE_HOUSEHOLDS = 1000

"""
    household_opacity(households, schoolclasses, offices)

Alternative to `household_gradient`: black for households, or blue for schools and green for
workplaces below the households, whose opacity grows with their number (up to `SM_MAX_OPACITY` from
`SM_OPAQUE_HOUSEHOLDS` households).
"""
function household_opacity(households, schoolclasses, offices)
    if schoolclasses > 0 || offices > 0
        opacity = households == 0 ? 0.0 : SM_MAX_OPACITY * min(1.0, log(households) / log(SM_OPAQUE_HOUSEHOLDS))
        return FIG_PIXEL(((1 - opacity) .* sm_setting_color(schoolclasses))...)
    end
    return households == 0 ? FIG_WHITE : FIG_PIXEL(0, 0, 0)
end


# 3. FIGURE

# legend (inches): distance from the upper edge, height of a row, positions of the markers and labels,
# and width
const SM_LEGEND = (top = 0.05, row = 0.18, marker = 0.06, label = 0.16, width = 1.1)
const SM_LEGEND_ENTRIES = (("Households", FIG_COLORS.darkgray), ("Schools", FIG_COLORS.blue),
    ("Workplaces", FIG_COLORS.green))

# pixels of the legend
function settings_legend()
    l = SM_LEGEND
    height = l.row * length(SM_LEGEND_ENTRIES)
    p = plot(size = (100l.width, 100height), xlims = (0, l.width), ylims = (0, height), framestyle = :none,
        ticks = nothing, legend = false, grid = false, widen = false, margin = -2mm, fontfamily = FIG_FONT)
    for (k, (label, color)) in enumerate(SM_LEGEND_ENTRIES)
        y = height - (k - 0.5) * l.row
        scatter!(p, [l.marker], [y], markercolor = color, markerstrokewidth = 0, markersize = 3.5, label = "")
        annotate!(p, l.label, y, Plots.text(label, FIG_SMALL_FONTSIZE, FIG_FONT, :left, :vcenter))
    end
    return figure_image(p)
end

"""
    settings_map_figure(data; coloring = household_gradient)

Pixels of the settings map of `data` (`settings_per_pixel`) with its legend in the upper left corner.
`coloring` returns the color of a pixel with the given numbers of households, school classes, and
offices (`household_gradient` or `household_opacity`).
"""
function settings_map_figure(data; coloring = household_gradient)
    image = coloring.(data.households, data.schoolclasses, data.offices)
    legend = settings_legend()
    image[CartesianIndices(legend) .+ CartesianIndex(fig_px(SM_LEGEND.top), 0)] .= legend
    return image
end

"""
    settings_map(;
        folder = nothing,
        folder_tag = nothing,
        sim = nothing,
        data = nothing,
        coloring = household_gradient)

Draws the map of all households, schools, and workplaces (`settings_map_figure`) and saves the figure
(named by `figure_name`) and its data (`settings_map.csv.gz`) to `folder`. Pass `data` (from
`settings_per_pixel` or `load_settings_per_pixel`) to redraw the figure without loading the
population. Otherwise the data come from `sim` (default: the German population).

If `folder` is `nothing`, a new result folder is created, with `folder_tag` appended to its
timestamp-based name (see `result_folder`).
"""
function settings_map(;
        folder = nothing,
        folder_tag = nothing,
        sim = nothing,
        data = nothing,
        coloring = household_gradient)

    isnothing(folder) ? folder = result_folder(folder_tag) : mkpath(folder)

    if isnothing(data)
        sim = isnothing(sim) ? Simulation(population = "DE", seed = 0) : sim
        data = settings_per_pixel(sim)
    end
    save_settings_per_pixel(data, folder)
    return save_figure(settings_map_figure(data, coloring = coloring), :settings_map, folder)
end
