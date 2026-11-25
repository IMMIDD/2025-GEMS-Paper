###
### FIGURE STYLE
###
# Sizes, fonts, colors, and file names shared by the figures of the paper, and the
# functions that title, combine, and save them. All sizes are print sizes: figures
# are drawn at 100 px per inch with font sizes in pt and saved at `FIG_DPI`.

export FIGURE_NUMBERS, figure_name
export figure_image, load_panel, scale_image
export title_strip, section_strip, combine_panels, save_figure


# 1. SIZES (PLOS)

# widths of the figures (the only two) and maximum height (inches), resolution
const FIG_FULL_WIDTH_IN = 7.5
const FIG_TEXT_WIDTH_IN = 5.2
const FIG_MAX_HEIGHT_IN = 8.75
const FIG_DPI = 300

# maximum file size (bytes)
const FIG_MAX_FILE_SIZE = 10 * 10^6

# pixels of a length in print
fig_px(inches) = round(Int, inches * FIG_DPI)


# 2. FONTS AND COLORS

const FIG_FONT = "Times Roman"
const FIG_FONT_BOLD = "Times Bold"

# section headings; panel titles; axis labels; tick labels, legends, and notes (pt)
const FIG_SECTION_FONTSIZE = 12
const FIG_HEADING_FONTSIZE = 10
const FIG_FONTSIZE = 9
const FIG_SMALL_FONTSIZE = 8

# thin lines next to section headings
const FIG_LINE_COLOR = "#9aa0a8"

# theme colors of the drawn figures (figures/drawn_figures.pptx): the first four are GEMS colors
const FIG_COLORS = (
    blue = "#009AFA",
    green = "#6C9F33",
    red = "#D62728",
    purple = "#9467BD",
    orange = "#FFA500",
    gray = "#9AA0A8",
    darkgray = "#444444",
    lightgray = "#D9D9D9")


# 3. FIGURE NUMBERS

"""
    FIGURE_NUMBERS

Numbers of the figures in the paper by name. Figures 13 and 14 are in the repository
of the RSV example. The drawn figures (1-3, 9, 10) carry their file names in the notes
of their slides.
"""
const FIGURE_NUMBERS = (
    setting_overview = 1,
    layered_architecture = 2,
    simulation_routine = 3,
    settings_map = 4,
    population_validation = 5,
    demographic_maps = 6,
    commuting = 7,
    contact_matrices = 8,
    disease_model = 9,
    trism_example = 10,
    first_wave = 11,
    first_wave_states = 12,
    household_isolation = 15,
    baseline_readouts = 16,
    parameter_sweep = 17,
    testing_scenarios = 18,
    benchmarks = 19)

"""
    figure_name(figure::Symbol; panel = "", name = figure)

File name (without extension) of a figure of the paper (`FIGURE_NUMBERS`), e.g.
`figure_name(:parameter_sweep)` is `"Fig_17_parameter_sweep"`, or of one of its panels, e.g.
`figure_name(:parameter_sweep, panel = "B", name = "sweep_heatmap")` is
`"Fig_17B_sweep_heatmap"`.
"""
function figure_name(figure::Symbol; panel = "", name = figure)
    haskey(FIGURE_NUMBERS, figure) ||
        throw(ArgumentError("Unknown figure $figure. Add it to FIGURE_NUMBERS."))
    return "Fig_$(lpad(FIGURE_NUMBERS[figure], 2, '0'))$(panel)_$name"
end


# 4. TEXT

# estimated character widths of Times, as shares of the font size (digits and capitals are wider
# than the average lowercase letter, bold text by about 10 %)
const FIG_CHAR_WIDTH = 0.42
char_width(c) = c == '%' ? 0.83 : isdigit(c) ? 0.5 : c in "MW" ? 0.89 : isuppercase(c) ? 0.66 : c == ' ' ? 0.25 :
    c in ".,:;()|/–-" ? 0.33 : c in "ijlftr" ? 0.3 : c in "mw" ? 0.72 : FIG_CHAR_WIDTH

text_width_in(s, size; bold = false) = (bold ? 1.1 : 1.0) * sum(char_width, s; init = 0.0) * size / 72
line_height_in(size) = 1.2 * size / 72

"""
    title_text!(p, x, y, title, width_in; size = FIG_HEADING_FONTSIZE, valign = :bottom)

Draws a panel title at the relative position (`x`, `y`) of `p`, left
aligned: its letter ("A)") in bold, as in the other figures of the paper,
followed by the rest in regular weight. `width_in` is the width of the plot
area, to place the rest after the letter.
"""
function title_text!(p, x, y, title, width_in; size = FIG_HEADING_FONTSIZE, valign = :bottom)
    m = match(r"^([A-Z]\))\s*(.*)$", title)
    letter, rest = isnothing(m) ? ("", title) : (m.captures[1], m.captures[2])
    dx = 0.0
    if !isempty(letter)
        annotate!(p, (x, y), Plots.text(letter, size, FIG_FONT_BOLD, :left, valign))
        dx = (text_width_in(letter, size, bold = true) + 0.3size / 72) / width_in
    end
    return annotate!(p, (x + dx, y), Plots.text(rest, size, FIG_FONT, :left, valign))
end

"""
    section_header(txt; width_in = FIG_FULL_WIDTH_IN - 4 / 25.4, gap_in = 0.12)

Empty subplot carrying the heading of a section of the figure: centered as
the section titles of the other figures of the paper (e.g.
`plot_paper_figure_1`), bold and larger than the panel titles, between two
thin lines that set it apart from them. `width_in` is the width of the
heading's plot area, `gap_in` the distance of the lines from the (estimated) text.
"""
function section_header(txt; width_in = FIG_FULL_WIDTH_IN - 4 / 25.4, gap_in = 0.12)
    p = plot(framestyle = :none, ticks = false, legend = false, grid = false,
        xlims = (0, 1), ylims = (0, 1), margin = 0mm, top_margin = -2mm, bottom_margin = -2mm,
        fontfamily = FIG_FONT)
    half = (0.5 * text_width_in(txt, FIG_SECTION_FONTSIZE, bold = true) + gap_in) / width_in
    plot!(p, [0, 0.5 - half, NaN, 0.5 + half, 1], fill(0.5, 5), linecolor = FIG_LINE_COLOR, linewidth = 0.6)
    return annotate!(p, 0.5, 0.5, Plots.text(txt, FIG_SECTION_FONTSIZE, FIG_FONT_BOLD, :center, :vcenter))
end


# 5. CANVAS
# Texts, legends, and colorbars around the plots of a figure are drawn on a transparent canvas in
# inches above them.

"""
    plot_corners(p, k; aspect = nothing)

Position (left, right, bottom, top; inches from the lower left corner) of the plot area of the
subplot `k` of the figure `p`, which must be laid out (`Plots.prepare_output`). `aspect` is its width
per height, if it is fixed (then it is centered in the area that Plots reserves).
"""
function plot_corners(p, k; aspect = nothing)
    inches(x) = x.value / 25.4
    area = p[k].plotarea
    w, h = inches(Plots.width(area)), inches(Plots.height(area))
    x, y = inches(Plots.left(area)) + w / 2, p.attr[:size][2] / 100 - inches(Plots.top(area)) - h / 2
    width = isnothing(aspect) ? w : min(w, h * aspect)
    height = isnothing(aspect) ? h : width / aspect
    return (left = x - width / 2, right = x + width / 2, bottom = y - height / 2, top = y + height / 2)
end

"""
    add_canvas!(p)

Adds a transparent canvas above the figure `p` and returns its subplot number. Its coordinates are
inches from the lower left corner.
"""
function add_canvas!(p)
    canvas = length(p.subplots) + 1
    width, height = p.attr[:size] ./ 100
    plot!(p, inset = bbox(0, 0, 1, 1), subplot = canvas, framestyle = :none, ticks = nothing, legend = false,
        grid = false, widen = false, margin = -2mm, background_color_subplot = :transparent,
        background_color_inside = :transparent, xlims = (0, width), ylims = (0, height))
    return canvas
end

"""
    canvas_text!(p, canvas, x, y, txt; size = FIG_FONTSIZE, font = FIG_FONT, halign = :center, rotation = 0)

Draws `txt` at the position (`x`, `y`) of the `canvas` of the figure `p`, vertically centered.
"""
canvas_text!(p, canvas, x, y, txt; size = FIG_FONTSIZE, font = FIG_FONT, halign = :center, rotation = 0) =
    annotate!(p, [(x, y, Plots.text(txt, size, font, halign, :vcenter, rotation = rotation))], subplot = canvas)

"""
    canvas_title!(p, canvas, x, y, letter, title)

Draws a panel title (bold letter and title) on the `canvas` of the figure `p`, starting at (`x`, `y`).
"""
function canvas_title!(p, canvas, x, y, letter, title)
    canvas_text!(p, canvas, x, y, "$letter)", size = FIG_HEADING_FONTSIZE, font = FIG_FONT_BOLD, halign = :left)
    x += text_width_in("$letter)", FIG_HEADING_FONTSIZE, bold = true) + 0.3FIG_HEADING_FONTSIZE / 72
    return canvas_text!(p, canvas, x, y, title, size = FIG_HEADING_FONTSIZE, halign = :left)
end


# sizes of a legend on the canvas (inches): distance of its frame from the axes and from its entries,
# length of the keys, and distance between keys and labels
const FIG_LEGEND_GAP_IN = 0.02
const FIG_LEGEND_PADDING_IN = 0.05
const FIG_LEGEND_KEY_IN = 0.16
const FIG_LEGEND_LABEL_GAP_IN = 0.05

"""
    canvas_legend!(p, canvas, k, labels)

Draws a legend of the series in the subplot `k` of the figure `p` on the `canvas`, in the upper right
corner of the plot area. `labels` are pairs of the labels of the series and those in the legend.
"""
function canvas_legend!(p, canvas, k, labels)
    area = plot_corners(p, k)
    row = line_height_in(FIG_SMALL_FONTSIZE)
    w = FIG_LEGEND_KEY_IN + FIG_LEGEND_LABEL_GAP_IN + 2FIG_LEGEND_PADDING_IN +
        1.05 * maximum(text_width_in.(last.(labels), FIG_SMALL_FONTSIZE))
    h = row * length(labels) + 2FIG_LEGEND_PADDING_IN
    x, top = area.right - FIG_LEGEND_GAP_IN - w, area.top - FIG_LEGEND_GAP_IN
    plot!(p, x .+ [0, w, w, 0], top .- [0, 0, h, h], subplot = canvas, label = "", seriestype = :shape,
        fillcolor = :white, linecolor = FIG_COLORS.lightgray, linewidth = 1)
    for (i, (name, label)) in enumerate(labels)
        s = first(s for s in p[k].series_list if s[:label] == name)
        y = top - FIG_LEGEND_PADDING_IN - (i - 0.5) * row
        key = x .+ FIG_LEGEND_PADDING_IN .+ [0, FIG_LEGEND_KEY_IN]
        if s[:seriestype] == :shape || !isnothing(s[:fillrange])
            plot!(p, key[[1, 2, 2, 1]], y .+ 0.3row .* [-1, -1, 1, 1], subplot = canvas, label = "",
                seriestype = :shape, fillcolor = s[:fillcolor], linewidth = 0, linecolor = s[:fillcolor])
        else
            plot!(p, key, [y, y], subplot = canvas, label = "", linecolor = s[:linecolor],
                linestyle = s[:linestyle], linewidth = 1.5)
        end
        canvas_text!(p, canvas, key[2] + FIG_LEGEND_LABEL_GAP_IN, y, label, size = FIG_SMALL_FONTSIZE, halign = :left)
    end
    return p
end


# 6. IMAGES
# Panels from different sources (plots, drawn panels, screenshots) are combined as pixels, so
# that none of them is resampled.

const FIG_PIXEL = Plots.Colors.RGB{Plots.Colors.N0f8}
const FIG_WHITE = FIG_PIXEL(1, 1, 1)

# pixel without transparency (on white)
flatten_pixel(c) = FIG_PIXEL(Plots.Colors.color(c))
function flatten_pixel(c::Plots.Colors.TransparentColor)
    a = Plots.Colors.alpha(c)
    rgb = Plots.Colors.RGB{Float64}(Plots.Colors.color(c))
    on_white(x) = clamp(a * x + (1 - a), 0, 1)
    return FIG_PIXEL(on_white(rgb.r), on_white(rgb.g), on_white(rgb.b))
end

white_image(height, width) = fill(FIG_WHITE, height, width)

# drawn panels and screenshots that are combined with plots
const FIG_PANEL_FOLDER = joinpath(dirname(@__DIR__), "figures", "panels")

"""
    load_panel(file)

Pixels of a PNG file without transparency (on white).
"""
function load_panel(file)
    isfile(file) || throw(ArgumentError("The panel $file does not exist."))
    return flatten_pixel.(PNGFiles.load(file))
end

# share of each of `n` pixels in each pixel of the same row scaled down to `m` pixels
function scaling_weights(n, m)
    overlap(i, j) = max(0.0, min(i * m, j * n) - max((i - 1) * m, (j - 1) * n)) / n
    return [overlap(i, j) for j in 1:m, i in 1:n]
end

"""
    scale_image(image, width)

Scales `image` down to `width` pixels, keeping its aspect ratio (means of the covered pixels).
"""
function scale_image(image, width)
    h, w = size(image)
    width <= w || throw(ArgumentError("The image has only $w pixels, it cannot be scaled to $width."))
    rows, columns = scaling_weights(h, round(Int, h * width / w)), scaling_weights(w, width)
    scaled(channel) = clamp.(rows * Float64.(channel.(image)) * columns', 0, 1)
    return FIG_PIXEL.(scaled(Plots.Colors.red), scaled(Plots.Colors.green), scaled(Plots.Colors.blue))
end

"""
    figure_image(p)

Pixels of the plot `p`, drawn at `FIG_DPI` (plots have 100 px per inch). GR rounds sizes to even
numbers of pixels; differences of up to two pixels are evened out with white (or cut off).
"""
function figure_image(p)
    width, height = fig_px.(p.attr[:size] ./ 100)
    p.attr[:dpi] = FIG_DPI
    image = mktempdir() do dir
        file = joinpath(dir, "plot.png")
        png(p, file)
        load_panel(file)
    end
    all(abs.(size(image) .- (height, width)) .<= 2) ||
        throw(ErrorException("The plot has $(reverse(size(image))) instead of ($width, $height) pixels."))
    out = white_image(height, width)
    h, w = min.(size(image), (height, width))
    out[1:h, 1:w] .= image[1:h, 1:w]
    return out
end

"""
    title_strip(title, width; height_in = 0.25, offset_in = 0.0)

Strip of `width` pixels with a panel title (e.g. `"A) Strategy"`, see `title_text!`) that starts
`offset_in` inches from its left edge.
"""
function title_strip(title, width; height_in = 0.25, offset_in = 0.0)
    width_in = width / FIG_DPI
    # Plots adds 2 mm around a plot area
    strip = plot(framestyle = :none, ticks = false, legend = false, grid = false,
        xlims = (0, 1), ylims = (0, 1), widen = false,
        left_margin = -2mm, right_margin = -2mm, top_margin = -2mm, bottom_margin = -2mm,
        size = (100width_in, 100height_in), fontfamily = FIG_FONT)
    title_text!(strip, offset_in / width_in, 0.5, title, width_in, valign = :vcenter)
    return figure_image(strip)
end

"""
    section_strip(txt, width; height_in = 0.3)

Strip of `width` pixels with the heading of a section (see `section_header`).
"""
function section_strip(txt, width; height_in = 0.3)
    width_in = width / FIG_DPI
    strip = section_header(txt, width_in = width_in, gap_in = 0.25)
    plot!(strip, left_margin = -2mm, right_margin = -2mm, size = (100width_in, 100height_in))
    return figure_image(strip)
end

"""
    combine_panels(images; direction = :horizontal, gap_in = 0.0, align = :start)

Places `images` next to each other (`:horizontal`) or below each other (`:vertical`) on white,
without resampling them. `gap_in` is the space between them (inches), and `align` the position of
smaller images across this direction (`:start`, `:center`, or `:end`).
"""
function combine_panels(images; direction = :horizontal, gap_in = 0.0, align = :start)
    direction in (:horizontal, :vertical) || throw(ArgumentError("Unknown direction $direction."))
    align in (:start, :center, :end) || throw(ArgumentError("Unknown alignment $align."))
    along, across = direction == :horizontal ? (2, 1) : (1, 2)
    gap = fig_px(gap_in)
    total = sum(size.(images, along)) + gap * (length(images) - 1)
    extent = maximum(size.(images, across))
    out = direction == :horizontal ? white_image(extent, total) : white_image(total, extent)

    pos = 0
    for image in images
        free = extent - size(image, across)
        offset = align == :start ? 0 : align == :center ? free ÷ 2 : free
        along_range = pos .+ (1:size(image, along))
        across_range = offset .+ (1:size(image, across))
        if direction == :horizontal
            out[across_range, along_range] .= image
        else
            out[along_range, across_range] .= image
        end
        pos += size(image, along) + gap
    end
    return out
end

"""
    save_figure(figure, name::Symbol, folder; panel = "", panel_name = name)

Saves `figure` (a plot or pixels) as PNG at `FIG_DPI` in `folder`, named by `figure_name`, and returns
its path. Figures must be as wide as the page or the text (`FIG_FULL_WIDTH_IN`, `FIG_TEXT_WIDTH_IN`)
and not higher than `FIG_MAX_HEIGHT_IN`; single panels (`panel`, e.g. `"B"`) may have any width.
"""
function save_figure(figure, name::Symbol, folder; panel = "", panel_name = name)
    image = figure isa AbstractMatrix ? figure : figure_image(figure)
    height, width = size(image)
    widths = fig_px.((FIG_FULL_WIDTH_IN, FIG_TEXT_WIDTH_IN))
    isempty(panel) && !(width in widths) &&
        throw(ArgumentError("Figures must be $(join(widths, " or ")) px wide. Got $width px."))
    (width <= maximum(widths) && height <= fig_px(FIG_MAX_HEIGHT_IN)) ||
        throw(ArgumentError("The figure is too large: $width x $height px."))

    file = joinpath(folder, figure_name(name, panel = panel, name = panel_name) * ".png")
    PNGFiles.save(file, image, dpi = FIG_DPI)
    filesize(file) <= FIG_MAX_FILE_SIZE || @warn "$file is larger than 10 MB."
    return file
end
