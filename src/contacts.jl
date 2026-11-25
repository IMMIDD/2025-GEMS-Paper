###
### CONTACT PATTERNS (FIGURE OF SECTION 2.4)
###

# Mean daily contacts per person between the age groups (ego and contact), sampled
# with GEMS. First the data (`contact_matrices`), then the layout and the drawing
# of the figure (`plot_contact_matrices`), then the analysis (`contacts`).

export ContactMatrixPanel
export contact_age_matrix
export contact_matrices
export contact_matrix_table, load_contact_matrices
export plot_contact_matrices
export contacts


# 1. DATA

# age groups of the figure: five years, the last one open (85+)
const CM_INTERVAL = 5
const CM_GROUPS = 18

cm_group(age) = min(Int(age) ÷ CM_INTERVAL, CM_GROUPS - 1) + 1
cm_lower_age(group) = (group - 1) * CM_INTERVAL

"""
    ContactMatrixPanel

One panel of the contact figure: the `section` (heading of the row), the
`title`, a `name` for the exported data, and the `matrix` of the mean daily
contacts of a person of age group i (row, ego) with persons of age group j
(column, contact).
"""
struct ContactMatrixPanel
    section::String
    title::String
    name::String
    matrix::Matrix{Float64}
end

# sections and titles of the panels by name
const CM_UNIFORM = "Uniform mixing"
const CM_MIXING = "Mixing assumptions in municipalities"
const CM_TITLES = Dict(
    "household" => (CM_UNIFORM, "Households"),
    "workplace" => (CM_UNIFORM, "Workplaces"),
    "education" => (CM_UNIFORM, "Educational settings"),
    "municipality_uniform" => (CM_MIXING, "Uniform mixing"),
    "municipality_polymod" => (CM_MIXING, "POLYMOD"),
    "municipality_strict" => (CM_MIXING, "Strict mixing"))

contact_panel(name, matrix) = ContactMatrixPanel(CM_TITLES[name]..., name, matrix)

# counts the contacts (ego × contact age group) that the members of `setting` draw in
# one time step with the contact sampling method of the setting; `stream` separates
# the random numbers of the setting types
function count_contacts!(counts, setting, sim, stream, present, sampled)
    empty!(present)
    present_individuals!(present, setting, sim)
    length(present) < 2 && return counts # no one to meet
    rng = Random.Xoshiro(1_000_000_000 * stream + Int(id(setting)))
    for k in eachindex(present)
        empty!(sampled)
        GEMS.sample_contacts!(sampled, contact_sampling_method(setting), setting, k, present,
            Int16(0), true, rng)
        ego = cm_group(age(present[k]))
        for contact in sampled
            counts[ego, cm_group(age(contact))] += 1
        end
    end
    return counts
end

"""
    contact_age_matrix(sim, settingtypes)

Mean daily contacts of a person (rows, ego) with persons of each age group
(columns), in age groups of `CM_INTERVAL` years, summed over `settingtypes`
(e.g. all levels of the workplace hierarchy).

Every member of every setting draws contacts once with the contact sampling
method of the setting in `sim`. The sums are divided by all persons of the ego's
age group, including those without a setting of these types. Each row therefore
sums to the mean daily contacts of a person of that age. The random numbers
depend only on the setting, not on the number of threads.
"""
function contact_age_matrix(sim, settingtypes)
    counts = zeros(Int, CM_GROUPS, CM_GROUPS)
    for (stream, T) in enumerate(settingtypes)
        all_settings = settings(sim, T)
        (isnothing(all_settings) || isempty(all_settings)) && continue
        chunks = Iterators.partition(eachindex(all_settings),
            cld(length(all_settings), 8 * Threads.nthreads()))
        tasks = map(collect(chunks)) do chunk
            Threads.@spawn begin
                c = zeros(Int, CM_GROUPS, CM_GROUPS)
                present, sampled = Individual[], Individual[]
                for k in chunk
                    count_contacts!(c, all_settings[k], sim, stream, present, sampled)
                end
                c
            end
        end
        for t in tasks
            counts .+= fetch(t)
        end
    end

    persons = zeros(Int, CM_GROUPS)
    foreach(i -> persons[cm_group(age(i))] += 1, individuals(sim))
    return counts ./ max.(persons, 1)
end

# true if all municipalities of `sim` use age-based contact sampling
has_age_based_municipalities(sim) =
    all(m -> contact_sampling_method(m) isa AgeBasedContactSampling, settings(sim, Municipality))

"""
    contact_matrices(sim)

Panels of the contact figure (`ContactMatrixPanel`):

- uniform mixing in households, workplaces, and educational settings, with the
  contact sampling of `sim` and summed over all levels of the setting hierarchy
- the three mixing assumptions in municipalities of the first wave analysis:
  uniform mixing, the age-based sampling of `sim` (POLYMOD in the reference
  model), and strict mixing (own age group only, see `use_municipality_mixing!`)

`sim` needs age-based sampling in the municipalities (reference model). Their
contact sampling is restored at the end.
"""
function contact_matrices(sim)
    has_age_based_municipalities(sim) ||
        throw(ArgumentError("The contact figure needs age-based contact sampling in the municipalities."))
    municipalities = settings(sim, Municipality)
    original = deepcopy(contact_sampling_method.(municipalities))

    panels = [
        contact_panel("household", contact_age_matrix(sim, [Household])),
        contact_panel("workplace", contact_age_matrix(sim, [Office, Department, Workplace, WorkplaceSite])),
        contact_panel("education", contact_age_matrix(sim, [SchoolClass, SchoolYear, School, SchoolComplex]))]

    municipality(name) = contact_panel(name, contact_age_matrix(sim, [Municipality]))
    try
        polymod = municipality("municipality_polymod")
        use_municipality_mixing!(sim, :uniform)
        push!(panels, municipality("municipality_uniform"), polymod)
        use_municipality_mixing!(sim, :same_age)
        push!(panels, municipality("municipality_strict"))
    finally
        contact_sampling_method!.(municipalities, original)
    end
    return panels
end

"""
    contact_matrix_table(panels)

The matrices of `contact_matrices` as a long table: one row per panel and pair
of age groups (lower age limits, the last group is open).
"""
contact_matrix_table(panels) = DataFrame([(
        panel = p.name,
        ego_age = cm_lower_age(i),
        contact_age = cm_lower_age(j),
        contacts_per_day = p.matrix[i, j])
    for p in panels for i in 1:CM_GROUPS for j in 1:CM_GROUPS])


const CM_DATA_FILE = "contact_matrices.csv"

"""
    load_contact_matrices(folder)

Reads the panels of `contact_matrices` from the table saved by `contacts`, to redraw the figure without
loading the population.
"""
function load_contact_matrices(folder)
    table = CSV.read(joinpath(folder, CM_DATA_FILE), DataFrame)
    return map(unique(table.panel)) do name
        matrix = zeros(CM_GROUPS, CM_GROUPS)
        for r in eachrow(table[table.panel .== name, :])
            matrix[cm_group(r.ego_age), cm_group(r.contact_age)] = r.contacts_per_day
        end
        contact_panel(name, matrix)
    end
end


# 2. LAYOUT

# The figure is drawn on one canvas in inches, so that all distances are fixed: size of a panel, gaps
# between the columns and rows, height of the section headings and panel titles, and room for the axis
# labels and the colorbar (together as wide as the text, `FIG_TEXT_WIDTH_IN`)
const CM_PANEL = 1.32
const CM_COLUMN_GAP = 0.1
const CM_ROW_GAP = 0.08
const CM_HEADING_HEIGHT = 0.28
const CM_TITLE_HEIGHT = 0.22
const CM_LEFT = 0.45
const CM_BOTTOM = 0.4
const CM_COLORBAR_GAP = 0.08
const CM_COLORBAR_WIDTH = 0.09
const CM_COLORBAR_LABELS = 0.42

# distances of the labels (inches): panel titles above a panel, length of the ticks, tick labels and axis
# labels from the axis, label of the colorbar
const CM_TITLE_OFFSET = 0.11
const CM_TICK = 0.03
const CM_MINOR_TICK = 0.02
const CM_TICKLABEL_OFFSET = 0.06
const CM_XTICKLABEL_OFFSET = 0.11
const CM_XLABEL_OFFSET = 0.28
const CM_YLABEL_OFFSET = 0.35
const CM_COLORBAR_LABEL_OFFSET = 0.34

# widths of the lines of the axes and of the grid
const CM_LINEWIDTH = 0.6
const CM_GRID_LINEWIDTH = 0.4

# logarithmic color scale: all values up to the lower limit get the lowest color;
# from a pale tint over the blue of the GEMS data to navy
const CM_LIMITS = (0.01, 5.0)
const CM_COLORBAR_TICKS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]
const CM_COLORS = cgrad([RGB(0.89, 0.93, 0.97), RGB(0.66, 0.82, 0.95), PV_MODEL_COLOR,
    RGB(0.0, 0.16, 0.34)], [0.0, 0.35, 0.68, 1.0])
const CM_BLACK = RGB(0, 0, 0)
const CM_WHITE = RGB(1, 1, 1)

# labeled ages (at the lower limit of their group), the open group in the middle of its cell
const CM_AGE_TICKS = 0:10:70

cm_log(v) = log10(clamp(v, CM_LIMITS...))
cm_position(age) = age / CM_INTERVAL * CM_PANEL / CM_GROUPS
cm_age_ticks() = [[(cm_position(a), string(a)) for a in CM_AGE_TICKS];
    (cm_position(cm_lower_age(CM_GROUPS) + CM_INTERVAL / 2), "$(cm_lower_age(CM_GROUPS))+")]


# 3. DRAWING

cm_text!(p, x, y, s, halign; size = FIG_SMALL_FONTSIZE, rotation = 0) =
    canvas_text!(p, 1, x, y, s, size = size, halign = halign, rotation = rotation)

cm_line!(p, x, y; color = CM_BLACK, linewidth = CM_LINEWIDTH) =
    plot!(p, x, y, color = color, linewidth = linewidth, label = "")

cm_frame!(p, x, y, w, h) = cm_line!(p, [x, x + w, x + w, x, x], [y, y, y + h, y + h, y])

# section heading between two lines (as in the population validation figure)
function cm_heading!(p, heading, x, y, width)
    half = 0.5 * heading_width_in(heading) + FIG_HEADING_GAP_IN
    center = x + width / 2
    canvas_text!(p, 1, center, y, heading, size = FIG_SECTION_FONTSIZE, font = FIG_FONT_BOLD)
    for side in ([x, center - half], [center + half, x + width])
        cm_line!(p, side, [y, y], color = FIG_LINE_COLOR)
    end
end

# heatmap of one panel with its letter and title; lower left corner at (x, y)
function cm_panel!(p, panel, letter, x, y)
    edges = CM_PANEL / CM_GROUPS .* (0:CM_GROUPS)
    heatmap!(p, x .+ edges, y .+ edges, permutedims(cm_log.(panel.matrix)),
        color = CM_COLORS, clims = cm_log.(CM_LIMITS), colorbar = false)
    for e in edges[2:end-1]
        cm_line!(p, [x + e, x + e], [y, y + CM_PANEL], color = CM_WHITE, linewidth = CM_GRID_LINEWIDTH)
        cm_line!(p, [x, x + CM_PANEL], [y + e, y + e], color = CM_WHITE, linewidth = CM_GRID_LINEWIDTH)
    end
    cm_frame!(p, x, y, CM_PANEL, CM_PANEL)
    canvas_title!(p, 1, x, y + CM_PANEL + CM_TITLE_OFFSET, letter, panel.title)
end

function cm_xaxis!(p, x, y)
    for (pos, label) in cm_age_ticks()
        cm_line!(p, [x + pos, x + pos], [y, y - CM_TICK])
        cm_text!(p, x + pos, y - CM_XTICKLABEL_OFFSET, label, :center)
    end
    cm_text!(p, x + CM_PANEL / 2, y - CM_XLABEL_OFFSET, "Ego age", :center, size = FIG_FONTSIZE)
end

function cm_yaxis!(p, x, y)
    for (pos, label) in cm_age_ticks()
        cm_line!(p, [x, x - CM_TICK], [y + pos, y + pos])
        cm_text!(p, x - CM_TICKLABEL_OFFSET, y + pos, label, :right)
    end
    cm_text!(p, x - CM_YLABEL_OFFSET, y + CM_PANEL / 2, "Contact age", :center, size = FIG_FONTSIZE,
        rotation = 90)
end

# vertical colorbar with labeled ticks and the minor ticks of a logarithmic axis
function cm_colorbar!(p, x, y, height)
    low, high = cm_log.(CM_LIMITS)
    position(v) = y + height * (log10(v) - low) / (high - low)
    edges = range(low, high, length = 257)
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    right = x + CM_COLORBAR_WIDTH
    heatmap!(p, [x, right], y .+ height .* (edges .- low) ./ (high - low),
        reshape(collect(centers), :, 1), color = CM_COLORS, clims = (low, high), colorbar = false)
    cm_frame!(p, x, y, CM_COLORBAR_WIDTH, height)
    for v in [k * 10.0^e for e in -3:1 for k in 1:9]
        CM_LIMITS[1] <= v <= CM_LIMITS[2] || continue
        cm_line!(p, [right, right + CM_MINOR_TICK], fill(position(v), 2), linewidth = CM_GRID_LINEWIDTH)
    end
    for v in CM_COLORBAR_TICKS
        label = isinteger(v) ? string(Int(v)) : string(v)
        cm_line!(p, [right, right + CM_TICK], fill(position(v), 2))
        cm_text!(p, right + CM_TICKLABEL_OFFSET, position(v),
            v == CM_LIMITS[1] ? "≤$label" : label, :left)
    end
    cm_text!(p, right + CM_COLORBAR_LABEL_OFFSET, y + height / 2,
        "Contacts per person and day (log scale)", :center, size = FIG_FONTSIZE, rotation = 90)
end

"""
    plot_contact_matrices(panels)

Contact figure of the paper: the matrices of `panels` (`contact_matrices`), one
row per section, with the age of the ego on the x-axis and the age of the contact
on the y-axis, and one logarithmic color scale for all panels.
"""
function plot_contact_matrices(panels)
    sections = unique(p.section for p in panels)
    rows = [[p for p in panels if p.section == s] for s in sections]
    columns = maximum(length, rows)
    panels_width = columns * CM_PANEL + (columns - 1) * CM_COLUMN_GAP
    row_height = CM_HEADING_HEIGHT + CM_TITLE_HEIGHT + CM_PANEL
    width = CM_LEFT + panels_width + CM_COLORBAR_GAP + CM_COLORBAR_WIDTH + CM_COLORBAR_LABELS
    height = length(rows) * row_height + (length(rows) - 1) * CM_ROW_GAP + CM_BOTTOM

    p = plot(
        size = (100width, 100height),
        xlims = (0, width),
        ylims = (0, height),
        framestyle = :none,
        ticks = false,
        grid = false,
        legend = false,
        colorbar = false,
        widen = false,
        margin = -2mm,
        fontfamily = FIG_FONT)

    letters = Iterators.Stateful('A':'Z')
    for (r, (section, row)) in enumerate(zip(sections, rows))
        top = height - (r - 1) * (row_height + CM_ROW_GAP)
        y = top - row_height
        cm_heading!(p, section, CM_LEFT, top - CM_HEADING_HEIGHT / 2, panels_width)
        for (c, panel) in enumerate(row)
            x = CM_LEFT + (c - 1) * (CM_PANEL + CM_COLUMN_GAP)
            cm_panel!(p, panel, popfirst!(letters), x, y)
            c == 1 && cm_yaxis!(p, x, y)
            r == length(rows) && cm_xaxis!(p, x, y)
        end
    end
    cm_colorbar!(p, CM_LEFT + panels_width + CM_COLORBAR_GAP, CM_BOTTOM,
        height - CM_BOTTOM - CM_HEADING_HEIGHT - CM_TITLE_HEIGHT)
    return p
end


# 4. ANALYSIS

"""
    contacts(;
        folder = nothing,
        folder_tag = nothing,
        panels = nothing,
        sim = isnothing(panels) ? DE_sim() : nothing)

Runs the contact analysis: draws the contact figure (`Fig_08_contact_matrices.png`) and
writes its data (`contact_matrices.csv`). By default with the reference model and
the population of Germany, as in the paper. Simulations without age-based contact
sampling in the municipalities are skipped with a warning. Pass `panels` (from
`load_contact_matrices`) to redraw the figure without loading the population.

If `folder` is `nothing`, a new result folder is created, with `folder_tag`
appended to its timestamp-based name (see `result_folder`).
"""
function contacts(;
        folder = nothing,
        folder_tag = nothing,
        panels = nothing,
        sim = isnothing(panels) ? DE_sim() : nothing)

    if isnothing(panels) && !has_age_based_municipalities(sim)
        @warn "The contact figure is skipped, as it needs age-based contact sampling in the municipalities (e.g. models/DE_model_first_wave.toml)."
        return
    end

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder(folder_tag) : nothing

    isnothing(panels) && (panels = contact_matrices(sim))
    CSV.write(joinpath(folder, CM_DATA_FILE), contact_matrix_table(panels))
    save_figure(plot_contact_matrices(panels), :contact_matrices, folder)
end
