###
### SETTING SIZES: GEMS SETTINGS VS. OFFICIAL STATISTICS
###
# Sizes of settings in the GEMS population of Germany and in official
# statistics, from the perspective of the persons in them: the share of
# workers and students in workplaces and universities of a size class, and
# the mean size of the groups (school classes, daycare groups) in which
# contacts take place. Built like the population validation
# (`population_validation.jl`, whose readers, style and panels are reused):
#   1. reading the reference tables in `inputs/population/` (the KMK
#      workbook is downloaded into `data/`, see `obtain_kmk_classes`)
#   2. measuring the same sizes in the simulation's population
#   3. joining both (saved and loaded with the population validation data)
#   4. drawing the panels (part of the figure of the population validation)

export setting_size_data, setting_size_panels
export SizeClass, size_class, WORKPLACE_SIZE_CLASSES, UNIVERSITY_SIZE_CLASSES
export GroupType, GROUP_TYPES
export official_establishments, official_university_sizes, official_group_sizes, kmk_pupils_and_classes
export model_workplace_sizes, model_university_sizes, model_group_sizes


# SIZE CLASSES AND GROUP TYPES

"""
    SizeClass(label, min, max, reference_label)

A size class from `min` to `max` (both inclusive; `typemax(Int)` for an
open-ended class), its `label` in figures, and the label of the
corresponding row or column in the reference table (`""` if the reference
is classified by its numeric sizes).
"""
struct SizeClass
    label::String
    min::Int
    max::Int
    reference_label::String
end

SizeClass(label, min, max) = SizeClass(label, min, max, "")

"""
    size_class(value, classes)

Label of the class in `classes` that contains `value`. Each value must
fall into exactly one class.
"""
function size_class(value, classes)
    hits = filter(c -> c.min <= value <= c.max, classes)
    length(hits) == 1 || error("Size $value falls into $(length(hits)) classes.")
    return only(hits).label
end

"""
    WORKPLACE_SIZE_CLASSES

Size classes of establishments by number of employees. Each class is a
union of the size classes of the employment statistics of the BA
(`official_establishments`); GEMS workplaces are classified by their number
of members.
"""
const WORKPLACE_SIZE_CLASSES = [
    SizeClass("1–9", 1, 9),
    SizeClass("10–49", 10, 49),
    SizeClass("50–249", 50, 249),
    SizeClass("250–499", 250, 499),
    SizeClass("500+", 500, typemax(Int)),
]

"""
    UNIVERSITY_SIZE_CLASSES

Size classes of universities by number of students. The reference
(21311-0002) lists students per university, so both sides are classified
by their numeric sizes.
"""
const UNIVERSITY_SIZE_CLASSES = [
    SizeClass("<0.5k", 0, 499),
    SizeClass("0.5–2k", 500, 1_999),
    SizeClass("2–5k", 2_000, 4_999),
    SizeClass("5–10k", 5_000, 9_999),
    SizeClass("10–20k", 10_000, 19_999),
    SizeClass("20k+", 20_000, typemax(Int)),
]

"""
    GroupType(key, label, setting, median_ages, reference, kmk_sheet)

A type of group in which children and pupils have their closest contacts:
its `key`, its `label` in figures, the education setting of GEMS whose
classes form these groups (`EDUCATION_SETTINGS`), the range of the median
age of the members of a GEMS class to count (`(min, max)`, `min` inclusive,
`max` exclusive), and the source of the reference: `:daycare` (daycare
statistics, `official_daycare_groups`), `:kmk` (the KMK workbook sheet
`kmk_sheet`, `kmk_pupils_and_classes`) or `:none` (no official group sizes).

The groups are compared by their mean size (members per group), the only
size measure that the official statistics publish for them.
"""
struct GroupType
    key::String
    label::String
    setting::String
    median_ages::Tuple{Float64, Float64}
    reference::Symbol
    kmk_sheet::String
end

"""
    GROUP_TYPES

Group types of `official_group_sizes` and `model_group_sizes`. Together
they cover all classes of GEMS except those of universities, which are
compared by university size instead (`UNIVERSITY_SIZE_CLASSES`):
- daycare: groups in daycare facilities (daycare statistics, 22541-0001)
- primary: classes of the primary level (KMK, Primarstufe)
- secondary (5–10): classes of grades 5 to 10 (KMK, Sekundarstufe I). GEMS
  has no school types: all secondary and vocational schools have the same
  class size. Gesyland forms the classes per grade (grade = age − 6), so
  they are split by the median age of their members: up to 16 (grade 10)
- secondary (11+): the other secondary classes (median age 17 or older):
  upper secondary (Gymnasium) and vocational classes, which GEMS does not
  distinguish (classes with a median age of 20 or more are vocational).
  Official statistics publish no class sizes for them (they are taught in
  courses), so they have no reference.
- special: classes of special schools (KMK, Förderschulen)
"""
const GROUP_TYPES = [
    GroupType("daycare", "Daycare", "daycare", (-Inf, Inf), :daycare, ""),
    GroupType("primary", "Primary", "primary", (-Inf, Inf), :kmk, "A.I.1.2a"),
    GroupType("secondary_5_10", "Sec. 5–10", "secondary", (-Inf, 17.0), :kmk, "A.I.1.3a"),
    GroupType("secondary_11", "Sec. 11+", "secondary", (17.0, Inf), :none, ""),
    GroupType("special", "Special", "special", (-Inf, Inf), :kmk, "A.I.1.5a"),
]


# 1. REFERENCE DATA

"""
    official_establishments(file = joinpath(REFERENCE_DIR, "ba_establishments_by_size_2019.csv"))

Establishments in Germany and their employees subject to social insurance
per size class (`WORKPLACE_SIZE_CLASSES`) on 30.06.2019, from the
employment statistics of the Bundesagentur für Arbeit (BA): category,
units (establishments), persons (employees). `file` is an extract of
sheet 2.3 of the BA table "Betriebe und sozialversicherungspflichtig
Beschäftigte" (see its header and `inputs/README.md`).

An establishment (Betrieb) is a local unit: a company with branches in
several municipalities consists of several establishments. The BA counts
establishments with at least one employee subject to social insurance and
classifies them by the number of these employees. GEMS workplaces include
all employed persons, so they are slightly larger by this measure.
"""
function official_establishments(file = joinpath(REFERENCE_DIR, "ba_establishments_by_size_2019.csv"))
    df = CSV.read(file, DataFrame, comment = "#")
    df.max = coalesce.(df.max, typemax(Int))
    within(c) = (df.min .>= c.min) .& (df.max .<= c.max)
    all(i -> count(c -> within(c)[i], WORKPLACE_SIZE_CLASSES) == 1, 1:nrow(df)) ||
        error("The size classes of $file do not fall into exactly one of WORKPLACE_SIZE_CLASSES each.")
    return DataFrame(
        category = getfield.(WORKPLACE_SIZE_CLASSES, :label),
        units = [sum(df.establishments[within(c)]) for c in WORKPLACE_SIZE_CLASSES],
        persons = [sum(df.employees[within(c)]) for c in WORKPLACE_SIZE_CLASSES])
end

"""
    official_university_sizes(file = joinpath(REFERENCE_DIR, "21311-0002_de.csv"); semester = "WS 2019/20")

Students per university (Hochschule) in Germany in `semester`, all
nationalities and both sexes, from the student statistics table
21311-0002: code, name, students. Universities without students (closed
or merged) are left out; locations of a university that the table lists
separately (e.g. "Universität Heidelberg in Mannheim") count as
universities of their own.
"""
function official_university_sizes(file = joinpath(REFERENCE_DIR, "21311-0002_de.csv"); semester = "WS 2019/20")
    df = read_genesis_csv(file)
    df = df[(df.row1 .== semester) .& (df.col1 .== "Insgesamt") .& (df.col2 .== "Insgesamt"), :]
    isempty(df) && error("No data for $semester in $file.")
    df = df[coalesce.(df.value, 0.0) .> 0, :]
    return DataFrame(code = df.row2, name = df.row3, students = Int.(df.value))
end

"""
    KMK_CLASSES_FILE

Location of the KMK workbook of pupils and classes 2010–2019, Teil A
(`SKL_Teil_A_Schueler_2019.xlsx`). The KMK terms of use do not permit
redistribution, so the workbook is not part of the repository: it is
downloaded from `KMK_CLASSES_URL` if it does not exist.
"""
const KMK_CLASSES_FILE = joinpath("data", "SKL_Teil_A_Schueler_2019.xlsx")

"""
    KMK_CLASSES_URL

KMK, Dokumentation 227: Schüler/-innen, Klassen, Lehrkräfte und
Absolventen/-innen der Schulen 2010–2019, Teil A (allgemeinbildende und
berufliche Schulen). © Sekretariat der Kultusministerkonferenz.
"""
const KMK_CLASSES_URL = "https://www.kmk.org/fileadmin/Dateien/pdf/Statistik/Dokumentationen/" *
    "SKL_Teil_A_Schueler_2019.xlsx"

"""
    KMK_CLASSES_SHA256

SHA-256 checksum of the file at `KMK_CLASSES_URL` (retrieved 21.11.2026).
"""
const KMK_CLASSES_SHA256 = "10d9fa5d4908587c08948d28d063aa476cf5ce8d00ad999f5ba7eacc4bf4236f"

"""
    obtain_kmk_classes(file = KMK_CLASSES_FILE)

Downloads the KMK workbook from `KMK_CLASSES_URL` if `file` does not exist,
and checks that it is the version used in the paper.
"""
function obtain_kmk_classes(file = KMK_CLASSES_FILE)
    if !isfile(file)
        printinfo("Downloading the KMK workbook of pupils and classes from $KMK_CLASSES_URL")
        mkpath(dirname(file))
        Downloads.download(KMK_CLASSES_URL, file)
    end
    checksum = bytes2hex(open(SHA.sha256, file))
    if checksum != KMK_CLASSES_SHA256
        @warn "$file is not the version of the KMK workbook used in the paper ($KMK_CLASSES_URL)."
    end
    return file
end

"""
    kmk_pupils_and_classes(sheet; file = obtain_kmk_classes(), year = 2019)

Pupils and classes in Germany ("BG") in `year` from a sheet of the KMK
workbook "Schüler/-innen, Klassen, Lehrkräfte und Absolvierende der
Schulen", Teil A (`KMK_CLASSES_URL`; e.g. `"A.I.1.2a"` for the primary
level). Each sheet holds blocks with one row per federal state and one
column per year, titled above the row "Land" (e.g. "Schüler/-innen",
"Klassen").
"""
function kmk_pupils_and_classes(sheet; file = obtain_kmk_classes(), year = 2019)
    d = XLSX.readxlsx(file)[sheet][:]
    text(v) = ismissing(v) ? "" : strip(string(v))
    blocks = Dict{String, Float64}()
    for r in 1:size(d, 1)
        text(d[r, 1]) == "Land" || continue
        title = last(filter(!isempty, text.(d[max(1, r - 6):(r - 1), 2])))
        col = findfirst(v -> text(v) == string(year), d[r, :])
        total = findnext(i -> text(d[i, 1]) == "BG", 1:size(d, 1), r)
        (isnothing(col) || isnothing(total)) && error("No total for $year in block \"$title\" of sheet $sheet.")
        blocks[title] = Float64(d[total, col])
    end
    all(haskey.(Ref(blocks), ("Schüler/-innen", "Klassen"))) ||
        error("Sheet $sheet has no blocks of pupils and classes (found: $(collect(keys(blocks)))).")
    return (pupils = blocks["Schüler/-innen"], classes = blocks["Klassen"])
end

"""
    official_daycare_groups(file = joinpath(REFERENCE_DIR, "22541-0001_de.csv"); date = "01.03.2019")

Children and groups in daycare facilities in Germany on `date`, from the
daycare statistics table 22541-0001. The children include those in
facilities without a group structure (about 13 % of the facilities) and
school children in after-school care, so their mean per group is an upper
bound for the groups of children before school age.
"""
function official_daycare_groups(file = joinpath(REFERENCE_DIR, "22541-0001_de.csv"); date = "01.03.2019")
    df = read_genesis_csv(file)
    df = df[df.row1 .== date, :]
    value(m) = only(df.value[df.col1 .== m])
    return (children = value("Kinder in Tageseinrichtungen"), groups = value("Gruppen in Kindertageseinrichtungen"))
end

"""
    official_group_sizes()

Members and groups per group type (`GROUP_TYPES`) in Germany: children and
groups in daycare (`official_daycare_groups`), pupils and classes of the
school levels (`kmk_pupils_and_classes`): category, members, groups
(`missing` for group types without a reference).
"""
function official_group_sizes()
    rows = map(GROUP_TYPES) do g
        if g.reference == :daycare
            d = official_daycare_groups()
            (category = g.key, members = d.children, groups = d.groups)
        elseif g.reference == :kmk
            k = kmk_pupils_and_classes(g.kmk_sheet)
            (category = g.key, members = k.pupils, groups = k.classes)
        else
            (category = g.key, members = missing, groups = missing)
        end
    end
    return DataFrame(rows)
end


# 2. MODEL DATA

"""
    model_workplace_sizes(sim)

Number of members of every workplace (`Workplace` setting, the sum of its
offices) in the simulation's population.

GEMS workplaces contain departments, which contain offices; members are
assigned to offices. In the GEMS population of Germany (v1.0.2), every
workplace site contains exactly one workplace, so sites and workplaces
have the same sizes.
"""
function model_workplace_sizes(sim)
    offices, departments = settings(sim, Office), settings(sim, Department)
    department_size = zeros(Int, length(departments))
    foreach(o -> department_size[contained(o)] += size(o), offices)
    workplace_size = zeros(Int, length(settings(sim, Workplace)))
    foreach(d -> workplace_size[contained(d)] += department_size[id(d)], departments)
    return workplace_size
end

"""
    model_university_sizes(institutions)

Number of students of every university (institution classified as
university) in the simulation's population; `institutions` from
`education_institutions(sim)`.
"""
function model_university_sizes(institutions)
    (; classes, institution, type) = institutions
    students = Dict{Int, Int}()
    for (k, c) in enumerate(classes)
        type[institution[k]] == "university" || continue
        students[institution[k]] = get(students, institution[k], 0) + size(c)
    end
    return collect(values(students))
end

"""
    model_group_sizes(institutions)

Members and groups per group type (`GROUP_TYPES`) in the simulation's
population: the classes of the institutions of the type's education
setting whose members' median age lies in the type's range; `institutions`
from `education_institutions(sim)`: category, members, groups.
"""
function model_group_sizes(institutions)
    (; classes, institution, type) = institutions
    median_age = [size(c) > 0 ? median(age.(individuals(c))) : NaN for c in classes]
    rows = map(GROUP_TYPES) do g
        lo, hi = g.median_ages
        selected = [k for k in eachindex(classes) if type[institution[k]] == g.setting &&
            size(classes[k]) > 0 && lo <= median_age[k] < hi]
        (category = g.key, members = sum(k -> size(classes[k]), selected; init = 0), groups = length(selected))
    end
    return DataFrame(rows)
end


# 3. DATA

"""
    class_totals(sizes, classes)

Number of units and their members per size class, in the order of
`classes`: category, units, persons.
"""
function class_totals(sizes, classes)
    labels = size_class.(sizes, Ref(classes))
    return DataFrame(
        category = getfield.(classes, :label),
        units = [count(==(c.label), labels) for c in classes],
        persons = [sum(sizes[labels .== c.label]; init = 0) for c in classes])
end

# model and reference side by side (columns prefixed with `model_` and `reference_`)
function side_by_side(model, reference)
    prefixed(df, p) = rename(df, [n => "$(p)_$n" for n in names(df) if n != "category"])
    return innerjoin(prefixed(model, "model"), prefixed(reference, "reference"), on = :category)
end

"""
    setting_size_data(sim)

Model and reference sizes of settings in Germany, as a named tuple of tables:
- `workplaces`: per size class (`WORKPLACE_SIZE_CLASSES`) the number of
  workplaces (GEMS) and establishments (reference), and their members or
  employees subject to social insurance (BA): category, model_units, model_persons,
  reference_units, reference_persons
- `universities`: per size class (`UNIVERSITY_SIZE_CLASSES`) the number of
  universities and their students, same columns
- `groups`: per group type (`GROUP_TYPES`) the members and groups:
  category, model_members, model_groups, reference_members, reference_groups
"""
function setting_size_data(sim)
    institutions = education_institutions(sim)
    groups = model_group_sizes(institutions)
    students = model_university_sizes(institutions)

    # every member of a school setting is in a group type or at a university
    members = sum(size, institutions.classes)
    covered = sum(groups.members) + sum(students; init = 0)
    covered == members || error("Group types and universities cover $covered of $members members of school settings.")

    return (
        workplaces = side_by_side(class_totals(filter(>(0), model_workplace_sizes(sim)), WORKPLACE_SIZE_CLASSES),
            official_establishments()),
        universities = side_by_side(class_totals(students, UNIVERSITY_SIZE_CLASSES),
            class_totals(official_university_sizes().students, UNIVERSITY_SIZE_CLASSES)),
        groups = side_by_side(groups, official_group_sizes()))
end

# 4. PANELS

# share of the persons per size class
person_shares(df, side) = df[!, "$(side)_persons"] ./ sum(df[!, "$(side)_persons"])

# mean size of the groups per group type
mean_group_size(df, side) = df[!, "$(side)_members"] ./ df[!, "$(side)_groups"]

"""
    setting_size_panels(data; letters = ('E', 'F', 'G'))

Setting sizes in GEMS and the reference (the tables of `setting_size_data`,
as part of `validation_data`), from the perspective of their members: the
share of workers per workplace size class, the mean size of daycare
groups and school classes, and the share of students per university size
class, titled with `letters`.
"""
function setting_size_panels(data; letters = ('E', 'F', 'G'))
    w, g, u = data.workplaces, data.groups, data.universities
    group_labels = [only(t.label for t in GROUP_TYPES if t.key == k) for k in g.category]
    return [
        plot_paired_bars(person_shares(w, "reference"), person_shares(w, "model"), w.category,
            title = "$(letters[1])) Workplaces", xlabel = "Employees", ylabel = "Share of workers",
            reference_label = "BA", xrotation = 30),
        plot_paired_bars(mean_group_size(g, "reference"), mean_group_size(g, "model"), group_labels,
            title = "$(letters[2])) Mean group size", ylabel = "Members per group",
            reference_label = "Destatis, KMK", shares = false, ytick_step = 10, xrotation = 30),
        plot_paired_bars(person_shares(u, "reference"), person_shares(u, "model"), u.category,
            title = "$(letters[3])) Universities", xlabel = "Students", ylabel = "Share of students",
            reference_label = "Destatis", xrotation = 30),
    ]
end
