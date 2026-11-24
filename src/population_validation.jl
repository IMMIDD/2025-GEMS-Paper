###
### POPULATION VALIDATION: GEMS POPULATION VS. OFFICIAL STATISTICS
###
# National panels (top) and agreement plots of the federal states (bottom).
# The code is split into four parts that can be reused separately:
#   1. reading the reference tables in `inputs/population/` (GENESIS-Online
#      and Zensus 2022 database exports)
#   2. aggregating the simulation's population in the same categories
#   3. drawing the panels (national bars, age pyramid, agreement plots)
#   4. composing and saving the figure
# All sizes are given in print units: the figure is drawn at its print size
# (100 px per inch) with font sizes in pt, and the resolution is only raised
# when saving.

export population_validation
export validation_data, save_validation_data, load_validation_data
export validation_figure
export plot_agreement, state_agreement_panels
export read_genesis_csv
export destatis_age_sex, mikrozensus_household_sizes
export zensus_household_types, zensus_senior_status
export HouseholdComposition, HouseholdCategory, HOUSEHOLD_TYPES, CORE_FAMILY_TYPES, FAMILY_TYPES, SENIOR_STATUS, classify
export is_other_family_type, family_type_counts
export PersonActivity, PersonCategory, ACTIVITY_STATUS, zensus_activity_status, ACTIVITY_INDICATORS, activity_indicators
export is_below_minimum_age, is_employed, is_in_education, is_not_working, is_retired
export SchoolProfile, EDUCATION_SETTINGS, COMPARED_EDUCATION_SETTINGS, is_daycare, is_primary_school, is_special_school, is_secondary_school, is_university
export official_education_settings, official_daycare_by_age
export is_single_senior, is_two_seniors_only, is_three_plus_seniors_only, is_seniors_and_younger, is_without_seniors
export is_single, is_couple_without_children, is_couple_with_children, is_single_mother, is_single_father,
    is_without_core_family
export plot_age_pyramid, plot_household_size_panel


# STYLE

const PV_FONT = "Times Roman"
const PV_FONTSIZE = 9
const PV_SMALL_FONTSIZE = 8
const PV_HEADING_FONTSIZE = 10
const PV_SECTION_FONTSIZE = 12
const PV_FONT_BOLD = "Times Bold"
# GEMS in blue, the reference data in orange (the first two GEMS colors)
const PV_MODEL_COLOR = gemscolors(2)[1]
const PV_REFERENCE_COLOR = gemscolors(2)[2]
const PV_REFERENCE_ALPHA = 0.45
const PV_LINE_COLOR = "#9aa0a8"

# print width (inches; the height follows from the layout) and export resolution
const PV_WIDTH_IN = 7.5
const PV_DPI = 300

# name of the model data: the default population of GEMS, generated with Gesyland
const PV_MODEL_LABEL = "GEMS (Gesyland)"

# English state names in the order of the state numbers in the AGS
const STATE_NAMES_EN = ["Schleswig-Holstein", "Hamburg", "Lower Saxony", "Bremen",
    "North Rhine-Westphalia", "Hesse", "Rhineland-Palatinate", "Baden-Württemberg",
    "Bavaria", "Saarland", "Berlin", "Brandenburg",
    "Mecklenburg-Vorpommern", "Saxony", "Saxony-Anhalt", "Thuringia"]

# German state names in the order of the state numbers in the AGS
const STATE_NAMES_DE = ["Schleswig-Holstein", "Hamburg", "Niedersachsen", "Bremen",
    "Nordrhein-Westfalen", "Hessen", "Rheinland-Pfalz", "Baden-Württemberg",
    "Bayern", "Saarland", "Berlin", "Brandenburg",
    "Mecklenburg-Vorpommern", "Sachsen", "Sachsen-Anhalt", "Thüringen"]


# 1. REFERENCE DATA

const REFERENCE_DIR = joinpath("inputs", "population")

"""
    read_genesis_csv(file)

Reads a table exported as CSV from GENESIS-Online or the Zensus 2022
database (semicolon-separated, decimal commas, title rows above and notes
below the table, a value flag after every value) and returns it in long
format:
- `row1`, `row2`, ...: the row labels (e.g. reference date, age)
- `col1`, `col2`, ...: the column labels, one per header row from the top
  (e.g. region, sex, unit)
- `value`: the value (`missing` for suppressed values such as `-`)

Value flags are dropped.
"""
function read_genesis_csv(file)
    rows = split.(readlines(file), ';')
    rows[1][1] = lstrip(rows[1][1], '﻿')

    # column header rows start with an empty field; the data follows below them
    first_header = findfirst(r -> length(r) > 1 && isempty(r[1]), rows)
    isnothing(first_header) && error("No table found in $file.")
    last_header = first_header
    while isempty(rows[last_header + 1][1])
        last_header += 1
    end
    headers = rows[first_header:last_header]
    nlabels = findfirst(!isempty, headers[1]) - 1

    # data rows: up to the first row without a label (notes start with "____")
    last_data = last_header
    while last_data < length(rows) && length(rows[last_data + 1]) > nlabels &&
            !isempty(rows[last_data + 1][1]) && !startswith(rows[last_data + 1][1], "_")
        last_data += 1
    end

    rowlabels = [String[] for _ in 1:nlabels]
    collabels = [String[] for _ in headers]
    value = Union{Float64, Missing}[]
    for r in rows[(last_header + 1):last_data], j in (nlabels + 1):2:length(r)
        isempty(headers[1][j]) && continue
        for k in 1:nlabels
            push!(rowlabels[k], r[k])
        end
        for (k, h) in enumerate(headers)
            push!(collabels[k], h[j])
        end
        v = tryparse(Float64, replace(r[j], "," => "."))
        push!(value, isnothing(v) ? missing : v)
    end

    return DataFrame(
        ["row$k" => l for (k, l) in enumerate(rowlabels)]...,
        ["col$k" => l for (k, l) in enumerate(collabels)]...,
        "value" => value)
end

"""
    state_number(region)

State number of a region label: `"01 Schleswig-Holstein"` (Zensus) or
`"Schleswig-Holstein"` (GENESIS) -> 1, `"Deutschland"` -> 0.
"""
function state_number(region)
    region == "Deutschland" && return 0
    m = match(r"^(\d+) ", region)
    !isnothing(m) && return parse(Int, m.captures[1])
    i = findfirst(==(region), STATE_NAMES_DE)
    isnothing(i) && error("Unknown region \"$region\".")
    return i
end

# "Unter 1 Jahr"/"unter 1 Jahr" -> 0, "1 Jahr"/"1-Jährige" -> 1, "90 Jahre und mehr" -> 90
age_from_label(label) = startswith(lowercase(label), "unter") ? 0 : parse(Int, match(r"\d+", label).match)

# "1 Person", "Haushalte mit 5 und mehr Personen" -> 1, 5
size_from_label(label) = parse(Int, match(r"\d+", label).match)

const SEX_LABELS = Dict("männlich" => "male", "weiblich" => "female")

"""
    destatis_age_sex(file = joinpath(REFERENCE_DIR, "12411-0013_de.csv");
        date = "31.12.2019")

Persons per federal state, sex and single year of age (0 to 90, where 90
means 90 and older) on `date` from GENESIS table 12411-0013 (Fortschreibung
des Bevölkerungsstandes).
"""
function destatis_age_sex(file = joinpath(REFERENCE_DIR, "12411-0013_de.csv"); date = "31.12.2019")
    df = read_genesis_csv(file)
    df = df[(df.row1 .== date) .& haskey.(Ref(SEX_LABELS), df.col2) .& (df.row2 .!= "Insgesamt"), :]
    isempty(df) && error("No data for $date in $file.")
    return DataFrame(
        state = state_number.(df.col1),
        sex = getindex.(Ref(SEX_LABELS), df.col2),
        age = age_from_label.(df.row2),
        count = df.value)
end

"""
    mikrozensus_household_sizes(file = joinpath(REFERENCE_DIR, "12211-9034_de.csv");
        year = "2019")

Private households per federal state and size (1 to 5, where 5 means five
or more persons) in `year` from GENESIS table 12211-9034 (Mikrozensus).
The table gives thousands of households; they are converted to counts.
"""
function mikrozensus_household_sizes(file = joinpath(REFERENCE_DIR, "12211-9034_de.csv"); year = "2019")
    df = read_genesis_csv(file)
    df = df[(df.row1 .== year) .& (df.col2 .!= "Insgesamt"), :]
    isempty(df) && error("No data for $year in $file.")
    return DataFrame(
        state = state_number.(df.col1),
        size = size_from_label.(df.col2),
        count = 1000 .* df.value)
end

"""
    zensus_household_categories(file, categories)

Private households per federal state and category from a Zensus 2022
household table (e.g. 5000H-1005, 5000H-1003). `categories` is a list of
`HouseholdCategory`s (see `HOUSEHOLD_TYPES`, `SENIOR_STATUS`); each one names
the row of the table it corresponds to (`census_label`). Totals and other
rows are dropped.
"""
function zensus_household_categories(file, categories)
    keys = Dict(c.census_label => c.key for c in categories)
    df = read_genesis_csv(file)
    df.state = state_number.(df.col1)
    df.label = strip.(df.row2)
    df = df[(df.col3 .== "Anzahl") .& haskey.(Ref(keys), df.label) .& (df.state .> 0), :]
    length(unique(df.label)) == length(categories) ||
        error("Not all categories found in $file: $(setdiff(keys |> Base.keys |> collect, df.label)).")
    return DataFrame(state = df.state, category = getindex.(Ref(keys), df.label), count = df.value)
end

"""
    zensus_household_types(file = joinpath(REFERENCE_DIR, "5000H-1005_de.csv"))

Private households per federal state and household type (`HOUSEHOLD_TYPES`)
on 15.05.2022 from Zensus table 5000H-1005.
"""
zensus_household_types(file = joinpath(REFERENCE_DIR, "5000H-1005_de.csv")) =
    zensus_household_categories(file, HOUSEHOLD_TYPES)

"""
    zensus_senior_status(file = joinpath(REFERENCE_DIR, "5000H-1003_de.csv"))

Private households per federal state and senior status (`SENIOR_STATUS`)
on 15.05.2022 from Zensus table 5000H-1003 (detailed senior status).
"""
zensus_senior_status(file = joinpath(REFERENCE_DIR, "5000H-1003_de.csv")) =
    zensus_household_categories(file, SENIOR_STATUS)

"""
    zensus_activity_status(file = joinpath(REFERENCE_DIR, "2000S-1002_de.csv"))

Persons per federal state and activity status (`ACTIVITY_STATUS`) on
15.05.2022 from Zensus table 2000S-1002 (Erwerbsstatus, ausführlich). Each
category sums the table rows it corresponds to (`census_labels`).
"""
function zensus_activity_status(file = joinpath(REFERENCE_DIR, "2000S-1002_de.csv"))
    keys = Dict(l => c.key for c in ACTIVITY_STATUS for l in c.census_labels)
    df = read_genesis_csv(file)
    df.state = state_number.(df.col1)
    df.label = strip.(df.row2)
    df = df[(df.col3 .== "Anzahl") .& haskey.(Ref(keys), df.label) .& (df.state .> 0), :]
    length(unique(df.label)) == length(keys) ||
        error("Not all categories found in $file: $(setdiff(collect(Base.keys(keys)), df.label)).")
    df.category = getindex.(Ref(keys), df.label)
    return combine(groupby(df, [:state, :category]), :value => sum => :count)
end

"""
    official_daycare_by_age(file = joinpath(REFERENCE_DIR, "22541-0112_de.csv"); date = "01.03.2019")

Children in daycare who do not attend school yet, per federal state (by
location of the facility) and single year of age (0 to 6) on `date`, from
GENESIS table 22541-0112. School children in after-school care (Hort) and
the open-ended group "7 Jahre und mehr" are left out.
"""
function official_daycare_by_age(file = joinpath(REFERENCE_DIR, "22541-0112_de.csv"); date = "01.03.2019")
    df = read_genesis_csv(file)
    df = df[(df.row1 .== date) .& (df.col1 .== "Alter der Nichtschulkinder") .&
            (df.col2 .!= "Insgesamt") .& .!startswith.(df.col2, "7"), :]
    isempty(df) && error("No data for $date in $file.")
    return DataFrame(state = state_number.(df.row2), age = age_from_label.(df.col2), count = df.value)
end

"""
    official_education_settings(; school_year = "2019/20", date = "01.03.2019")

Participants per federal state (by location of the institution) and
education setting (`COMPARED_EDUCATION_SETTINGS`), combined from three
GENESIS tables:
- daycare: children not attending school yet (22541-0112, see `official_daycare_by_age`)
- primary: Grundschulen, Vorklassen and Schulkindergärten (21111-0011)
- special: Förderschulen (21111-0011)
- secondary: all other general schools (21111-0011) and all vocational
  schools (21121-0010)
"""
function official_education_settings(; school_year = "2019/20", date = "01.03.2019")
    rows = DataFrame(state = Int[], category = String[], count = Float64[])
    add!(states, category, counts) = append!(rows, DataFrame(state = states, category = category, count = counts))

    daycare = combine(groupby(official_daycare_by_age(date = date), :state), :count => (c -> sum(skipmissing(c))) => :count)
    add!(daycare.state, "daycare", daycare.count)

    general = read_genesis_csv(joinpath(REFERENCE_DIR, "21111-0011_de.csv"))
    general = general[(general.col2 .== school_year) .& (general.col4 .== "Insgesamt") .& (general.row3 .== "Insgesamt"), :]
    isempty(general) && error("No general school data for $school_year.")
    general.state = state_number.(general.row1)
    count_of(types) = combine(groupby(general[in.(general.row2, Ref(types)), :], :state),
        :value => (v -> sum(skipmissing(v))) => :count)
    primary_types = ["Grundschulen", "Vorklassen", "Schulkindergärten"]
    primary, special = count_of(primary_types), count_of(["Förderschulen"])
    all_general = count_of(["Insgesamt"])
    add!(primary.state, "primary", primary.count)
    add!(special.state, "special", special.count)

    vocational = read_genesis_csv(joinpath(REFERENCE_DIR, "21121-0010_de.csv"))
    vocational = vocational[(vocational.row1 .== school_year) .& (vocational.row2 .== "Insgesamt"), :]
    isempty(vocational) && error("No vocational school data for $school_year.")
    vocational.state = state_number.(vocational.col2)
    secondary = innerjoin(innerjoin(all_general, primary, on = :state, makeunique = true),
        special, on = :state, makeunique = true)
    secondary = innerjoin(secondary, DataFrame(state = vocational.state, vocational = coalesce.(vocational.value, 0.0)), on = :state)
    add!(secondary.state, "secondary", secondary.count .- secondary.count_1 .- secondary.count_2 .+ secondary.vocational)
    return rows
end

"""
    VALIDATION_REFERENCES

Reference data sets for the validation: functions returning the reference
tables, and the names of their sources for figures (without reference
dates, which belong in the caption).
- `:official2019`: population on 31.12.2019 (Destatis, 12411-0013) and
  households 2019 (Microcensus, 12211-9034); household types, senior status
  and activity status from the Zensus 2022, which reports them for all
  federal states in categories that can be derived from the GEMS population;
  education settings from the school and daycare statistics 2019

Further data sets (e.g. for a later version of the population) are added
as further entries with the same fields.
"""
const VALIDATION_REFERENCES = Dict(
    :official2019 => (
        age_sex = destatis_age_sex,
        household_sizes = mikrozensus_household_sizes,
        household_types = zensus_household_types,
        senior_status = zensus_senior_status,
        activity_status = zensus_activity_status,
        education = official_education_settings,
        labels = (age_sex = "Destatis", household_sizes = "Microcensus",
            household_types = "Census", senior_status = "Census", activity_status = "Census",
            education = "Destatis")))

# source names of the reference data of `validation_data` (see `VALIDATION_REFERENCES`)
reference_labels(data) = VALIDATION_REFERENCES[Symbol(data.reference)].labels


# 2. MODEL DATA

# Gesyland models care homes as households of 90+ residents (sizes 10-89
# do not occur); official household statistics cover private households only
const PRIVATE_HH_MAX = 9

# federal state number (first two AGS digits) of an AGS
state_code(a) = id(a) ÷ 10^6

"""
    GESYLAND_SEX

Sex codes of the Gesyland population. Gesyland follows the microcensus
(1 = male, 2 = female), while GEMS defines `FEMALE = 1` and `MALE = 2`.
The Gesyland convention is confirmed by the population: code 2 dominates
at ages 85 and older, as women do in the official statistics.
"""
const GESYLAND_SEX = Dict(1 => "male", 2 => "female")

"""
    model_age_sex(individual_data_df; sex_codes = GESYLAND_SEX)

Persons per federal state, sex and single year of age in the simulation's
population (`individual_data_df = individual_data(sim)`).
"""
function model_age_sex(individual_data_df; sex_codes = GESYLAND_SEX)
    df = DataFrame(
        state = state_code.(individual_data_df.ags),
        sex = getindex.(Ref(sex_codes), Int.(individual_data_df.sex)),
        age = Int.(individual_data_df.age))
    return combine(groupby(df, [:state, :sex, :age]), nrow => :count)
end

"""
    model_household_sizes(sim)

Private households (up to `PRIVATE_HH_MAX` members) per federal state and
size in the simulation's population. Larger households are care homes.
"""
function model_household_sizes(sim)
    hhs = filter(h -> size(h) <= PRIVATE_HH_MAX, households(sim))
    df = DataFrame(state = state_code.(ags.(hhs)), size = size.(hhs))
    return combine(groupby(df, [:state, :size]), nrow => :count)
end


# HOUSEHOLD COMPOSITION
#
# The Zensus classifies households by the relationships of their members
# (couples, parents, children). GEMS only knows the members' ages and sexes,
# so each Zensus category is approximated by a predicate on these. Every
# predicate documents the Zensus definition and the approximation.

"""
    HouseholdComposition(ages, sexes)

Ages and sexes (`"male"`, `"female"`) of the members of one household.
"""
struct HouseholdComposition
    ages::Vector{Int}
    sexes::Vector{String}
end

members(h::HouseholdComposition) = length(h.ages)

"""
    SENIOR_AGE

Age from which a person counts as a senior (Zensus: 65 years and older).
"""
const SENIOR_AGE = 65

"""
    GENERATION_GAP

Minimum age difference between parents and children (years). Members at
least this much younger than the oldest member count as children; all
others (the oldest member and anyone less than this much younger) form
the parent generation.
"""
const GENERATION_GAP = 18

is_senior(age) = age >= SENIOR_AGE

"""
    parent_generation(h)

Indices of the members within `GENERATION_GAP` years of the oldest member.
"""
parent_generation(h::HouseholdComposition) = findall(a -> maximum(h.ages) - a < GENERATION_GAP, h.ages)

"""
    children(h)

Indices of the members at least `GENERATION_GAP` years younger than the
oldest member. Like the Zensus, this includes adult children.
"""
children(h::HouseholdComposition) = findall(a -> maximum(h.ages) - a >= GENERATION_GAP, h.ages)


# Senior status (Zensus 5000H-1003): purely age-based, reproduced exactly.

"""
    is_single_senior(h)

Zensus: "Singlehaushalte mit Senior/-in", a single-person household whose
member is a senior. GEMS: one member, aged `SENIOR_AGE` or older. Exact.
"""
is_single_senior(h::HouseholdComposition) = members(h) == 1 && is_senior(h.ages[1])

"""
    is_two_seniors_only(h)

Zensus: "Zweipersonenhaushalte mit ausschließlich Senioren/-innen". GEMS:
two members, both seniors. Exact.
"""
is_two_seniors_only(h::HouseholdComposition) = members(h) == 2 && all(is_senior, h.ages)

"""
    is_three_plus_seniors_only(h)

Zensus: "Mehrpersonenhaushalte mit ausschließlich Senioren/-innen" (three or
more members). GEMS: three or more members, all seniors. Exact.
"""
is_three_plus_seniors_only(h::HouseholdComposition) = members(h) >= 3 && all(is_senior, h.ages)

"""
    is_seniors_and_younger(h)

Zensus: "Haushalte mit Senioren/-innen und Jüngeren". GEMS: at least one
senior and at least one younger member. Exact.
"""
is_seniors_and_younger(h::HouseholdComposition) = any(is_senior, h.ages) && !all(is_senior, h.ages)

"""
    is_without_seniors(h)

Zensus: "Haushalte ohne Senioren/-innen". GEMS: no member is a senior. Exact.
"""
is_without_seniors(h::HouseholdComposition) = !any(is_senior, h.ages)


# Household type (Zensus 5000H-1005): based on relationships, approximated
# by ages (children = members at least `GENERATION_GAP` years younger than
# the oldest member) and sexes.

"""
    is_single(h)

Zensus: "Einpersonenhaushalte (Singlehaushalte)". GEMS: one member. Exact.
"""
is_single(h::HouseholdComposition) = members(h) == 1

"""
    is_couple_without_children(h)

Zensus: "Paare ohne Kind(er)", a married or unmarried couple (of any sex)
without children in the household. GEMS: two members of the parent
generation and no children.

Deviations: two adults of similar age who are not a couple (e.g. a flat
share of two, two siblings) count as a couple.
"""
is_couple_without_children(h::HouseholdComposition) =
    members(h) >= 2 && isempty(children(h)) && length(parent_generation(h)) == 2

"""
    is_couple_with_children(h)

Zensus: "Paare mit Kind(ern)", a couple with at least one unmarried child
of any age in the household. GEMS: two members of the parent generation
and at least one child.

Deviations: a single parent living with a sibling or a flatmate of similar
age counts as a couple; a couple whose youngest member is less than
`GENERATION_GAP` years younger than the oldest (e.g. a young stepchild of
a much younger partner) may be misclassified.
"""
is_couple_with_children(h::HouseholdComposition) =
    !isempty(children(h)) && length(parent_generation(h)) == 2

"""
    is_single_mother(h)

Zensus: "Alleinerziehende Mütter", a mother without a partner in the
household, living with at least one unmarried child of any age. GEMS: one
member of the parent generation, who is female, and at least one child.

Deviations: in three-generation households (e.g. grandmother, mother,
child), the grandmother is the only member of the parent generation, and
the household counts as a single mother.
"""
is_single_mother(h::HouseholdComposition) =
    !isempty(children(h)) && length(parent_generation(h)) == 1 && h.sexes[only(parent_generation(h))] == "female"

"""
    is_single_father(h)

Zensus: "Alleinerziehende Väter". GEMS: as `is_single_mother`, with a male
member of the parent generation. Same deviations.
"""
is_single_father(h::HouseholdComposition) =
    !isempty(children(h)) && length(parent_generation(h)) == 1 && h.sexes[only(parent_generation(h))] == "male"

"""
    is_without_core_family(h)

Zensus: "Mehrpersonenhaushalte ohne Kernfamilie", several members but no
couple and no parent-child relationship (e.g. flat shares, siblings).
GEMS: several members and three or more members in the parent generation.

Deviations: flat shares of two count as couples (see
`is_couple_without_children`); couples living with a third adult of similar
age (e.g. a parent of one partner within `GENERATION_GAP` years) count here.
"""
is_without_core_family(h::HouseholdComposition) = members(h) >= 2 && length(parent_generation(h)) >= 3

"""
    HouseholdCategory(key, label, predicate, census_label)

One category of a household classification: an identifier, a short label
for figures, the predicate deciding whether a household belongs to it, and
the label of the corresponding row in the Zensus table.
"""
struct HouseholdCategory
    key::String
    label::String
    predicate::Function
    census_label::String
end

"""
    HOUSEHOLD_TYPES

Household types of Zensus table 5000H-1005 and their GEMS predicates.
"""
const HOUSEHOLD_TYPES = [
    HouseholdCategory("single", "Single", is_single, "Einpersonenhaushalte (Singlehaushalte)"),
    HouseholdCategory("couple", "Couple", is_couple_without_children, "Paare ohne Kind(er)"),
    HouseholdCategory("couple_children", "Couple, kids", is_couple_with_children, "Paare mit Kind(ern)"),
    HouseholdCategory("single_mother", "Lone mother", is_single_mother, "Alleinerziehende Mütter"),
    HouseholdCategory("single_father", "Lone father", is_single_father, "Alleinerziehende Väter"),
    HouseholdCategory("no_core_family", "Other", is_without_core_family, "Mehrpersonenhaushalte ohne Kernfamilie"),
]

"""
    SENIOR_STATUS

Senior status of Zensus table 5000H-1003 (detailed) and the GEMS predicates.
"""
const SENIOR_STATUS = [
    HouseholdCategory("single_senior", "1 senior", is_single_senior, "Singlehaushalte mit Senior/-in"),
    HouseholdCategory("two_seniors", "2 seniors", is_two_seniors_only,
        "Zweipersonenhaushalte mit ausschl. Senioren/-innen"),
    HouseholdCategory("three_plus_seniors", "3+ seniors", is_three_plus_seniors_only,
        "Mehrpersonenhaushalte mit ausschl. Senioren/-innen"),
    HouseholdCategory("seniors_and_younger", "Mixed", is_seniors_and_younger,
        "Haushalte mit Senioren/-innen und Jüngeren"),
    HouseholdCategory("no_seniors", "No seniors", is_without_seniors, "Haushalte ohne Senioren/-innen"),
]

"""
    classify(h, categories)

Key of the category of household `h`. Each household must match exactly one
predicate of `categories`; an error names the household otherwise, so that
gaps and overlaps of the predicates are found.
"""
function classify(h::HouseholdComposition, categories)
    hits = filter(c -> c.predicate(h), categories)
    length(hits) == 1 ||
        error("Household with ages $(h.ages) and sexes $(h.sexes) matches $(length(hits)) categories: $(getfield.(hits, :key)).")
    return only(hits).key
end

"""
    household_compositions(sim; sex_codes = GESYLAND_SEX)

Federal state and `HouseholdComposition` of every private household (up to
`PRIVATE_HH_MAX` members) in the simulation's population.
"""
function household_compositions(sim; sex_codes = GESYLAND_SEX)
    hhs = filter(h -> size(h) <= PRIVATE_HH_MAX, households(sim))
    comps = map(hhs) do h
        inds = individuals(h)
        HouseholdComposition(Int.(age.(inds)), getindex.(Ref(sex_codes), Int.(sex.(inds))))
    end
    return state_code.(ags.(hhs)), comps
end

"""
    model_household_categories(states, compositions, categories)

Private households per federal state and category (see `classify`).
"""
function model_household_categories(states, compositions, categories)
    df = DataFrame(state = states, category = classify.(compositions, Ref(categories)))
    return combine(groupby(df, [:state, :category]), nrow => :count)
end


# ACTIVITY STATUS
#
# The Zensus classifies persons by the labour-force concept: employment
# comes first (anyone with paid work in the reference week, including
# working students and apprentices), then unemployment, then the groups of
# persons not in the labour force. GEMS knows whether a person is assigned
# to a workplace or to a school setting (daycare, school, university), and
# no person has both.

"""
    PersonActivity(age, workplace, school)

Age of a person and whether they are assigned to a workplace (an office)
and to a school setting (a school class).
"""
struct PersonActivity
    age::Int
    workplace::Bool
    school::Bool
end

"""
    MINIMUM_WORKING_AGE

Age from which persons are part of the labour-force classification
(Zensus: "Personen unterhalb des Mindestalters" are those under 15).
"""
const MINIMUM_WORKING_AGE = 15

"""
    is_below_minimum_age(p)

Zensus: "Personen unterhalb des Mindestalters", persons under 15 (whether
in daycare, school or neither). GEMS: age under `MINIMUM_WORKING_AGE`. Exact.
"""
is_below_minimum_age(p::PersonActivity) = p.age < MINIMUM_WORKING_AGE

"""
    is_employed(p)

Zensus: "Erwerbstätige", persons aged 15 and older with paid work in the
reference week, including working students and apprentices in dual
vocational training. GEMS: aged 15 or older and assigned to a workplace.

Deviations: GEMS assigns nobody to both a workplace and a school setting,
so working students and apprentices (who attend vocational school) are in
education in GEMS (see `is_in_education`).
"""
is_employed(p::PersonActivity) = p.age >= MINIMUM_WORKING_AGE && p.workplace

"""
    is_in_education(p)

Zensus: "Schüler/-innen u. Studierende (nicht erwerbsaktiv)", pupils and
students aged 15 and older who are neither employed nor unemployed. GEMS:
aged 15 or older, assigned to a school setting and not to a workplace.

Deviations: GEMS counts all pupils and students aged 15 or older here,
including those who would work in reality (working students, apprentices).
"""
is_in_education(p::PersonActivity) = p.age >= MINIMUM_WORKING_AGE && p.school && !p.workplace

"""
    RETIREMENT_AGE

Age from which persons without a workplace and without an education
setting count as retired in GEMS (statutory retirement age).
"""
const RETIREMENT_AGE = 65

"""
    is_not_working(p)

Zensus: "Erwerbslose" (unemployed: not employed, looking for work and
available) and "Sonstige" (others not in the labour force, e.g. homemakers,
carers, persons unable to work). GEMS: aged 15 to under `RETIREMENT_AGE`,
assigned to neither a workplace nor an education setting.

Deviations: GEMS has no attribute for retirement, unemployment or other
reasons for not working, so age separates this category from the retired.
Early retirees under 65, who are retired in the Zensus, are here in GEMS;
"Sonstige" aged 65 and older are retired in GEMS.
"""
is_not_working(p::PersonActivity) =
    MINIMUM_WORKING_AGE <= p.age < RETIREMENT_AGE && !p.workplace && !p.school

"""
    is_retired(p)

Zensus: "Empfänger/-innen von Ruhegehalt/Kapitalerträgen", persons living
on pensions or capital income (at any age). GEMS: aged `RETIREMENT_AGE`
or older, assigned to neither a workplace nor an education setting.

Deviations: see `is_not_working` (early retirees under 65 are not here).
"""
is_retired(p::PersonActivity) = p.age >= RETIREMENT_AGE && !p.workplace && !p.school

"""
    PersonCategory(key, label, predicate, census_labels)

One category of a classification of persons: an identifier, a short label
for figures, the predicate deciding whether a person belongs to it, and
the labels of the rows of the Zensus table it corresponds to (their counts
are summed).
"""
struct PersonCategory
    key::String
    label::String
    predicate::Function
    census_labels::Vector{String}
end

"""
    ACTIVITY_STATUS

Activity status of Zensus table 2000S-1002 (Erwerbsstatus, ausführlich)
and the GEMS predicates; every person matches exactly one.
"""
const ACTIVITY_STATUS = [
    PersonCategory("below_15", "Under 15", is_below_minimum_age, ["Personen unterhalb des Mindestalters"]),
    PersonCategory("employed", "Employed", is_employed, ["Erwerbstätige"]),
    PersonCategory("in_education", "In education", is_in_education,
        ["Schüler/-innen u. Studierende (nicht erwerbsaktiv)"]),
    PersonCategory("not_working", "Not working (15–64)", is_not_working, ["Erwerbslose", "Sonstige"]),
    PersonCategory("retired", "Retired (65+)", is_retired, ["Empfänger/-innen von Ruhegehalt/Kapitalerträgen"]),
]

"""
    classify(p, categories)

Key of the category of person `p`. Each person must match exactly one
predicate of `categories`; an error names the person otherwise.
"""
function classify(p::PersonActivity, categories)
    hits = filter(c -> c.predicate(p), categories)
    length(hits) == 1 ||
        error("Person $(p) matches $(length(hits)) categories: $(getfield.(hits, :key)).")
    return only(hits).key
end

"""
    model_activity_status(sim)

Persons in private households (up to `PRIVATE_HH_MAX` members, as in the
Zensus table) per federal state and activity status (`ACTIVITY_STATUS`).
"""
function model_activity_status(sim)
    unset = GEMS.DEFAULT_SETTING_ID
    states, keys = Int[], String[]
    for h in households(sim)
        size(h) <= PRIVATE_HH_MAX || continue
        st = state_code(ags(h))
        for i in individuals(h)
            p = PersonActivity(Int(age(i)), office_id(i) != unset, class_id(i) != unset)
            push!(states, st)
            push!(keys, classify(p, ACTIVITY_STATUS))
        end
    end
    return combine(groupby(DataFrame(state = states, category = keys), [:state, :category]), nrow => :count)
end


# EDUCATION SETTINGS
#
# Gesyland does not label the type of an education institution (GEMS
# "School"). Each type has a fixed class capacity, though (crèche 12,
# daycare 20, primary school 20, special school 15, orientation stage 25,
# secondary and vocational schools 30, universities of applied sciences 50,
# universities 90), so the most common class size of an institution and the
# median age of its attendees identify its type. Official statistics count
# participants at the location of the institution, so GEMS counts them at
# the location of their school class, too.

"""
    SchoolProfile(class_size, median_age)

Most common class size of an education institution and the median age of
its attendees.
"""
struct SchoolProfile
    class_size::Int
    median_age::Float64
end

# class sizes that identify an institution type
const KNOWN_CLASS_SIZES = (12, 15, 20, 25, 30, 50, 90)
known_class_size(s::SchoolProfile) = s.class_size in KNOWN_CLASS_SIZES

"""
    is_daycare(s)

Official: children in daycare who do not attend school yet (crèche and
Kita). GEMS: class size 12 (crèche) or 20 with a median age under 6 (Kita);
institutions without a known class size (all classes not full) with a
median age under 6.

Deviations: after-school care (Hort) is not assigned in GEMS and not
counted in the reference.
"""
is_daycare(s::SchoolProfile) = s.class_size == 12 || (s.class_size == 20 && s.median_age < 6) ||
    (!known_class_size(s) && s.median_age < 6)

"""
    is_primary_school(s)

Official: Grundschulen, Vorklassen and Schulkindergärten. GEMS: class size
20 with a median age of 6 or older; institutions without a known class
size with a median age from 6 to under 10.
"""
is_primary_school(s::SchoolProfile) = (s.class_size == 20 && s.median_age >= 6) ||
    (!known_class_size(s) && 6 <= s.median_age < 10)

"""
    is_special_school(s)

Official: Förderschulen. GEMS: class size 15. Pupils with special needs
taught in regular schools count in their school type in both.
"""
is_special_school(s::SchoolProfile) = s.class_size == 15

"""
    is_secondary_school(s)

Official: all general schools except primary and special schools
(Gymnasium, Realschule, Hauptschule, schools with several tracks, Waldorf
schools, evening schools) and all vocational schools. GEMS: class size 30
(all secondary and vocational school types) or 25 (orientation stage);
institutions without a known class size with a median age from 10 to under 19.

Deviations: GEMS cannot tell the secondary and vocational school types
apart, as they share the class size; Waldorf schools also teach primary
grades but count here entirely.
"""
is_secondary_school(s::SchoolProfile) = s.class_size in (25, 30) ||
    (!known_class_size(s) && 10 <= s.median_age < 19)

"""
    is_university(s)

GEMS: class size 90 (universities) or 50 (universities of applied
sciences); institutions without a known class size with a median age of 19
or older.

Universities are classified so that every institution has exactly one
category, but not compared with official statistics (see
`COMPARED_EDUCATION_SETTINGS`).
"""
is_university(s::SchoolProfile) = s.class_size in (50, 90) ||
    (!known_class_size(s) && s.median_age >= 19)

"""
    EDUCATION_SETTINGS

Education settings and their GEMS predicates (on `SchoolProfile`s); the
official counts are combined in `official_education_settings`.
"""
const EDUCATION_SETTINGS = [
    PersonCategory("daycare", "Daycare", is_daycare, String[]),
    PersonCategory("primary", "Primary", is_primary_school, String[]),
    PersonCategory("special", "Special", is_special_school, String[]),
    PersonCategory("secondary", "Secondary, vocational", is_secondary_school, String[]),
    PersonCategory("university", "University", is_university, String[]),
]

"""
    COMPARED_EDUCATION_SETTINGS

Keys of the education settings compared with official statistics. Not
universities: the official statistics count students at the location of
the university, but GEMS assigns students to universities regardless of
distance, so shares per state would reflect this assignment rather than
enrollment.
"""
const COMPARED_EDUCATION_SETTINGS = ("daycare", "primary", "special", "secondary")

"""
    classify(s, categories)

Key of the education setting of an institution with profile `s`. Each
profile must match exactly one predicate of `categories`.
"""
function classify(s::SchoolProfile, categories)
    hits = filter(c -> c.predicate(s), categories)
    length(hits) == 1 ||
        error("Institution $(s) matches $(length(hits)) categories: $(getfield.(hits, :key)).")
    return only(hits).key
end

"""
    education_institutions(sim)

Education institutions (`School` settings) of the simulation's population
and their classes, as a named tuple:
- `classes`: all `SchoolClass` settings
- `institution`: index of the institution of each class
- `type`: key of the education setting of each institution
  (`EDUCATION_SETTINGS`, classified by its `SchoolProfile`); `""` for
  institutions without attendees

GEMS has no institution types, so the type follows from the class sizes
and ages (see the predicates of `EDUCATION_SETTINGS`).
"""
function education_institutions(sim)
    classes, years, schools = settings(sim, SchoolClass), settings(sim, SchoolYear), settings(sim, School)
    institution = [contained(years[contained(c)]) for c in classes]

    # profile and type of each institution
    sizes = [Int[] for _ in schools]
    ages = [Int[] for _ in schools]
    for (k, c) in enumerate(classes)
        push!(sizes[institution[k]], size(c))
        append!(ages[institution[k]], Int.(age.(individuals(c))))
    end
    type = [isempty(ages[s]) ? "" : classify(SchoolProfile(mode(sizes[s]), median(ages[s])), EDUCATION_SETTINGS)
        for s in eachindex(schools)]
    return (classes = classes, institution = institution, type = type)
end

"""
    model_education(sim)

Participants in the education settings `COMPARED_EDUCATION_SETTINGS` per
federal state (by location of their school class) and setting: state,
category, count. Every institution is classified (`education_institutions`);
participants of the other settings are not counted.
"""
function model_education(sim)
    (; classes, institution, type) = education_institutions(sim)
    education = Dict{Tuple{Int, String}, Int}()
    for (k, c) in enumerate(classes)
        t = type[institution[k]]
        t in COMPARED_EDUCATION_SETTINGS || continue
        key = (state_code(ags(c)), t)
        education[key] = get(education, key, 0) + size(c)
    end
    return DataFrame(state = first.(keys(education)), category = last.(keys(education)),
        count = collect(values(education)))
end

"""
    topcode(df, col, top)

Merges all values of `col` above `top` into `top` and sums their counts.
"""
topcode(df, col, top) = combine(
    groupby(transform(df, col => (x -> min.(x, top)) => col), setdiff(names(df), ["count"])),
    :count => (c -> sum(skipmissing(c))) => :count)

"""
    join_counts(model, reference, keys)

Joins model and reference counts on `keys` into one table with the
columns `model` and `reference` (missing combinations count as 0).
"""
function join_counts(model, reference, keys)
    df = outerjoin(
        rename(model, :count => :model),
        rename(reference, :count => :reference),
        on = keys)
    df.model = coalesce.(df.model, 0)
    df.reference = coalesce.(df.reference, 0)
    return sort!(df, keys)
end

"""
    validation_data(sim;
        reference = :official2019,
        individual_data_df = individual_data(sim))

Model and reference counts in common categories (see
`VALIDATION_REFERENCES`), as a named tuple:
- `age_sex`: state, sex, age, model, reference
- `household_sizes`: state, size, model, reference
- `household_types`: state, category (`HOUSEHOLD_TYPES`), model, reference
- `senior_status`: state, category (`SENIOR_STATUS`), model, reference
- `activity_status`: state, category (`ACTIVITY_STATUS`), model, reference (persons)
- `education`: state, category (`COMPARED_EDUCATION_SETTINGS`), model, reference (participants, by location)
- `workplaces`, `universities`, `groups`: setting sizes in Germany (`setting_size_data`)
- `reference`: name of the reference data set (a key of `VALIDATION_REFERENCES`)

The oldest age and the largest household size are the top groups that
both the model and the reference provide (e.g. 90 and older, 5 and more).
"""
function validation_data(sim;
        reference = :official2019,
        individual_data_df = individual_data(sim))

    ref = VALIDATION_REFERENCES[reference]
    model_ages, ref_ages = model_age_sex(individual_data_df), ref.age_sex()
    model_hh, ref_hh = model_household_sizes(sim), ref.household_sizes()

    top_age = min(maximum(model_ages.age), maximum(ref_ages.age))
    top_size = min(maximum(model_hh.size), maximum(ref_hh.size))

    states, comps = household_compositions(sim)

    return (
        age_sex = join_counts(topcode(model_ages, :age, top_age), topcode(ref_ages, :age, top_age),
            [:state, :sex, :age]),
        household_sizes = join_counts(topcode(model_hh, :size, top_size), topcode(ref_hh, :size, top_size),
            [:state, :size]),
        household_types = join_counts(model_household_categories(states, comps, HOUSEHOLD_TYPES),
            ref.household_types(), [:state, :category]),
        senior_status = join_counts(model_household_categories(states, comps, SENIOR_STATUS),
            ref.senior_status(), [:state, :category]),
        activity_status = join_counts(model_activity_status(sim), ref.activity_status(), [:state, :category]),
        education = join_counts(model_education(sim), ref.education(), [:state, :category]),
        setting_size_data(sim)...,
        reference = string(reference))
end

# tables of `validation_data` and their file names
const VALIDATION_TABLES = (:age_sex, :household_sizes, :household_types, :senior_status, :activity_status,
    :education, :workplaces, :universities, :groups)
validation_file(folder, table) = joinpath(folder, "validation_$(table).csv")

"""
    save_validation_data(data, folder)

Writes the tables of `validation_data` as CSV files and the name of the
reference data set to `folder`.
"""
function save_validation_data(data, folder)
    mkpath(folder)
    foreach(t -> CSV.write(validation_file(folder, t), getproperty(data, t)), VALIDATION_TABLES)
    open(io -> TOML.print(io, Dict("reference" => data.reference)),
        joinpath(folder, "validation_reference.toml"), "w")
end

"""
    load_validation_data(folder)

Reads the data written by `save_validation_data`, so that the figure can
be redrawn without loading the population again. Categories are read as
strings, as written.
"""
function load_validation_data(folder)
    tables = map(VALIDATION_TABLES) do t
        t => CSV.read(validation_file(folder, t), DataFrame,
            types = Dict(:category => String, :sex => String), validate = false)
    end
    reference = TOML.parsefile(joinpath(folder, "validation_reference.toml"))["reference"]
    return (; tables..., reference = reference)
end


# 3. PANELS

"""
    rectangles(x0, x1, y0, y1)

Coordinates of axis-aligned rectangles, separated by `NaN`, so that all of
them can be drawn as one `:shape` series (one legend entry, no seams).
"""
function rectangles(x0, x1, y0, y1)
    xs, ys = Float64[], Float64[]
    for i in eachindex(x0)
        append!(xs, [x0[i], x1[i], x1[i], x0[i], x0[i], NaN])
        append!(ys, [y0[i], y0[i], y1[i], y1[i], y0[i], NaN])
    end
    return xs, ys
end

"""
    staircase(values, positions; closed = false)

Outline of horizontal bars of length `values` at `positions` (spaced by 1),
as one line from bottom to top. With `closed = true`, the outline starts and
ends on the axis, so that it can be filled as one area.
"""
function staircase(values, positions; closed = false)
    xs = repeat(values, inner = 2)
    ys = vec([positions' .- 0.5; positions' .+ 0.5])
    closed || return xs, ys
    return [0.0; xs; 0.0], [first(ys); ys; last(ys)]
end

pv_text(s, size = PV_FONTSIZE, args...) = Plots.text(s, size, PV_FONT, args...)

percent_label(x; digits = 0) = string(round(100x, digits = digits) |> (v -> digits == 0 ? Int(v) : v), "%")

# label of a top-coded category (e.g. "5+")
category_label(v, top) = v == top ? "$(v)+" : string(v)

"""
    national_shares(df, keys)

Sums the counts of all states per `keys` and converts model and reference
counts into shares of their totals.
"""
function national_shares(df, keys)
    nat = combine(groupby(df, keys), :model => sum => :model, :reference => sum => :reference)
    nat.model = nat.model ./ sum(nat.model)
    nat.reference = nat.reference ./ sum(nat.reference)
    return sort!(nat, keys)
end

"""
    plot_age_pyramid(age_sex; title = "A) Age and sex", reference_label = "Reference",
        show_top_group = true, kwargs...)

Population pyramid of Germany: shares of the total population per sex and
single year of age. Reference as area, GEMS as a line. The oldest age is a
top-coded group (e.g. 90 and older), which appears as a wide bar; with
`show_top_group = false` it is not drawn (shares still refer to the whole
population).
"""
function plot_age_pyramid(age_sex; title = "A) Age and sex", reference_label = "Reference",
        show_top_group = true, age_label = "Age", xtick_step = 0.002, legend = false, kwargs...)

    nat = national_shares(age_sex, [:sex, :age])
    show_top_group || (nat = nat[nat.age .< maximum(nat.age), :])
    top_age = maximum(nat.age)
    side(sex) = sex == "male" ? -1 : 1
    xmax = 1.1 * maximum([nat.model; nat.reference])
    xticks = collect(xtick_step:xtick_step:xmax)
    yticks = 0:10:top_age

    p = plot(
        title = title,
        titlelocation = :left,
        titlefontsize = PV_HEADING_FONTSIZE,
        xlims = (-xmax, xmax),
        ylims = (-0.5, top_age + 0.5),
        xticks = ([-reverse(xticks); 0; xticks], percent_label.(abs.([-reverse(xticks); 0; xticks]), digits = 1)),
        yticks = (yticks, category_label.(yticks, top_age)),
        ylabel = age_label,
        grid = :x,
        gridalpha = 0.2,
        legend = legend,
        legendfontsize = PV_SMALL_FONTSIZE,
        foreground_color_legend = nothing,
        background_color_legend = nothing,
        tickfontsize = PV_SMALL_FONTSIZE,
        guidefontsize = PV_FONTSIZE,
        fontfamily = PV_FONT;
        kwargs...)

    for sex in ("male", "female")
        sub = nat[nat.sex .== sex, :]
        plot!(p, staircase(side(sex) .* sub.reference, sub.age, closed = true),
            seriestype = :shape, fillcolor = PV_REFERENCE_COLOR, fillalpha = PV_REFERENCE_ALPHA,
            linewidth = 0, label = sex == "male" ? reference_label : "")
        plot!(p, staircase(side(sex) .* sub.model, sub.age),
            linecolor = PV_MODEL_COLOR, linewidth = 1, label = sex == "male" ? PV_MODEL_LABEL : "")
    end

    vline!(p, [0], linecolor = :black, linewidth = 0.5, label = "")
    # sex labels in the lower corners, which the youngest ages leave free
    annotate!(p, -0.95xmax, 3, pv_text("Men", PV_FONTSIZE, :left, :bottom))
    annotate!(p, 0.95xmax, 3, pv_text("Women", PV_FONTSIZE, :right, :bottom))
    return p
end

"""
    plot_paired_bars(reference, model, labels;
        title = "", xlabel = "", reference_label = "Reference",
        legend = :auto, xrotation = 0, shares = true, ytick_step = nothing, kwargs...)

Values per category as pairs of bars: reference (left) and GEMS (right).
`reference` and `model` are given in the order of `labels`: shares (y axis
in percent, at most 100 %) or, with `shares = false`, other values such as
group sizes (y axis with whole numbers in steps of `ytick_step`). A
`missing` reference value is marked "n/a" instead of a bar, so that it
is not read as zero.

With `legend = :auto`, the legend goes into the upper corner above the
lower bars (left or right, compared by the two outermost categories on each
side); its background is transparent, so it never hides a bar, and the
y axis leaves room above the highest bar.
"""
function plot_paired_bars(reference, model, labels;
        title = "", xlabel = "", reference_label = "Reference",
        legend = :auto, xrotation = 0, shares = true, ytick_step = nothing, kwargs...)

    n = length(labels)
    x = collect(1.0:n)
    # headroom for the legend; shares never above 100 %
    ymax = 1.3 * maximum(skipmissing([model; reference]))
    shares && (ymax = min(ymax, 1.0))
    step = something(ytick_step, shares ? (ymax > 0.5 ? 0.2 : 0.1) : 5)
    yticks = 0:step:ymax
    ylabels = shares ? percent_label.(yticks) : string.(round.(Int, yticks))
    if legend == :auto
        height(i) = max(coalesce(reference[i], 0), model[i])
        edge = min(2, n)
        left, right = maximum(height.(1:edge)), maximum(height.((n - edge + 1):n))
        legend = left <= right ? :topleft : :topright
    end

    p = plot(
        title = title,
        titlelocation = :left,
        titlefontsize = PV_HEADING_FONTSIZE,
        xlims = (0.4, n + 0.6),
        ylims = (0, ymax),
        xticks = (x, labels),
        xrotation = xrotation,
        yticks = (yticks, ylabels),
        xlabel = xlabel,
        grid = :y,
        gridalpha = 0.2,
        legend = legend,
        legendfontsize = PV_SMALL_FONTSIZE,
        foreground_color_legend = nothing,
        background_color_legend = nothing,
        tickfontsize = PV_SMALL_FONTSIZE,
        guidefontsize = PV_FONTSIZE,
        fontfamily = PV_FONT;
        kwargs...)

    w = 0.38
    zero = zeros(n)
    plot!(p, rectangles(x .- w, x, zero, coalesce.(reference, 0.0)),
        seriestype = :shape, fillcolor = PV_REFERENCE_COLOR, fillalpha = PV_REFERENCE_ALPHA,
        linewidth = 0, label = reference_label)
    plot!(p, rectangles(x, x .+ w, zero, model),
        seriestype = :shape, fillcolor = PV_MODEL_COLOR, fillalpha = 0.85,
        linewidth = 0, label = PV_MODEL_LABEL)
    for i in findall(ismissing, reference)
        annotate!(p, x[i] - w / 2, 0.02ymax,
            Plots.text("n/a", PV_SMALL_FONTSIZE, PV_FONT, :left, :vcenter, rotation = 90, color = PV_LINE_COLOR))
    end
    return p
end

"""
    plot_household_size_panel(household_sizes; title = "B) Household size",
        reference_label = "Reference", kwargs...)

Shares of private households per size in Germany: reference and GEMS as
pairs of bars.
"""
function plot_household_size_panel(household_sizes; title = "B) Household size",
        reference_label = "Reference", kwargs...)
    nat = national_shares(household_sizes, [:size])
    return plot_paired_bars(nat.reference, nat.model, category_label.(nat.size, maximum(nat.size));
        title = title, xlabel = "Persons per household", reference_label = reference_label, kwargs...)
end

"""
    only_categories(df, categories)

Rows of `df` that belong to one of `categories`, so that shares computed
from the result refer to these categories only.
"""
only_categories(df, categories) = df[in.(df.category, Ref(getfield.(categories, :key))), :]

"""
    CORE_FAMILY_TYPES

The household types (`HOUSEHOLD_TYPES`) of households with a core family:
couples without children, couples with children, single mothers and single
fathers.
"""
const CORE_FAMILY_TYPES = filter(c -> !(c.key in ("single", "no_core_family")), HOUSEHOLD_TYPES)

"""
    is_other_family_type(h)

All households that match none of the `CORE_FAMILY_TYPES` predicates.
Since every household matches exactly one predicate of `HOUSEHOLD_TYPES`
(see `classify`), these are single-person households (`is_single`) and
multi-person households without a core family (`is_without_core_family`).

Zensus: the total number of households minus the four core family types of
table 5000H-1005, i.e. "Einpersonenhaushalte (Singlehaushalte)" and
"Mehrpersonenhaushalte ohne Kernfamilie", whose rows add up with the core
family types to the total.

Deviations: two-person flat shares count as couples in GEMS (see
`is_couple_without_children`) but belong here in the Zensus.
"""
is_other_family_type(h::HouseholdComposition) = !any(c -> c.predicate(h), CORE_FAMILY_TYPES)

"""
    FAMILY_TYPES

Family types of all households: the `CORE_FAMILY_TYPES` and "Other" for all
remaining households (`is_other_family_type`). Their shares refer to all
private households. The counts are derived from the household types with
`family_type_counts`.
"""
const FAMILY_TYPES = [CORE_FAMILY_TYPES;
    HouseholdCategory("other", "Other", is_other_family_type,
        "Einpersonenhaushalte (Singlehaushalte) + Mehrpersonenhaushalte ohne Kernfamilie")]

"""
    family_type_counts(household_types)

Counts per federal state and family type (`FAMILY_TYPES`) from the counts
per household type (the `household_types` table of `validation_data`): all
household types that are not a core family type are summed up as "other",
for the model and the reference alike.
"""
function family_type_counts(household_types)
    core = getfield.(CORE_FAMILY_TYPES, :key)
    df = transform(household_types, :category => ByRow(c -> c in core ? c : "other") => :category)
    return combine(groupby(df, [:state, :category]), :model => sum => :model, :reference => sum => :reference)
end

"""
    plot_household_category_panel(df, categories; title, reference_label = "Reference", kwargs...)

Shares of private households per category (`HOUSEHOLD_TYPES`,
`FAMILY_TYPES`, `SENIOR_STATUS`) in Germany: reference and GEMS as pairs of
bars. Shares refer to the households in `categories`.
"""
function plot_household_category_panel(df, categories; title, reference_label = "Reference", kwargs...)
    nat = national_shares(only_categories(df, categories), [:category])
    order = [findfirst(==(c.key), nat.category) for c in categories]
    return plot_paired_bars(nat.reference[order], nat.model[order], getfield.(categories, :label);
        title = title, reference_label = reference_label, kwargs...)
end


"""
    state_shares(df, key, categories)

Matrix of shares (16 states × categories) for model and reference.
"""
function state_shares(df, key, categories)
    shares(col) = [begin
            sub = df[df.state .== s, :]
            tot = sum(sub[!, col])
            [sum(sub[sub[!, key] .== c, col]) / tot for c in categories]
        end for s in 1:16] |> rows -> permutedims(reduce(hcat, rows))
    return shares(:model), shares(:reference)
end

# AGREEMENT PLOTS

"""
    plot_agreement(reference, model;
        labels = String[],
        nlabels = 3,
        deviation = abs.(model .- reference),
        groups = nothing,
        summary = "",
        log_axes = false,
        ticks = nothing,
        tick_step = nothing,
        tick_label = nothing,
        reference_label = "Reference",
        show_ylabel = true,
        title = "",
        panel_width_pt = 145,
        panel_height_pt = 150,
        markersize = 4,
        markeralpha = nothing,
        reverse_groups = false,
        kwargs...)

Scatter plot of model values (y) against reference values (x), with the
diagonal of perfect agreement (grey). The `nlabels` points with the largest
`deviation` are labeled with their `labels` (see `place_labels!`); `summary`
is printed in the lower right corner. With `log_axes = true`, both axes are
logarithmic. Values are shares and labeled as percentages, unless
`tick_label` formats the ticks otherwise (e.g. for counts, with `tick_step`).

`groups` (one integer from 1 to n per point) colors the points by group
(`group_colors`); without groups, all points have the model color. The
legend of the groups is drawn separately (`legend_strip`), so that panels
with and without groups keep the same size. `markersize` and `markeralpha`
(default 0.8 without groups, 0.9 with groups) set the points, e.g. small and
transparent for many points. The groups are drawn in their order, the last
on top, or with `reverse_groups = true` the first on top.
"""
function plot_agreement(reference, model;
        labels = String[],
        nlabels = 3,
        deviation = abs.(model .- reference),
        groups = nothing,
        summary = "",
        log_axes = false,
        ticks = nothing,
        tick_step = nothing,
        tick_label = nothing,
        reference_label = "Reference",
        show_ylabel = true,
        title = "",
        panel_width_pt = 145,
        panel_height_pt = 150,
        markersize = 4,
        markeralpha = nothing,
        reverse_groups = false,
        kwargs...)

    vals = [reference; model]
    # linear axes start at 0/0
    lo = log_axes ? 0.8 * minimum(vals) : 0.0
    hi = log_axes ? 1.25 * maximum(vals) : 1.05 * maximum(vals)
    step = something(tick_step, hi - lo > 0.6 ? 0.2 : hi - lo > 0.3 ? 0.1 : 0.05)
    ticks = isnothing(ticks) ? collect((step * ceil(lo / step)):step:hi) : filter(t -> lo <= t <= hi, ticks)
    ticklabels = isnothing(tick_label) ? percent_label.(ticks, digits = any(t -> 0 < t < 0.01, ticks) ? 1 : 0) :
        tick_label.(ticks)

    # positions on the axes as fractions (0 to 1), for label placement
    frac(v) = log_axes ? log(v / lo) / log(hi / lo) : (v - lo) / (hi - lo)
    unfrac(f) = log_axes ? lo * (hi / lo)^f : lo + f * (hi - lo)

    p = plot(
        title = title,
        titlelocation = :left,
        titlefontsize = PV_HEADING_FONTSIZE,
        xlims = (lo, hi),
        ylims = (lo, hi),
        xscale = log_axes ? :log10 : :identity,
        yscale = log_axes ? :log10 : :identity,
        xticks = (ticks, ticklabels),
        yticks = (ticks, ticklabels),
        xlabel = reference_label,
        ylabel = show_ylabel ? PV_MODEL_LABEL : "",
        aspect_ratio = :equal,               # same limits on both axes: always square
        grid = true,
        gridalpha = 0.15,
        tick_direction = :in,
        legend = false,
        tickfontsize = PV_SMALL_FONTSIZE,
        guidefontsize = PV_FONTSIZE,
        fontfamily = PV_FONT;
        kwargs...)

    plot!(p, [lo, hi], [lo, hi], linecolor = PV_LINE_COLOR, linewidth = 1.0, label = "")
    if isnothing(groups)
        scatter!(p, reference, model,
            markercolor = PV_MODEL_COLOR, markerstrokewidth = 0, markersize = markersize,
            markeralpha = something(markeralpha, 0.8), label = "")
    else
        colors = group_colors(maximum(groups))
        for g in (reverse_groups ? (maximum(groups):-1:1) : (1:maximum(groups)))
            sel = groups .== g
            scatter!(p, reference[sel], model[sel],
                markercolor = colors[g], markerstrokewidth = 0, markersize = markersize,
                markeralpha = something(markeralpha, 0.9), label = "")
        end
    end

    # the summary text (lower right corner) must not be covered by labels
    summary_width = 0.6 * PV_SMALL_FONTSIZE * length(summary) / panel_width_pt
    blocked = isempty(summary) ? [] : [(0.95 - summary_width, 1.0, 0.0, 0.03 + 1.3 * PV_SMALL_FONTSIZE / panel_height_pt)]
    place_labels!(p, frac.(reference), frac.(model), labels, deviation, unfrac;
        nlabels = nlabels, blocked = blocked, panel_width_pt = panel_width_pt, panel_height_pt = panel_height_pt)

    isempty(summary) || annotate!(p, unfrac(0.97), unfrac(0.03),
        pv_text(summary, PV_SMALL_FONTSIZE, :right, :bottom))
    return p
end

"""
    place_labels!(p, xf, yf, labels, deviation, unfrac;
        nlabels = 3, min_visible = 0.03, blocked = [], panel_width_pt = 145,
        panel_height_pt = 150, fontsize = PV_SMALL_FONTSIZE)

Labels the `nlabels` points with the largest `deviation` in an agreement
plot, but only those whose distance from the diagonal is visible (at least
`min_visible` of the axis length). Positions `xf`, `yf` are fractions of the axes (0 to 1), `unfrac`
converts them back to data coordinates. `blocked` lists further areas
(x0, x1, y0, y1 as fractions) that labels must not cover, such as the
legend or the summary text.

Text sizes are estimated from the number of characters and the font size
relative to the panel size (`panel_width_pt`, `panel_height_pt`). Each label
tries several positions: next to its point (to the right in the left half
of the panel, to the left in the right half), then moved away from the
diagonal in steps, then on the other side, then moved toward the diagonal.
The first position that stays inside the panel, on the same side of the
diagonal as its point, and does not cover a point, another label or a
blocked area is used; if there is none, the point stays unlabeled. Labels
that are not directly next to their point get a faint line to it, which
cannot cross the diagonal either, since label and point are on the same side.
"""
function place_labels!(p, xf, yf, labels, deviation, unfrac;
        nlabels = 3, min_visible = 0.03, blocked = [], panel_width_pt = 145, panel_height_pt = 150,
        fontsize = PV_SMALL_FONTSIZE)

    pad, inner, step = 0.025, 0.02, 0.07
    width(s) = 0.48 * fontsize * length(s) / panel_width_pt
    height = 1.6 * fontsize / panel_height_pt      # line height with some spacing
    dot = 0.015                                 # radius of a point
    taken = Tuple{Float64, Float64, Float64, Float64}[blocked...]

    overlaps(a, b) = a[1] < b[2] && b[1] < a[2] && a[3] < b[4] && b[3] < a[4]
    # a box lies entirely above the diagonal if its lower right corner does,
    # entirely below it if its upper left corner does
    same_side(box, own) = yf[own] > xf[own] ? box[3] - box[2] > dot : box[1] - box[4] > dot
    function free(box, own)
        box[1] >= inner && box[2] <= 1 - inner && box[3] >= inner && box[4] <= 1 - inner || return false
        same_side(box, own) || return false
        any(t -> overlaps(box, t), taken) && return false
        for k in eachindex(xf)
            k == own && continue
            overlaps(box, (xf[k] - dot, xf[k] + dot, yf[k] - dot, yf[k] + dot)) && return false
        end
        return true
    end

    # the largest deviations, but only points whose distance from the diagonal is visible
    idx = filter(i -> abs(yf[i] - xf[i]) >= min_visible,
        collect(partialsortperm(deviation, 1:min(nlabels, length(deviation)), rev = true)))
    for i in idx
        s, w = labels[i], width(labels[i])
        away = yf[i] < xf[i] ? -1 : 1              # away from the diagonal
        preferred = xf[i] < 0.5                      # text to the right of the point
        choice = nothing
        # away from the diagonal first, then toward it
        for dir in (away, -away), right in (preferred, !preferred), k in (dir == away ? (0:9) : (1:9))
            y = yf[i] + dir * k * step
            x = xf[i] + (right ? pad : -pad)
            box = right ? (x, x + w, y - height / 2, y + height / 2) : (x - w, x, y - height / 2, y + height / 2)
            if free(box, i)
                choice = (x = x, y = y, right = right, box = box, leader = k > 0 || right != preferred)
                break
            end
        end
        # no free position: leave the point unlabeled rather than covering something
        isnothing(choice) && continue
        push!(taken, choice.box)

        if choice.leader
            # from the near edge of the text to just before the point
            x0 = clamp(xf[i], choice.box[1], choice.box[2])
            y0 = clamp(yf[i], choice.box[3], choice.box[4])
            dx, dy = xf[i] - x0, yf[i] - y0
            r = hypot(dx, dy)
            f = r > 0 ? max(0.0, (r - 1.2dot) / r) : 0.0
            plot!(p, unfrac.([x0, x0 + f * dx]), unfrac.([y0, y0 + f * dy]),
                linecolor = :gray, linealpha = 0.5, linewidth = 0.6, label = "")
        end
        annotate!(p, unfrac(choice.x), unfrac(choice.y), pv_text(s, fontsize, choice.right ? :left : :right))
    end
    return p
end

"""
    group_colors(n)

`n` colors for groups such as age groups: the GEMS default palette.
"""
group_colors(n) = gemscolors(n)

"""
    ACTIVITY_INDICATORS

Indicators of activity, each a share of a federal state's population:
- employed: Zensus "Erwerbstätige" (including working students and
  apprentices) / GEMS persons with a workplace; by place of residence,
  shares of the persons in the Zensus table (private households)
- daycare: children in daycare (22541-0112, without after-school care) /
  GEMS persons in daycare settings (`is_daycare`)
- pupils: pupils at general and vocational schools (school statistics
  21111-0011 and 21121-0010; primary, special and secondary/vocational
  education settings) / GEMS persons in these settings
- daycare and pupils are counted by location of the institution,
  as shares of the population (Destatis 31.12.2019, GEMS population); the
  education settings of GEMS are identified by `EDUCATION_SETTINGS`
- not working (15–64): Zensus "Erwerbslose" and "Sonstige" / GEMS
  `is_not_working`; by place of residence, as for the employed
- retired (65+): Zensus "Empfänger/-innen von Ruhegehalt/Kapitalerträgen" /
  GEMS `is_retired`; by place of residence, as for the employed

University students are not an indicator: the official statistics count
them at the location of the university, but GEMS assigns students to
universities regardless of distance, so the shares per state reflect this
assignment rather than enrollment.

The indicators are not mutually exclusive in the official data (apprentices
and working students are both employed and pupils), so they do not add up,
and no remainder is derived from them. In GEMS, apprentices and working
students are pupils or students only.
"""
const ACTIVITY_INDICATORS = [
    (key = "employed", label = "Employed"),
    (key = "daycare", label = "Daycare"),
    (key = "pupils", label = "Pupils"),
    (key = "not_working", label = "Not working (15–64)"),
    (key = "retired", label = "Retired (65+)"),
]

"""
    activity_indicators(data)

Shares of `ACTIVITY_INDICATORS` per federal state for the model and the
reference, from the `activity_status`, `education` and `age_sex` tables of
`validation_data`, as dictionaries indexed by (state, key).
"""
function activity_indicators(data)
    out = (model = Dict{Tuple{Int, String}, Float64}(), reference = Dict{Tuple{Int, String}, Float64}())
    # education indicators and the education settings (`EDUCATION_SETTINGS`) they count
    education_indicators = ("daycare" => ("daycare",), "pupils" => ("primary", "special", "secondary"))
    for s in 1:16, (col, d) in ((:model, out.model), (:reference, out.reference))
        act = data.activity_status[data.activity_status.state .== s, :]
        persons = sum(act[!, col])
        share(k) = sum(act[act.category .== k, col]) / persons
        d[(s, "employed")] = share("employed")
        d[(s, "not_working")] = share("not_working")
        d[(s, "retired")] = share("retired")
        population = sum(data.age_sex[data.age_sex.state .== s, col])
        for (key, categories) in education_indicators
            edu = data.education[(data.education.state .== s) .& in.(data.education.category, Ref(categories)), :]
            d[(s, key)] = sum(edu[!, col]) / population
        end
    end
    return out
end

"""
    state_agreement_panels(data; first_letter = 'C')

Agreement plots of the federal states: share of the national population,
age structure (10-year groups), household sizes, household family types,
household senior status and activity indicators (`ACTIVITY_INDICATORS`),
each GEMS against the reference (axes labeled with the source names
of `reference_labels(data)`). Returns the panels and, separately, the
contents of their color legends (`nothing` for the population panel): the
group labels and a title.
"""
function state_agreement_panels(data; first_letter = 'C')
    letter(k) = string(first_letter + k)
    sources = reference_labels(data)

    # share of the national population
    pop = combine(groupby(data.age_sex, :state), :model => sum => :model, :reference => sum => :reference)
    sort!(pop, :state)
    pop_m, pop_r = pop.model ./ sum(pop.model), pop.reference ./ sum(pop.reference)
    rel = pop_m ./ pop_r .- 1
    p_pop = plot_agreement(pop_r, pop_m,
        labels = STATE_NAMES_EN,
        deviation = abs.(rel),
        summary = "mean |Δ| = $(round(100 * mean(abs.(rel)), digits = 1))%",
        reference_label = sources.age_sex,
        title = "$(letter(0))) Population")

    # shares of categories per state, flattened to one point per state and
    # category, colored by category (labels name the state; the color names the category)
    function category_panel(df, key, categories, reference_label, title; nlabels = 3)
        m, r = state_shares(df, key, categories)
        lab = [STATE_NAMES_EN[s] for s in 1:16, _ in eachindex(categories)]
        grp = [c for _ in 1:16, c in eachindex(categories)]
        return plot_agreement(vec(r), vec(m),
            labels = vec(lab),
            nlabels = nlabels,
            groups = vec(grp),
            summary = "mean |Δ| = $(round(100 * mean(abs.(m .- r)), digits = 2)) pp",
            reference_label = reference_label,
            show_ylabel = false,
            title = title)
    end

    age = transform(data.age_sex, :age => (a -> min.(a .÷ 10, 8)) => :group)
    # decades (the oldest group is open-ended); the caption names the 10-year groups
    age_labels = ["<10"; [string(10g, "s") for g in 1:7]; "80+"]
    p_age = category_panel(age, :group, 0:8, sources.age_sex, "$(letter(1))) Age groups", nlabels = 4)

    hh_sizes = 1:maximum(data.household_sizes.size)
    hh_labels = category_label.(hh_sizes, last(hh_sizes))
    p_hh = category_panel(data.household_sizes, :size, hh_sizes, sources.household_sizes,
        "$(letter(3))) Household sizes")

    # family types: shares of all households, "other" = no core family type (see `FAMILY_TYPES`)
    p_types = category_panel(family_type_counts(data.household_types), :category,
        getfield.(FAMILY_TYPES, :key), sources.household_types, "$(letter(4))) Household family types")

    seniors = getfield.(SENIOR_STATUS, :key)
    p_seniors = category_panel(data.senior_status, :category, seniors, sources.senior_status,
        "$(letter(5))) Household senior status")

    # activity indicators: separate shares of the population, each against its own source
    # (not mutually exclusive in the official data, see `activity_indicators`)
    ind = activity_indicators(data)
    keys_ind = getfield.(ACTIVITY_INDICATORS, :key)
    im = [ind.model[s, k] for s in 1:16, k in keys_ind]
    ir = [ind.reference[s, k] for s in 1:16, k in keys_ind]
    p_activity = plot_agreement(vec(ir), vec(im),
        labels = vec([STATE_NAMES_EN[s] for s in 1:16, _ in keys_ind]),
        groups = vec([j for _ in 1:16, j in eachindex(keys_ind)]),
        summary = "mean |Δ| = $(round(100 * mean(abs.(im .- ir)), digits = 2)) pp",
        reference_label = "Census / Destatis",
        show_ylabel = false,
        title = "$(letter(2))) Employment and education")

    # persons (population, age, activity), then households (sizes, family types, senior status)
    panels = [p_pop, p_age, p_activity, p_hh, p_types, p_seniors]
    legends = [
        nothing,
        (labels = age_labels, title = "Age"),
        (labels = getfield.(ACTIVITY_INDICATORS, :label), title = ""),
        (labels = hh_labels, title = "Persons"),
        (labels = getfield.(FAMILY_TYPES, :label), title = ""),
        (labels = getfield.(SENIOR_STATUS, :label), title = ""),
    ]
    return panels, legends
end


# 4. FIGURE

"""
    section_header(txt; width_in = PV_WIDTH_IN - 4 / 25.4)

Empty subplot carrying the heading of a section of the figure: centered as
the section titles of the other figures of the paper (e.g.
`plot_paper_figure_1`), bold and larger than the panel titles, between two
thin lines that set it apart from them. `width_in` is the width of the
heading's plot area.
"""
function section_header(txt; width_in = PV_WIDTH_IN - 4 / 25.4)
    p = plot(framestyle = :none, ticks = false, legend = false, grid = false,
        xlims = (0, 1), ylims = (0, 1), margin = 0mm, top_margin = -2mm, bottom_margin = -2mm,
        fontfamily = PV_FONT)
    half = (0.5 * text_width_in(txt, PV_SECTION_FONTSIZE, bold = true) + 0.12) / width_in
    plot!(p, [0, 0.5 - half, NaN, 0.5 + half, 1], fill(0.5, 5), linecolor = PV_LINE_COLOR, linewidth = 0.6)
    return annotate!(p, 0.5, 0.5, Plots.text(txt, PV_SECTION_FONTSIZE, PV_FONT_BOLD, :center, :vcenter))
end


# AXIS TEXT, LEGENDS AND LAYOUT OF THE FIGURE
#
# GR places tick labels and axis labels at distances relative to the largest
# dimension of the figure, which leaves wide gaps in a tall figure, and
# centers rotated tick labels on their length. The figure therefore draws
# them itself (`axis_text!`), at fixed distances in inches, using the sizes
# of the plot areas that `fit_validation_layout` measures.

# common angle of rotated tick labels (degrees)
const PV_TICK_ROTATION = 35
# distance of tick labels from the axis, and of axis labels from the tick labels (inches)
const PV_TICK_GAP_IN = 0.03
const PV_GUIDE_GAP_IN = 0.04
# distance below a panel title (and below the color legend of an agreement plot, inches)
const PV_TITLE_GAP_IN = 0.07
# legend entries: distance from the marker to its label, and between entries (inches)
const PV_MARKER_GAP_IN = 0.1
const PV_ENTRY_GAP_IN = 0.14
# estimated character widths of Times, as shares of the font size (digits and capitals are wider
# than the average lowercase letter, bold text by about 10 %)
const PV_CHAR_WIDTH = 0.42
char_width(c) = c == '%' ? 0.83 : isdigit(c) ? 0.5 : c in "MW" ? 0.89 : isuppercase(c) ? 0.66 : c == ' ' ? 0.25 :
    c in ".,:;()|/–-" ? 0.33 : c in "ijlftr" ? 0.3 : c in "mw" ? 0.72 : PV_CHAR_WIDTH
# name of the reference data in the figure (the sources are named in the caption)
const PV_REFERENCE_NAME = "Official statistics"

text_width_in(s, size; bold = false) = (bold ? 1.1 : 1.0) * sum(char_width, s; init = 0.0) * size / 72
line_height_in(size) = 1.2 * size / 72
"""
    title_text!(p, x, y, title, width_in; size = PV_HEADING_FONTSIZE, valign = :bottom)

Draws a panel title at the relative position (`x`, `y`) of `p`, left
aligned: its letter ("A)") in bold, as in the other figures of the paper,
followed by the rest in regular weight. `width_in` is the width of the plot
area, to place the rest after the letter.
"""
function title_text!(p, x, y, title, width_in; size = PV_HEADING_FONTSIZE, valign = :bottom)
    m = match(r"^([A-Z]\))\s*(.*)$", title)
    letter, rest = isnothing(m) ? ("", title) : (m.captures[1], m.captures[2])
    dx = 0.0
    if !isempty(letter)
        annotate!(p, (x, y), Plots.text(letter, size, PV_FONT_BOLD, :left, valign))
        dx = (text_width_in(letter, size, bold = true) + 0.3size / 72) / width_in
    end
    return annotate!(p, (x + dx, y), Plots.text(rest, size, PV_FONT, :left, valign))
end

"""
    panel_title!(p, width_in, height_in)

Replaces the title of the panel `p` by `title_text!` above its plot area
(`width_in` × `height_in` inches), at the distance `PV_TITLE_GAP_IN`.
Returns the room the title needs above the plot area (inches).
"""
function panel_title!(p, width_in, height_in)
    title = p[1][:title]
    plot!(p, title = "")
    title_text!(p, 0.0, 1 + PV_TITLE_GAP_IN / height_in, title, width_in)
    return PV_TITLE_GAP_IN + line_height_in(PV_HEADING_FONTSIZE)
end

"""
    axis_text!(p, width_in, height_in; guide_gap_in = PV_GUIDE_GAP_IN)

Replaces the tick labels and axis labels that GR would draw for the panel
`p` (a plot with one subplot whose ticks are given as values and labels) by
annotations at fixed distances from a plot area of `width_in` × `height_in`
inches: y tick labels right-aligned next to the axis, the y axis label next
to them (`guide_gap_in` beyond the estimated width of the tick labels), x
tick labels centered below their ticks or, if rotated, ending at their
ticks, and the x axis label below them. Returns the room the text needs
left of and below the plot area (inches).
"""
function axis_text!(p, width_in, height_in; guide_gap_in = PV_GUIDE_GAP_IN)
    sp = p[1]
    xaxis, yaxis = sp[:xaxis], sp[:yaxis]
    (xt, xl), (yt, yl) = xaxis[:ticks], yaxis[:ticks]
    xguide, yguide, rot = xaxis[:guide], yaxis[:guide], xaxis[:rotation]
    (x0, x1), (y0, y1) = xlims(p), ylims(p)
    # positions as fractions of the axis, linear or logarithmic
    scaled(axis, v) = axis[:scale] == :log10 ? log10(v) : v
    relx(v) = (scaled(xaxis, v) - scaled(xaxis, x0)) / (scaled(xaxis, x1) - scaled(xaxis, x0))
    rely(v) = (scaled(yaxis, v) - scaled(yaxis, y0)) / (scaled(yaxis, y1) - scaled(yaxis, y0))
    s, g = PV_SMALL_FONTSIZE, PV_FONTSIZE
    plot!(p, xticks = (xt, fill("", length(xt))), yticks = (yt, fill("", length(yt))),
        xguide = "", yguide = "")

    shown(ticks, labels, lo, hi) = [(v, l) for (v, l) in zip(ticks, labels) if lo <= v <= hi]
    ys, xs = shown(yt, yl, y0, y1), shown(xt, xl, x0, x1)

    # y axis: tick labels, then the axis label (rotated, bottom toward the axis)
    for (v, l) in ys
        annotate!(p, (-PV_TICK_GAP_IN / width_in, rely(v)), Plots.text(l, s, PV_FONT, :right, :vcenter))
    end
    left = PV_TICK_GAP_IN + maximum(((_, l),) -> text_width_in(l, s), ys; init = 0.0)
    if !isempty(yguide)
        x = -(left + guide_gap_in) / width_in
        annotate!(p, (x, 0.5), Plots.text(yguide, g, PV_FONT, :center, :bottom, rotation = 90))
        left += guide_gap_in + line_height_in(g)
    end

    # x axis: tick labels (rotated ones end at their tick), then the axis label
    h = line_height_in(s)
    if rot == 0
        for (v, l) in xs
            annotate!(p, (relx(v), -PV_TICK_GAP_IN / height_in), Plots.text(l, s, PV_FONT, :center, :top))
        end
        bottom = PV_TICK_GAP_IN + h
    else
        dy = PV_TICK_GAP_IN + 0.5h * cosd(rot)
        for (v, l) in xs
            annotate!(p, (relx(v), -dy / height_in), Plots.text(l, s, PV_FONT, :right, :vcenter, rotation = rot))
        end
        bottom = dy + maximum(((_, l),) -> text_width_in(l, s) * sind(rot) + 0.5h * cosd(rot), xs; init = 0.0)
    end
    if !isempty(xguide)
        y = -(bottom + PV_GUIDE_GAP_IN) / height_in
        annotate!(p, (0.5, y), Plots.text(xguide, g, PV_FONT, :center, :top))
        bottom += PV_GUIDE_GAP_IN + line_height_in(g)
    end
    return (left = left, bottom = bottom)
end

"""
    reserve!(p, room, figure_width_in, figure_height_in; cushion_mm = 0.5, bottom_cushion_mm = cushion_mm)

Sets the left and bottom margins of the panel `p` so that the room taken by
the text of `axis_text!` is kept free, plus a cushion (`cushion_mm` left,
`bottom_cushion_mm` below). Plots adds its own padding for ticks (2 mm,
plus 1 mm and 1 % of the figure size per axis with ticks), which is
subtracted.
"""
function reserve!(p, room, figure_width_in, figure_height_in; cushion_mm = 0.5, bottom_cushion_mm = cushion_mm)
    pad(fig_in) = 3.0 + 0.254fig_in
    plot!(p,
        left_margin = (25.4room.left - pad(figure_width_in) + cushion_mm) * mm,
        bottom_margin = (25.4room.bottom - pad(figure_height_in) + bottom_cushion_mm) * mm)
end

"""
    legend_flow(legend, width_in)

Positions of the entries of a color legend (`legend` from
`state_agreement_panels`: labels and an optional title) that flow like text
over a width of `width_in` inches: each entry (marker and label) follows the
previous one and moves to the next line when the line is full. Lines after
the first start where the entries of the first line start (after the
title). Returns the start of the entries (`x0`), the line and position
(inches) of each entry, and the number of lines.
"""
function legend_flow(legend, width_in)
    s = PV_SMALL_FONTSIZE
    x0 = isempty(legend.title) ? 0.0 : text_width_in(legend.title, s) + 0.08
    positions = Tuple{Int, Float64}[]
    line, x = 1, x0
    for l in legend.labels
        w = PV_MARKER_GAP_IN + text_width_in(l, s)
        if x > x0 && x + w > width_in
            line, x = line + 1, x0
        end
        push!(positions, (line, x))
        x += w + PV_ENTRY_GAP_IN
    end
    return (x0 = x0, positions = positions, lines = isempty(positions) ? 0 : last(positions)[1])
end

legend_lines(legend, width_in) = isnothing(legend) ? 0 : legend_flow(legend, width_in).lines

"""
    legend_strip(title, legend, width_in, height_in; flow_width_in = width_in)

Strip above an agreement plot: the panel `title` and, below it, the color
legend of its groups (`legend` from `state_agreement_panels`, or `nothing`):
an optional legend title followed by the group labels, flowing like text
(`legend_flow`) over `flow_width_in` inches, which may reach into the space
between the columns (the strip's plot area is `width_in` wide).
"""
function legend_strip(title, legend, width_in, height_in; flow_width_in = width_in)
    p = plot(framestyle = :none, ticks = false, legend = false, grid = false,
        xlims = (0, 1), ylims = (0, 1), margin = 0mm, top_margin = -2mm, bottom_margin = -2mm,
        fontfamily = PV_FONT)
    title_text!(p, 0.0, 1.0, title, width_in, valign = :top)
    isnothing(legend) && return p

    s = PV_SMALL_FONTSIZE
    line = line_height_in(s) / height_in
    top = 1 - (line_height_in(PV_HEADING_FONTSIZE) + PV_TITLE_GAP_IN) / height_in
    line_y(r) = top - (r - 0.5) * line
    flow = legend_flow(legend, flow_width_in)
    colors = group_colors(length(legend.labels))
    isempty(legend.title) || annotate!(p, 0, line_y(1), Plots.text(legend.title, s, PV_FONT, :left, :vcenter))
    for (k, (l, (r, x))) in enumerate(zip(legend.labels, flow.positions))
        # a marker, not a text symbol: GR draws text in the nearest color of its color table
        scatter!(p, [(x + 0.03) / width_in], [line_y(r)], markercolor = colors[k], markerstrokewidth = 0,
            markersize = 4)
        annotate!(p, (x + PV_MARKER_GAP_IN) / width_in, line_y(r), Plots.text(l, s, PV_FONT, :left, :vcenter))
    end
    return p
end

# height of a legend strip with `rows` rows of entries (inches)
strip_height_in(rows) = line_height_in(PV_HEADING_FONTSIZE) + PV_TITLE_GAP_IN +
    rows * line_height_in(PV_SMALL_FONTSIZE) + PV_TITLE_GAP_IN

"""
    reference_legend!(p, width_in)

Legend of the colors of the age pyramid and the bar charts, the reference
data (`PV_REFERENCE_NAME`) and GEMS, drawn next to each other from the left
edge of the subplot `p` (plot area `width_in` inches wide, x from 0 to 1).
Returns the width of the legend (inches).
"""
function reference_legend!(p, width_in)
    entries = ((PV_REFERENCE_NAME, PV_REFERENCE_COLOR, PV_REFERENCE_ALPHA), (PV_MODEL_LABEL, PV_MODEL_COLOR, 0.85))
    entry_width(label) = 0.17 + text_width_in(label, PV_FONTSIZE)
    total = sum(e -> entry_width(e[1]), entries) + 0.15 * (length(entries) - 1)
    x = 0.0
    for (label, color, alpha) in entries
        scatter!(p, [(x + 0.05) / width_in], [0.5], markershape = :rect, markersize = 6, markercolor = color,
            markeralpha = alpha, markerstrokewidth = 0)
        annotate!(p, (x + 0.17) / width_in, 0.5, Plots.text(label, PV_FONTSIZE, PV_FONT, :left, :vcenter))
        x += entry_width(label) + 0.15
    end
    return total
end

"""
    ValidationLayout

Sizes of the figure of the population validation (`validation_figure`), in
inches unless noted:
- `header_in`: section headers; `bar_row_in`: each of the two rows of bar
  charts;
  `state_row_in`: each row of agreement plots (without its legend strip)
- `pyramid_width`: share of the figure width taken by the age pyramid;
  `state_widths`: widths of the three columns of agreement plots (shares)
- `pyramid_left_mm`, `bars_right_mm`: extra margins (mm) that align the
  pyramid and the rightmost bar charts with the agreement plots;
  `last_bar_bottom_mm`: extra bottom margin of the lower bar charts, which
  puts their x axis at the x axis of the pyramid
- `bar_plot_in`: target height of the plot areas of the bar charts
- `state_gap_mm`, `state_row_gap_mm`: extra space (mm) between the columns
  and between the rows of agreement plots
- `state_scale`: size of the agreement plots relative to the largest size
  the width allows; the width they leave free goes into the space between
  the columns (`state_extra_gap_mm`, fitted)
- measured plot areas (`fit_validation_layout`), used to place the axis
  text: `pyramid_size`, `bar_width_in`, `state_side_in`
"""
Base.@kwdef struct ValidationLayout
    header_in::Float64 = 0.2
    bar_row_in::Float64 = 1.7
    state_row_in::Float64 = 2.4
    pyramid_width::Float64 = 0.3
    state_widths::Vector{Float64} = fill(1 / 3, 3)
    pyramid_left_mm::Float64 = 0.0
    bars_right_mm::Float64 = 1.0
    last_bar_bottom_mm::Float64 = 0.0
    bar_plot_in::Float64 = 0.7
    pyramid_size::Tuple{Float64, Float64} = (1.8, 2.8)
    bar_width_in::Float64 = 1.0
    state_side_in::Float64 = 1.9
    state_gap_mm::Float64 = 5.0
    state_scale::Float64 = 0.95
    state_extra_gap_mm::Float64 = 0.0
    state_row_gap_mm::Float64 = 1.5
end

with_layout(l::ValidationLayout; changes...) = ValidationLayout(;
    (f => get(changes, f, getfield(l, f)) for f in fieldnames(ValidationLayout))...)

"""
    build_validation_figure(data, layout::ValidationLayout; dpi = PV_DPI)

Composes the figure of `validation_figure` with the sizes of `layout`.
Each panel carries an empty series named after its letter (`panel_mark`),
so that `plot_areas` finds it in the figure.
"""
function build_validation_figure(data, layout::ValidationLayout; dpi = PV_DPI)
    l = layout

    # the agreement plots first, as their axis text and legend strips set the height of the figure;
    # the first column labels the y axis, the last row the x axis
    panels, legends = state_agreement_panels(data, first_letter = 'H')
    panel_right_mm(k) = mod1(k, 3) < 3 ? l.state_gap_mm + l.state_extra_gap_mm : 1.0
    titles = [p[1][:title] for p in panels]
    rooms = map(enumerate(panels)) do (k, p)
        plot!(p, title = "", xlabel = k > 3 ? PV_REFERENCE_NAME : "", ylabel = mod1(k, 3) == 1 ? PV_MODEL_LABEL : "",
            # square by the layout (`fit_validation_layout`); a fixed aspect ratio would shrink
            # the plot areas while the layout is fitted
            aspect_ratio = :auto, top_margin = -2mm,
            right_margin = panel_right_mm(k) * mm)
        axis_text!(p, l.state_side_in, l.state_side_in)
    end
    # legend strips in rows of their own, so that their plot areas reach across the space between
    # the columns (markers outside a plot area are clipped). Plots divides a row by the widths of
    # the plot areas, so each strip gets the width of its panel's plot area plus the space right
    # of it (Plots adds 2 mm to each side of a plot area): the cells of both rows coincide. Each
    # strip starts where its panel's plot area starts (`reserve!` leaves `room.left` plus its
    # cushion of 0.5 mm left of it).
    strip_left_in(k) = rooms[k].left + 0.5 / 25.4
    strip_width_in(k) = l.state_side_in + (2 + panel_right_mm(k)) / 25.4
    strip_rows = [maximum(k -> legend_lines(legends[k], strip_width_in(k) - 0.05), r .+ (1:3)) for r in (0, 3)]
    strips_in = strip_height_in.(strip_rows)
    rows = [strips_in[1], l.state_row_in, strips_in[2], l.state_row_in]
    heights = [l.header_in, 2l.bar_row_in, l.header_in, sum(rows)]
    width_in, height_in = PV_WIDTH_IN, sum(heights)

    # Germany: age pyramid with the reference legend below it, bar charts of households and settings
    pyramid = plot_age_pyramid(data.age_sex, title = "A) Age and sex", show_top_group = true, age_label = "Age", xlabel = "Share of population",
        xtick_step = 0.004, legend = false)
    households = [
        plot_household_size_panel(data.household_sizes, title = "B) Household size"),
        plot_household_category_panel(family_type_counts(data.household_types), FAMILY_TYPES,
            title = "C) Family type"),
        plot_household_category_panel(data.senior_status, SENIOR_STATUS, title = "D) Senior status"),
    ]
    foreach(p -> plot!(p, ylabel = "Share of households"), households)
    bars = [households; setting_size_panels(data, letters = ('E', 'F', 'G'))]
    for p in bars
        rotated = p[1][:xaxis][:rotation] != 0 || p === households[2] || p === households[3]
        plot!(p, legend = false, xrotation = rotated ? PV_TICK_ROTATION : 0)
    end
    # titles and axis text drawn by the figure; Plots adds 2 mm above a plot area without title
    title_margin(p, w, h) = plot!(p, top_margin = (25.4panel_title!(p, w, h) - 2 + 0.5) * mm)
    title_margin(pyramid, l.pyramid_size...)
    # below the axis label of the pyramid: the legend of the reference data and GEMS (an inset, see
    # below); the room is reserved here, so the pyramid keeps the x axis of the lower bar charts
    pyramid_room = axis_text!(pyramid, l.pyramid_size...)
    legend_top_in = pyramid_room.bottom + PV_TITLE_GAP_IN
    legend_height_in = line_height_in(PV_FONTSIZE) + 0.04
    reserve!(pyramid, (left = pyramid_room.left, bottom = legend_top_in + legend_height_in), width_in, height_in)
    plot!(pyramid, left_margin = pyramid[1][:left_margin] + l.pyramid_left_mm * mm)
    for (k, p) in enumerate(bars)
        title_margin(p, l.bar_width_in, l.bar_plot_in)
        reserve!(p, axis_text!(p, l.bar_width_in, l.bar_plot_in), width_in, height_in)
        k > 3 && plot!(p, bottom_margin = p[1][:bottom_margin] + l.last_bar_bottom_mm * mm)
        k in (3, 6) && plot!(p, right_margin = l.bars_right_mm * mm)
    end
    mark_panel!(pyramid, "A")
    foreach(((k, p),) -> mark_panel!(p, string('B' + k - 1)), enumerate(bars))

    # the pyramid spans both rows of bar charts (Plots gives it the bottom padding of the lower row)
    top = plot(pyramid, plot(bars..., layout = grid(2, 3)),
        layout = grid(1, 2, widths = [l.pyramid_width, 1 - l.pyramid_width]))

    # State level: agreement plots in two rows of three, each with a strip of title and color
    # legend above it
    strips = map(enumerate(panels)) do (k, p)
        reserve!(p, rooms[k], width_in, height_in)
        k <= 3 && plot!(p, bottom_margin = p[1][:bottom_margin] + l.state_row_gap_mm * mm)
        mark_panel!(p, string('H' + k - 1))
        s = legend_strip(titles[k], legends[k], strip_width_in(k), strips_in[k <= 3 ? 1 : 2],
            flow_width_in = strip_width_in(k) - 0.05)
        # Plots adds 2 mm on each side of a plot area
        plot!(s, left_margin = (25.4strip_left_in(k) - 2) * mm, right_margin = -2mm)
    end
    row(v) = plot(v..., layout = grid(1, 3, widths = l.state_widths))
    strip_widths = strip_width_in.(1:3) ./ sum(strip_width_in.(1:3))
    strip_row(v) = plot(v..., layout = grid(1, 3, widths = strip_widths))
    bottom = plot(strip_row(strips[1:3]), row(panels[1:3]), strip_row(strips[4:6]), row(panels[4:6]),
        layout = grid(4, 1, heights = rows ./ sum(rows)))

    figure = plot(section_header("National level"), top, section_header("State level"), bottom,
        layout = grid(4, 1, heights = heights ./ sum(heights)),
        size = (100width_in, 100height_in),
        dpi = dpi,
        fontfamily = PV_FONT)

    # legend of the reference data and GEMS below the pyramid: an inset, positioned relative to the
    # plot area of the pyramid and therefore outside the alignment of the grid
    a = findfirst(sp -> any(sr -> sr[:label] == panel_mark("A"), sp.series_list), figure.subplots)
    k = length(figure.subplots) + 1
    h = l.pyramid_size[2]
    plot!(figure, inset = (a, bbox(0.0, 1 + legend_top_in / h, 1.0, legend_height_in / h)), subplot = k,
        framestyle = :none, ticks = false, legend = false, grid = false, xlims = (0, 1), ylims = (0, 1))
    reference_legend!(figure[k], l.pyramid_size[1])
    return figure
end

# empty series that names a panel by its letter (subplots are copied when plots are combined)
panel_mark(letter) = "panel $letter"
mark_panel!(p, letter) = plot!(p, Float64[], Float64[], label = panel_mark(letter), primary = false)

"""
    plot_areas(figure)

Plot areas (in inches) of the panels of `figure` (marked by `mark_panel!`),
by panel letter: left, right, top, width, height.
"""
function plot_areas(figure)
    Plots.prepare_output(figure)
    inch(x) = x.value / 25.4
    areas = Dict{String, NamedTuple}()
    for sp in figure.subplots, series in sp.series_list
        m = match(r"^panel ([A-Z])$", string(series[:label]))
        isnothing(m) && continue
        a = Plots.plotarea(sp)
        areas[m.captures[1]] = (left = inch(Plots.left(a)), right = inch(Plots.right(a)), top = inch(Plots.top(a)),
            width = inch(Plots.width(a)), height = inch(Plots.height(a)))
    end
    return areas
end

"""
    fit_validation_layout(data; layout = ValidationLayout(), iterations = 8, tolerance_in = 0.005)

Adjusts `layout` to the plot areas that Plots computes for the figure
(`plot_areas`), until all deviations are below `tolerance_in`:
- the three columns of agreement plots have plot areas of the same width,
  and their rows are as high as these are wide (square plots: the axes have
  the same limits)
- the plot area of the age pyramid starts where those of the agreement
  plots in the first column start
- the plot areas of the rightmost bar charts end where those of the
  agreement plots in the last column end
- the plot areas of the bar charts are `layout.bar_plot_in` high
- the x axis of the lower bar charts is at the height of the x axis of the
  age pyramid
- the sizes used to place the axis text are the measured ones
"""
function fit_validation_layout(data; layout = ValidationLayout(), iterations = 8, tolerance_in = 0.005)
    l = layout
    for _ in 1:iterations
        a = plot_areas(build_validation_figure(data, l))
        widths = [a[k].width for k in ("H", "I", "J")]
        # padding of the columns without the fitted extra gaps; the plots take `state_scale` of the
        # largest side the width allows, the rest goes into the two gaps between the columns
        extra = [1, 1, 0] .* l.state_extra_gap_mm ./ 25.4
        padding = l.state_widths .* PV_WIDTH_IN .- widths .- extra
        side = l.state_scale * (PV_WIDTH_IN - sum(padding)) / 3
        extra_gap_in = (PV_WIDTH_IN - sum(padding) - 3side) / 2
        shift_left = a["H"].left - a["A"].left
        shift_right = a["J"].right - max(a["D"].right, a["G"].right)
        # x axis of the lower bar charts at the x axis of the pyramid
        shift_bottom = (a["E"].top + a["E"].height) - (a["A"].top + a["A"].height)
        pyramid = (a["A"].width, a["A"].height)
        deviation = maximum(abs, [widths .- side; a["H"].height - side; a["B"].height - l.bar_plot_in;
            extra_gap_in - l.state_extra_gap_mm / 25.4;
            shift_left; shift_right; shift_bottom; collect(pyramid .- l.pyramid_size); a["B"].width - l.bar_width_in;
            side - l.state_side_in])
        deviation < tolerance_in && return l
        l = with_layout(l,
            bar_row_in = l.bar_row_in - a["B"].height + l.bar_plot_in,
            state_widths = (padding .+ [1, 1, 0] .* extra_gap_in .+ side) ./ PV_WIDTH_IN,
            state_extra_gap_mm = 25.4extra_gap_in,
            state_row_in = l.state_row_in - a["H"].height + side,
            pyramid_left_mm = l.pyramid_left_mm + 25.4shift_left,
            bars_right_mm = l.bars_right_mm - 25.4shift_right,
            last_bar_bottom_mm = l.last_bar_bottom_mm + 25.4shift_bottom,
            pyramid_size = pyramid,
            bar_width_in = a["B"].width,
            state_side_in = side)
    end
    @warn "The layout of the validation figure did not converge in $iterations iterations."
    return l
end

"""
    validation_figure(data; layout = fit_validation_layout(data), dpi = PV_DPI)

Figure of the population validation at the print width of a PLOS figure
(`PV_WIDTH_IN`):
- "National level": the age pyramid (A) with the legend of the reference data and
  GEMS below it; to its right the household panels (B–D: size, family type,
  senior status; shares of households) and the setting panels (E–G,
  `setting_size_panels`)
- "State level": the agreement plots (`state_agreement_panels`, H–M) in
  two rows of three, each with its title and color legend above it; the
  first column labels the y axis (GEMS), the last row the x axis

The reference data are named `PV_REFERENCE_NAME`; the sources belong in the
caption. The sizes come from `fit_validation_layout`: square agreement plots
that use the full width, a pyramid that starts at their left edge and bar
charts that end at their right edge. The height of the figure follows.
"""
validation_figure(data; layout = fit_validation_layout(data), dpi = PV_DPI) =
    build_validation_figure(data, layout, dpi = dpi)

"""
    population_validation(;
        folder = nothing,
        folder_tag = nothing,
        sim = nothing,
        reference = :official2019,
        data = nothing,
        individual_data_df = nothing)

Compares the GEMS population with official statistics (`reference`, see
`VALIDATION_REFERENCES`) and saves the figure (PNG) and its data (CSV) to
`folder`. Pass `data` (from `validation_data` or
`load_validation_data`) to redraw the figure without loading the
population; otherwise `sim` (default: the German population) is used.
Pass `individual_data_df` to reuse an already computed
`individual_data(sim)`.

If `folder` is `nothing`, a new result folder is created, with
`folder_tag` appended to its timestamp-based name (see `result_folder`).
"""
function population_validation(;
        folder = nothing,
        folder_tag = nothing,
        sim = nothing,
        reference = :official2019,
        data = nothing,
        individual_data_df = nothing)

    isnothing(folder) ? folder = result_folder(folder_tag) : mkpath(folder)

    if isnothing(data)
        sim = isnothing(sim) ? Simulation(population = "DE", seed = 0) : sim
        data = validation_data(sim, reference = reference,
            individual_data_df = something(individual_data_df, individual_data(sim)))
    end
    save_validation_data(data, folder)

    p = validation_figure(data)
    png(p, joinpath(folder, "population_validation.png"))
    return p
end
