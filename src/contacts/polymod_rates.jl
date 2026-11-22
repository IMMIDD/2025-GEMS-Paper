###
### CONTACT DATA: DAILY CONTACTS PER GEMS SETTING LAYER FROM POLYMOD
###

export contact_layer_rates
export polymod_rate_table
export contact_analyses


# Each POLYMOD contact is assigned to one GEMS setting layer, and the contacts
# of a layer are averaged over the participants who are members of that layer
# in GEMS (GEMS samples `contactparameter` contacts per member and day). The
# result is a daily average: GEMS does not distinguish days of the week, so
# weekdays count 5/7 and weekend days 2/7. Within each day type, participants
# are weighted by age to the members of the layer in the GEMS population.


"""
    CONTACT_LAYERS

GEMS setting layers that POLYMOD contacts are assigned to.
"""
const CONTACT_LAYERS = [:household, :school, :work, :municipality]

"""
    LAYER_SETTING_TYPES

GEMS setting types of each layer. The contact parameters of a layer's setting
types add up to the contacts per member and day of that layer.
"""
const LAYER_SETTING_TYPES = Dict(
    :household => ["Household"],
    :school => ["SchoolClass", "SchoolYear", "School", "SchoolComplex"],
    :work => ["Office", "Department", "Workplace", "WorkplaceSite"],
    :municipality => ["Municipality"])

"""
    CONTACT_AGE_LIMITS

Lower limits of the age groups used for weighting (as in Table 1 of
Mossong et al. 2008).
"""
const CONTACT_AGE_LIMITS = [0, 5, 10, 15, 20, 30, 40, 50, 60, 70]

"""
    WEEKDAY_SHARE

Share of weekdays (Monday to Friday) in the daily average.
"""
const WEEKDAY_SHARE = 5 / 7

# index of the age group of `age` in `CONTACT_AGE_LIMITS`
age_group(age) = ismissing(age) ? missing : searchsortedlast(CONTACT_AGE_LIMITS, age)


## PARTICIPANTS

"""
    polymod_status(participants)

Status of each POLYMOD participant (rows as in `participants`):

- `in_education`: member of a GEMS education setting (daycare to university).
  Children's diaries (`type` 3): attends day care or school (`child_care` "Y"),
  except day care at home (`child_care_detail` 1, a childminder, no GEMS
  setting). The detail is missing in some countries (all German diaries), then
  `child_care` alone decides. Own diaries (`type` 1 or 2): in full-time
  education (`part_occupation` 5).
- `working`: own diary and working (`part_occupation` 1). In children's diaries,
  the occupation is the parent's.
- `age_group` (see `CONTACT_AGE_LIMITS`), `weekend` (`missing` if the day is
  unknown) and `hh_size`.
"""
function polymod_status(participants)
    child = coalesce.(participants.type .== 3, false)
    own = coalesce.(in.(participants.type, Ref((1, 2))), false)
    occupation = participants.part_occupation
    attends = coalesce.(participants.child_care .== "Y", false) .&
        .!coalesce.(participants.child_care_detail .== 1, false)
    return DataFrame(
        part_id = participants.part_id,
        in_education = (child .& attends) .| (own .& coalesce.(occupation .== 5, false)),
        working = own .& coalesce.(occupation .== 1, false),
        age_group = age_group.(participants.part_age),
        weekend = [ismissing(d) ? missing : d in (0, 6) for d in participants.dayofweek],
        hh_size = participants.hh_size)
end


## CONTACTS

"""
    contact_layer(home, school, work, in_education, working)

GEMS layer of one contact from its location flags and the participant's status.
Home contacts go to the household. Contacts at school or work go to the layer the
participant is a member of (school contacts of a working, not studying,
participant, e.g. a teacher, go to work). All other contacts (transport, leisure,
other places, no location, or school/work without membership) go to the
municipality.
"""
function contact_layer(home::Bool, school::Bool, work::Bool, in_education::Bool, working::Bool)
    home && return :household
    school && in_education && return :school
    work && working && return :work
    school && working && return :work
    return :municipality
end

"""
    contact_layers(contacts, status)

GEMS layer of each contact (see `contact_layer`).
"""
function contact_layers(contacts, status)
    st = Dict(r.part_id => (r.in_education, r.working) for r in eachrow(status))
    flag(c) = coalesce.(contacts[!, c], false)
    home, school, work = flag(:cnt_home), flag(:cnt_school), flag(:cnt_work)
    return [contact_layer(home[i], school[i], work[i], st[contacts.part_id[i]]...)
        for i in 1:nrow(contacts)]
end

"""
    layer_assignment_summary(contacts, layers)

Number and share of contacts per layer and origin (location flags), to report how
contacts without a matching setting were assigned.
"""
function layer_assignment_summary(contacts, layers)
    flags = POLYMOD_LOCATIONS
    origin = map(1:nrow(contacts)) do i
        set = [string(f)[5:end] for f in flags if coalesce(contacts[i, f], false)]
        isempty(set) ? "none" : join(set, "+")
    end
    df = combine(groupby(DataFrame(layer = layers, origin = origin), [:layer, :origin]),
        nrow => :contacts)
    df.share = df.contacts ./ nrow(contacts)
    return sort(df, [:layer, order(:contacts, rev = true)])
end

"""
    participant_layer_counts(status, contacts, layers; household_cap = false)

Contacts per participant (rows as in `status`, zeros included) and layer.
With `household_cap`, household contacts beyond the other household members
(`hh_size - 1`) count as municipality contacts (visitors at home).
"""
function participant_layer_counts(status, contacts, layers; household_cap = false)
    row = Dict(id => k for (k, id) in enumerate(status.part_id))
    counts = DataFrame([l => zeros(Int, nrow(status)) for l in CONTACT_LAYERS])
    for (id, l) in zip(contacts.part_id, layers)
        counts[row[id], l] += 1
    end
    if household_cap
        for k in 1:nrow(status)
            ismissing(status.hh_size[k]) && continue
            excess = max(0, counts.household[k] - (status.hh_size[k] - 1))
            counts.household[k] -= excess
            counts.municipality[k] += excess
        end
    end
    return counts
end


## WEIGHTING

"""
    age_weighted_mean(x, groups, target)

Mean of `x`, with each age group weighted to its share in `target` (shares per
age group of `CONTACT_AGE_LIMITS`). Values with a missing age group are left out;
target shares of groups without values are spread over the other groups.
"""
function age_weighted_mean(x, groups, target)
    keep = .!ismissing.(groups)
    x, groups = x[keep], collect(skipmissing(groups[keep]))
    present = unique(groups)
    total = sum(target[g] for g in present)
    return sum(target[g] / total * mean(x[groups .== g]) for g in present)
end

"""
    gems_age_shares(sim)

Shares per age group (`CONTACT_AGE_LIMITS`) of the members of each GEMS layer:
all individuals (household, municipality), individuals in an education setting
(school) and in a workplace (work).
"""
function gems_age_shares(sim)
    n = length(CONTACT_AGE_LIMITS)
    counts = Dict(l => zeros(n) for l in CONTACT_LAYERS)
    for i in individuals(sim)
        g = age_group(Int(age(i)))
        counts[:household][g] += 1
        counts[:municipality][g] += 1
        is_student(i) && (counts[:school][g] += 1)
        is_working(i) && (counts[:work][g] += 1)
    end
    return Dict(l => c ./ sum(c) for (l, c) in counts)
end


## RATES

"""
    contact_layer_rates(polymod, age_shares; country = "DE", household_cap = false)

Contacts per member and day of each GEMS layer from the POLYMOD participants of
`country` (all countries if `nothing`). Members: all participants (household,
municipality), participants in education (school), working participants (work).

Columns: `layer`, `members` (participants), `contacts`, `raw` (unweighted mean),
`weekday` and `weekend` (age-weighted to `age_shares`, see `gems_age_shares`),
and `daily` (`WEEKDAY_SHARE` weekday + rest weekend, the value for GEMS).
"""
function contact_layer_rates(polymod, age_shares; country = "DE", household_cap = false)
    participants = isnothing(country) ? polymod.participants :
        polymod.participants[coalesce.(polymod.participants.country .== country, false), :]
    status = polymod_status(participants)
    contacts = polymod.contacts[in.(polymod.contacts.part_id, Ref(Set(status.part_id))), :]
    layers = contact_layers(contacts, status)
    counts = participant_layer_counts(status, contacts, layers, household_cap = household_cap)

    members = Dict(
        :household => trues(nrow(status)), :municipality => trues(nrow(status)),
        :school => status.in_education, :work => status.working)

    rows = map(CONTACT_LAYERS) do l
        m = members[l]
        x, groups, weekend = counts[m, l], status.age_group[m], status.weekend[m]
        day(w) = coalesce.(weekend .== w, false)
        weekday = age_weighted_mean(x[day(false)], groups[day(false)], age_shares[l])
        weekend_mean = age_weighted_mean(x[day(true)], groups[day(true)], age_shares[l])
        (layer = l, members = count(m), contacts = sum(x), raw = mean(x),
            weekday = weekday, weekend = weekend_mean,
            daily = WEEKDAY_SHARE * weekday + (1 - WEEKDAY_SHARE) * weekend_mean)
    end
    return DataFrame(rows)
end

"""
    household_contact_probabilities(polymod; country = "DE", max_size = 5)

Share of the other household members met per day (input of
`HouseholdSizeSampling`) by household size, from the POLYMOD participants of
`country` (all countries if `nothing`); sizes from `max_size` up are pooled.
Home contacts are capped at the other household members (`hh_size - 1`), the
rest are visitors. Daily average: weekdays `WEEKDAY_SHARE`, weekend the rest.

Columns: `hh_size` (smallest size of the group), `participants`, `contacts`
(daily household contacts), `others` (mean number of other members) and
`probability` (`contacts / others`).
"""
function household_contact_probabilities(polymod; country = "DE", max_size = 5)
    participants = isnothing(country) ? polymod.participants :
        polymod.participants[coalesce.(polymod.participants.country .== country, false), :]
    status = polymod_status(participants)
    contacts = polymod.contacts[in.(polymod.contacts.part_id, Ref(Set(status.part_id))), :]
    counts = participant_layer_counts(status, contacts, contact_layers(contacts, status),
        household_cap = true)

    rows = map(2:max_size) do n
        m = coalesce.(n == max_size ? status.hh_size .>= n : status.hh_size .== n, false)
        day(w) = m .& coalesce.(status.weekend .== w, false)
        daily(x) = WEEKDAY_SHARE * mean(x[day(false)]) + (1 - WEEKDAY_SHARE) * mean(x[day(true)])
        c = daily(counts.household)
        others = daily(status.hh_size .- 1)
        (hh_size = n, participants = count(m), contacts = c, others = others,
            probability = c / others)
    end
    return DataFrame(rows)
end

"""
    config_layer_contacts(configfile)

Contacts per member and day of each GEMS layer in `configfile` (sum of the
contact parameters of the layer's setting types, see `LAYER_SETTING_TYPES`).
"""
function config_layer_contacts(configfile)
    settings = TOML.parsefile(configfile)["Settings"]
    param(t) = haskey(settings, t) ?
        settings[t]["contact_sampling_method"]["parameters"]["contactparameter"] : 0.0
    # rounded: sums of decimal parameters are not exact in floating point
    return Dict(l => round(sum(param, types), digits = 6) for (l, types) in LAYER_SETTING_TYPES)
end

"""
    polymod_rate_table(polymod, age_shares; configfile)

`contact_layer_rates` for Germany and all countries, each with all home contacts
and with household contacts capped at the other household members, next to the
contacts per layer in `configfile`.
"""
function polymod_rate_table(polymod, age_shares; configfile)
    config = config_layer_contacts(configfile)
    tables = map(Iterators.product(["DE", nothing], [false, true])) do (country, cap)
        df = contact_layer_rates(polymod, age_shares, country = country, household_cap = cap)
        insertcols!(df, 1, :country => something(country, "all"),
            :household_rule => cap ? "capped" : "all home")
    end
    df = vcat(vec(tables)...)
    df.config = [config[l] for l in df.layer]
    return df
end


## ANALYSES

"""
    contact_analyses(;
        folder = nothing,
        folder_tag = nothing,
        configfile = "models/DE_model_first_wave.toml",
        sim = nothing,
        population = "DE")

Derives the daily contacts per GEMS layer from POLYMOD (`polymod_rate_table`,
weighted to the age structure of `sim`, built from `configfile` and `population`
if not passed) and writes them, together with the assignment of the contacts to
the layers, to `folder`.

If `folder` is `nothing`, a new result folder is created, with `folder_tag`
appended to its timestamp-based name (see `result_folder`).
"""
function contact_analyses(;
        folder = nothing,
        folder_tag = nothing,
        configfile = "models/DE_model_first_wave.toml",
        sim = nothing,
        population = "DE")

    isnothing(folder) ? folder = result_folder(folder_tag) : mkpath(folder)
    sim = isnothing(sim) ? Simulation(configfile = configfile, population = population, seed = 0) : sim

    polymod = load_polymod()
    rates = polymod_rate_table(polymod, gems_age_shares(sim), configfile = configfile)
    CSV.write(joinpath(folder, "polymod_contact_rates.csv"), rates)

    # household contact probabilities by household size (HouseholdSizeSampling)
    households = vcat([insertcols!(household_contact_probabilities(polymod, country = c), 1,
        :country => something(c, "all")) for c in ["DE", nothing]]...)
    CSV.write(joinpath(folder, "polymod_household_probabilities.csv"), households)

    # assignment of the German contacts to the layers
    de = polymod.participants[coalesce.(polymod.participants.country .== "DE", false), :]
    status = polymod_status(de)
    contacts = polymod.contacts[in.(polymod.contacts.part_id, Ref(Set(de.part_id))), :]
    CSV.write(joinpath(folder, "polymod_layer_assignment_DE.csv"),
        layer_assignment_summary(contacts, contact_layers(contacts, status)))

    printinfo("POLYMOD contact rates written to $folder")
    return rates
end
