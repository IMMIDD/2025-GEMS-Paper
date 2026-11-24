###
### CONTACT DATA: AGE-AGE MATRIX OF THE MUNICIPALITY CONTACTS FROM POLYMOD
###

# Input of `AgeBasedContactSampling` for municipalities: for each age group of an
# individual, the probability distribution over the age groups of its contacts,
# relative to the population (GEMS multiplies each row by the age structure of
# the municipality before sampling, see `AgeBasedContactSampling`).


"""
    MATRIX_INTERVAL

Width of the age groups of the municipality contact matrix (years).
"""
const MATRIX_INTERVAL = 5

"""
    MATRIX_GROUPS

Number of age groups of the municipality contact matrix: 0-4, ..., 65-69, 70+
(in GEMS, the last group collects all older individuals).
"""
const MATRIX_GROUPS = 15

# age group (1 to MATRIX_GROUPS) of an age in years
matrix_group(age) = ismissing(age) ? missing : min(Int(age) ÷ MATRIX_INTERVAL, MATRIX_GROUPS - 1) + 1

"""
    contact_age(contacts)

Age of each contact: the exact age, otherwise the midpoint of the estimated age
range, otherwise `missing`.
"""
function contact_age(contacts)
    return map(eachrow(contacts)) do c
        !ismissing(c.cnt_age_exact) && return c.cnt_age_exact
        (ismissing(c.cnt_age_est_min) || ismissing(c.cnt_age_est_max)) && return missing
        return (c.cnt_age_est_min + c.cnt_age_est_max) / 2
    end
end

"""
    population_age_counts(sim)

Number of individuals of `sim` per age group of the municipality contact matrix.
"""
function population_age_counts(sim)
    counts = zeros(MATRIX_GROUPS)
    for i in individuals(sim)
        counts[matrix_group(age(i))] += 1
    end
    return counts
end


"""
    municipality_contact_matrix(polymod, population; country = "DE")

Age-age matrix of the municipality contacts of the POLYMOD participants of
`country` (contacts assigned to the municipality layer by `contact_layers`
except visitors at home, whose contacts are not identifiable), with
`population` the number of individuals per age group (`population_age_counts`).

1. `contacts`: mean daily contacts of a participant of age group i with persons
   of age group j (weekdays 5/7, weekend 2/7; participants without contacts
   count with 0).
2. `reciprocal`: made reciprocal, so that the contacts of group i with group j
   equal those of group j with group i in the population:
   (c_ij N_i + c_ji N_j) / (2 N_i).
3. `per_capita`: divided by the population of the contacts' age group, N_j.
4. `probabilities`: rows normalized to 1 (input of `AgeBasedContactSampling`).

Also returns `participants` per age group and the share of contacts without an
age (`missing_age`), which are left out.
"""
function municipality_contact_matrix(polymod, population; country = "DE")
    length(population) == MATRIX_GROUPS ||
        throw(ArgumentError("Population counts must have $MATRIX_GROUPS age groups."))
    participants = polymod.participants[coalesce.(polymod.participants.country .== country, false), :]
    status = polymod_status(participants)
    contacts = polymod.contacts[in.(polymod.contacts.part_id, Ref(Set(status.part_id))), :]
    contacts = contacts[contact_layers(contacts, status) .== :municipality, :]
    contact_groups = matrix_group.(contact_age(contacts))
    missing_age = count(ismissing, contact_groups) / length(contact_groups)

    # participants with age and day of the week
    group = Dict(r.part_id => (matrix_group(a), r.weekend)
        for (r, a) in zip(eachrow(status), participants.part_age))
    n = Dict(w => zeros(MATRIX_GROUPS) for w in (false, true))
    for (g, w) in values(group)
        (ismissing(g) || ismissing(w)) && continue
        n[w][g] += 1
    end
    sum(n[false]) + sum(n[true]) > 0 ||
        throw(ErrorException("No participants with age and day of the week."))

    # contacts per day type: sums from group i to group j
    s = Dict(w => zeros(MATRIX_GROUPS, MATRIX_GROUPS) for w in (false, true))
    for (id, j) in zip(contacts.part_id, contact_groups)
        ismissing(j) && continue
        i, w = group[id]
        (ismissing(i) || ismissing(w)) && continue
        s[w][i, j] += 1
    end

    # daily mean per participant of group i; a day type without participants of a
    # group is replaced by the other one
    c = zeros(MATRIX_GROUPS, MATRIX_GROUPS)
    for i in 1:MATRIX_GROUPS
        weekday = n[false][i] > 0 ? s[false][i, :] ./ n[false][i] : nothing
        weekend = n[true][i] > 0 ? s[true][i, :] ./ n[true][i] : nothing
        isnothing(weekday) && isnothing(weekend) &&
            throw(ErrorException("No participants in age group $i."))
        c[i, :] = isnothing(weekday) ? weekend : isnothing(weekend) ? weekday :
            WEEKDAY_SHARE .* weekday .+ (1 - WEEKDAY_SHARE) .* weekend
    end

    N = population
    reciprocal = [(c[i, j] * N[i] + c[j, i] * N[j]) / (2 * N[i])
        for i in 1:MATRIX_GROUPS, j in 1:MATRIX_GROUPS]
    per_capita = reciprocal ./ N'
    probabilities = per_capita ./ sum(per_capita, dims = 2)

    return (contacts = c, reciprocal = reciprocal, per_capita = per_capita,
        probabilities = probabilities, participants = n[false] .+ n[true],
        missing_age = missing_age)
end

"""
    matrix_table(matrix)

The matrices of `municipality_contact_matrix` as a long table: one row per
pair of age groups (`age_group`, `contact_age_group`, lower age limits).
"""
function matrix_table(matrix)
    limits = [(k - 1) * MATRIX_INTERVAL for k in 1:MATRIX_GROUPS]
    rows = [(age_group = limits[i], contact_age_group = limits[j],
        contacts = matrix.contacts[i, j], reciprocal = matrix.reciprocal[i, j],
        per_capita = matrix.per_capita[i, j], probability = matrix.probabilities[i, j])
        for i in 1:MATRIX_GROUPS for j in 1:MATRIX_GROUPS]
    return DataFrame(rows)
end

"""
    write_contact_matrix(file, probabilities)

Writes `probabilities` in the format read by `AgeBasedContactSampling`: one row
per line, values separated by spaces, no header.
"""
function write_contact_matrix(file, probabilities)
    open(file, "w") do io
        for i in axes(probabilities, 1)
            println(io, join(string.(probabilities[i, :]), " "))
        end
    end
    return file
end
