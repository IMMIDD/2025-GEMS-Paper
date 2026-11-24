###
### CONTACT SAMPLING: HOUSEHOLD-SIZE SPECIFIC CONTACTS
###

export HouseholdSizeSampling


"""
    HouseholdSizeSampling <: ContactSamplingMethod

Contact sampling for households: each other present member is contacted once
per day with a probability that depends on the household size. POLYMOD counts
the distinct persons met per day; in households, participants met about 90% of
the other household members, almost independent of the household size (see
`household_contact_probabilities`).

Households with at least `large_size` members (care homes) instead draw a
Poisson-distributed number of contacts with mean `large_contacts` from the
present members (as `ContactparameterSampling`): POLYMOD has no care home
residents, so they get the mean daily home contacts of all participants.

# Parameters

- `probabilities::NTuple{N, Float64}`: probability to contact each other member,
  for household sizes 2, 3, ...; the last entry applies to all larger households
  below `large_size` (e.g. `[q2, q3, q4, q5plus]`, given as a vector or tuple).
- `large_size::Int = typemax(Int)`: smallest household size treated as a care
  home (none by default).
- `large_contacts::Float64 = 0.0`: mean daily contacts in care homes.

Single-person households have no household contacts.

The method holds only plain values (`isbitstype`): GEMS gives every setting its
own copy of the sampling method (`deepcopy`), which for such a struct returns
the method itself, so all households share one object.
"""
struct HouseholdSizeSampling{N} <: ContactSamplingMethod
    probabilities::NTuple{N, Float64}
    large_size::Int
    large_contacts::Float64

    function HouseholdSizeSampling{N}(probabilities::NTuple{N, Float64},
            large_size::Int, large_contacts::Float64) where N
        N >= 1 ||
            throw(ArgumentError("At least one contact probability is required."))
        all(p -> 0 <= p <= 1, probabilities) ||
            throw(ArgumentError("Contact probabilities must be in [0, 1]. Got $probabilities."))
        large_size >= 2 ||
            throw(ArgumentError("'large_size' must be at least 2. Got $large_size."))
        large_contacts >= 0 ||
            throw(ArgumentError("'large_contacts' must be non-negative. Got $large_contacts."))
        return new{N}(probabilities, large_size, large_contacts)
    end
end

function HouseholdSizeSampling(probabilities::Union{AbstractVector{<:Real}, Tuple{Vararg{Real}}},
        large_size::Integer = typemax(Int), large_contacts::Real = 0.0)
    p = Tuple(Float64.(collect(probabilities)))
    return HouseholdSizeSampling{length(p)}(p, Int(large_size), Float64(large_contacts))
end

# keyword constructor, as GEMS builds contact sampling methods from config files
HouseholdSizeSampling(; probabilities, large_size = typemax(Int), large_contacts = 0.0) =
    HouseholdSizeSampling(probabilities, large_size, large_contacts)

"""
    contact_probability(csm::HouseholdSizeSampling, household_size)

Probability to contact each other member in a household of `household_size`
(below `csm.large_size`).
"""
contact_probability(csm::HouseholdSizeSampling, household_size::Integer) =
    csm.probabilities[clamp(household_size - 1, 1, length(csm.probabilities))]


"""
    sample_contacts!(indivs, csm::HouseholdSizeSampling, setting, individual_index,
        present_inds, tick, replace, rng)

Fills `indivs` with the contacts of the individual at `individual_index` in
`present_inds`. Households below `csm.large_size` members (the setting's size,
not only the present members): each other present member once, with
`contact_probability` (`replace` is ignored). Larger households: Poisson
(`csm.large_contacts`) contacts drawn with replacement from the other present
members.
"""
function GEMS.sample_contacts!(indivs::Vector{Individual}, csm::HouseholdSizeSampling,
        setting::Setting, individual_index::Int, present_inds::Vector{Individual},
        tick::Int16, replace::Bool, rng::Random.Xoshiro)

    setting isa IndividualSetting ||
        throw(ArgumentError("HouseholdSizeSampling needs a setting with individuals, got $(typeof(setting))."))
    if isempty(present_inds)
        throw(ArgumentError("No Individual is present in $setting. Please provide a Setting, where at least 1 Individual is present!"))
    end

    empty!(indivs)
    n_present = length(present_inds)
    n_present == 1 && return
    n = size(setting)

    if n >= csm.large_size
        k = GEMS.gems_rand(rng, GEMS.Distributions.Poisson(csm.large_contacts))
        for _ in 1:k
            # any other present member (skipping the individual itself)
            offset = GEMS.gems_rand(rng, 1:(n_present - 1))
            push!(indivs, present_inds[mod(individual_index + offset - 1, n_present) + 1])
        end
        return
    end

    q = contact_probability(csm, n)
    for (k, other) in enumerate(present_inds)
        k == individual_index && continue
        GEMS.gems_rand(rng) < q && push!(indivs, other)
    end
    return
end
