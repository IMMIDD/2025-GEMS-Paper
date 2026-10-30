###
### FIRST WAVE 2020: CONTACT PHASES AND TESTING
###

"""
    ASCERTAINMENT_2020

Share of the infections that are detected by tests.
"""
const ASCERTAINMENT_2020 = 0.22

"""
    PHASE_DATES

Start dates of the contact phases 2 to 4 in 2020: 09.03. (large events
cancelled), 16.03. (schools closed), and 23.03. (contact restrictions).
Phase 1 is the time before.
"""
const PHASE_DATES = [Date(2020, 3, 9), Date(2020, 3, 16), Date(2020, 3, 23)]

"""
    PHASE_TRANSITIONS

Days of the linear transitions between the phases (Dehning et al. 2020, Science 369, eabb9789).
"""
const PHASE_TRANSITIONS = [3, 3, 4]

"""
    SCHOOL_CLOSURE_2020

Date from which all schools are closed (16.03.2020).
"""
const SCHOOL_CLOSURE_2020 = Date(2020, 3, 16)

"""
    SCHOOL_TYPES

Setting types closed by the school closure.
"""
const SCHOOL_TYPES = (SchoolClass, SchoolYear, School, SchoolComplex)

"""
    CARE_HOME_SIZE

Households with at least this many members are treated as care homes: their
contacts are scaled with the phases like those of the other settings.
"""
const CARE_HOME_SIZE = 11


# days since s (dates or ticks)
elapsed(x::Date, s::Date) = Dates.value(x - s)
elapsed(x::Real, s::Real) = x - s


"""
    phase_weights(x, switches, transitions)

Weights of the phases on day `x` (date or tick): change k moves phase k to
phase k + 1 linearly over `transitions[k]` days from `switches[k]` on (0: step).
Each day gets the mean of the ramp over the day.
"""
function phase_weights(x, switches, transitions)
    v = zeros(length(switches) + 1)
    v[1] = 1.0
    for (k, (s, days)) in enumerate(zip(switches, transitions))
        e = elapsed(x, s)
        a = days == 0 ? Float64(e >= 0) : clamp((e + 0.5) / days, 0.0, 1.0)
        v[k + 1] += a * v[k]
        v[k] -= a * v[k]
    end
    return v
end


"""
    phase_matrix(xs, switches, transitions)

Days x phases matrix of `phase_weights`.
"""
function phase_matrix(xs, switches, transitions)
    rows = [permutedims(phase_weights(x, switches, transitions)) for x in xs]
    return reduce(vcat, rows; init = zeros(0, length(switches) + 1))
end


"""
    check_transitions(dates, transitions)

Checks that there is one transition (in days, at least 0) per sorted phase date.
"""
function check_transitions(dates, transitions)
    if length(transitions) != length(dates)
        throw(ArgumentError("Need one transition per phase date ($(length(dates))), " *
            "got $(length(transitions))."))
    end
    all(>=(0), transitions) || throw(ArgumentError("Transitions must be at least 0 days."))
    issorted(dates) || throw(ArgumentError("Phase dates must be sorted."))
end



###
### CONTACT PHASES AND TESTING
###

"""
    progression_shares(sim)

Shares of symptomatic and hospitalized (`Hospitalized`, `Critical`) infections
in the population of `sim`.
"""
function progression_shares(sim::Simulation)
    pa = progression_assignment(pathogen(sim))
    pa isa GEMS.AgeBasedProgressionAssignment ||
        throw(ArgumentError("Needs an AgeBasedProgressionAssignment, got $(typeof(pa))."))
    categories = nameof.(pa.progression_categories)
    n = zeros(length(pa.age_groups))
    for i in individuals(sim)
        n[findfirst(g -> in_group(i.age, g), pa.age_groups)] += 1
    end
    M = pa.stratification_matrix # categories x age groups
    share(c) = sum(M[in.(categories, Ref(c)), :] * n) / sum(n)
    return (
        symptomatic = 1 - share((:Asymptomatic,)),
        hospitalized = share((:Hospitalized, :Critical)))
end


"""
    symptom_test_probability(ascertainment, symptomatic, hospitalized)

Probability that a symptomatic case is tested, such that a share `ascertainment`
of all infections is detected if all hospitalized cases are tested as well.
"""
function symptom_test_probability(ascertainment::Real, symptomatic::Real, hospitalized::Real)
    p = (ascertainment - hospitalized) / (symptomatic - hospitalized)
    0 <= p <= 1 || throw(ArgumentError("Ascertainment $ascertainment is not reachable " *
        "(symptomatic $symptomatic, hospitalized $hospitalized)."))
    return p
end


# setting types whose contacts are scaled (all except households) and a function that sets their
# contact parameters to a factor times the config value; care homes (households with at least
# `care_homes` members) are scaled as well
function contact_scaling(sim::Simulation, care_homes::Integer)
    scaled(T) = T != Household && !isempty(settings(sim, T))
    types = [T for T in settingtypes(settingscontainer(sim)) if scaled(T)]
    base = Dict(map(types) do T
        csm = contact_sampling_method(first(settings(sim, T)))
        if !(csm isa ContactparameterSampling)
            throw(ArgumentError("$T does not use ContactparameterSampling."))
        end
        T => csm.contactparameter
    end)
    care = filter(h -> length(individuals(h)) >= care_homes, settings(sim, Household))
    c_care = isempty(care) ? 0.0 : contact_sampling_method(first(care)).contactparameter
    residents = sum(h -> length(individuals(h)), care; init = 0)
    printinfo("Care homes (households with at least $care_homes members): $(length(care)), " *
        "$residents residents, contact parameter $c_care scaled")

    function scale(f)
        for (T, c) in base
            csm = ContactparameterSampling(f * c)
            foreach(s -> contact_sampling_method!(s, csm), settings(sim, T))
        end
        if !isempty(care)
            csm = ContactparameterSampling(f * c_care)
            foreach(s -> contact_sampling_method!(s, csm), care)
        end
    end
    return (types = types, scale = scale)
end


"""
    add_contact_phases!(sim, factors;
        dates,
        transitions,
        school_closure,
        school_reopening = nothing,
        care_homes = CARE_HOME_SIZE)

Multiplies the contact parameters of all settings except households with
`factors[k]` in phase k (phase 1 before `dates[1]`, linear transitions over
`transitions` days). Care homes (households with at least `care_homes` members)
are scaled as well. Schools are closed from `school_closure` on and reopened
from `school_reopening` on (`nothing`: never).
"""
function add_contact_phases!(sim::Simulation, factors::AbstractVector{<:Real};
    dates,
    transitions,
    school_closure,
    school_reopening = nothing,
    care_homes = CARE_HOME_SIZE)

    if length(factors) != length(dates) + 1
        throw(ArgumentError("Need one factor per phase ($(length(dates) + 1)), " *
            "got $(length(factors))."))
    end
    check_transitions(dates, transitions)
    reopens_early = !isnothing(school_reopening) &&
        (isnothing(school_closure) || school_reopening <= school_closure)
    if reopens_early
        throw(ArgumentError("Schools can only reopen after their closure ($school_closure), " *
            "got $school_reopening."))
    end

    # stepmod runs at the end of a tick, so a change made at this tick applies from the next tick on
    tick_before(d) = isnothing(d) ? typemax(Int) : Dates.value(d - sim.startdate) - 1
    switches = [Dates.value(d - sim.startdate) for d in dates]
    factor_at(t) = sum(phase_weights(t, switches, transitions) .* factors)
    contacts = contact_scaling(sim, care_homes)
    schools = [T for T in SCHOOL_TYPES if T in contacts.types]
    set_schools!(f) = foreach(T -> foreach(f, settings(sim, T)), schools)

    closing, reopening = tick_before(school_closure), tick_before(school_reopening)
    current = Ref(factor_at(0))
    closed = Ref(closing < 0)
    reopened = Ref(false)
    contacts.scale(current[])
    if closed[]
        set_schools!(close!)
    end

    inner = sim.stepmod
    sim.stepmod = s -> begin
        inner(s)
        t = Int(tick(s))
        f = factor_at(t + 1)
        if f != current[]
            contacts.scale(f)
            current[] = f
        end
        if !closed[] && closing <= t
            set_schools!(close!)
            closed[] = true
        end
        if closed[] && !reopened[] && reopening <= t
            set_schools!(open!)
            reopened[] = true
        end
        s
    end
    return sim
end


# uniform draw per individual and seed (thread-safe, the same individuals in runs with the same
# seed)
person_draw(seed, i) = rand(Random.Xoshiro(hash((seed, id(i)))))


"""
    add_pcr_testing!(sim; p_symptomatic, seed = 0)

Adds a perfect PCR test at symptom onset (with probability `p_symptomatic`) and
at hospital admission. Positives are not isolated.
"""
function add_pcr_testing!(sim::Simulation; p_symptomatic::Real, seed = 0)
    pcr_symptoms = TestType("PCR symptoms", pathogen(sim), sim)
    pcr_hospital = TestType("PCR hospital", pathogen(sim), sim)

    symptom_testing = IStrategy("symptom PCR", sim,
        condition = i -> person_draw(seed, i) < p_symptomatic)
    add_measure!(symptom_testing, Test("symptom PCR", pcr_symptoms))
    add_symptom_trigger!(sim, SymptomTrigger(symptom_testing))

    # only cases that were not detected before
    hospital_testing = IStrategy("hospital PCR", sim,
        condition = i -> last_reported_at(i) < exposure(i))
    add_measure!(hospital_testing, Test("hospital PCR", pcr_hospital))
    add_hospitalization_trigger!(sim, HospitalizationTrigger(hospital_testing))
    return sim
end
