###
### CONTACT SAMPLING: AGE-BASED MUNICIPALITY CONTACTS WITH THE NATIONAL AGE STRUCTURE
###

export national_age_based_sampling
export use_age_based_municipalities!
export with_contactparameter


"""
    MUNICIPALITY_MATRIX_FILE

Age-age matrix of the municipality contacts derived from POLYMOD (written by
`derive_contact_parameters`, see `municipality_contact_matrix`).
"""
const MUNICIPALITY_MATRIX_FILE = joinpath(CONTACT_INPUT_DIR, "polymod_municipality_matrix_DE.txt")


"""
    national_age_based_sampling(sim;
        contactparameter,
        matrix_file = MUNICIPALITY_MATRIX_FILE,
        interval = MATRIX_INTERVAL)

`AgeBasedContactSampling` with the age-age matrix in `matrix_file` (rows sum to
1, age groups of `interval` years, the last one open) and the age structure of
the whole population of `sim`.

GEMS v1.0.2 computes the age structure on the first sampling call from the
members of that setting and keeps it for all settings sharing the method. With
the national age structure passed here, that computation never runs: all
settings use the same structure, and runs are reproducible.

With the national age structure, an individual of age group i in a setting with
age structure q has on average `contactparameter * (q ⋅ m_i) / (Q ⋅ m_i)` contacts
(m_i the matrix row, Q the national age structure): the national mean is
`contactparameter`, locally the number follows the availability of contacts of
the preferred ages.
"""
function national_age_based_sampling(sim;
        contactparameter::Real,
        matrix_file = MUNICIPALITY_MATRIX_FILE,
        interval::Integer = MATRIX_INTERVAL)

    data = Matrix{Float64}(GEMS.DelimitedFiles.readdlm(matrix_file))
    matrix = ContactMatrix{Float64}(data, Int(interval), size(data, 1) * interval)

    # age structure with the binning of GEMS: min(age, max_age) ÷ interval + 1
    max_age = matrix.aggregation_bound - 1
    pyramid = zeros(size(data, 1))
    for i in individuals(sim)
        pyramid[min(Int(age(i)), max_age) ÷ interval + 1] += 1
    end
    pyramid ./= sum(pyramid)

    # the keyword constructor of GEMS checks the rows; this one does not
    all(s -> isapprox(s, 1, atol = 1e-10), sum(data, dims = 2)) ||
        throw(ArgumentError("Rows of $matrix_file must sum to 1."))
    return AgeBasedContactSampling(Float64(contactparameter), Int(interval), matrix, pyramid)
end

"""
    use_age_based_municipalities!(sim;
        contactparameter = nothing,
        matrix_file = MUNICIPALITY_MATRIX_FILE)

Replaces the contact sampling of all municipalities of `sim` with
`national_age_based_sampling`, keeping the mean number of contacts of the
config (or `contactparameter`). Returns the new sampling method.
"""
function use_age_based_municipalities!(sim;
        contactparameter = nothing,
        matrix_file = MUNICIPALITY_MATRIX_FILE)

    municipalities = settings(sim, Municipality)
    c = something(contactparameter, contact_sampling_method(first(municipalities)).contactparameter)
    csm = national_age_based_sampling(sim, contactparameter = c, matrix_file = matrix_file)
    foreach(s -> contact_sampling_method!(s, csm), municipalities)
    return csm
end


"""
    with_contactparameter(csm, c)

Copy of the contact sampling method `csm` with the mean number of contacts `c`.
`AgeBasedContactSampling` keeps its matrix and age structure.
"""
with_contactparameter(csm::ContactparameterSampling, c::Real) = ContactparameterSampling(Float64(c))
with_contactparameter(csm::AgeBasedContactSampling, c::Real) =
    AgeBasedContactSampling(Float64(c), csm.interval, csm.contact_matrix, csm.age_pyramid)
with_contactparameter(csm, c::Real) =
    throw(ArgumentError("Cannot change the contact parameter of $(typeof(csm))."))
