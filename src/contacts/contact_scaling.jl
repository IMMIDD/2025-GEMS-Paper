###
### CONTACT SAMPLING: SCALING THE MEAN NUMBER OF CONTACTS
###

export with_contactparameter


"""
    with_contactparameter(csm, c)

Copy of the contact sampling method `csm` with the mean number of contacts `c`
(e.g. for the contact phases of the first wave, `contact_scaling`, or reduced
municipality contacts). `AgeBasedContactSampling` keeps its matrix and starts
with an empty age structure, so each setting it is assigned to takes the age
structure of its own members (instead of the one `csm` may already hold).
"""
with_contactparameter(csm::ContactparameterSampling, c::Real) = ContactparameterSampling(Float64(c))
with_contactparameter(csm::AgeBasedContactSampling, c::Real) =
    AgeBasedContactSampling(Float64(c), csm.interval, csm.contact_matrix, Float64[])
with_contactparameter(csm, c::Real) =
    throw(ArgumentError("Cannot change the contact parameter of $(typeof(csm))."))
