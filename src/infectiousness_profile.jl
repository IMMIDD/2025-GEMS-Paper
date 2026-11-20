###
### INFECTIOUSNESS PROFILE (TRANSMISSION FUNCTION)
###

export InfectiousnessProfile
export with_transmission_scale

"""
    InfectiousnessProfile <: TransmissionFunction

A `TransmissionFunction` that scales the transmission probability by an infectiousness
profile over the days since the infectiousness onset of the infecting individual.

The per-contact transmission probability on day `τ` (`τ = 0` is the tick of the
infectiousness onset) is `transmission_scale * w[τ]`, with `w[τ]` the probability mass
of the profile in `[τ, τ + 1)`. The profile is `distribution(parameters...)` (as for the
progression parameters of a config), shifted by `offset` days: `τ = X - offset` for a draw
`X` of the distribution. Its mass before `τ = 0` is dropped and the rest normalized to 1,
so `transmission_scale` is the sum of the daily per-contact transmission probabilities
over the infection. Recovered individuals are immune, as with `ConstantTransmissionRate`.

# Config

```toml
[Pathogens.Covid19.transmission_function]
    type = "InfectiousnessProfile"
    [Pathogens.Covid19.transmission_function.parameters]
        transmission_scale = 0.46
        distribution = "Gamma"
        parameters = [97.1875, 0.26891]
        offset = 21.625
```
"""
mutable struct InfectiousnessProfile <: TransmissionFunction
    transmission_scale::Float64
    distribution::String
    parameters::Vector{Float64}
    offset::Float64
    weights::Vector{Float64}

    function InfectiousnessProfile(;
            transmission_scale::Real,
            distribution::AbstractString,
            parameters::AbstractVector,
            offset::Real = 0.0,
            max_days::Integer = 60)

        transmission_scale < 0 && throw(ArgumentError("Transmission scale must be non-negative."))

        d = GEMS.create_distribution(Dict("distribution" => distribution, "parameters" => parameters))
        F(x) = GEMS.Distributions.cdf(d, x + offset)
        w = [F(t + 1) - F(t) for t in 0:max_days]
        sum(w) > 0 || throw(ArgumentError("The profile has no mass after the infectiousness onset."))

        return new(transmission_scale, distribution, Float64.(parameters), offset, w ./ sum(w))
    end
end

Base.show(io::IO, ip::InfectiousnessProfile) = write(io,
    "InfectiousnessProfile(scale=$(ip.transmission_scale), $(ip.distribution)$(Tuple(ip.parameters)), offset=$(ip.offset))")


"""
    transmission_probability(tf::InfectiousnessProfile, infecter::Individual, infectee::Individual,
        setting::Setting, tick::Int16, rng::Xoshiro)

Returns `transmission_scale * w[τ]` for the day `τ` since the infectiousness onset of the
`infecter`, and `0.0` if the `infectee` has already recovered (natural immunity).
"""
function GEMS.transmission_probability(tf::InfectiousnessProfile, infecter::Individual,
        infectee::Individual, setting::Setting, tick::Int16, rng::Random.Xoshiro)::Float64
    !infected(infecter) && throw(ArgumentError("Infecting individual must be infected to calculate transmission probability."))

    # natural immunity
    if -1 < recovery(infectee) <= tick
        return 0.0
    end

    τ = tick - infectiousness_onset(infecter)
    (τ < 0 || τ >= length(tf.weights)) && return 0.0
    return tf.transmission_scale * tf.weights[τ + 1]
end


"""
    with_transmission_scale(tf::TransmissionFunction, x)

The transmission function `tf` with another scale `x`: the rate of a
`ConstantTransmissionRate`, the `transmission_scale` of an `InfectiousnessProfile`
(same profile).
"""
with_transmission_scale(::ConstantTransmissionRate, x) = ConstantTransmissionRate(transmission_rate = x)
with_transmission_scale(tf::InfectiousnessProfile, x) = InfectiousnessProfile(
    transmission_scale = x,
    distribution = tf.distribution,
    parameters = tf.parameters,
    offset = tf.offset)
