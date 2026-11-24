using GEMSPaper

# Derivation of the contact parameters of the config from the POLYMOD survey
# (not part of runall()): writes the parameter values per config entry and the
# tables they are derived from to inputs/contacts/ and lists the config entries
# that differ from the derived values (see `derive_contact_parameters`).

configfile = "models/DE_model_first_wave.toml"

# age structure of the German population model for weighting the participants
population = "DE"

# POLYMOD participants the parameters are derived from
country = "DE"

derive_contact_parameters(
    configfile = configfile,
    population = population,
    country = country)
