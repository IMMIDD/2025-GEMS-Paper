###
### CONTACT DATA: CONTACT PARAMETERS OF THE CONFIG FROM POLYMOD
###

export derive_contact_parameters


"""
    CONTACT_INPUT_DIR

Folder of the contact parameters derived from POLYMOD and the tables they are
derived from (written by `derive_contact_parameters`).
"""
const CONTACT_INPUT_DIR = joinpath("inputs", "contacts")

"""
    LEVEL_SHARES

Shares of the school and work contacts per member at each level of the setting
hierarchy: most contacts at the most immediate level (own class or office), each
higher level half of the one below, the two outermost levels the same.
"""
const LEVEL_SHARES = Dict(
    "SchoolClass" => 1 / 2, "SchoolYear" => 1 / 4, "School" => 1 / 8, "SchoolComplex" => 1 / 8,
    "Office" => 1 / 2, "Department" => 1 / 4, "Workplace" => 1 / 8, "WorkplaceSite" => 1 / 8)

"""
    POLYMOD_SOURCE_NOTE

Header of the derived files: source, license and changes.
"""
const POLYMOD_SOURCE_NOTE = """
    # Source: Mossong J, Hens N, Jit M, et al. (2020). POLYMOD social contact data (Version 2)
    # [Data set]. Zenodo. https://doi.org/10.5281/zenodo.3874557. Survey: Mossong J, et al. (2008).
    # PLoS Med 5(3):e74. https://doi.org/10.1371/journal.pmed.0050074
    # License: CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/)
    # Changes: own calculation (derive_contact_parameters.jl, src/contacts/): contacts assigned to
    # GEMS setting layers, averaged per layer member, daily average (weekdays 5/7, weekend 2/7)."""


"""
    contact_parameters(rates, households; care_home_size = CARE_HOME_SIZE)

Values of the contact sampling parameters in the config, one row per entry
(`setting_type`, `parameter`, `index` for vector parameters, `value`,
`derivation`), from the layer rates of one country with capped household
contacts (`contact_layer_rates`) and the household contact probabilities
(`household_contact_probabilities`) of the same country.
"""
function contact_parameters(rates, households; care_home_size = CARE_HOME_SIZE)
    daily(l) = only(rates.daily[rates.layer .== l])
    rows = NamedTuple[]
    add!(t, p, i, v, d) = push!(rows,
        (setting_type = t, parameter = p, index = i, value = v, derivation = d))

    # households: probability per other member by household size, care homes
    for (k, r) in enumerate(eachrow(households))
        label = k == nrow(households) ? "$(r.hh_size)+" : string(r.hh_size)
        add!("Household", "probabilities", k, r.probability,
            "share of the other household members met per day, household size $label")
    end
    add!("Household", "large_size", missing, care_home_size,
        "smallest household size treated as a care home (not from POLYMOD)")
    add!("Household", "large_contacts", missing, daily(:household_all_home),
        "mean daily home contacts of all participants (care homes)")

    # education and work: contacts per member split over the setting hierarchy
    for (layer, name) in [(:school, "an education setting"), (:work, "a workplace")]
        for t in LAYER_SETTING_TYPES[layer]
            add!(t, "contactparameter", missing, LEVEL_SHARES[t] * daily(layer),
                "$(LEVEL_SHARES[t]) of the daily contacts per member of $name")
        end
    end

    add!("Municipality", "contactparameter", missing, daily(:municipality),
        "daily contacts outside home, education and work, incl. visitors at home")
    return DataFrame(rows)
end

"""
    config_contact_parameters(configfile)

Contact sampling parameters in `configfile`, as rows like `contact_parameters`
(without `derivation`).
"""
function config_contact_parameters(configfile)
    settings = TOML.parsefile(configfile)["Settings"]
    rows = NamedTuple[]
    for (t, s) in settings
        haskey(s, "contact_sampling_method") || continue
        for (p, v) in get(s["contact_sampling_method"], "parameters", Dict())
            if v isa AbstractVector
                for (i, x) in enumerate(v)
                    push!(rows, (setting_type = t, parameter = p, index = i, config = x))
                end
            else
                push!(rows, (setting_type = t, parameter = p, index = missing, config = v))
            end
        end
    end
    return DataFrame(rows)
end

"""
    compare_contact_parameters(parameters, configfile; digits = 3)

`parameters` (see `contact_parameters`) next to the values in `configfile`, with
`matches` if the config holds the value rounded to `digits`.
"""
function compare_contact_parameters(parameters, configfile; digits = 3)
    df = leftjoin(parameters, config_contact_parameters(configfile),
        on = [:setting_type, :parameter, :index], matchmissing = :equal)
    df.matches = [!ismissing(c) && isapprox(c, round(v, digits = digits), atol = 10.0^-digits / 2)
        for (v, c) in zip(df.value, df.config)]
    return df
end

# CSV file with the POLYMOD source note above the table
function write_polymod_csv(file, df, title)
    open(file, "w") do io
        println(io, "# ", title)
        println(io, POLYMOD_SOURCE_NOTE)
        CSV.write(io, df, append = true, writeheader = true)
    end
    return file
end


"""
    derive_contact_parameters(;
        configfile = "models/DE_model_first_wave.toml",
        sim = nothing,
        population = "DE",
        country = "DE",
        dir = CONTACT_INPUT_DIR)

Derives the contact parameters of the config from the POLYMOD participants of
`country`, weighted to the age structure of `sim` (built from `configfile` and
`population` if not passed), and writes to `dir`:

- `contact_parameters_<country>.csv`: one row per config entry with its value
  and derivation (`contact_parameters`), the values to enter in the config.
- `polymod_contact_rates.csv`: daily contacts per member of each layer, for
  Germany and all countries, with all and with capped home contacts
  (`polymod_rate_table`).
- `polymod_household_probabilities.csv`: share of the other household members
  met per day by household size (`household_contact_probabilities`).
- `polymod_municipality_matrix_<country>.txt`: age-age matrix of the
  municipality contacts in the format of `AgeBasedContactSampling`, and
  `polymod_municipality_matrix_<country>.csv` with its steps
  (`municipality_contact_matrix`).

Returns the parameters next to the values in `configfile` and warns about the
entries that differ (rounded to three decimals).
"""
function derive_contact_parameters(;
        configfile = "models/DE_model_first_wave.toml",
        sim = nothing,
        population = "DE",
        country = "DE",
        dir = CONTACT_INPUT_DIR)

    mkpath(dir)
    sim = isnothing(sim) ? Simulation(configfile = configfile, population = population, seed = 0) : sim
    shares = gems_age_shares(sim)
    polymod = load_polymod()

    rates = polymod_rate_table(polymod, shares, configfile = configfile)
    households = vcat([insertcols!(household_contact_probabilities(polymod, country = c), 1,
        :country => something(c, "all")) for c in ["DE", nothing]]...)

    # layer rates of `country`: capped household contacts, plus all home contacts (care homes)
    capped = contact_layer_rates(polymod, shares, country = country, household_cap = true)
    all_home = contact_layer_rates(polymod, shares, country = country, household_cap = false)
    push!(capped, merge(NamedTuple(only(eachrow(all_home[all_home.layer .== :household, :]))),
        (layer = :household_all_home,)))
    parameters = contact_parameters(capped,
        household_contact_probabilities(polymod, country = country))

    write_polymod_csv(joinpath(dir, "contact_parameters_$country.csv"), parameters,
        "Contact sampling parameters of the GEMS config derived from POLYMOD ($country)")
    write_polymod_csv(joinpath(dir, "polymod_contact_rates.csv"), rates,
        "Daily contacts per member of each GEMS setting layer (POLYMOD; config: $configfile)")
    write_polymod_csv(joinpath(dir, "polymod_household_probabilities.csv"), households,
        "Share of the other household members met per day by household size (POLYMOD)")

    # age-age matrix of the municipality contacts (AgeBasedContactSampling)
    matrix = municipality_contact_matrix(polymod, population_age_counts(sim), country = country)
    write_contact_matrix(joinpath(dir, "polymod_municipality_matrix_$country.txt"),
        matrix.probabilities)
    write_polymod_csv(joinpath(dir, "polymod_municipality_matrix_$country.csv"), matrix_table(matrix),
        "Age-age matrix of the municipality contacts (POLYMOD $country; contacts without age: " *
        "$(round(100 * matrix.missing_age, digits = 1))%, left out)")
    printinfo("Contact parameters written to $dir")

    comparison = compare_contact_parameters(parameters, configfile)
    differing = comparison[.!comparison.matches, :]
    if !isempty(differing)
        @warn "$(nrow(differing)) contact parameters differ from $configfile. Enter (rounded):\n" *
            join(["$(r.setting_type).$(r.parameter)" * (ismissing(r.index) ? "" : "[$(r.index)]") *
                " = $(round(r.value, digits = 3)) (config: $(r.config))" for r in eachrow(differing)], "\n")
    end
    return comparison
end
