module GEMSPaper

# PACKAGES
using CSV
using DataFrames
using Dates
using Measures
using Plots
using Proj
using Statistics
using StatsBase
using Revise
using Reexport

@reexport using GEMS

# EXPORTS

export result_folder
export DE_sim
export map_dims
export map_titlefontsize
export map_legendfontsize
export runall

# COMMON FUNCTIONS
"""
result_folder()

Creates the path to a new result folder.
"""
function result_folder() 
    folder = joinpath("results", "$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS_sss"))")
    mkpath(folder)
    return folder
end

"""
    DE_sim(;configfile = "models/DE_model_R0_3.28.toml")

Wrapper for German simulation model.
"""
DE_sim(;configfile = "models/DE_model_R0_3.28.toml") = Simulation(
    configfile=configfile,
    population="DE",
    seed = 0)

"""
    map_dims()

Returns standard map dimensions (for Germany).
"""
map_dims() = (1600, 2000)
map_titlefontsize() = 34
map_legendfontsize() = 24

# INCLUDE STYLES
include("styles/rd_contacts_only.jl")

# INCLUDE ANALYSES
include("demographics.jl")
include("settings.jl")
include("mobility.jl")
include("contacts.jl")
include("benchmarking.jl")
include("rki_data.jl")
include("rki_seeding.jl")
include("rki_imports.jl")

# INCLUDE SCENARIOS
include("sc_outbreak.jl")
include("sc_self_isolation.jl")
include("sc_sweeping.jl")
include("sc_testing.jl")
include("sc_rki_validation.jl")



# RUN ALL
"""
    runall(;
        sim_obj = nothing,
        configfile = "models/DE_model_R0_3.28.toml",
        population = "DE")

Executes all model analyses.
"""
function runall(;
        sim_obj = nothing,
        configfile = "models/DE_model_R0_3.28.toml",
        population = "DE")

    # create result folder
    folder = result_folder()

    # sim object
    sim = isnothing(sim_obj) ? Simulation(configfile=configfile, population=population, seed = 0) : sim_obj

    
    printinfo("### RUNNING DEMOGRAPHICS ANALYSES")
    # run mobility analyses
    demographics(
        folder = folder,
        sim = sim
    )

    printinfo("### RUNNING SETTINGS ANALYSES")
    # run setting analyses
    settinganalyses(
        folder = folder,
        sim = sim
    )

    printinfo("### RUNNING MOBILITY ANALYSES")
    # run mobility analyses
    mobility(
        folder = folder,
        sim = sim
    )

    printinfo("### RUNNING CONTACT ANALYSES")
    # run contact analyses
    contacts(
        folder = folder,
        sim = sim
    )

    # break here, if sim-object was provided
    if !isnothing(sim_obj)
        @warn "If you pass a sim-object to runall(), it will not perform the simulation scenarios."
        return;
    end

    printinfo("### RUNNING OUTBREAK SCENARIO")
    scenario_outbreak(
        population = population,
        folder = folder
    )

    printinfo("### RUNNING SIMPLE INTERVENTION SCENARIO")
    scenario_self_isolation(
        population = population,
        folder = folder
    )

    printinfo("### RUNNING TESTING SCENARIOS")
    scenario_testing(
        population = population,
        folder = folder
    )

    printinfo("### RUNNING BENCHMARK")
    benchmark(
        folder = folder
    )

    printinfo("### RUNNING PARAMETER SWEEPS")
    scenario_sweeping(
        folder = folder
    )

end

end # module GEMSPaper
