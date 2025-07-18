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
    configfile,
    "data/people.jld2",
    "data/settings.jld2")

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

# INCLUDE SCCENARIOS
include("sc_outbreak.jl")



# RUN ALL
"""
    runall(;
        sim_obj = nothing,
        configfile = "models/DE_model_R0_3.28.toml",
        popfile = "data/people.jld2",
        settingsfile = "data/settings.jld2")

Executes all model analyses.
"""
function runall(;
        sim_obj = nothing,
        configfile = "models/DE_model_R0_3.28.toml",
        popfile = "data/people.jld2",
        settingsfile = "data/settings.jld2",
        sc_outbreak_config = "models/DE_model_R0_3.28_Heinsberg_Seed.toml")

    # create result folder
    folder = result_folder()

    # sim object
    sim = isnothing(sim_obj) ? Simulation(configfile, popfile, settingsfile) : sim_obj

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
    if !is_nothing(sim_obj)
        @warn "If you pass a sim-object to runall(), it will not perform the simulation scenarios."
        return;
    end

    printinfo("### RUNNING OUTBREAK SCENARIO")

    
end

end # module GEMSPaper
