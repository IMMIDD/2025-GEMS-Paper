module GEMSPaper

# PACKAGES
using DataFrames
using Dates
using Measures
using Plots
using Statistics
using Revise
using Reexport

@reexport using GEMS

# EXPORTS

export result_folder
export DE_sim
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
    DE_sim()

Wrapper for German simulation model.
"""
DE_sim() = Simulation(
        "models/DE_model_R0_3.28.toml",
        "data/people.jld2",
        "data/settings.jld2")

# INCLUDE ANALYSES
include("demographics.jl")
include("settings.jl")
include("mobility.jl")



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
        settingsfile = "data/settings.jld2")

    # create result folder
    folder = result_folder()

    # sim object
    sim = isnothing(sim_obj) ? Simulation(configfile, popfile, settingsfile) : sim_obj

    # run mobility analyses
    demographics(
        folder = folder,
        sim = sim
    )

    # run setting analyses
    settinganalyses(
        folder = folder,
        sim = sim
    )

    # run mobility analyses
    mobility(
        folder = folder,
        sim = sim
    )
end

end # module GEMSPaper
