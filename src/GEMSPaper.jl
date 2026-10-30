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
import LinearAlgebra
import Random
import TOML

@reexport using GEMS

# EXPORTS

export result_folder
export DE_sim
export map_dims
export map_titlefontsize
export map_legendfontsize
export runall
export run_first_wave

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

# INCLUDE RKI DATA AND MODEL (for the first wave scenario)
include("first_wave/rki_data.jl")
include("first_wave/rki_imports.jl")
include("first_wave/rki_seeding.jl")
include("first_wave/rki_phases.jl")
include("first_wave/rki_reproduction.jl")
include("first_wave/rki_simulation.jl")
include("first_wave/rki_comparison.jl")
include("first_wave/rki_calibration.jl")
include("first_wave/rki_figures.jl")

# INCLUDE SCENARIOS
include("sc_outbreak.jl")
include("sc_self_isolation.jl")
include("sc_sweeping.jl")
include("sc_testing.jl")
include("sc_first_wave.jl")



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


# RUN FIRST WAVE SCENARIO
"""
    run_first_wave(;
        configfile = "models/DE_model_RKI_2020_short_gt.toml",
        population = "DE",
        folder = nothing,
        parts = FIRST_WAVE_PARTS)

Executes the first wave 2020 scenario, compared with the RKI data (not part
of `runall`). To continue an earlier run or to run single `parts`, pass its `folder`.
"""
function run_first_wave(;
        configfile = "models/DE_model_RKI_2020_short_gt.toml",
        population = "DE",
        folder = nothing,
        parts = FIRST_WAVE_PARTS)

    # create result folder
    isnothing(folder) ? folder = result_folder() : nothing

    printinfo("### RUNNING FIRST WAVE SCENARIO")
    scenario_first_wave(
        configfile = configfile,
        population = population,
        folder = folder,
        parts = parts
    )

end

end # module GEMSPaper
