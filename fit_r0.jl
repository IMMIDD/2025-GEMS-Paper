using GEMSPaper, DataFrames, Statistics

# Global settings
people_file = "localdata/DE/people.jld2"
stngs_file = "localdata/DE/settings.jld2"
target_r0_value = 3.28

# Define custom ResultDataStyle for effective R extraction
mutable struct CustomRDfit <: GEMS.ResultDataStyle
    data::Dict{String, Any}
    function CustomRDfit(pP::PostProcessor)
        funcs = Dict(
            "dataframes" =>
                Dict(
                    "effectiveR" => () -> pP |> effectiveR,
                )
        )
        return(
            new(GEMS.process_funcs(funcs))
        )
    end
end

# Run simulation and calculate average R0
function get_avg_r0(configfile, tr, runs)
    r0_sum = 0.0
    for y in 1:runs
        sim = Simulation(configfile = configfile, population = people_file, settingsfile = stngs_file)
        sim.label = "TR: $tr"
        pathogen(sim).transmission_function = ConstantTransmissionRate(transmission_rate = tr)
        sim.stop_criterion = TimesUp(limit = 70)
        run!(sim)
        
        rd = ResultData(sim, style = "CustomRDfit")
        r0_sum += effectiveR(rd).effective_R[1:10] |> mean
    end
    return r0_sum / runs
end

# Bisection search algorithm to find the target transmission rate
function fit_transmission_rate_bin(configfile, target_r0; tr_low=0.01, tr_high=0.15, tolerance=0.02, max_iterations=15, nr_runs=10)
    println("\nStarting fit for: $configfile")
    println("Target R0: $target_r0")
    
    best_tr = tr_low

    for iter in 1:max_iterations
        tr_mid = (tr_low + tr_high) / 2.0
        
        println("Iteration $iter: Testing TR = $(round(tr_mid, digits=5))")
        current_r0 = get_avg_r0(configfile, tr_mid, nr_runs)
        
        println("Result: R0 = $(round(current_r0, digits=3))")

        if abs(current_r0 - target_r0) <= tolerance
            println("Converged. Transmission rate $tr_mid yields R0 ≈ $current_r0")
            best_tr = tr_mid
            break
        end

        if current_r0 < target_r0
            tr_low = tr_mid 
        else
            tr_high = tr_mid
        end
        
        if iter == max_iterations
            println("Max iterations reached. Best estimate is TR = $tr_mid")
            best_tr = tr_mid
        end
    end
    
    return best_tr
end

function fit_transmission_rate(configfile, target_r0; tr_start=0.05, tolerance=0.02, max_iterations=15, nr_runs=10)
    println("\nStarting linear fit for: $configfile")
    println("Target R0: $target_r0")
    
    current_tr = tr_start

    for iter in 1:max_iterations
        println("Iteration $iter: Testing TR = $(round(current_tr, digits=5))")
        
        current_r0 = get_avg_r0(configfile, current_tr, nr_runs)
        
        println("Result: R0 = $(round(current_r0, digits=3))")

        if abs(current_r0 - target_r0) <= tolerance
            println("Converged. Transmission rate $current_tr yields R0 ≈ $current_r0")
            return current_tr
        end

        current_tr = current_tr * (target_r0 / current_r0)
        
        if iter == max_iterations
            println("Max iterations reached. Best estimate is TR = $current_tr")
        end
    end
    
    return current_tr
end

# Main execution

# Group 1: Standard Municipality Contacts
config_group1 = "models/DE_model_R0_3.28.toml"
best_tr_group1 = fit_transmission_rate(config_group1, target_r0_value, tr_start=0.03)

# Group 2: Fewer Municipality Contacts
config_group2 = "models/DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml"
best_tr_group2 = fit_transmission_rate(config_group2, target_r0_value, tr_start=0.07)

println("\nFitting complete.")
println("1. Update 'DE_model_R0_3.28.toml' and 'DE_model_R0_3.28_Heinsberg_Seed.toml' with:")
println("   transmission_rate = $(round(best_tr_group1, digits=5))")
println("2. Update 'DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml' and 'DE_model_R0_3.28_Beckingen_Seed.toml' with:")
println("   transmission_rate = $(round(best_tr_group2, digits=5))")


output_summary = """
Fitting complete.
1. Update 'DE_model_R0_3.28.toml' and 'DE_model_R0_3.28_Heinsberg_Seed.toml' with:
   transmission_rate = $(round(best_tr_group1, digits=5))
2. Update 'DE_model_R0_3.28_Heinsberg_Seed_fewer_munContacts.toml' and 'DE_model_R0_3.28_Beckingen_Seed.toml' with:
   transmission_rate = $(round(best_tr_group2, digits=5))
"""

println("\n" * output_summary)
 
output_filename = "transmission_rate_results.txt"
open(output_filename, "w") do file
    write(file, output_summary)
end

println("Results have also been saved to: $output_filename")