using GEMSPaper, DataFrames, Statistics, TOML

# Calibration of the first wave scenario (not part of runall()):
# 1. transmission rate (or scale) for R0 = 3.28, as in fit_r0.jl
# 2. generation time of the model and contact factors of the phases (written to inputs/)

# Global settings
configfile = "models/DE_model_first_wave.toml"
target_r0_value = 3.28

# population of the transmission rate fit: contact structure closest to Germany
population = "RP"

# days per run: the secondary infections of the infections up to tick 10 are complete
days = 50


# Run simulation and calculate average R0
function get_avg_r0(configfile, tr, runs)
    # config copy with 0.1 % random initial infections (the first wave config has none)
    config = TOML.parsefile(configfile)
    config["Simulation"]["StartCondition"]["parameters"]["fraction"] = 0.001
    fit_configfile = tempname() * ".toml"
    open(fit_configfile, "w") do io
        TOML.print(io, config)
    end

    r0_sum = 0.0
    for _ in 1:runs
        sim = Simulation(configfile = fit_configfile, population = population)
        sim.label = "TR: $tr"
        pathogen(sim).transmission_function = with_transmission_scale(pathogen(sim).transmission_function, tr)
        sim.stop_criterion = TimesUp(limit = days)
        run!(sim)

        effective_r = effectiveR(PostProcessor(sim))
        r0_sum += mean(effective_r.effective_R[1:10])
    end
    rm(fit_configfile)
    return r0_sum / runs
end

# Linear search for the transmission rate (as in fit_r0.jl)
function fit_transmission_rate(configfile, target_r0;
        tr_start,
        tolerance = 0.02,
        max_iterations = 15,
        nr_runs = 10)

    current_tr = tr_start
    for iter in 1:max_iterations
        println("Iteration $iter: Testing TR = $(round(current_tr, digits = 5))")
        current_r0 = get_avg_r0(configfile, current_tr, nr_runs)
        println("Result: R0 = $(round(current_r0, digits = 3))")

        if abs(current_r0 - target_r0) <= tolerance
            return current_tr
        end
        current_tr = current_tr * (target_r0 / current_r0)
    end
    return current_tr
end


# Main execution

# 1. transmission rate (transmission_scale of an InfectiousnessProfile)
tf_parameters = TOML.parsefile(configfile)["Pathogens"]["Covid19"]["transmission_function"]["parameters"]
tr_name = haskey(tf_parameters, "transmission_scale") ? "transmission_scale" : "transmission_rate"
config_tr = tf_parameters[tr_name]
best_tr = fit_transmission_rate(configfile, target_r0_value, tr_start = config_tr)
println("\n$tr_name for R0 = $target_r0_value: $(round(best_tr, digits = 5)) " *
    "(config: $config_tr)")

# the calibration needs the fitted rate in the config
if abs(best_tr / config_tr - 1) > 0.01
    error("Enter $tr_name = $(round(best_tr, digits = 5)) in $configfile and run again.")
end

# 2. generation time and contact factors (on Germany)
calibrate_first_wave(configfile = configfile)
