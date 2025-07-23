    using GEMS # load GEMS package
    
    sim_baseline = Simulation(label = "Baseline")   # initialize baseline simulation
    run!(sim_baseline)                              # run baseline simulation
    rd_baseline = ResultData(sim_baseline)          # process baseline simulation results

    sim_scenario = Simulation(label = "Scenario")               # initialize scenario simulation
    self_isolation = IStrategy("Self Isolation", sim_scenario)  # define self-isolation strategy
    add_measure!(self_isolation, SelfIsolation(14))             # add 14-day self-isolation measure to strategy
    trigger = SymptomTrigger(self_isolation)                    # define trigger for self-isolation based on symptoms
    add_symptom_trigger!(sim_scenario, trigger)                 # add trigger to scenario simulation

    run!(sim_scenario)                      # run scenario simulation
    rd_scenario = ResultData(sim_scenario)  # process scenario simulation results

    gemsplot([rd_baseline, rd_scenario])    # plot results of both simulations



    