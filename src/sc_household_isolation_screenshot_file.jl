    using GEMS # load GEMS package
    
    sim_baseline = Simulation(label = "Baseline")   # initialize baseline simulation
    run!(sim_baseline)                              # run baseline simulation
    rd_baseline = ResultData(sim_baseline)          # process baseline simulation results

    sim_scenario = Simulation(label = "Scenario")                      # initialize scenario simulation
    member_isolation = IStrategy("Isolation", sim_scenario)            # define isolation strategy
    add_measure!(member_isolation, SelfIsolation(14))                  # add 14-day isolation measure
    member_identification = IStrategy("Find Household", sim_scenario)  # define household strategy
    add_measure!(member_identification,
        FindSettingMembers(Household, member_isolation))               # isolate all household members
    trigger = SymptomTrigger(member_identification)                    # trigger upon symptom onset
    add_symptom_trigger!(sim_scenario, trigger)                        # add trigger to scenario simulation

    run!(sim_scenario)                      # run scenario simulation
    rd_scenario = ResultData(sim_scenario)  # process scenario simulation results

    gemsplot([rd_baseline, rd_scenario])    # plot results of both simulations



    