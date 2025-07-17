export ContactsOnlyStyle

mutable struct ContactsOnlyStyle <: GEMS.ResultDataStyle
    data::Dict{String, Any}
    function ContactsOnlyStyle(pP::PostProcessor)
        funcs = Dict(
            "meta_data" =>
                Dict(),
            "sim_data" =>
                Dict(
                    "setting_data" => () -> pP |> settingdata,
                    "final_tick" => () -> pP |> simulation |> tick,
                ),
            
            # system data
            "system_data" => 
                Dict(),

            "setting_age_contacts" =>
                Dict(),

            "aggregated_setting_age_contacts" =>
                Dict(
                    # TODO: interval_steps shouldn't be hard coded. They rather should be defined in the config file.
                    # TODO: This list should be determined dynamically depending on what settings are present in the simulation
                    "Household" => () -> mean_contacts_per_age_group(pP, Household, 5),
                    "SchoolClass" => () -> mean_contacts_per_age_group(pP, SchoolClass, 5),
                    "SchoolYear" => () -> mean_contacts_per_age_group(pP, SchoolYear, 5),
                    "School" => () -> mean_contacts_per_age_group(pP, School, 5),
                    "SchoolComplex" => () -> mean_contacts_per_age_group(pP, SchoolComplex, 5),
                    "Office" => () -> mean_contacts_per_age_group(pP, Office, 5), 
                    "Department" => () -> mean_contacts_per_age_group(pP, Department, 5), 
                    "Workplace" => () -> mean_contacts_per_age_group(pP, Workplace, 5), 
                    "WorkplaceSite" => () -> mean_contacts_per_age_group(pP, WorkplaceSite, 5), 
                    "Municipality" => () -> mean_contacts_per_age_group(pP, Municipality, 5)
                ),

            "dataframes" =>
                Dict()      
        )

        # call all provided functions and replace
        # the dicts with their return values
        return(
            new(GEMS.process_funcs(funcs))
        )
    end
end