export office_data
export plot_mun_office_sizes
export plot_county_office_sizes

export settinganalyses

# DATA PREPARATION

function office_data(sim)
    return DataFrame(
        id = id.(offices(sim)),
        ags = ags.(offices(sim)),
        size = size.(offices(sim))
    )
end


# PLOTTING


function plot_mun_office_sizes(office_data_df)
    return office_data_df |>
        df -> groupby(df, :ags) |>
        df -> combine(df, :size => mean => :mean_size) |>
        df -> agsmap(df,
            linecolor = nothing,
            fontfamily = "Times Roman",
            clims = (0, maximum(df.mean_size)),
            title = "Mean office size per municipality")
end


function plot_county_office_sizes(office_data_df)
    return office_data_df |>
        df -> transform(df, :ags => (a -> county.(a)) => :ags) |>
        df -> groupby(df, :ags) |>
        df -> combine(df, :size => mean => :mean_size) |>
        df -> agsmap(df,
            linecolor = nothing,
            fontfamily = "Times Roman",
            clims = (0, maximum(df.mean_size)),
            title = "Mean office size per county")
end



# RUN ALL ANALYSES

"""
    settinganalyses(;
        folder = nothing,
        sim = Simulation(population = "SL"))

Runs all setting-related analyses.
"""
function settinganalyses(;
        folder = nothing,
        sim = Simulation(population = "SL"))

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder() : nothing

    # calculate muncipality & county setting sizes
    od_m = office_data(sim)

    # plot average office sizes
    p_od_m = plot_mun_office_sizes(od_m)
    png(p_od_m, joinpath(folder, "mun_mean_office_sizes.png"))
    p_od_c = plot_county_office_sizes(od_m)
    png(p_od_c, joinpath(folder, "county_mean_office_sizes.png"))


end
