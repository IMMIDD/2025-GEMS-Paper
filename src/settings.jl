export office_data
export convert_points

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


function convert_points(df)
    # projection transformation
    wgs84 = "+proj=longlat +datum=WGS84 +no_defs"
    dhdn_gk3 = """
    +proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel
    +datum=potsdam +units=m +no_defs
    """
    trans = Proj.Transformation(wgs84, dhdn_gk3)

    # convert points
    return trans.(df.lon, df.lat) |>
        tpl -> DataFrame(
            lon = first.(tpl),
            lat = last.(tpl))

end


# RUN ALL ANALYSES

"""
    settinganalyses(;
        folder = nothing,
        folder_tag = nothing,
        sim = Simulation(population = "DE"))

Runs all setting-related analyses.

If `folder` is `nothing`, a new result folder is created, with `folder_tag`
appended to its timestamp-based name (see `result_folder`).
"""
function settinganalyses(;
        folder = nothing,
        folder_tag = nothing,
        sim = Simulation(population = "DE", seed = 0))

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder(folder_tag) : nothing

    # calculate muncipality & county setting sizes
    od_m = office_data(sim)

    # plot average office sizes
    p_od_m = plot_mun_office_sizes(od_m)
    png(p_od_m, joinpath(folder, "mun_mean_office_sizes.png"))
    p_od_c = plot_county_office_sizes(od_m)
    png(p_od_c, joinpath(folder, "county_mean_office_sizes.png"))


end
