export office_data
export convert_points

export plot_settings
export plot_all_settings
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



"""
    plot_settings(hhlds, sclasses, offcs; max_households = 2_000_000)

Plots setting locations as dots on a map.
Input must be setting-vectors of households, school classes, and offices.
The `max_households` parameter provides an upper limit for household
dots. The number of school classes and offices will be scaled accordingly. 
"""
function plot_settings(hhlds, sclasses, offcs; max_households = 2_000_000)

    # load counties from households
    df = county.(ags.([hhlds; sclasses; offcs])) |>
        vec -> unique(vec) |>
        vec -> DataFrame(
            ags = vec,
            cnt = fill(1, length(vec)))

    # calculate number of points that go onto the map
    hh_samples = min(length(hhlds), max_households)
    sample_ratio = hh_samples / length(hhlds)
    sc_samples = Int(ceil(length(sclasses) * sample_ratio))
    off_samples = Int(ceil(length(offcs) * sample_ratio))

    # geolocations
    hh_pnts = DataFrame(
            lat = lat.(hhlds),
            lon = lon.(hhlds)) |>
        df -> df[sample(1:nrow(df), hh_samples; replace=false), :] |>
        convert_points

    sc_pnts = DataFrame(
            lat = lat.(sclasses),
            lon = lon.(sclasses)) |>
        df -> df[sample(1:nrow(df), sc_samples; replace=false), :] |>
        convert_points

    off_pnts = DataFrame(
            lat = lat.(offcs),
            lon = lon.(offcs)) |>
        df -> df[sample(1:nrow(df), off_samples; replace=false), :] |>
        convert_points

    # build base map
    p = agsmap(df,
        fillcolor = :greys,
        clims = (0, 1),
        colorbar = false,
        fontfamily = "Times Roman",
        size = 5 .* map_dims())

    # add households
    scatter!(p, hh_pnts.lon, hh_pnts.lat,
        label = "Households",
        legendfontsize = 100,
        markerstrokewidth=0,
        markersize=1,
        markercolor = :black)

    # add schools
    scatter!(p, sc_pnts.lon, sc_pnts.lat,
        label = "Schools",
        legendfontsize = 100,
        markerstrokewidth=0,
        markersize=1,
        markercolor = :blue)

    # add workplaces
    scatter!(p, off_pnts.lon, off_pnts.lat,
        label = "Workplaces",
        legendfontsize = 100,
        markerstrokewidth=0,
        markersize=1,
        markercolor = :red)

    return p
end


"""
    plot_all_settings(sim; max_households = 2_000_000)

Plots all households, school classes, and offices as dots on a map 
for the provided simulation object.
The `max_households` parameter provides an upper limit for household
dots. The number of school classes and offices will be scaled accordingly. 
"""
function plot_all_settings(sim; max_households = 2_000_000)

    # load regions
    regions = county.(ags.(municipalities(sim))) |> unique

    # load settings (filer for settings in
    # any of the ags in municipality list)
    hhlds = households(sim)
    # filter out schools that are not in the current model scope
    sclasses = schoolclasses(sim) |>
        sc -> sc[(c -> ags(c) in regions).(sc)]
    # filter out offices that are not in the current model scope
    offcs = offices(sim) |>
        off -> off[(c -> ags(c) in regions).(off)]

    return plot_settings(hhlds, sclasses, offcs, max_households = max_households)

end


# RUN ALL ANALYSES

"""
    settinganalyses(;
        folder = nothing,
        folder_tag = nothing,
        sim = Simulation(population = "SL"))

Runs all setting-related analyses.

If `folder` is `nothing`, a new result folder is created, with `folder_tag`
appended to its timestamp-based name (see `result_folder`).
"""
function settinganalyses(;
        folder = nothing,
        folder_tag = nothing,
        sim = Simulation(population = "SL", seed = 0))

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
