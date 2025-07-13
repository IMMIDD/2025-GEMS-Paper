export mun_commutes
export county_commutes
export plot_mun_commutes_fraction
export plot_county_commutes_fraction

"""
    mun_commutes(sim)

Returns a dataframe with within- and cross-municipalitiy commutes.
Input must be `Simulation` object.
"""
function mun_commutes(sim)

    # households
    hhlds = DataFrame(
        id = id.(households(sim)),
        hh_ags = ags.(households(sim))
    )

    # offices
    offs = DataFrame(
        id = id.(offices(sim)),
        off_ags = ags.(offices(sim))
    )

    # individuals
    inds = DataFrame(
        hh_id = household_id.(individuals(sim)),
        off_id = office_id.(individuals(sim))
    )

    return inds |>
        df -> leftjoin(df, hhlds, on = :hh_id => :id) |>
        df -> leftjoin(df, offs, on = :off_id => :id) |>
        # filter jobless people
        df -> df[.!ismissing.(df.off_ags), :] |> 
        df -> groupby(df, [:hh_ags, :off_ags]) |>
        df -> combine(df, nrow => :commutes)

end


"""
    county_communtes(mun_commutes_df)

Returns a dataframe with within- and cross-county commutes.
Input must be cross-municipality dataframe (output of `mun_commutes()` function).
"""
function county_commutes(mun_commutes_df)
    return mun_commutes_df |>
        # get county part of all municipality AGSs
        df -> transform(df,
            :hh_ags => (h -> county.(h)) => :hh_ags,
            :off_ags => (o -> county.(o)) => :off_ags) |>
        df -> groupby(df, [:hh_ags, :off_ags]) |>
        df -> combine(df, :commutes => sum => :commutes)
end


"""
    plot_mun_commutes_fraction(mun_commutes_df)

Returns a map with the fraction of (working) people whose
workplace (office) is in a different municipality from their
household. Input must be cross-municipality dataframe (output of `mun_commutes()` function).
"""
function plot_mun_commutes_fraction(mun_commutes_df)
    return mun_commutes_df |>
        # calculate fraction of outgoing commutes
        df -> groupby(df, :hh_ags) |>
        df -> combine(df, [:hh_ags, :off_ags, :commutes] => ((h, o, c) -> sum(c[o .!= h[1]]) / sum(c)) => :out_fraction) |>
        df -> rename!(df, :hh_ags => :ags) |>
        # mapping
        df -> transform(df, :ags => (a -> convert(Vector{AGS} , a)) => :ags) |>
        df -> prepare_map_df!(df, level = 3) |>
        df -> agsmap(df, linecolor = nothing, clims = (0, 1))
end

"""
    plot_county_commutes_fraction(county_commutes_df)

Returns a map with the fraction of (working) people whose
workplace (office) is in a different county from their
household. Input must be cross-municipality dataframe (output of `county_commutes()` function).
"""
function plot_county_commutes_fraction(county_commutes_df)
    return county_commutes_df |>
        # calculate fraction of outgoing commutes
        df -> groupby(df, :hh_ags) |>
        df -> combine(df, [:hh_ags, :off_ags, :commutes] => ((h, o, c) -> sum(c[o .!= h[1]]) / sum(c)) => :out_fraction) |>
        df -> rename!(df, :hh_ags => :ags) |>
        # mapping
        df -> transform(df, :ags => (a -> convert(Vector{AGS} , a)) => :ags) |>
        df -> prepare_map_df!(df, level = 2) |>
        df -> agsmap(df, linecolor = nothing, clims = (0, 1))
end

"""
    mobility(;
        folder = nothing,
        sim = Simulation(population = "SL"))

Runs all mobilty-related analyses.
"""
function mobility(;
        folder = nothing,
        sim = Simulation(population = "SL"))

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder() : nothing

    # calculate municipality & county commutes
    mc = mun_commutes(sim)
    cc = county_commutes(mc)
    
    # plot fraction of outgoing commuters per municipality
    p_mc = plot_mun_commutes_fraction(mc)
    png(p_mc, joinpath(folder, "outgoing_commuters_fraction_per_municipality.png"))

    # plot fraction of outgoing commuters per county
    p_cc = plot_county_commutes_fraction(cc)
    png(p_cc, joinpath(folder, "outgoing_commuters_fraction_per_county.png"))

end