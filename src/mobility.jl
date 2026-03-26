export mun_commutes
export county_commutes
export mun_commutes_out_fraction
export county_commutes_out_fraction
export mun_commutes_in_fraction
export county_commutes_in_fraction

export plot_mun_commutes_out_fraction
export plot_county_commutes_out_fraction
export plot_mun_commutes_in_fraction
export plot_county_commutes_in_fraction
export plot_county_commutes_combined

export mobility

# DATA PREPARATION

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
    mun_commutes_out_fraction(mun_commutes_df)

Returns a dataframe with the fraction of (working) people whose
workplace (office) is in a different municipality from their
household. Input must be cross-municipality dataframe (output of `mun_commutes()` function).
"""
function mun_commutes_out_fraction(mun_commutes_df)
    return mun_commutes_df |>
        # calculate fraction of outgoing commutes
        df -> groupby(df, :hh_ags) |>
        df -> combine(df, [:hh_ags, :off_ags, :commutes] => ((h, o, c) -> 100 * sum(c[o .!= h[1]]) / sum(c)) => :out_fraction)
end


"""
    county_commutes_out_fraction(county_commutes_df)

Returns a dataframe with the fraction of (working) people whose
workplace (office) is in a different county from their
household. Input must be cross-municipality dataframe (output of `county_commutes()` function).
"""
function county_commutes_out_fraction(county_commutes_df)
    return county_commutes_df |>
        # calculate fraction of outgoing commutes
        df -> groupby(df, :hh_ags) |>
        df -> combine(df, [:hh_ags, :off_ags, :commutes] => ((h, o, c) -> 100 * sum(c[o .!= h[1]]) / sum(c)) => :out_fraction) |>
        df -> leftjoin(df, county_data(), on = [:hh_ags => :ags])
end


"""
    mun_commutes_out_fraction(mun_commutes_df)

Returns a dataframe with the fraction of people who work in this
municipality and live in another.
Input must be cross-municipality dataframe (output of `mun_commutes()` function).
"""
function mun_commutes_in_fraction(mun_commutes_df)
    return mun_commutes_df |>
        # calculate fraction of incoming commutes
        df -> groupby(df, :off_ags) |>
        df -> combine(df, [:hh_ags, :off_ags, :commutes] => ((h, o, c) -> 100 * sum(c[h .!= o[1]]) / sum(c)) => :in_fraction)
end


"""
    county_commutes_in_fraction(county_commutes_df)

Returns a dataframe with the fraction of people who work in this
county and live in another.
Input must be cross-municipality dataframe (output of `county_commutes()` function).
"""
function county_commutes_in_fraction(county_commutes_df)
    return county_commutes_df |>
        # calculate fraction of outgoing commutes
        df -> groupby(df, :off_ags) |>
        df -> combine(df, [:hh_ags, :off_ags, :commutes] => ((h, o, c) -> 100 * sum(c[h .!= o[1]]) / sum(c)) => :in_fraction) |>
        df -> leftjoin(df, county_data(), on = [:off_ags => :ags])
end


# PLOTTING

"""
    plot_mun_commutes_out_fraction(mun_commutes_out_fraction_df)

Returns a map with the fraction of (working) people whose
workplace (office) is in a different municipality from their
household. Input must be municipality out fraction dataframe
(output of `mun_commutes_out_fraction()` function).
"""
function plot_mun_commutes_out_fraction(mun_commutes_out_fraction_df)
    return mun_commutes_out_fraction_df |>
        df -> rename(df, :hh_ags => :ags) |>
        # mapping
        df -> transform(df, :ags => (a -> convert(Vector{AGS} , a)) => :ags) |>
        df -> prepare_map_df(df, level = 3) |>
        df -> agsmap(df,
            size = map_dims(),
            titlefontsize = map_titlefontsize(),
            ytickfontsize = map_legendfontsize(),
            colorbar_titlefontsize = map_legendfontsize(),
            right_margin=5mm,
            linecolor = nothing,
            clims = (0, 100),
            colorbar_title = "Percentage (%)", 
            fontfamily = "Times Roman",
            title = "Fraction of working people commuting \n to a different municipality")
end


"""
    plot_county_commutes_out_fraction(county_commutes_out_fraction_df)

Returns a map with the fraction of (working) people whose
workplace (office) is in a different county from their
household. Input must be county out fraction dataframe
(output of `county_commutes_out_fraction()` function).
"""
function plot_county_commutes_out_fraction(county_commutes_out_fraction_df)
    return county_commutes_out_fraction_df |>
        df -> rename(df, :hh_ags => :ags) |>
        df -> DataFrames.select(df, :ags, :out_fraction) |>
        # mapping
        df -> transform(df, :ags => (a -> convert(Vector{AGS} , a)) => :ags) |>
        df -> prepare_map_df(df, level = 2) |>
        df -> agsmap(df,
            size = map_dims(),
            titlefontsize = map_titlefontsize(),
            ytickfontsize = map_legendfontsize(),
            colorbar_titlefontsize = map_legendfontsize(),
            right_margin=5mm,
            linecolor = nothing,
            clims = (0, 75),
            colorbar_title = "Percentage (%)",
            fontfamily = "Times Roman",
            title = "Fraction of working people commuting \n to a different county")
end


"""
    plot_mun_commutes_out_fraction(mun_commutes_in_fraction_df)

Returns a map with the fraction of people who work in this
municipality and live in another.
Input must be municipality in fraction dataframe
(output of `mun_commutes_in_fraction()` function).
"""
function plot_mun_commutes_in_fraction(mun_commutes_in_fraction_df)
    return mun_commutes_in_fraction_df |>
        df -> rename(df, :off_ags => :ags) |>
        # mapping
        df -> transform(df, :ags => (a -> convert(Vector{AGS} , a)) => :ags) |>
        df -> prepare_map_df(df, level = 3) |>
        df -> agsmap(df,
            size = map_dims(),
            titlefontsize = map_titlefontsize(),
            ytickfontsize = map_legendfontsize(),
            colorbar_titlefontsize = map_legendfontsize(),
            right_margin=5mm,
            linecolor = nothing,
            clims = (0, 100),
            colorbar_title = "Percentage (%)",
            fontfamily = "Times Roman",
            title = "Fraction of workforce commuting \n from a different municipality")
end


"""
    plot_county_commutes_in_fraction(county_commutes_in_fraction_df)

Returns a map with the fraction of people who work in this
county and live in another.
Input must be county in fraction dataframe 
(output of `county_commutes_in_fraction()` function).
"""
function plot_county_commutes_in_fraction(county_commutes_in_fraction_df)
    return county_commutes_in_fraction_df |>
        df -> rename(df, :off_ags => :ags) |>
        df -> DataFrames.select(df, :ags, :in_fraction) |>
        # mapping
        df -> transform(df, :ags => (a -> convert(Vector{AGS} , a)) => :ags) |>
        df -> prepare_map_df(df, level = 2) |>
        df -> agsmap(df,
            size = map_dims(),
            titlefontsize = map_titlefontsize(),
            ytickfontsize = map_legendfontsize(),
            colorbar_titlefontsize = map_legendfontsize(),
            right_margin=5mm,
            linecolor = nothing,
            clims = (0, 75),
            colorbar_title = "Percentage (%)", 
            fontfamily = "Times Roman",
            title = "Fraction of workforce commuting \n from a different county")
end


"""
    plot_county_commutes_combined(p_in_cc, p_out_cc)

Returns the combined plot of county in- and out commutes.
Input must be the in- and out commute plots (output of
`plot_county_commutes_in_fraction()` and
`plot_county_commutes_out_fraction` function).
"""
function plot_county_commutes_combined(p_in_cc, p_out_cc)
    l = @layout [a{0.41w} b{0.59w}]
    
    return plot(plot!(p_in_cc, colorbar = false), p_out_cc,
        titlefontsize = map_titlefontsize(),
        ytickfontsize = map_legendfontsize(),
        colorbar_titlefontsize = map_legendfontsize(),
        right_margin=5mm,
        layout = l,
        size = (2 * map_dims()[1], map_dims()[2]))
end

# RUN ALL ANALYSES

"""
    mobility(;
        folder = nothing,
        sim = Simulation(population = "SL"))

Runs all mobilty-related analyses.
"""
function mobility(;
        folder = nothing,
        sim = Simulation(population = "SL", seed = 0))

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder() : nothing

    # calculate municipality & county commutes
    mc = mun_commutes(sim)
    cc = county_commutes(mc)

    mcif = mun_commutes_in_fraction(mc)
    mcof = mun_commutes_out_fraction(mc)
    ccif = county_commutes_in_fraction(cc)
    ccof = county_commutes_out_fraction(cc)
    
    # store data
    CSV.write(joinpath(folder, "county_commutes_in_fraction.csv"), ccif)
    CSV.write(joinpath(folder, "county_commutes_out_fraction.csv"), ccof)

    # plot fraction of outgoing commuters per municipality
    # p_out_mc = plot_mun_commutes_out_fraction(mcof)
    # png(p_out_mc, joinpath(folder, "outgoing_commuters_fraction_per_municipality.png"))

    # plot fraction of outgoing commuters per county
    p_out_cc = plot_county_commutes_out_fraction(ccof)
    png(p_out_cc, joinpath(folder, "outgoing_commuters_fraction_per_county.png"))

    # plot fraction of workforce as incoming commuters per municipality
    # p_in_mc = plot_mun_commutes_in_fraction(mcif)
    # png(p_in_mc, joinpath(folder, "incoming_commuters_fraction_per_municipality.png"))

    # plot fraction of workforce as incoming commuters per county
    p_in_cc = plot_county_commutes_in_fraction(ccif)
    png(p_in_cc, joinpath(folder, "incoming_commuters_fraction_per_county.png"))

    p_comb = plot_county_commutes_combined(p_in_cc, p_out_cc)
    png(p_comb, joinpath(folder, "in_and_out_commuter_fraction_per_county.png"))
end