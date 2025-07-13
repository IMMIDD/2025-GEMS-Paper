export individual_data
export plot_age_per_state

export demographics


# DATA PREPARATION

function individual_data(sim)
    return individuals(sim) |>
        inds -> DataFrame(
            id = id.(inds),
            age = age.(inds),
            sex = sex.(inds),
            ags = (i -> ags(household(i, sim))).(inds)
        )
end


# PLOTTING

function plot_age_per_state(individual_data_df)

    function buildplot(ages, plts)
        p = histogram(ages,
            fillcolor = :blue,        # Fill color of bars
            linecolor = :black,       # Outline color of bars
            linewidth = 1.5,          # Width of the bar outline
            fillalpha = 0.0,
            xlabel = "",
            ylabel = "",
            legend = false,
            xlims = (0, 100),
            yticks = nothing,
            #linewidth = 0,
            fontfamily = "Times Roman",
            bins = 0:10:100,
        )

        push!(plts, p)
    end

    state_names = germanshapes(1) |>
        shps -> DataFrame(ags = AGS.(shps.AGS_0), gen = shps.GEN) |>
        # some shapefiles have mulitple entries
        df -> groupby(df, :ags) |>
        df -> combine(df, :gen => first => :gen)

    plts = []

    individual_data_df |>
        df -> transform(df, :ags => (a -> state.(a)) => :ags) |>
        df -> groupby(df, :ags) |>
        df -> combine(df, :age => (ages -> buildplot(ages, plts)) => :age_plot)# |>
        # add state names
        #df -> leftjoin(df, state_names, on = :ags) |>
        #df -> transform(df, [:age_plot, :gen] => ByRow((plt, name) -> plot!(plt, title = name)) => :age_plot) |>
        # df -> plot((df.age_plot)...,
        #     layout = (4, 4),
        #     size = (1000, 1000))

        return plot(plts...,
            layout = (4, 4),
            size = (1000, 1000))

end



# RUN ALL ANALYSES

"""
    demographics(;
        folder = nothing,
        sim = Simulation(population = "SL"))

Runs all mobilty-related analyses.
"""
function demographics(;
        folder = nothing,
        sim = Simulation(population = "SL"))

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder() : nothing

    # calculate individual data
    id = individual_data(sim)
    
    # plot age structure per state
    p_age_state = plot_age_per_state(id)
    png(p_age_state, joinpath(folder, "age_structure_per_state.png"))
    
end