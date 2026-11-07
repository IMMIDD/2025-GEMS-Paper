export plot_contact_matrix
export plot_contact_matrices

export contacts



"""
    plot_contact_matrix(rd, type)

Returns a contact matrix for a specified setting type (e.g. Household).
Needs an ResultData object that contains contact sampling data.
"""
function plot_contact_matrix(rd, type)
    return gemsplot(rd, type = :AggregatedSettingAgeContacts, settingtype = type,
        title = "$type Contacts",    
        plot_title = "", 
        xlabel = "Age",
        ylabel = "Contact Age",
        show_values = false,
        #margin = 5mm,
        #top_margin = 10mm,
        #bottom_margin = 25mm,
        left_margin = 10mm,
        right_margin = 15mm,
        titlefontsize = 30,
        xtickfontsize = 18,
        ytickfontsize = 18,
        colorbar_titlefontsize = 18,
        xlabelfontsize = 22,
        ylabelfontsize = 22,
        size = (1000, 900),
        xrotation = 45)
end


"""
    plot_contact_matrices(rd)

Returns a multiplot of contact matrices.
Needs an ResultData object that contains contact sampling data.
"""
function plot_contact_matrices(rd)
    plts = [plot_contact_matrix(rd, t) for t in
        [
            Household, Municipality,
            # Office, SchoolClass,
            # Department, SchoolYear,
            Workplace, School,
            # WorkplaceSite, SchoolComplex
        ]]

    return plot(plts...,
        layout = (2, 2),
        size = (2400, 2000))
end


"""
    contacts(;
        folder = nothing,
        folder_tag = nothing,
        sim = Simulation(population = "SL"))

Runs all contact-related analyses.

If `folder` is `nothing`, a new result folder is created, with `folder_tag`
appended to its timestamp-based name (see `result_folder`).
"""
function contacts(;
        folder = nothing,
        folder_tag = nothing,
        sim = Simulation(population = "SL", seed = 0))

    # create result folder if not passed
    isnothing(folder) ? folder = result_folder(folder_tag) : nothing

    # build ResultData Object for contacts only
    rd = ResultData(sim, style = "ContactsOnlyStyle")
    p_cm = plot_contact_matrices(rd)
    png(p_cm, joinpath(folder, "contact_matrices.png"))

end