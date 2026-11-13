# plotting from different folders and combing that data to make one data fike
# then averaged the data and plot the mean ± SD of the data
# useful for one parameter change analysis
using CSV, DataFrames, Plots, Statistics

function combine_and_plot(parent_folder::String, subfolder_paths::Vector{String})
    all_dfs = DataFrame[]

    for folder in subfolder_paths
        csv_path = joinpath(folder, "output_prob_no_curvature.csv")

        if isfile(csv_path)
            println("Reading: $csv_path")
            df = CSV.read(csv_path, DataFrame)
            push!(all_dfs, df)
        else
            @warn "File not found, skipping: $csv_path"
        end
    end

    if isempty(all_dfs)
        error("None of the requested output_prob_no_curvature.csv files were found.")
    end

    # Concatenate all selected files
    final_df = vcat(all_dfs..., cols=:union)
    final_path = joinpath(parent_folder, "final_all_no_curvature_output_R1.5_v5.0_pf0.05.csv")
    CSV.write(final_path, final_df)
    println("Wrote combined file → $final_path")
    println("Total rows: $(nrow(final_df))")
    xcol= :aspect
    ycol= :prob_no_curv
    # Mean ± SD of Curvature_Effect vs pf
    if !(xcol in propertynames(final_df)) || !(ycol in propertynames(final_df))
        error("Required columns $xcol and $ycol not found")
    end

    clean = dropmissing(final_df, [xcol, ycol])
    grouped = groupby(clean, xcol)
    summary = combine(grouped,
                      ycol => mean   => :mean_curv,
                      ycol => std    => :sd_curv,
                      ycol => length => :n)
    sort!(summary, xcol)

    println("Summary (mean ± SD):")
    println(summary)

    # Plot
    p = plot(summary.aspect, summary.mean_curv,
             yerror = summary.sd_curv,
             xlabel = "aspect ratio (b/a)",
             ylabel = "Probability of No Curvature",
             title  = "R=1.5, v = 5.0, pf= 0.05",
             marker = :circle,
             markersize = 7,
             linewidth = 2,
             legend = false,
             seriestype = :scatter)

    plot!(p, summary.aspect, summary.mean_curv, linewidth=1.5, color=:black, alpha=0.5)

    plot_path = joinpath(parent_folder, "nocurvatureprob_vs_aspect_R1.5_v5.0.png")
    savefig(p, plot_path)
    println("Saved plot → $plot_path")

    return final_df, summary, p
end

# ======================  RUN THIS  ======================
parent_folder = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.3\\"

analyzed_folders = [
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.3\aspect 0.1\\",
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.3\aspect 0.25\\",
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.3\aspect 0.4\\",
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.3\aspect 0.55\\"

    # add the other folders you want below (copy the line and change the path)
    #raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.2\aspect 0.25",
    #raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.3\aspect 0.25",
    #raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.4\aspect 0.25"
]

combine_and_plot(parent_folder, analyzed_folders)
