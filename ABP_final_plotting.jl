using CSV, DataFrames, Plots, Statistics

function combine_and_plot(parent_folder::String, subfolder_paths::Vector{String})
    all_dfs = DataFrame[]

    for folder in subfolder_paths
        csv_path = joinpath(folder, "output_combined.csv")

        if isfile(csv_path)
            println("Reading: $csv_path")
            df = CSV.read(csv_path, DataFrame)
            push!(all_dfs, df)
        else
            @warn "File not found, skipping: $csv_path"
        end
    end

    if isempty(all_dfs)
        error("None of the requested output_combined.csv files were found.")
    end

    # Concatenate all selected files
    final_df = vcat(all_dfs..., cols=:union)
    final_path = joinpath(parent_folder, "final_all_output_R2.0.csv")
    CSV.write(final_path, final_df)
    println("Wrote combined file → $final_path")
    println("Total rows: $(nrow(final_df))")

    # Mean ± SD of Curvature_Effect vs pf
    if !(:pf in propertynames(final_df)) || !(:Curvature_Effect in propertynames(final_df))
        error("Required columns :pf and :Curvature_Effect not found")
    end

    clean = dropmissing(final_df, [:pf, :Curvature_Effect])
    grouped = groupby(clean, :v)
    summary = combine(grouped,
                      :Curvature_Effect => mean   => :mean_curv,
                      :Curvature_Effect => std    => :sd_curv,
                      :Curvature_Effect => length => :n)
    sort!(summary, :v)

    println("Summary (mean ± SD):")
    println(summary)

    # Plot
    p = plot(summary.v, summary.mean_curv,
             yerror = summary.sd_curv,
             xlabel = "v",
             ylabel = "Curvature_Effect",
             title  = "R=1.5, a=50.0, b=12.5, pf=0.1",
             marker = :circle,
             markersize = 7,
             linewidth = 2,
             legend = false,
             seriestype = :scatter)

    plot!(p, summary.v, summary.mean_curv, linewidth=1.5, color=:black, alpha=0.5)

    plot_path = joinpath(parent_folder, "curvature_vs_v_R2.0_a50.0_b12.5_pf0.1.png")
    savefig(p, plot_path)
    println("Saved plot → $plot_path")

    return final_df, summary, p
end

# ======================  RUN THIS  ======================
parent_folder = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data"

analyzed_folders = [
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R2.0 pf 0.1\aspect 0.25\v 5.0\\",
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R2.0 pf 0.1\aspect 0.25\v 10.0\\"
    # add the other folders you want below (copy the line and change the path)
    #raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.2\aspect 0.25",
    #raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.3\aspect 0.25",
    #raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.4\aspect 0.25"
]

combine_and_plot(parent_folder, analyzed_folders)
