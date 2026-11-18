using CSV, DataFrames, Plots, StatsPlots, LaTeXStrings, CategoricalArrays,  Statistics

function combine_and_violin(parent_folder::String, subfolder_paths::Vector{String})
    all_dfs = DataFrame[]

    for folder in subfolder_paths
        csv_path = joinpath(folder, "output_zero_intervals_1000-10000.csv")
        if isfile(csv_path)
            println("Reading: $csv_path")
            push!(all_dfs, CSV.read(csv_path, DataFrame))
        else
            @warn "File not found, skipping: $csv_path"
        end
    end

    if isempty(all_dfs)
        error("No CSV files found.")
    end

    final_df = vcat(all_dfs..., cols=:union)
    final_path = joinpath(parent_folder, "final_all_zero_intervals_R1.5_v5.0_pf0.1_1000-10000.csv")
    CSV.write(final_path, final_df)
    println("Wrote combined file → $final_path")

    clean = dropmissing(final_df, [:aspect, :interval])
    clean = filter(row -> row.interval > 20, clean)

    # round first, then make categorical
    clean.aspect_rounded = round.(clean.aspect, digits=2)
    clean.aspect_cat = categorical(clean.aspect_rounded, ordered=true)

    # four shades of blue
    blue_shades = [colorant"#a6cee3", colorant"#1f78b4", colorant"#08519c", colorant"#08306b"]

    p = violin(clean.aspect_cat, clean.interval,
               xlabel = "b/a",
               ylabel = L"$T_{cross(\kappa = 0})$",
               title  = "R=1.5, v = 5.0, pf= 0.1",
               legend = false,
               linewidth = 1,
               trim = true,
               fillcolor = blue_shades,
               linecolor = :black,
               xticks = ([0.1, 0.25, 0.4, 0.55], ["0.1", "0.25", "0.4", "0.55"]))

    # overlay medians
    medians = combine(groupby(clean, :aspect_cat), :interval => median => :med)
    scatter!(p, medians.aspect_cat, medians.med,
             marker = (:circle, 6, :red), label = "")

    plot_path = joinpath(parent_folder, "zero_intervals_violin_vs_aspect_1000_10000.png")
    savefig(p, plot_path)
    println("Saved violin plot → $plot_path")

    return final_df, p
end

parent_folder = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.1\\"

analyzed_folders = [
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.1\aspect 0.1\\",
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.1\aspect 0.25\\",
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.1\aspect 0.4\\",
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.1\aspect 0.55\\"
]

combine_and_violin(parent_folder, analyzed_folders)





