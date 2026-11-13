using CSV, DataFrames, Plots, Statistics

function plot_curvature_vs_aspect_multi_pf(parent_folder::String, pf_files::Vector)
    colors  = [:blue, :red, :green, :orange, :purple]
    markers = [:circle, :diamond, :utriangle, :square, :star5]

    p = plot(xlabel = "aspect ratio (b/a)",
             ylabel = "No Curvature Probability",
             title  = "R=1.5, v=5.0",
             legend = :topright,
             grid   = true)

    for (idx, (pf_val, csv_path)) in enumerate(pf_files)
        if !isfile(csv_path)
            @warn "File not found, skipping pf=$pf_val: $csv_path"
            continue
        end

        println("Reading (pf=$pf_val): $csv_path")
        df = CSV.read(csv_path, DataFrame)
         xcol= :aspect
         ycol= :prob_no_curv
        if !(xcol in propertynames(df)) || !(ycol in propertynames(df))
            @warn "Missing required columns in $csv_path"
            continue
        end

        clean = dropmissing(df, [xcol, ycol])
        if nrow(clean) == 0
            @warn "No valid data for pf = $pf_val"
            continue
        end

        grouped = groupby(clean, xcol)
        summary = combine(grouped,
                          ycol => mean   => :mean_curv,
                          ycol => std    => :sd_curv)
        sort!(summary, xcol)

        color  = colors[mod1(idx, length(colors))]
        marker = markers[mod1(idx, length(markers))]

        plot!(p, summary.aspect, summary.mean_curv,
              yerror = summary.sd_curv,
              label  = "Pf = $pf_val",
              legend = :topleft,
              color  = color,
              marker = marker,
              markersize = 7,
              linewidth = 2,
              markerstrokewidth = 1)
    end

    plot_path = joinpath(parent_folder, "nocurvature_vs_aspect_multi_pf.png")
    savefig(p, plot_path)
    println("Saved plot → $plot_path")

    return p
end

# ======================  RUN THIS  ======================
parent = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\\"

pf_files = [
    (0.05, raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.05\final_all_no_curvature_output_R1.5_v5.0_pf0.05.csv"),
    (0.1, raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.1\final_all_no_curvature_output_R1.5_v5.0_pf0.1.csv"),
    (0.2, raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.2\final_all_no_curvature_output_R1.5_v5.0_pf0.2.csv"),
    (0.3, raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.3\final_all_no_curvature_output_R1.5_v5.0_pf0.3.csv"),
]

plot_curvature_vs_aspect_multi_pf(parent, pf_files)