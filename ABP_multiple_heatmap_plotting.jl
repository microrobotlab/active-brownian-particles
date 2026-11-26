using CSV, DataFrames, Plots, Statistics

function plot_curvature_scatter(parent_folder::String, csv_paths::Vector{String};
                                pf_col::Symbol = :pf_critical,
                                aspect_col::Symbol = :aspect,
                                curv_col::Symbol = :Curvature_Effect)

    all_data = DataFrame(eccentricity=Float64[], pf_critical=Float64[], mean_curv=Float64[])

    for csv_path in csv_paths
        if !isfile(csv_path)
            @warn "File not found, skipping: $csv_path"
            continue
        end

        println("Reading $csv_path")
        df = CSV.read(csv_path, DataFrame)

        required = (aspect_col, pf_col, curv_col)
        if !all(c -> c in propertynames(df), required)
            @warn "Missing required columns $(required) in $csv_path"
            continue
        end

        clean = dropmissing(df, collect(required))
        if nrow(clean) == 0
            @warn "No valid data in $csv_path"
            continue
        end

        # Compute eccentricity from aspect = b/a
        clean.eccentricity = sqrt.(1.0 .- clean[!, aspect_col].^2)

        # Group by eccentricity + pf_critical and take mean curvature
        grouped = groupby(clean, [:eccentricity, pf_col])
        summary = combine(grouped, curv_col => mean => :mean_curv)

        for row in eachrow(summary)
            push!(all_data, (
                eccentricity = row.eccentricity,
                pf_critical  = row[pf_col],
                mean_curv    = row.mean_curv
            ))
        end
    end

    if nrow(all_data) == 0
        @warn "No data collected — nothing to plot"
        return nothing
    end

    p = Plots.scatter(all_data.pf_critical, all_data.mean_curv;
                marker_z = all_data.eccentricity,
                marker   = :circle,
                markersize = 8,
                markerstrokewidth = 0.5,
                xlabel = L"$Pf_{critical}$",
                ylabel = L"$\kappa$",
                colorbar = true,
                colorbar_title = "Eccentricity",
                c = :viridis,
                legend = false,
                size = (800, 500),
                grid = true,
                clims = (minimum(all_data.eccentricity), maximum(all_data.eccentricity)))

    plot_path = joinpath(parent_folder, "curvature_vs_pf_colored_by_eccentricity.png")
    Plots.savefig(p, plot_path)
    println("Saved plot → $plot_path")

    return p
end

# ======================  RUN THIS  ======================
parent = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\\"

csv_paths = [
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.05\final_all_output_R1.5_v5.0_pf0.05.csv",
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.1\final_all_output_R1.5_v5.0_pf0.1.csv",
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.2\final_all_output_R1.5_v5.0_pf0.2.csv",
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.3\final_all_output_R1.5_v5.0_pf0.3.csv",
    raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.4\final_all_output_R1.5_v5.0_pf0.4.csv",
]

plot_curvature_scatter(parent, csv_paths)