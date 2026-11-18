using CSV, DataFrames, Statistics, LaTeXStrings, CairoMakie

"""
plot_interval_and_curvature(csv_path::String, output_png::String;
                            curvature_x = ["0.1", "0.25", "0.4", "0.55"],
                            curvature_y = [1024, 64, 16, 6])

Reads the zero-intervals CSV, computes the median interval for each aspect,
and plots it together with the curvature values on a dual y-axis.
"""
function plot_interval_and_curvature(csv_path::String, output_png::String;
                                     curvature_x = ["0.1", "0.25", "0.4", "0.55"],
                                     curvature_y = [1024, 64, 16, 6])

    println("Reading $csv_path")
    df = CSV.read(csv_path, DataFrame)

    # required columns
    for col in (:aspect, :interval)
        col in propertynames(df) || error("Column $col not found")
    end

    # optional filter (remove very small intervals if desired)
    df = filter(row -> row.interval > 20, df)

    # median interval per aspect
    
    summary = combine(groupby(df, :aspect), :interval => mean => :median_interval)
    sort!(summary, :aspect)

    println("Median intervals:")
    println(summary)

    # convert aspect to string labels for matching
    aspect_labels = string.(round.(summary.aspect, digits=2))

    #=

    # ---------- left axis: median interval ----------
    p1 = plot(aspect_labels, summary.median_interval,
              xlabel = "aspect ratio (b/a)",
              ylabel = L"Median time period to cross $\kappa=0$",
              label = "Median interval",
              linewidth = 2,
              marker = :circle,
              markersize = 7,
              color = :blue,
              legend = :topright,
              grid = true)

    # ---------- right axis: curvature (overlaid) ----------
    p2 = plot(curvature_x, curvature_y,
              ylabel = L"$K$",
              label = "Curvature",
              linewidth = 2,
              marker = :diamond,
              markersize = 7,
              color = :red,
              legend = :topleft,
              yaxis = :right,          # puts this y-axis on the right
              grid = false)

    # combine them into one plot
    p = plot(p1, p2, layout = (1, 2), size = (900, 400))
    =#

    #savefig(p, output_png)
    # convert string labels to numbers for plotting (or keep them as categorical)
    x_vals = parse.(Float64, curvature_x)

    # ---------- Makie dual-axis plot ----------
    f = CairoMakie.Figure()

    ax1 = Axis(f[1, 1],
               xlabel = "aspect ratio (b/a)",
               ylabel = L"Mean time period to cross $\kappa=0$",
               #ylimits= (0,7000),
               
               xticklabelcolor = :black,
               yticklabelcolor = :blue)

    ax2 = Axis(f[1, 1],
               ylabel = L"$K$",
               xticklabelcolor = :black,
               yticklabelcolor = :red,
               #xaxisposition = :top,
               
               yaxisposition = :right,
               yticklabelsvisible = true)

    CairoMakie.ylims!(ax1, 0, 1700)
    CairoMakie.ylims!(ax2, 0, 1700)

    # left axis (blue)
    scatterlines!(ax1, x_vals, summary.median_interval, color = :blue, marker = :circle,
               markersize = 7, linewidth = 2, label = "Time interval")


    # right axis (red)
    scatterlines!(ax2, x_vals, curvature_y, color = :red, marker = :circle,
               markersize = 7,linewidth = 2, label = "Curvature")

    # optional legends
    axislegend(ax1, position = :rt)
    axislegend(ax2, position = :rc)

    display(f)
    save(output_png, f)

    println("Saved plot → $output_png")

    return summary, p
end

csv_path   = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.05\final_all_zero_intervals_R1.5_v5.0_pf0.05_1000-10000.csv"
output_png = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.05\interval_and_curvature_vs_aspect.png"

plot_interval_and_curvature(csv_path, output_png;
                            curvature_x = ["0.1", "0.25", "0.4", "0.55"],
                            curvature_y = [1024, 64, 16, 6])
