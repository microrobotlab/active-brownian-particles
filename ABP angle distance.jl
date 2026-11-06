using CSV, DataFrames, Plots, Statistics

"""
export_angle_distance(input_csv::String, output_csv::String, a, b;
                      plot_range::Union{Nothing,Tuple{Int,Int}}=nothing)

plot_range = (start_step, end_step)   e.g. (200000, 400000)
for one frame plotting put start and end as same value e.g. (200000, 200000)
All particles with StepN in [start_step, end_step] are used for the histograms.
"""
function export_angle_distance(input_csv::String, output_csv::String, a, b;
                               plot_range::Union{Nothing,Tuple{Int,Int}}=nothing)

    println("Reading $input_csv")
    df = CSV.read(input_csv, DataFrame)

    for col in (:Pn, :StepN, :xpos, :ypos)
        col in propertynames(df) || error("Column $col not found in $input_csv")
    end

    df.angle    = atan.(df.ypos, df.xpos)
    df.distance = sqrt.(df.xpos.^2 .+ df.ypos.^2)
    df.r_norm   = sqrt.((df.xpos ./ a).^2 .+ (df.ypos ./ b).^2)

    out = df[:, [:Pn, :StepN, :angle, :distance, :r_norm]]
    CSV.write(output_csv, out)
    println("Wrote $(nrow(out)) rows → $output_csv")

    # ---------- plots over a range of StepN ----------
    if plot_range !== nothing
        start_step, end_step = plot_range
        sub = filter(row -> start_step <= row.StepN <= end_step, out)

        if nrow(sub) == 0
            @warn "No particles found in StepN range $start_step:$end_step"
        else
            println("Plotting $(nrow(sub)) particles from StepN $start_step to $end_step")
            base = replace(output_csv, ".csv" => "")
            range_str = "$(start_step)_to_$(end_step)"

            # 1. Angle histogram
            p1 = histogram(sub.angle, bins=60,
                           xlabel = "angle (rad)",
                           ylabel = "count",
                           title  = "Angle distribution\nStepN $start_step – $end_step",
                           legend = false,
                           xguidefontsize = 16,
                           tickfontsize = 13,
                           xticks = ([-π, -π/2, 0, π/2, π],
                                     [L"-\pi", L"-\pi/2", "0", L"\pi/2", L"\pi"]))
            savefig(p1, base * "_angle_$(range_str).png")
            println("Saved angle plot")

            # 2. Distance histogram
            p2 = histogram(sub.distance, bins=60,
                           xlabel = "radial distance",
                           ylabel = "count",
                           title  = "Distance distribution\nStepN $start_step – $end_step",
                           legend = false,
                           xguidefontsize = 16,
                           tickfontsize = 13)
            savefig(p2, base * "_distance_$(range_str).png")
            println("Saved distance plot")

            # 3. r_norm histogram
            p3 = histogram(sub.r_norm, bins=60,
                           xlabel = "r_norm",
                           ylabel = "count",
                           title  = "r_norm distribution\nStepN $start_step – $end_step",
                           legend = false,
                           xguidefontsize = 16,
                           tickfontsize = 13)
            savefig(p3, base * "_rnorm_$(range_str).png")
            println("Saved r_norm plot")

            display(p1); display(p2); display(p3)
        end
    end

    return out
end

parent_folder = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2025\03.March\20250325-165047\R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1\run1\\"
filename = "20250325-165047 R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1 run1"

input  = joinpath(parent_folder, filename * ".csv")
output = joinpath(parent_folder, filename * "_angle_distance.csv")

# Plot all particles between StepN 200000 and 400000
export_angle_distance(input, output, 50.0, 12.5; plot_range=(100000, 500000))



