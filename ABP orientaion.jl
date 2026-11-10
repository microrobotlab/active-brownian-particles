using CSV, DataFrames, Plots, Statistics
22

"""
zone_velocity_orientation(input_csv, output_csv, a, b, R, eθ, δt)

Reads xpos, ypos and orient from the CSV.
Computes speed from successive positions.
Averages speed & orientation per zone.
Saves full zone plots + separate EqL/EqR plots.
"""
function zone_velocity_orientation(input_csv::String, output_csv::String, a, b, R, eθ, δt)

    # ---- change if your orientation column has a different name ----
    orient_col = :orientation
    # -----------------------------------------------------------------

    println("Reading $input_csv")
    df = CSV.read(input_csv, DataFrame)

    for col in (:Pn, :StepN, :xpos, :ypos, orient_col)
        col in propertynames(df) || error("Missing column $col in the CSV")
    end

    # Compute speed from positions
    sort!(df, [:Pn, :StepN])
    df.speed = fill(NaN, nrow(df))

    for sub in groupby(df, :Pn)
        if nrow(sub) < 2
            continue
        end
        dx = diff(sub.xpos)
        dy = diff(sub.ypos)
        spd = sqrt.(dx.^2 .+ dy.^2) ./ δt
        df.speed[parentindices(sub)[1][2:end]] = spd
    end

    # Zone averages
    gdf = groupby(df, :StepN, sort=true)
    steps = unique(df.StepN)
    records = NamedTuple[]

    for (i, step) in enumerate(steps)
        group = gdf[i]
        zones = Dict(
            :EqL   => (speeds=Float64[], orients=Float64[]),
            :EqR   => (speeds=Float64[], orients=Float64[]),
            :PoleU => (speeds=Float64[], orients=Float64[]),
            :PoleD => (speeds=Float64[], orients=Float64[])
        )

        for j in 1:nrow(group)
            x = group.xpos[j]
            y = group.ypos[j]
            r = sqrt((x/a)^2 + (y/b)^2)

            if 0 < abs(r - 1) < 2R && !isnan(group.speed[j])
                θ = atan(y, x)
                spd = group.speed[j]
                ori = group[!, orient_col][j]

                if (π - eθ) <= θ <= π || -π <= θ <= -(π - eθ)
                    push!(zones[:EqL].speeds, spd);  push!(zones[:EqL].orients, ori)
                elseif -eθ <= θ <= eθ
                    push!(zones[:EqR].speeds, spd);  push!(zones[:EqR].orients, ori)
                elseif eθ < θ < (π - eθ)
                    push!(zones[:PoleU].speeds, spd); push!(zones[:PoleU].orients, ori)
                elseif -(π - eθ) < θ < -eθ
                    push!(zones[:PoleD].speeds, spd); push!(zones[:PoleD].orients, ori)
                end
            end
        end

        avg(v) = isempty(v) ? NaN : mean(v)

        push!(records, (
            StepN            = step,
            avg_v_EqL        = avg(zones[:EqL].speeds),
            avg_orient_EqL   = avg(zones[:EqL].orients),
            avg_v_EqR        = avg(zones[:EqR].speeds),
            avg_orient_EqR   = avg(zones[:EqR].orients),
            avg_v_PoleU      = avg(zones[:PoleU].speeds),
            avg_orient_PoleU = avg(zones[:PoleU].orients),
            avg_v_PoleD      = avg(zones[:PoleD].speeds),
            avg_orient_PoleD = avg(zones[:PoleD].orients),
            n_EqL   = length(zones[:EqL].speeds),
            n_EqR   = length(zones[:EqR].speeds),
            n_PoleU = length(zones[:PoleU].speeds),
            n_PoleD = length(zones[:PoleD].speeds)
        ))
    end

    result = DataFrame(records)
    CSV.write(output_csv, result)
    println("Wrote zone averages → $output_csv")

    base = replace(output_csv, ".csv" => "")

    # ---------- Full zone plots ----------
    p_speed = plot(result.StepN, result.avg_v_EqL,   label="EqL",   xlabel="StepN", ylabel="mean speed")
    plot!(p_speed, result.StepN, result.avg_v_EqR,   label="EqR")
    plot!(p_speed, result.StepN, result.avg_v_PoleU, label="PoleU")
    plot!(p_speed, result.StepN, result.avg_v_PoleD, label="PoleD")
    savefig(p_speed, base * "_avg_speed.png")

    p_ori = plot(result.StepN, result.avg_orient_EqL,   label="EqL",   xlabel="StepN", ylabel="mean orientation (rad)")
    plot!(p_ori, result.StepN, result.avg_orient_EqR,   label="EqR")
    plot!(p_ori, result.StepN, result.avg_orient_PoleU, label="PoleU")
    plot!(p_ori, result.StepN, result.avg_orient_PoleD, label="PoleD")
    savefig(p_ori, base * "_avg_orientation.png")

    # ---------- Separate EqL & EqR plots ----------
    # Orientation EqL + EqR
    p_ori_eq = plot(result.StepN, result.avg_orient_EqL, label="EqL", xlabel="StepN", ylabel="mean orientation (rad)",
                    title="Equator orientation", linewidth=2)
    plot!(p_ori_eq, result.StepN, result.avg_orient_EqR, label="EqR", linewidth=2)
    plot!( df.orient_col, linestyle=:dash, color=:black, label="0 rad")
    savefig(p_ori_eq, base * "_avg_orientation_Eq.png")

    # Velocity EqL + EqR
    p_v_eq = plot(result.StepN, result.avg_v_EqL, label="EqL", xlabel="StepN", ylabel="mean speed",
                  title="Equator velocity", linewidth=2)
    plot!(p_v_eq, result.StepN, result.avg_v_EqR, label="EqR", linewidth=2)
    savefig(p_v_eq, base * "_avg_speed_Eq.png")

    println("Saved all plots")
    display(p_speed)
    display(p_ori)
    display(p_ori_eq)
    display(p_v_eq)

    return result
end

parent = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2025\03.March\20250325-165047\R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1\run1\\"
filename = "20250325-165047 R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1 run1"

input  = joinpath(parent, filename * ".csv")
output = joinpath(parent, filename * "_zone_avg_v_orient.csv")

a  = 50.0
b  = 12.5
R  = 1.5
δt = 0.001
eθ = perimeter_angle_find(a, b)   # or a numeric value

zone_velocity_orientation(input, output, a, b, R, eθ, δt)



