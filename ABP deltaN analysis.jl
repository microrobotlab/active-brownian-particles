using CSV, DataFrames, Statistics, Plots

function multianalysis_prob_no_curvature(folder_names::String, parent_folder::String)

    mainfolders = String[]
    open(folder_names) do file
        for line in eachline(file)
            folder = strip(line)
            if !isempty(folder)
                full_path = joinpath(parent_folder, folder)
                if isdir(full_path)
                    push!(mainfolders, full_path)
                else
                    @warn "Folder $full_path does not exist, skipping."
                end
            end
        end
    end

    if isempty(mainfolders)
        error("No valid folders found in $folder_names")
    end

    summary_results = NamedTuple[]
    interval_results = NamedTuple[]

    for mainfolder in mainfolders
        println("Processing mainfolder: $mainfolder")
        mainfolder_name = basename(mainfolder)

        subfolders = [joinpath(mainfolder, d) for d in readdir(mainfolder) if isdir(joinpath(mainfolder, d))]
        param_folders = filter(f -> occursin("R=", basename(f)) && occursin("pf=", basename(f)), subfolders)

        if isempty(param_folders)
            @warn "No parameter sub-folder found in $mainfolder, skipping."
            continue
        end
        sub_subfolder = param_folders[1]
        params = parse_parameters(basename(sub_subfolder))

        run_folders = [joinpath(sub_subfolder, d) for d in readdir(sub_subfolder) if isdir(joinpath(sub_subfolder, d))]
        ICS = length(run_folders)
        if ICS == 0
            @warn "No run folders in $sub_subfolder, skipping."
            continue
        end

        # collect intervals for this condition only
        condition_intervals = Float64[]

        for i in 1:ICS
            run_folder = joinpath(sub_subfolder, "run$i")
            if !isdir(run_folder)
                @warn "Run folder $run_folder missing, skipping."
                continue
            end

            filename = "$(mainfolder_name) R=$(params.R) v=$(params.v) a=$(params.a) b=$(params.b) pf=$(params.pf) run$i"
            param_csv = joinpath(run_folder, filename * "_parameter_p.csv")

            if !isfile(param_csv)
                @warn "Missing $param_csv, skipping run $i"
                continue
            end

            df = CSV.read(param_csv, DataFrame)

            for col in (:t, :NeqL, :NeqR, :NpoleU, :NpoleD)
                col in propertynames(df) || error("Column $col not found in $param_csv")
            end

            delta = df.NeqL .+ df.NeqR .- df.NpoleU .- df.NpoleD
            t_zeros = sort(df.t[delta .== 0])

            n_zero = length(t_zeros)
            t_max = maximum(df.t)
            prob_no_curv = n_zero / t_max

            push!(summary_results, (
                Run           = i,
                Main_Folder   = mainfolder_name,
                R             = params.R,
                v             = params.v,
                a             = params.a,
                b             = params.b,
                aspect        = params.b / params.a,
                pf            = params.pf,
                n_zero        = n_zero,
                t_max         = t_max,
                prob_no_curv  = prob_no_curv
            ))

            if n_zero >= 2
                intervals = diff(t_zeros)
                append!(condition_intervals, intervals)

                for Δt in intervals
                    push!(interval_results, (
                        Run           = i,
                        Main_Folder   = mainfolder_name,
                        R             = params.R,
                        v             = params.v,
                        a             = params.a,
                        b             = params.b,
                        aspect        = params.b / params.a,
                        pf            = params.pf,
                        interval      = Δt
                    ))
                end
            end
        end
     # ---------- two histograms for this condition ----------
        if !isempty(condition_intervals)
            intervals_1s   = filter(x -> abs(x - 1.0) < 1e-6, condition_intervals)
            intervals_other = filter(x -> x > 1.0 + 1e-6, condition_intervals)

            # 1 s histogram
            if !isempty(intervals_1s)
                p1 = histogram(intervals_1s,
                               bins=10,
                               xlabel="Interval (s)",
                               ylabel="Count",
                               title="Intervals = 1 s\n$(mainfolder_name)",
                               legend=false)
                savefig(p1, joinpath(mainfolder, "intervals_1s.png"))
                println("Saved 1 s histogram")
            end

            # other intervals histogram
            if !isempty(intervals_other)
                p2 = histogram(intervals_other,
                               bins=40,
                               xlabel="Interval (s)",
                               ylabel="Count",
                               title="Intervals > 1 s\n$(mainfolder_name)",
                               legend=false)
                savefig(p2, joinpath(mainfolder, "intervals_other.png"))
                println("Saved other intervals histogram")
            end
        end
    end

    # ---------- write the two CSVs ----------
    summary_df = DataFrame(summary_results)
    intervals_df = DataFrame(interval_results)

    if !isempty(summary_df)
        CSV.write(joinpath(parent_folder, "output_prob_no_curvature.csv"), summary_df)
        println("Summary CSV written")
    end

    if !isempty(intervals_df)
        CSV.write(joinpath(parent_folder, "output_zero_intervals.csv"), intervals_df)
        println("Intervals CSV written")
    end

    return summary_df, intervals_df
end

# ---------- helper ----------
function parse_parameters(folder_name::String)
    d = Dict{Symbol, Float64}()
    for part in split(folder_name)
        if occursin('=', part)
            k, v = split(part, '=')
            try
                d[Symbol(k)] = parse(Float64, v)
            catch
            end
        end
    end
    return (
        R  = get(d, :R,  NaN),
        v  = get(d, :v,  NaN),
        a  = get(d, :a,  NaN),
        b  = get(d, :b,  NaN),
        pf = get(d, :pf, NaN)
    )
end

parent_folder = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.05\aspect 0.1\\"
folder_names  = joinpath(parent_folder, "folder_names.txt")

summary_df, intervals_df = multianalysis_prob_no_curvature(folder_names, parent_folder)