using CSV, DataFrames, Statistics

"""
multianalysis_prob_no_curvature(folder_names::String, parent_folder::String)

For each run, reads the _parameter_p.csv and computes:
  P_no_curv = (number of steps where NeqL+NeqR-NpoleU-NpoleD == 0) / maximum(t)
"""
function multianalysis_prob_no_curvature(folder_names::String, parent_folder::String)

    # ---------- read list of main folders ----------
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

    results = NamedTuple[]

    for mainfolder in mainfolders
        println("Processing mainfolder: $mainfolder")
        mainfolder_name = basename(mainfolder)

        # find the parameter sub-folder
        subfolders = [joinpath(mainfolder, d) for d in readdir(mainfolder) if isdir(joinpath(mainfolder, d))]
        param_folders = filter(f -> occursin("R=", basename(f)) && occursin("pf=", basename(f)), subfolders)

        if isempty(param_folders)
            @warn "No parameter sub-folder found in $mainfolder, skipping."
            continue
        end
        sub_subfolder = param_folders[1]
        params = parse_parameters(basename(sub_subfolder))

        # run folders
        run_folders = [joinpath(sub_subfolder, d) for d in readdir(sub_subfolder) if isdir(joinpath(sub_subfolder, d))]
        ICS = length(run_folders)
        if ICS == 0
            @warn "No run folders in $sub_subfolder, skipping."
            continue
        end

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

            # required columns
            for col in (:t, :NeqL, :NeqR, :NpoleU, :NpoleD)
                col in propertynames(df) || error("Column $col not found in $param_csv")
            end

            # delta N
            delta = df.NeqL .+ df.NeqR .- df.NpoleU .- df.NpoleD

            # count how many times delta == 0
            n_zero = count(==(0), delta)

            # maximum t in the file
            t_max = maximum(df.t)

            # probability
            prob_no_curv = n_zero / t_max

            push!(results, (
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
        end
    end

    if isempty(results)
        @warn "No results collected."
        return nothing
    end

    combined_df = DataFrame(results)
    output_file = joinpath(parent_folder, "output_prob_no_curvature.csv")
    CSV.write(output_file, combined_df)
    println("Combined results written to $output_file")

    return combined_df
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

parent_folder = "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07Coding\\workstationMRL\\2026\\10.October\\final paper data\\R1.5 pf 0.05\\aspect 0.55\\" # destination folder path
folder_names= parent_folder*"folder_names.txt" # file containing the names of the folders to be analyzed
multianalysis_prob_no_curvature(folder_names, parent_folder)