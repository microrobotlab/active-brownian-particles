# Here the analysis for multiple files is done, the data is read from multiple folders 
# Analysis is done on each file and 
# One can alsso plot average of all data from multiple files
# INPUT: mainfolder path of the main folder which has all the run folders
# OUTPUT: Analysis of each file 
function multianalysis(folder_names::String, parent_folder::String)
    # Read main folder names from .txt file
    mainfolders = String[]
    try
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
    catch e
        error("Failed to read $folder_names: $e")
    end

    if isempty(mainfolders)
        error("No valid folders found in $folder_names")
    end

    # Collect results as NamedTuples
    results = NamedTuple[]

    # Process each main folder
    for mainfolder in mainfolders
        println("Processing mainfolder: $mainfolder")
        if !isdir(mainfolder)
            @warn "Mainfolder $mainfolder does not exist, skipping."
            continue
        end

        mainfolder_name = basename(mainfolder)

        # ---------- Find the parameter sub-folder dynamically ----------
        subfolders = [joinpath(mainfolder, d) for d in readdir(mainfolder) if isdir(joinpath(mainfolder, d))]
        # Keep only folders whose name looks like "R=... v=... a=... b=... pf=..."
        param_folders = filter(f -> occursin("R=", basename(f)) && occursin("pf=", basename(f)), subfolders)

        if isempty(param_folders)
            @warn "No parameter sub-folder found in $mainfolder, skipping."
            continue
        end

        # Use the first (and normally only) matching folder
        sub_subfolder = param_folders[1]
        println("Using sub-subfolder: $(basename(sub_subfolder))")

        # Now parse parameters from the real folder name
        params = parse_parameters(basename(sub_subfolder))

        # Dynamically determine number of run folders
        run_folders = [joinpath(sub_subfolder, d) for d in readdir(sub_subfolder) if isdir(joinpath(sub_subfolder, d))]
        ICS = length(run_folders)
        if ICS == 0
            @warn "No run folders found in $sub_subfolder, skipping."
            continue
        end
        println("Found $ICS run folders in $sub_subfolder")

        # Storage for this mainfolder (for plots)
        peak_frequencies = Float64[]
        curvature_vals   = Float64[]

        # Process each run
        for i in 1:ICS
            println("Analysis for run $i in $sub_subfolder")
            run_folder = joinpath(sub_subfolder, "run$i")
            if !isdir(run_folder)
                @warn "Run folder $run_folder does not exist, skipping."
                continue
            end

            filename = "$(mainfolder_name) R=$(params.R) v=$(params.v) a=$(params.a) b=$(params.b) pf=$(params.pf) run$i"
            t = joinpath(run_folder, filename)

            # Statistical / curvature analysis
            curvature = stat_analysis_perimeter(params.a, params.b, params.R, t, δt, 2, 1000, 10000)

            # FFT analysis
            peak_freq_mHz, peak_period, freqs_mHz, psd = FFT_analysis_window(t, run_folder, δt, resample_data, 3)

            # Store one complete result row
            push!(results, (
                Run                = i,
                Main_Folder        = mainfolder_name,
                R                  = params.R,
                v                  = params.v,
                a                  = params.a,
                b                  = params.b,
                aspect             = params.b / params.a,
                pf                 = params.pf,
                Peak_Frequency_mHz = peak_freq_mHz,
                Curvature_Effect   = curvature
            ))

            push!(peak_frequencies, peak_freq_mHz)
            push!(curvature_vals, curvature)
        end

        # Plots for this mainfolder
        if !isempty(peak_frequencies)
            plot(1:length(peak_frequencies), peak_frequencies,
                 title="Peak Frequencies for $mainfolder_name",
                 xlabel="Run", ylabel="Frequency (mHz)", legend=false)
            savefig(joinpath(mainfolder, "peak_frequencies_$(mainfolder_name).png"))

            plot(1:length(curvature_vals), curvature_vals,
                 title="Curvature Analysis Results for $mainfolder_name",
                 xlabel="Run", ylabel="Curvature effect", legend=false)
            savefig(joinpath(mainfolder, "curvature_analysis_$(mainfolder_name).png"))
        end
    end

    # One big combined table
    if isempty(results)
        @warn "No results collected."
        return nothing
    end

    combined_df = DataFrame(results)
    output_file = joinpath(parent_folder, "output_combined.csv")
    CSV.write(output_file, combined_df)
    println("Combined results written to $output_file")

    return combined_df
end

# Helper: parse "R=1.5 v=5.0 a=50.0 b=12.5 pf=0.4" into a NamedTuple
function parse_parameters(folder_name::String)
    d = Dict{Symbol, Float64}()
    for part in split(folder_name)
        if occursin('=', part)
            k, v = split(part, '=')
            try
                d[Symbol(k)] = parse(Float64, v)
            catch
                # ignore unparseable parts
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
