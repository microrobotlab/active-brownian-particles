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
        error("Failed to read $txt_file: $e")
    end

    if isempty(mainfolders)
        error("No valid folders found in $txt_file")
    end

    # Initialize storage for all results
    all_peak_frequencies = Dict{String, Vector{Float64}}()
# Initialize storage for results in long format
results = Tuple{Int,String,Float64}[]  # (Run, Main_Folder, Peak_Frequency_mHz)
    # Process each main folder
    for mainfolder in mainfolders
        println("Processing mainfolder: $mainfolder")
        if !isdir(mainfolder)
            @warn "Mainfolder $mainfolder does not exist, skipping."
            continue
        end
# Define mainfolder_name at the start of the loop
mainfolder_name = basename(mainfolder)          # e.g., 20250415-561314
        # Access the single sub-subfolder
        sub_subfolder = joinpath(mainfolder, "R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1")
        sub_subfolder_name = basename(sub_subfolder)  # R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1
        println("Processing sub-subfolder: $sub_subfolder_name")
        if !isdir(sub_subfolder)
            @warn "Sub-subfolder $sub_subfolder does not exist, skipping."
            continue
        end
                # Dynamically determine ICS from the number of run folders
                run_folders = [joinpath(sub_subfolder, d) for d in readdir(sub_subfolder) if isdir(joinpath(sub_subfolder, d))]
                ICS = length(run_folders)
                if ICS == 0
                    @warn "No run folders found in $sub_subfolder, skipping."
                    continue
                end
                println("Found $ICS run folders in $sub_subfolder")

        # Initialize storage for this mainfolder
        peak_frequencies = Float64[]
        psds = []

        # Process runs (run1, run2, ..., runICS)
        for i in 1:ICS
            println("Analysis for run $i in $sub_subfolder")
            run_folder = joinpath(sub_subfolder, "run$i")
            if !isdir(run_folder)
                @warn "Run folder $run_folder does not exist, skipping."
                continue
            end

            # Construct filename using main folder name and constant parameters
            mainfolder_name = basename(mainfolder)  # e.g., 20250415-561314
            filename = "$mainfolder_name R=$R v=5.0 a=$a b=$b pf=0.1 run$i"
            t = joinpath(run_folder, filename)

            # Perform statistical analysis
            inside_Np = stat_analysis_perimeter(a, b, R, t, δt, 2)

            # FFT analysis
            peak_freq_mHz, peak_period, freqs_mHz, psd = FFT_analysis_window(t, run_folder, δt, resample, 3)

            # Store results
            push!(results, (i, mainfolder_name, peak_freq_mHz))
            push!(peak_frequencies, peak_freq_mHz)
            push!(psds, psd)
        end

        # Store peak frequencies for this mainfolder
        all_peak_frequencies[mainfolder_name] = peak_frequencies

        # Save plot for this mainfolder
        plot(1:length(peak_frequencies), peak_frequencies, 
             title="Peak Frequencies for $mainfolder_name", 
             xlabel="Run", ylabel="Frequency (mHz)", legend=false)
        savefig(joinpath(mainfolder, "peak_frequencies_$(mainfolder_name).png"))
    end

    # Create combined CSV for all mainfolders
    combined_df = DataFrame(
      Run = [r[1] for r in results],
      Main_Folder = [r[2] for r in results],
      Peak_Frequency_mHz = [r[3] for r in results]
  )
    output_file = joinpath(parent_folder, "peak_frequencies_combined.csv")
    CSV.write(output_file, combined_df)
end




# function multianalysis(mainfolder,ICS)
#   # Initialize storage for peak frequencies and PSDs
#   peak_frequencies = Float64[]
#   psds = []
#     for i in 1:ICS # number of runs
#       println("analysis for run $i")
#     mainfolder1= mainfolder*"run$i\\"
#     filename="20250430-174406 R=1.5 v=10.0 a=50.0 b=12.5 pf=0.1 run$i" # dont put \\ after the filename
#     t= mainfolder1 * filename 
   
#   #f1000= joinpath(mainfolder1,filename * ".csv") 
 
#  # #    f1= "D:\\j.sharma\\P07\\workstationMRL\\20241104-121620\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.2\\run1\\20241104-121620 R=2.0 v=10.0 a=50.0 b=25.0 pf=0.2 run1_p.csv\\"
# # df= CSV.read(f1000, DataFrame) 
#   inside_Np=stat_analysis_perimeter(a,b,R,t,δt,2) # 0 for pole, equator, 1 for only right left, 2 for entire
#  #FFT_analysis(t,δt,resample,3) # FFT analysis of the data, number at the end is for different type of analyis

#  ########################## FFT analysis of the data using Welch's method ##########################
#  peak_freq_mHz, peak_period, freqs_mHz, psd = FFT_analysis_window(t,δt,resample,3)

#  push!(peak_frequencies, peak_freq_mHz)
#  push!(psds, psd)

# end
#    # Save peak frequencies to a CSV file
#    output_file = joinpath(mainfolder, "peak_frequencies.csv")
#    peak_df = DataFrame(Run = 1:ICS, Peak_Frequency_mHz = peak_frequencies)
#    CSV.write(output_file, peak_df)
#    plot(1:ICS, peak_frequencies, title="Peak Frequencies for $ICS runs", xlabel="Run", ylabel="Frequency (mHz)", legend=false)
#    savefig(mainfolder*"peak_frequencies.png")
  
    
# end