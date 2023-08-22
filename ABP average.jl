# PURPOSE: average the data from multiple folders/simulation runs and make a single averaged plot
# MethodL Data will be read from multiple folders and averaged data will be plotted and stored
# Why ? I am using this code for the fitting the average abd cheching for the transients
using CSV, FileIO, DataFrames

main_directory="C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\07.July\\parameter scan\\20230731-094339\\R=2.0 v=5.0 a=50.0 b=25.0 pf=0.1\\" # destination folder path


function aggregate_csv_data(main_directory::String)
    # This function will return the aggregated DataFrame

    # Initialize an empty DataFrame
    aggregated_data = DataFrame()

    # List all sub-folders inside the main directory
    subfolders = filter(isdir, readdir(main_directory, join=true))

    for folder in subfolders
        # List all .csv data files inside the sub-folder
        data_files = filter(filename -> occursin("_p.csv", filename) && isfile(filename), readdir(folder, join=true))

        for file in data_files
            # Read the data
            data = CSV.File(file) |> DataFrame #using chaining

            # Check if aggregated_data is empty
            if isempty(aggregated_data)
                aggregated_data = copy(data)
            else
                # Append data to the aggregated data
                aggregated_data = vcat(aggregated_data, data)
            end
        end
    end
    return aggregated_data
end

# Use the function to get the aggregated data
main_directory="C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\07.July\\parameter scan\\20230731-094339\\R=2.0 v=5.0 a=50.0 b=25.0 pf=0.1\\"
aggregated_data = aggregate_csv_data(main_directory)

averages = Dict()
for col in names(aggregated_data)
    if eltype(aggregated_data[!, col]) <: Number
        averages[col] = mean(skipmissing(aggregated_data[!, col]))
    end
end

mean(aggregated_data[19995:20000,:p1])