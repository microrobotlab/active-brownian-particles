# PURPOSE: average the data from multiple folders/simulation runs and make a single averaged plot
# MethodL Data will be read from multiple folders and averaged data will be plotted and stored
# Why ? I am using this code for the fitting the average abd cheching for the transients
using CSV, FileIO, DataFrames, Plots, LaTeXStrings


println(replace("aggregated_data" , "aggregated" => "all"))

function all_csv_data(mainfolder)
    # This function will return the aggregated DataFrame

    # Initialize an empty DataFrame
    all_data = DataFrame()
    t1= plot()
    f= "average_pf.png"

    # List all sub-folders inside the main directory
    subfolders = filter(isdir, readdir(mainfolder, join=true))

    for folder in subfolders
        # List all .csv data files inside the sub-folder, this command choose the specific data files with _p.csv extension.
        # if you have one .csv file per sunfolder this can be omitted 
        data_files = filter(filename -> occursin("_p.csv", filename) && isfile(filename), readdir(folder, join=true))

        for file in data_files
            # Read the data
            data = CSV.File(file) |> DataFrame #using chaining

          
            if isempty(all_data)
                all_data = copy(data)
            else
                # Append data to the aggregated data
                all_data = vcat(all_data, data)
            end
        end
    end

    #mean(all_data[!,:p1]), mean(all_data[!,:p2]) will give average off all packing fraction at equators and poles over time
    # but I would need the average of all runs over each time step or second, a different method is required, following steps

    #all_data[!,:t] = categorical(all_data[!,:t],compress=true) 
    gdf = groupby(all_data,:t,sort=true)
    
   
    time= [unique(g[!,:t]) for g in gdf]  # unique avoids repeation of time in two files and keeps the length of vector same as avg_eq and avg_p 
    avg_eq= [mean(g[!,:p1]) for g in gdf]

    for i in 1:length(gdf)

        #avg_eq= [mean(gdf[i][!,:p1])]

        avg_p=  [mean(gdf[i][!,:p2])]
    end
    scatter!(t1,[time],[avg_eq], ylimit=(0,35),legend=false) 
    avg_p=  [mean(g[!,:p2]) for g in gdf]
   
    scatter!(t1,time,avg_eq, ylimit=(0,0.3),legend=false) 
    xlabel!("Time (s)", xguidefont=font(16), xtickfont=font(11))
    plot!(ylabel=L"\mathrm{pf^{-}_{eqs}}",yguidefont=font(16), ytickfont=font(11))
    title!(" Equators ")
    plot(t1)
    savefig(f)
    #=
    scatter!(t2,[i],[pfp], ylimit=(0,0.3),legend=false) 
    xlabel!("Time (s)", xguidefont=font(16), xtickfont=font(11))
    plot!(ylabel=L"\mathrm{pf_{poles}}",yguidefont=font(16), ytickfont=font(11))
    title!(" Poles ")
    =#
    return nothing #avg_eq, time # mean(all_data[!,:p1]), mean(all_data[!,:p2]) all_data, gdf,
end

# Use the function to get the aggregated data
mainfolder="C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230823-102538\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\"
all_data = all_csv_data(mainfolder)
 
#averages = Dict()
#=for col in names(aggregated_data)
    if eltype(aggregated_data[!, col]) <: Number
        averages[col] = mean(skipmissing(aggregated_data[!, col]))
    end
end
=#
#mean(aggregated_data[19995:20000,:p1])