# PURPOSE: average the data from multiple folders/simulation runs and make a single averaged plot
# MethodL Data will be read from multiple folders and averaged data will be plotted and stored
# Why ? I am using this code for the fitting the average and checking for the transients
using CSV, FileIO, DataFrames, Plots, LaTeXStrings, Statistics 

gr()

function average(mainfolder)
  
    all_data = DataFrame()
    t1= plot()
    t2= plot()
    f= mainfolder*"V 10.0.png" # file of average of FFT of all individual FFTS
    # f1= mainfolder*"V 5.0_data.csv"

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
                # Append data
                all_data = vcat(all_data, data)
            end
        end
    end

    #mean(all_data[!,:p1]), mean(all_data[!,:p2]) will give average off all packing fraction at equators and poles over time
    # but I would need the average of all runs over each time step or second, a different method is required, following steps

   
    gdf = groupby(all_data,:t,sort=true)
    avg_eq= [mean(g[!,:p1]) for g in gdf]  # frequency
    avg_p=  [mean(g[!,:p2]) for g in gdf]   #power
    effect= avg_eq - avg_p
    time= all_data[1:10000,:t] # total time steps are now 10000, although in actual simulations they are 100 times more, because now i skipped 100 time steps in saving the data. keep this in mind
    # Data= DataFrame(t= time, p1= avg_eq, p2= avg_p)
    # CSV.write(f1,Data)
 

   t1= scatter(time./100.0, effect,legend=false)  
    xlabel!("Time (s)", xguidefont=font(16),xlimit=(0,10000),ylimit=(0,10), xtickfont=font(11))
    plot!(ylabel=L"\mathrm{N_{eqs}-{N_{poles}}}",yguidefont=font(16), ytickfont=font(11))
    title!(" V 10.0 100 avg")
 
#    t2= scatter(time./100.0, avg_p,legend=false) 
#     xlabel!("Time (s)", xguidefont=font(16),ylimit=(15,30), xtickfont=font(11))
#     plot!(ylabel=L"\mathrm{pf_{poles}}",xlimit=(0,10000),yguidefont=font(16), ytickfont=font(11))
#     title!(" Poles ")
    
    p= plot(t1)
    savefig(p,f)
    
    #k= plot(avg_f,avg_real, xlimit=(-0.5,0.5), ylimit=(0.02,200),seriestype=:stem, xlabel="Frequency(Hz)", ylabel="Power",legend=false)
    #savefig(k,f)
    #display(k)
    return nothing 
end

mainfolder= "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-205011\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\"
#all_data = average(mainfolder)

(average(mainfolder)) 

   #println(first(all_data,5))
    #println("Time:", time[1:5])
    #println("Avg Equators:", avg_eq[1:5])
    #println("Avg Poles:", avg_p[1:5])
