#function to create plots, histogram and box plot for a data file 

using Plots, CSV,DataFrames, StatsPlots

path = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2025\04.April\V5.0_peakfreq\\"

filename= "v5_perimeter_95runs.png"

f= joinpath(path, filename)
f1= joinpath(path, "boxplot_"*filename)
fv= path*"peak_frequencies_combined_data45_90runs.csv"

df= CSV.read(fv,DataFrame)
#x= ["30s" "60s"  "120s"]
x= df[:,:Peak_Frequency_mHz]
#y= [45.24 52.77 41.54]
#  y= df[:,:FFTcurvature]
# x= df[:,:Run]
 p1=histogram(x,xlabel="Frequency(mHz)",ylabel="Count",title="V5",ylims= (0,40), xtickfont=font(12), ytickfont=font(12),color= :lightblue, legend= false)

 println("Mean of frequencies: ", mean(x))
    println("Standard deviation of frequencies: ", std(x))

# Create a boxplot
p = @df df boxplot(:Peak_Frequency_mHz, legend=false,ylabel="Frequency(mHz)", ylims=(0, 5),title="V5", bar_width=0.02)

#plot!(p,  kind="box")
savefig(p1, f)
savefig(p, f1)

