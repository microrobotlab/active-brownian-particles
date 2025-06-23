#function to create plots, histogram and box plot for a data file 

using Plots, CSV,DataFrames, StatsPlots

path= raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2025\05.May\\"

filename= "v5_perimeterR2.0_diffequators_50runs.png"
#filename= "velocity.png"
f= joinpath(path, filename)
f1= joinpath(path, "boxplot_"*filename)
fv= path*"peak_frequencies_combined_pf0.2_1000_20000_perimeterR2.0.csv"

df= CSV.read(fv,DataFrame)
#x= ["30s" "60s"  "120s"]
x= df[:,:Peak_Frequency_mHz]
#y= [45.24 52.77 41.54]
#  y= df[:,:FFTcurvature]
# x= df[:,:Run]
 p1=histogram(x,xlabel="Frequency(mHz)",ylabel="Count",ylims= (0,10), xtickfont=font(12), ytickfont=font(12),color= :lightblue, legend= false,bins=0:0.05:2)
#p1=histogram(x,xlabel="velocity (μm/s)",ylabel="Number of particles",ylims= (0,10), xtickfont=font(12), ytickfont=font(12),color= :lightblue, legend= false,bins=0:0.5:15)
 #println("Mean of frequencies: ", mean(x))
  #  println("Standard deviation of frequencies: ", std(x))

# Create a boxplot
#p = @df df boxplot(:Peak_Frequency_mHz, legend=false,ylabel="Frequency(mHz)", ylims=(0, 5),title="V5", bins=1:0.1:1, bar_width=0.001)

#plot!(p,  kind="box")
savefig(p1, f)
#savefig(p, f1)

