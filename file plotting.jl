#function to create plots, histogram and box plot for a data file 

using Plots, CSV,DataFrames, StatsPlots

path= raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\\"

#filename= "v5_perimeterR2.0_diffequators_50runs.png"
filename= "curvature.png"
f= joinpath(path, filename)
#f1= joinpath(path, "boxplot_"*filename)
#fv= path*"peak_frequencies_combined_pf0.2_1000_20000_perimeterR2.0.csv"

#df= CSV.read(fv,DataFrame)
x= ["0.1", "0.25",  "0.4", "0.55"]
#x= df[:,:Peak_Frequency_mHz]
y= [1024, 64, 16, 6]
#  y= df[:,:FFTcurvature]
# x= df[:,:Run]
# p1=histogram(x,xlabel="Frequency(mHz)",ylabel="Count",ylims= (0,10), xtickfont=font(12), ytickfont=font(12),color= :lightblue, legend= false,bins=0:0.05:2)
#p1=histogram(x,xlabel="velocity (μm/s)",ylabel="Number of particles",ylims= (0,10), xtickfont=font(12), ytickfont=font(12),color= :lightblue, legend= false,bins=0:0.5:15)
 #println("Mean of frequencies: ", mean(x))
  #  println("Standard deviation of frequencies: ", std(x))

# Create a boxplot
#p = @df df boxplot(:Peak_Frequency_mHz, legend=false,ylabel="Frequency(mHz)", ylims=(0, 5),title="V5", bins=1:0.1:1, bar_width=0.001)
p = plot(x, y,
         xlabel = "b/a",
         ylabel = L"$K$",
         labelsize = 18,
         seriestype = :path,     # ← forces a line
         linewidth = 2,
         legend = false,
         marker = :circle,       # optional markers on the line
         markersize = 6)
#plot!(p,  kind="box")
savefig(p, f)
#savefig(p, f1)

