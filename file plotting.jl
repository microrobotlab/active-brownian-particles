#function to create plots, histogram and box plot for a data file 

using Plots, CSV,DataFrames, StatsPlots

path= raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\\"

csv_path   = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.05\final_all_zero_intervals_R1.5_v5.0_pf0.05_1000-10000.csv"
df = CSV.read(csv_path, DataFrame)

    # required columns
    for col in (:aspect, :interval)
        col in propertynames(df) || error("Column $col not found")
    end

    # optional filter (remove very small intervals if desired)
    df = filter(row -> row.interval > 20, df)

    # median interval per aspect
    
    summary = combine(groupby(df, :aspect), :interval => mean => :median_interval)
    sort!(summary, :aspect)

    println("Median intervals:")
    println(summary)

#filename= "v5_perimeterR2.0_diffequators_50runs.png"
filename= "curvature1.png"
f= joinpath(path, filename)
#f1= joinpath(path, "boxplot_"*filename)
#fv= path*"peak_frequencies_combined_pf0.2_1000_20000_perimeterR2.0.csv"

#df= CSV.read(fv,DataFrame)
x= ["0.1", "0.25",  "0.4", "0.55"]
#x = string.(round.(summary.aspect, digits=2))
x = parse.(Float64, x)
#x= df[:,:Peak_Frequency_mHz]
y= [1024, 64, 16, 6]
logK   = log10.(y)
logT   = log10.(summary.median_interval)
# Local slopes between consecutive points
slopes = diff(logT) ./ diff(logK)

# Mid-points of each segment (where we will place the slope label)
mid_logK = (logK[1:end-1] + logK[2:end]) ./ 2
mid_logT = (logT[1:end-1] + logT[2:end]) ./ 2

# Length of the small parallel line that indicates the slope (in data units)
line_len = 0.08

p = plot(logK, logT,
    xlabel = L"$K$",
    ylabel = L"Mean time period to cross $\kappa=0$ (s)",
    labelsize = 18,
    seriestype = :path,
    linewidth = 2,
    marker = :circle,
    markersize = 6,
    legend = false,
    # optional: make the axes look more like a typical log–log figure
    # xticks = (logK, string.(round.(K_vals, digits=0))),
    # yticks = ...
)

# Draw a short parallel segment + slope value for each regime
for i in 1:length(slopes)
    # direction of the local slope
    dx = line_len / sqrt(1 + slopes[i]^2)
    dy = slopes[i] * dx

    # short line centered on the mid-point
    plot!(p,
        [mid_logK[i] - dx/2, mid_logK[i] + dx/2],
        [mid_logT[i] - dy/2, mid_logT[i] + dy/2],
        color = :black, linewidth = 1.5, label = "")

    # slope value written next to the line
    annotate!(p, mid_logK[i] + 0.04, mid_logT[i] + 0.04,
              Plots.text("$(round(slopes[i], digits=2))", 10, :black, :left))
end

savefig(p, f)
println("Saved plot → $f")
#  y= df[:,:FFTcurvature]
# x= df[:,:Run]
# p1=histogram(x,xlabel="Frequency(mHz)",ylabel="Count",ylims= (0,10), xtickfont=font(12), ytickfont=font(12),color= :lightblue, legend= false,bins=0:0.05:2)
#p1=histogram(x,xlabel="velocity (μm/s)",ylabel="Number of particles",ylims= (0,10), xtickfont=font(12), ytickfont=font(12),color= :lightblue, legend= false,bins=0:0.5:15)
 #println("Mean of frequencies: ", mean(x))
  #  println("Standard deviation of frequencies: ", std(x))

# Create a boxplot
#p = @df df boxplot(:Peak_Frequency_mHz, legend=false,ylabel="Frequency(mHz)", ylims=(0, 5),title="V5", bins=1:0.1:1, bar_width=0.001)
#=
p = plot( log10.(summary.median_interval),log10.(y),
         ylabel = L"Mean time period to cross $\kappa=0$ (s)",
         xlabel = L"$K$",
         labelsize = 18,
         seriestype = :path,     # ← forces a line
         linewidth = 2,
         legend = false,
         marker = :circle,       # optional markers on the line
         markersize = 6)
#plot!(p,  kind="box")
savefig(p, f)
#savefig(p, f1)

=#