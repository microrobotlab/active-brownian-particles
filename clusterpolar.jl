using Plots, CSV, DataFrames, Clustering
using StatsBase
include("ABP orderparameters.jl")

path = "C:\\Users\\nikko\\OneDrive\\Documents\\Uni\\magistrale\\tesi\\simulations"

datestamp = "20250117-123605"

df_list = []
for p in readdir(joinpath(path,datestamp,"data"), join = true)
    if isdir(p)
        num = p[91:end]
        df = CSV.read(joinpath(p, datestamp*"_run$num.txt"), DataFrame)
        push!(df_list, df)
    end
end
df = df_list[1]

time = 1000

df1 = df[df.Time .== time, :]

xy = [df1.xpos df1.ypos]



res = dbscan(xy', 5, min_cluster_size = 2)

p = scatter(xy[:,1], xy[:,2], ratio = 1)
for i in 1:maximum(df.N)
    annotate!(p, xy[i,1], xy[i,2], text(string(df.N[i]), 8, :black, :left))
end
display(p)

df_cl = transform(groupby(df, :Time), [:xpos, :ypos] => ((x,y) -> assignments(dbscan([x y]', 5, min_cluster_size = 2))) => :dbscan)

df_cl[df_cl.Time .== time,:]

df_or = combine(groupby(df_cl, [:Time, :dbscan]), :orientation => mean_polarization => :polar, nrow => :partcount)
