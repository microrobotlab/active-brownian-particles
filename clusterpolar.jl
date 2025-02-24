using Plots, CSV, DataFrames, Clustering
using StatsBase
include("ABP orderparameters.jl")
include("ABP plot_animation.jl")
include("ABP radialdistribution.jl")
include("ABP main interactions heun.jl")

path = "C:\\Users\\nikko\\OneDrive\\Documents\\Uni\\magistrale\\tesi\\simulations"

datestamp = "20250123-123032"
dt = 5e-3
R = 2.
L = 175.
nbins = 100

function periodic_clustering!(cluster_df)
    df_sx, df_sy = copy(cluster_df), copy(cluster_df)
    df_sx[!,:xpos] = df_sx.xpos .- L/2
    df_sy[!,:ypos] = df_sy.ypos .- L/2

    for d in [df_sx, df_sy]
        xy = [d.xpos d.ypos]
        periodic_BC_array!(xy, L, R)
        d[!,:xpos] = xy[:,1]
        d[!,:ypos] = xy[:,2]
    end
    df_clx = transform(groupby(df_sx, :Time), [:xpos, :ypos] => ((x,y) -> assignments(dbscan([x y]', 4.5, min_cluster_size = 2))) => :dbscan)
    df_cly = transform(groupby(df_sy, :Time), [:xpos, :ypos] => ((x,y) -> assignments(dbscan([x y]', 4.5, min_cluster_size = 2))) => :dbscan)

    Threads.@threads for tocheck in groupby(cluster_df, :Time)
        tocheck.Time[1]%1000 == 0 && println(tocheck.Time[1])
        partcount = countmap(tocheck.dbscan)

        for df_cls in [df_clx[df_clx.Time .== tocheck.Time[1],:], df_cly[df_cly.Time .== tocheck.Time[1],:]]
            df_ors = combine(groupby(df_cls, [:Time, :dbscan]), :orientation => mean_polarization => :polar, nrow => :partcount)
            for row in eachrow(tocheck)
                pnum = row.N
                cl = row.dbscan
                shifted_row = df_cls[df_cls.N .== pnum,:]
                cl_shifted = shifted_row.dbscan[1]
            
                if (cl_shifted != 0) && (df_ors[df_ors.dbscan .== cl_shifted ,:partcount][1] > partcount[cl])
                    row.dbscan = -cl_shifted
                    # println(shifted_row)
                end
            end
            partcount = countmap(tocheck.dbscan)
        end
    end
end

df_list = []
df_cl_list = []
df_or_list = []
df_meanor_list = []
df_t_list = []
for p in readdir(joinpath(path,datestamp,"data"), join = true)
    if isdir(p)
        num = p[91:end]
        df = CSV.read(joinpath(p, datestamp*"_run$num.txt"), DataFrame)
        df_ds = df[df.Time .% 100 .== 0,:]
        push!(df_list, df_ds)
        df_cl = transform(groupby(df_ds, :Time), [:xpos, :ypos] => ((x,y) -> assignments(dbscan([x y]', 4.5, min_cluster_size = 2))) => :dbscan)

        periodic_clustering!(df_cl)

        local df_or = combine(groupby(df_cl, [:Time, :dbscan]), :orientation => mean_polarization => :polar, nrow => :partcount)
        local df_meanor = combine(groupby(df_or[df_or.dbscan .!== nothing,:], [:Time]), [:polar, :partcount] => ((x,y) -> mean(x,weights(y))) => :mean_polar)
        local df_meanor = combine(groupby(df_or[df_or.dbscan .!== nothing,:], [:Time]), :polar => mean => :mean_polar)
        local df_t = combine(groupby(df, :Time), :orientation => mean_polarization => :polar,)

        push!(df_cl_list, df_cl)
        push!(df_or_list, df_or)
        push!(df_meanor_list, df_meanor)
        push!(df_t_list, df_t)
    end
end

df_meanor_list

df_meanor = combine(groupby(vcat(df_meanor_list...), :Time), :mean_polar =>mean => :mean_polar, :mean_polar => std => :std_polar)
df_t = combine(groupby(vcat(df_t_list...), :Time), :polar =>mean => :mean_polar, :polar => std => :std_polar)

# p1 = Plots.plot()
# for d in df_meanor_list
#     Plots.plot!(p1, d.Time*dt, d.mean_polar, linewidth = 2.5)
# end
# for d in df_t_list
#     Plots.plot!(p1, d.Time*dt, d.polar, linewidth = 2.5)
# end

# begin
#     p2 = Plots.plot()
#     Plots.plot!(p2, df_t.Time*dt, df_t.mean_polar, linewidth = 1., ribbon = df_t.std_polar, label = "global")
#     Plots.plot!(p2,df_meanor.Time*dt, df_meanor.mean_polar, linewidth = 1., ribbon = df_meanor.std_polar, label = "local")
#     Plots.xlabel!(p2, "Time [s]")
#     Plots.ylabel!(p2, "Polarization")
#     Plots.xlims!(p2, (0,100))
#     # display(p1)
#     display(p2)
#     # savefig(p2, "C:\\Users\\nikko\\OneDrive\\Documents\\Uni\\magistrale\\tesi\\update\\20250124\\local_polarization.png")
# end

include("ABP radialdistribution.jl")
plot_one_timestep(df_list[1], 2, 100, 1)

# for ts in [1 25 50 100 200]
#     rdfplot = Plots.plot()
#     plot_one_timestep(df_list[1], 2, 100, ts,savepath = "situation$(dt*ts*100)s.png", title = "$(dt*ts*100)s")

#     rdf = combine(groupby(df_list[1], :Time), [:xpos, :ypos] => ((x,y) -> (radialdistributionfunction(copy(x),copy(y),R,L,nbins),)) => :RadialDistributionFunction)
#     transform!(rdf, :RadialDistributionFunction => ByRow(x -> x[1]) => :RadialDistributionFunction)

#     rs,bins = radialbinssquare(L,nbins)
#     bindata = DataFrame(Radius = rs, BinArea = bins)

#     pkindices, properties = findpeaks1d(rdf[ts,:RadialDistributionFunction]; height = 1.,)

#     Plots.plot!(rdfplot, bindata.Radius, rdf.RadialDistributionFunction[ts], xlabel = "Distance [μm]", ylabel = "g(r)", title = "$(dt*ts*100) s")
#     Plots.scatter!(rdfplot, bindata.Radius[pkindices], rdf.RadialDistributionFunction[ts][pkindices], label = false)
#     Plots.annotate!(bindata.Radius[pkindices[1]] + 5, rdf.RadialDistributionFunction[ts][pkindices[1]], ("peak position $(round(bindata.Radius[pkindices[1]], digits = 2))", :left, 7))
#     display(rdfplot)
#     savefig("rdf_$(dt*ts*100)s.png")
# end




# p3 = Plots.plot(df_meanor.Time, combine(groupby(df_or_list[1][df_or_list[1].dbscan .!== 0,:], [:Time]), :partcount => maximum => :partcount).partcount, label = "max")
# Plots.plot!(p3, df_meanor.Time, combine(groupby(df_or_list[1][df_or_list[1].dbscan .!== 0,:], [:Time]), :partcount => mean => :partcount).partcount, label = "mean")
