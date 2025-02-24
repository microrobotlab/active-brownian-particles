using Plots, CSV, DataFrames, Clustering, JLD2
using StatsBase
include("ABP orderparameters.jl")
include("ABP plot_animation.jl")
include("ABP radialdistribution.jl")
include("ABP main interactions heun.jl")

path = "D:\\nic_simulations\\lj_offcenter"
readdir(path)
datestamps = readdir(path)
dt = 5e-3
R = 2.
L = 175.
nbins = 100
mindist_cluster = 5.

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
    df_clx = transform(groupby(df_sx, :Time), [:xpos, :ypos] => ((x,y) -> assignments(dbscan([x y]', mindist_cluster, min_cluster_size = 2))) => :dbscan)
    df_cly = transform(groupby(df_sy, :Time), [:xpos, :ypos] => ((x,y) -> assignments(dbscan([x y]', mindist_cluster, min_cluster_size = 2))) => :dbscan)

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


for (i,datestamp ) in enumerate(datestamps)
    datestamp == "imgs" && continue

    info_dict = JLD2.load(joinpath(path, datestamp, "siminfo_dict.jld2"))
    velocity = info_dict["v"]
    offcenter = info_dict["offcenter"]
    # if i >3
    #     offcenter = -offcenter
    # end

    df_list = []
    df_cl_list = []
    df_or_list = []
    df_meanor_list = []
    df_t_list = []
    for p in readdir(joinpath(path,datestamp,"data"), join = true)
        if isdir(p)
            numloc = findfirst("run", p)[end] + 1
            num = p[numloc:end]
            df = CSV.read(joinpath(p, datestamp*"_run$num.txt"), DataFrame, ntasks = 16)
            df_ds = df[df.Time .% 100 .== 0,:]
            println("done")
            push!(df_list, df_ds)
            df_cl = transform(groupby(df_ds, :Time), [:xpos, :ypos] => ((x,y) -> assignments(dbscan([x y]', mindist_cluster, min_cluster_size = 2))) => :dbscan)

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


    p2 = Plots.plot(df_meanor.Time*dt, df_meanor.mean_polar, linewidth = 1., ribbon = df_meanor.std_polar, label = "local")
    Plots.plot!(p2, df_t.Time*dt, df_t.mean_polar, linewidth = 2., ribbon = df_t.std_polar, label = "global")
    Plots.xlabel!(p2, "Time [s]")
    Plots.ylabel!(p2, "Polarization")
    Plots.title!(p2, "v = $velocity μm s⁻¹, α = $offcenter")
    Plots.xlims!(p2, (0,250))

    display(p2)
    # savefig(p2, joinpath(path,"imgs", "polar_oc$(offcenter)_v$(velocity).png"))

    plot_one_timestep(df_list[2], 2, 175, 200)

    # rdf = combine(groupby(df_list[2], :Time), [:xpos, :ypos] => ((x,y) -> (radialdistributionfunction(copy(x),copy(y),R,L,nbins),)) => :RadialDistributionFunction)
    # transform!(rdf, :RadialDistributionFunction => ByRow(x -> x[1]) => :RadialDistributionFunction)
    # rs,bins = radialbinssquare(L,nbins)
    # bindata = DataFrame(Radius = rs, BinArea = bins)
    # Plots.plot(bindata.Radius, rdf.RadialDistributionFunction[200], xlabel = "Distance [μm]")

    p3 = Plots.plot(df_meanor.Time, combine(groupby(df_or_list[1][df_or_list[1].dbscan .!== 0,:], [:Time]), :partcount => maximum => :partcount).partcount, label = "max", linewidth = 2.)
    # Plots.plot!(p3, df_meanor.Time, combine(groupby(df_or_list[1][df_or_list[1].dbscan .!== 0,:], [:Time]), :partcount => mean => :partcount).partcount, label = "mean")
    Plots.title!(p3, "v = $velocity μm s⁻¹, α = $offcenter")
    display(p3)
    # savefig(p3, joinpath(path,"imgs", "cluster_oc$(offcenter)_v$(velocity).png"))

end
