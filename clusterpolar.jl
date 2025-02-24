using CairoMakie, CSV, DataFrames, Clustering, JLD2, Distributions
using StatsBase
include("ABP orderparameters.jl")
include("ABP plot_animation.jl")
include("ABP radialdistribution.jl")
include("ABP main interactions heun.jl")

path = "D:\\nic_simulations\\lj_offcenter"
datestamps = readdir(path)

nbins = 1000
mindist_cluster = 5.
L = 175.
R = 2.

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
        # tocheck.Time[1]%1000 == 0 && println(tocheck.Time[1])
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

maximum_polarization = Float64[]
maxpoltime = Float64[]
maxpolloc = Float64[]
maxpolloc_time = Float64[]
firstpeak_height = Float64[]
fpk_pos = Float64[]
oc = Float64[]

for (i,datestamp ) in enumerate(datestamps)
    datestamp == "imgs" && continue
    info_dict = JLD2.load(joinpath(path, datestamp, "siminfo_dict.jld2"))
    push!(oc, info_dict["offcenter"])

    # distr_str = "$(info_dict["v"])"
    # sigmaloc = findfirst("σ=", distr_str)[end] + 1
    # sigmastop = findfirst(")", distr_str)[end] - 1
    # σ = parse(Float64, distr_str[sigmaloc:sigmastop])
    velocity = info_dict["v"]

    offcenter = info_dict["offcenter"]
    dt = info_dict["dt"]
    # println(info_dict["dt"])
    if i >3
        offcenter = -offcenter
    end
    # for p in readdir(joinpath(path,datestamp,"data"), join = true)
    #     if isdir(p)
    #         numloc = findfirst("run", p)[end] + 1
    #         num = p[numloc:end]
    #         df = CSV.read(joinpath(p, datestamp*"_run$num.txt"), DataFrame, ntasks = 16)
    #         df_ds = df[df.Time .% 1000 .== 0,:]
    #         CSV.write(joinpath(p, datestamp*"_run$num"*"ds"*".txt"), df_ds)
    #         println("$p done")
    #     end
    # end
    # continue
    df_list = []
    df_cl_list = []
    df_or_list = []
    df_meanor_list = []
    df_t_list = []
    df_nclust_list = []


    for p in readdir(joinpath(path,datestamp,"data"), join = true)
        if isdir(p)
            numloc = findfirst("run", p)[end] + 1
            num = p[numloc:end]
            df = CSV.read(joinpath(p, datestamp*"_run$num"*"ds"*".txt"), DataFrame, ntasks = 16)    
            # println("done")
            push!(df_list, df)
            df_cl = transform(groupby(df, :Time), [:xpos, :ypos] => ((x,y) -> assignments(dbscan([x y]', mindist_cluster, min_cluster_size = 2))) => :dbscan)

            periodic_clustering!(df_cl)

            df_nclust = combine(groupby(df_cl[df_cl.dbscan .!== 0,:], [:Time, :dbscan]), nrow=>:partcount)
            df_nclust = combine(groupby(df_nclust, :Time), :partcount => maximum => :max_cl_size, nrow => :nclust)
            local df_or = combine(groupby(df_cl, [:Time, :dbscan]), :orientation => mean_polarization => :polar, nrow => :partcount)
            local df_meanor = combine(groupby(df_or[df_or.dbscan .!== nothing,:], [:Time]), [:polar, :partcount] => ((x,y) -> mean(x,weights(y))) => :mean_polar)
            local df_meanor = combine(groupby(df_or[df_or.dbscan .!== nothing,:], [:Time]), :polar => mean => :mean_polar)
            local df_t = combine(groupby(df, :Time), :orientation => mean_polarization => :polar,)

            push!(df_cl_list, df_cl)
            push!(df_or_list, df_or)
            push!(df_meanor_list, df_meanor)
            push!(df_t_list, df_t)
            push!(df_nclust_list, df_nclust)
        end
    end


    plot_one_timestep(df_list[1], 2, 175, Int(size(df_list[2], 1)/250), savepath = joinpath(path, "imgs", "situa$(velocity)_oc$(offcenter).png"))
    

    df_meannclust = combine(groupby(vcat(df_nclust_list...), :Time), :max_cl_size => mean => :max_cl_size, :nclust => mean => :nclust, :max_cl_size => std => :max_cl_size_std, :nclust => std => :nclust_std)

    df_meanor = combine(groupby(vcat(df_meanor_list...), :Time), :mean_polar =>mean => :mean_polar, :mean_polar => std => :std_polar)
    df_t = combine(groupby(vcat(df_t_list...), :Time), :polar =>mean => :mean_polar, :polar => std => :std_polar)

    CairoMakie.activate!(inline = true)

    maxlocpol= maximum(df_meanor.mean_polar)
    println("local polar $maxlocpol")
    println("at time $(df_meanor.Time[findfirst(x->x>0.95*maxlocpol, df_meanor.mean_polar)]*dt)")
    push!(maxpolloc, maxlocpol)
    push!(maxpolloc_time, df_meanor.Time[findfirst(x->x>0.95*maxlocpol, df_meanor.mean_polar)]*dt)
    maxpol = maximum(df_t.mean_polar)
    println("maximum polarizarion $maxpol")
    println("at time $(df_t.Time[findfirst(x->x>0.95*maxpol, df_t.mean_polar)]*dt)")
    push!(maximum_polarization, maxpol)
    push!(maxpoltime, df_t.Time[findfirst(x->x>0.95*maxpol, df_t.mean_polar)]*dt)

    p2 = Figure(fontsize = 20)
    ax = Axis(p2[1,1], limits = ((0,500), nothing),xlabel ="Time [s]", ylabel = "Polarization", xgridvisible = false, ygridvisible = false)
    lines!(ax, df_meanor.Time*dt, df_meanor.mean_polar, linewidth = 1., color = :blue, label = "local") #ribbon = df_meanor.std_polar, 
    band!(ax, df_meanor.Time*dt, df_meanor.mean_polar.-df_meanor.std_polar, df_meanor.mean_polar.+df_meanor.std_polar, color = :blue, alpha = 0.3)
    lines!(ax, df_t.Time*dt, df_t.mean_polar, linewidth = 1., color = :green, label = "global") #ribbon = df_meanor.std_polar, 
    band!(ax, df_t.Time*dt, df_t.mean_polar.-df_t.std_polar, df_t.mean_polar.+df_t.std_polar, color = :green, alpha = 0.3)
    # axislegend(ax)

    display(p2)

    p3 = Figure(fontsize = 20)
    ax1 = Axis(p3[1,1], limits = ((0,500), nothing), xlabel ="Time [s]", ylabel = "Max cluster size", xgridvisible = false, ygridvisible = false)
    ax2 = Axis(p3[2,1], limits = ((0,500), nothing), xlabel ="Time [s]", ylabel = "Number of clusters", xgridvisible = false, ygridvisible = false)
    linkxaxes!(ax1,ax2)

    lines!(ax1, df_meannclust.Time*dt, df_meannclust.max_cl_size, linewidth = 2., color = :black, label = "Maximum cluster size")
    band!(ax1, df_meannclust.Time*dt, df_meannclust.max_cl_size.-df_meannclust.max_cl_size_std, df_meannclust.max_cl_size.+df_meannclust.max_cl_size_std, color = :black, alpha = 0.3)
    lines!(ax2, df_meannclust.Time*dt, df_meannclust.nclust, linewidth = 2., color =:orange, label = "Number of clusters")
    band!(ax2, df_meannclust.Time*dt, df_meannclust.nclust.-df_meannclust.nclust_std, df_meannclust.nclust.+df_meannclust.nclust_std, color = :orange, alpha = 0.3)
    # axislegend(ax2)
    # axislegend(ax1)

    display(p3)

    rdf = combine(groupby(df_list[2], :Time), [:xpos, :ypos] => ((x,y) -> (radialdistributionfunction(copy(x),copy(y),R,L,nbins),)) => :RadialDistributionFunction)
    transform!(rdf, :RadialDistributionFunction => ByRow(x -> x[1]) => :RadialDistributionFunction)
    rs,bins = radialbinssquare(L,nbins)
    bindata = DataFrame(Radius = rs, BinArea = bins)
    p4 = Figure(fontsize = 20)
    ax = Axis(p4[1,1], xlabel = "Distance [μm]", ylabel = "Radial Distribution Function", limits = (nothing, (-1, 75)), xgridvisible = false, ygridvisible = false)
    lines!(bindata.Radius, rdf.RadialDistributionFunction[end], linewidth = 2.)

    firstpeak = maximum(rdf.RadialDistributionFunction[end])
    firstpeak_pos = rs[argmax(rdf.RadialDistributionFunction[end])]
    println("first peak $firstpeak")
    println("fp position $firstpeak_pos")
    push!(firstpeak_height, firstpeak)
    push!(fpk_pos, firstpeak_pos)

    display(p4)

    # save(joinpath(path,"imgs", "polar$(velocity)_oc$(offcenter).png"),p2)
    # save(joinpath(path,"imgs", "cluster$(velocity)_oc$(offcenter).png"),p3)
    # save(joinpath(path,"imgs", "rdf$(velocity)_oc$(offcenter).png"),p4)
end


begin
    fig = Figure(fontsize = 20)
    ax = Axis(fig[1,1], xticks = [0.25,0.5,0.75], xgridvisible=false,ygridvisible=false, xlabel = "α", ylabel = "Transient Time [s]")
    CairoMakie.scatter!(ax, oc[2:4], maxpoltime[2:4], markersize = 20, color = :black)

    display(fig)
    # save(joinpath(path,"imgs", "transient_time_alpha.pdf"),fig)
end
