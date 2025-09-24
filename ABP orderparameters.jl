using Clustering
using Markdown
using StatsBase
# include("ABP main interactions heun.jl")

"""
    mean_polarizantion(θ::Vector{Float64})

Calculates the mean polarization of a set of orientations.
Returns 1 for perfect alignment and 0 for no alignment.

# Examples
```julia-repl

julia> a = [π/2, π/2]
2-element Vector{Float64}:
 1.5707963267948966
 1.5707963267948966

julia> mean_polarizantion(a)
1.0

julia> b = [π/2, -π/2]
2-element Vector{Float64}:
  1.5707963267948966
 -1.5707963267948966

julia> mean_polarizantion(b)
6.123233995736766e-17

julia> mean_polarizantion(rand(1000)*2pi)
0.010532168326974655
```
"""
mean_polarization(θ) = abs(mean(exp.(im*θ)))

"""s
    orient_corr_func(x::Vector{Float64}, y::Vector{Float64}, θ::Vector{Float64})
In development


# Examples
```julia-repl
```
"""
function orient_corr_func(x::Vector{Float64}, y::Vector{Float64}, θ::Vector{Float64}, k::Int, thickness::Float64)
    rmin, rmax = k*thickness.*(1-1/2k, 1+1/2k)

end

"""
    cluster_size_distribution(xy::Array{Float64,2}, intrange::Float64)
Takes a position matrix and a range as arguments and returns the cluster size distribution using DBSCAN clustering algorithm.
"""
function cluster_size_distribution(xy::Array{Float64,2}, intrange::Float64)
    res = dbscan(xy', intrange, min_cluster_size = 2)
    return res.counts
end 

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