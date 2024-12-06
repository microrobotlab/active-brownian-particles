using Clustering
using Markdown
include("ABP main interactions opt.jl")

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
mean_polarizantion(θ::Vector{Float64}) = abs(mean(exp.(im*θ)))

"""
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
    cs_distribution = countmap(res.assignments)
    return cs_distribution
end 