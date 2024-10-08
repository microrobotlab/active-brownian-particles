using CSV, DataFrames, Distances, FFTW, LinearAlgebra, StatsBase, Plots
include("geo_toolbox.jl")

function radialbinssquare(L::Float64, nbins::Int)
    maxdist = L*sqrt(2.)/2
    rs = LinRange(maxdist/nbins,maxdist,nbins)
    areas = [intersection_area_sq(r,L) for r in rs] 

    bins = areas-circshift(areas,1) 
    bins[1] = areas[1]
    return rs, bins
end

function radialdistributionfunction(x::Array{Float64}, y::Array{Float64}, R::Float64, L::Float64, nbins::Int) # Very optimized
    rs,bins = RadialBinsSquare(L,nbins)
    Np = length(x)
    dens = Np/(L^2)

    xdist = pwdist(x)
    ydist = pwdist(y)

    threshold = L/2+R
    xind = abs.(xdist) .> threshold
    yind = abs.(ydist) .> threshold
    @. xdist[xind] = -sign(xdist[xind]) * (L - abs(xdist[xind]))
    @. ydist[yind] = -sign(ydist[yind]) * (L - abs(ydist[yind]))
    dist = @. sqrt(xdist^2 + ydist^2)

    r = vec(dist)
    nonzero_r = r[r .!= 0]
    cumulative = [sum(nonzero_r .< rb) for rb in rs]
    parts_in_bins = cumulative .- circshift(cumulative, 1)
    parts_in_bins[1] = cumulative[1]
    rdf = parts_in_bins./(Np*dens*bins)
    return rdf
end
