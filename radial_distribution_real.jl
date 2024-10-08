using CSV, DataFrames, Distances, FFTW, LinearAlgebra, Peaks, Plots, Random, Statistics
include("geo_toolbox.jl")

function periodic_BC_array!(xy::Array{Float64,2},L::Float64, R)   #when a particle crosses an edge it reappears on the opposite side
	# Boundary conditions: horizontal edge
	idx = abs.(xy[:,1]) .> L/2 + R #I create vector idx in which I have 1 where the absolute value of the x coordinate of the particles is outside the observation area
	if any(idx)
		xy[idx,1] .-= sign.(xy[idx,1]).*L   #where I have uni in idx I make the particle reappear on the opposite side of x with respect to 0
	end
	# Boundary conditions: vertical edge
	idy = abs.(xy[:,2]) .> L/2 + R
	if any(idy)
		xy[idy,2] .-= sign.(xy[idy,2]).*L
	end
	return nothing
end

rdfp = plot()
p2 = plot()

L = 300.
R = 2.
maxdist = L*sqrt(2.)/2
nbins = 1000
rs = LinRange(maxdist/nbins,maxdist,nbins)
areas = [intersection_area_sq(r,L) for r in rs]

bins = areas-circshift(areas,1)
bins[1] = areas[1]

for fname in ["..\\sims_to_analyze\\test2.csv"]
    df = CSV.read(fname, DataFrame)
    Np = maximum(unique(df[!,:N]))
    dens = Np/L^2
    bins_ensemble = Matrix{Int}(undef, nbins,0)
    for time in unique(df[end,:Time])
        inst_df = filter(:Time => t->t==time, df)
        xy = Matrix(inst_df[!, [:xpos, :ypos]])
        for xyc in eachrow(xy)
            xy_shifted = xy.-xyc'
            periodic_BC_array!(xy_shifted, L, R)
            r = sqrt.((xy_shifted[:,1]).*(xy_shifted[:,1]) + (xy_shifted[:,2]).*(xy_shifted[:,2]))
            cumulative = [sum(r[r.!=0].<rb) for rb in rs]
            parts_in_bins = cumulative-circshift(cumulative,1)
            parts_in_bins[1] = cumulative[1]
            bins_ensemble = hcat(bins_ensemble, parts_in_bins)
        end
    end
    rdf = mean(bins_ensemble, dims = 2)./(bins*dens)
    pos_id = vec(rdf.>0)
    # plot!(rdfp, rs[pos_id], rdf[pos_id])
    plot!(rdfp, rs, rdf)
    k = fftshift(fftfreq(length(rs), 1))
    ssf = 1 .+ dens.*abs.(fftshift(fft(rdf.-1.)))
    plot!(p2, k[k.>=0], ssf[k.>=0])#, xaxis=:log, yaxis=:log)
end
# yaxis!(rdfp, :log)
# xlims!(rdfp,0,20)
display(rdfp)
display(p2)
