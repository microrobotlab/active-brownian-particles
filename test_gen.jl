using LinearAlgebra, Plots, Random, StatsBase
include("test ellipsector.jl")

a = 10.
b = 5.
R = 0.

Np = 100000

hist_radius = histogram()
hist_ellipse = histogram2d()

rr = range(0.,9.5,20)
ar = Float64[]
for ri in rr
    push!(ar, intersection_area_el(a,b,ri,ri.+0.5))
end

θ = LinRange(0,2π, 1000)
rmax = (a-R)*(b-R)./(sqrt.((((a-R)*sin.(θ)).^2) .+ ((b-R)*cos.((θ))).^2))
xt = rmax.*cos.(θ)
yt = rmax.*sin.(θ)
p_sc = plot(xt, yt)

## First way. Problem: center biased. Advantage: fast, simple and correct
# α = rand(Np).*2π

# r_max= (a-R)*(b-R)./(sqrt.((((a-R)*sin.(α)).^2) .+ ((b-R)*cos.((α))).^2))

# r = rand(Np).*r_max

# x = r.*cos.(α)
# y = r.*sin.(α)

# histogram!(hist_radius, r, bins = 20)
# histogram2d!(hist_ellipse,x,y,bins=(50, 50), show_empty_bins=false, color = :plasma, aspect_ratio=:equal)
# scatter!(p_sc, x,y, legend = false)

# hline!([2π], linestyle=:dash)

## Second way. Presumably slower

##First generation

xy = (2 .*rand(Np, 2).-1).*repeat([a b], Np)

r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
rₚ = sqrt.(r)   
θ =atan.(xy[:,2], xy[:,1]) 

r_max= (a-R)*(b-R)./(sqrt.((((a-R)*sin.(θ)).^2) .+ ((b-R)*cos.((θ))).^2))
id = (rₚ .< (r_max))
xyacc = [xy[id,1] xy[id,2]]
Nacc = size(xyacc,1)
println(Nacc)

# Nacc = 0
# xyacc = Array{Float64}(undef, 0, 2)

while Nacc < Np
    xyt = (2 .*rand(1, 2).-1).*([a b])
    rt = sqrt((xyt[1])*(xyt[1]) + (xyt[2])*(xyt[2]))
    θt = atan(xyt[2], xyt[1]) 
    local rmax = (a-R)*(b-R)/(sqrt((((a-R)*sin(θt))^2) + ((b-R)*cos((θt)))^2))
    if rt < rmax
        global xyacc = vcat(xyacc, xyt)
        global Nacc = size(xyacc,1)
    end
end

x = xyacc[:,1]
y = xyacc[:,2]
r = sqrt.(x.^2 + y.^2)

histogram2d!(hist_ellipse,x,y,binds=(50, 50), show_empty_bins=false, color = :plasma, aspect_ratio=:equal)
hist = fit(Histogram, r, range(0,10,21))
hist.weights = round.(hist.weights./ar)
plot!(hist_radius, hist)
# scatter!(hist_radius, rr, ar*maximum(hist.weights)/maximum(ar))
# scatter!(p_sc, x,y, legend = false)

display(hist_ellipse)
display(hist_radius)
# display(p_sc)

