using LinearAlgebra
using Plots
using Random

a = 10.
b = 5.
R = 0.

Np = 100000
## First way. Problem: center biased. Advantage: fast, simple and correct
# α = rand(Np).*2π

# r_max= (a-R)*(b-R)./(sqrt.((((a-R)*sin.(α)).^2) .+ ((b-R)*cos.((α))).^2))

# r = rand(Np).*r_max

# x = r.*cos.(α)
# y = r.*sin.(α)


# # display(histogram(r, bins = 20))
# display(histogram2d(x,y,bins=(50, 50), show_empty_bins=false, color = :plasma, aspect_ratio=:equal))
# # display(scatter(x,y, legend = false))

# # p = hline!([2π], linestyle=:dash)


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
    rmax = (a-R)*(b-R)/(sqrt((((a-R)*sin(θt))^2) + ((b-R)*cos((θt)))^2))
    if rt < rmax
        global xyacc = vcat(xyacc, xyt)
        global Nacc = size(xyacc,1)
    end
end
println(Nacc)
x = xyacc[:,1]
y = xyacc[:,2]
display(histogram2d(x,y,bins=(50, 50), show_empty_bins=false, color = :plasma, aspect_ratio=:equal))

# θ = LinRange(0,2π, 1000)
# rmax = (a-R)*(b-R)./(sqrt.((((a-R)*sin.(θ)).^2) .+ ((b-R)*cos.((θ))).^2))
# xt = rmax.*cos.(θ)
# yt = rmax.*sin.(θ)

# display(plot(xt,yt),)

