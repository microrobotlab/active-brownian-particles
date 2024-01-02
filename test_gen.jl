using LinearAlgebra
using Plots
using Random

a = 10.
b = 5.
R = 0.

Np = 10000

α = rand(Np).*2π

r_max= (a-R)*(b-R)./(sqrt.((((a-R)*sin.(α)).^2) .+ ((b-R)*cos.((α))).^2))

r = rand(Np).*r_max

x = r.*cos.(α)
y = r.*sin.(α)


display(histogram(r, bins = 20))
# display(histogram2d(x,y,bins=(50, 50), show_empty_bins=false, color = :plasma, aspect_ratio=:equal))
# display(scatter(x,y, legend = false))

# p = hline!([2π], linestyle=:dash)


