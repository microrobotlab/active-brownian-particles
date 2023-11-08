using Distances
using LinearAlgebra
using Plots
using Random

function radial_density_sq(xy::Array{Float64,2}, nshells::Int64, L::Float64)
    thickness = .5*L/nshells
    lims = [i*thickness for i in 0:nshells-1]
    nums = zeros(nshells)
    for i in 1:nshells
        id = (abs.(xy[:,1]) .> lims[i]) .| (abs.(xy[:,2]) .> lims[i])
        nums[i] = sum(id)
    end
    area = (lims.+thickness).^2 .- lims.^2
    nc = [nums[i]-nums[i+1] for i in 1:length(nums)-1]
    nc = vcat(nc, nums[end])

    return nc./area
end

function radial_density_cr(xy::Array{Float64,2}, nshells::Int64, R::Float64)
    thickness = 0.5R/nshells
    lims = [i*thickness for i in 0:nshells-1]

    r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
  
    rₚ = sqrt.(r)   
    # α =atan.(xy[:,2], xy[:,1]) 
    nums = zeros(nshells)

    for i in 1:nshells
        id = rₚ .> lims[i]
        nums[i] = sum(id)
    end
    area = π*((lims.+thickness).^2 .- lims.^2)
    nc = [nums[i]-nums[i+1] for i in 1:length(nums)-1]
    nc = vcat(nc, nums[end])

    return nc./area
end

# function radial_density_el(xy::Array{Float64,2}, nshells::Int64, L::Float64)
#     r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
#     rₚ = sqrt.(r)   
#     α =atan.(xy[:,2], xy[:,1]) 

    
#     a = L/2 #These are the proportions used in the main code
#     b = L/4
#     rₑ = (a-R)*(b-R)./(sqrt.((((a-R)*sin.(α)).^2) .+ ((b-R)*cos.((α))).^2))
#     nums = zeros(nshells)

#     for i in 1:nshells
#         id = rₚ .> lims[i]
#         nums[i] = sum(id)
#     end
#     area = π*((lims.+thickness).^2 .- lims.^2)
#     nc = [nums[i]-nums[i+1] for i in 1:length(nums)-1]
#     nc = vcat(nc, nums[end])

#     return nc
# end



# rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
function circle(x,y,r)
    θ = 0:1:360
    Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
end



Random.seed!(42)
L = 10.
xy = L.*(rand(Float64, (5,2)).-0.5)



# xy = hcat(1:5, 1:5)
# println(xy)

dens = radial_density_cr(convert(Array{Float64}, xy),5,L)

print(dens)

p = scatter(xy[:,1], xy[:,2], legend = false, aspect_ratio = :equal)
for l in 1:5
plot!(circle(0,0,l), opacity=.2)
end
p
