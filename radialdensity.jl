using Distances
using LinearAlgebra
using Plots
using Random

function intersection_area(a,b,R₁,R₂) #a larger semiaxis, b smaller semiaxis, R₁ smaller circle radius

    if R₂ <= b
        area = pi*(R₂^2-R₁^2)
    end

    if R₁ < b && R₂ > b
        area = ellipse_sector(a,b,R₂)-(pi*R₁^2)
    end

    if R₁ >= b
        area = ellipse_sector(a,b,R₂) - ellipse_sector(a,b,R₁)
    end

    return area
end


function ellipse_sector(a,b,r) # area of intersection between ellipse and circle with r >b
    x = a*sqrt((r^2-b^2)/(a^2-b^2))
    y = b*sqrt((a^2-r^2)/(a^2-b^2))

    θ₁ = atan(y,x)
    θ₂ = atan(y,-x)
    Aₑ = 0.5a*b*atan(a*tan(θ₁)/b)
    Aᵢₙₜ = pi*a*b - 4*Aₑ
    # Aₑ = a*b*(θ₂-θ₁-atan((b-a)*sin(2θ₂)/(a+b+(b-a)*cos(2*θ₂)))+atan((b-a)*sin(2θ₁)/(a+b+(b-a)*cos(2*θ₁))))

    Aₒ = 2*θ₁*r^2

    return Aᵢₙₜ+Aₒ
end

function radial_density_sqsq(xy::Array{Float64,2}, nshells::Int64, L::Float64) #square norm!
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
    thickness = R/nshells
    add(x) = x+thickness
    lims = [i*thickness for i in 0:nshells-1]
    lims = hcat(lims, map(add, lims))

    r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
    rₚ = sqrt.(r)   
    nums = zeros(nshells)

    for row in eachrow(lims)
        r1 = row[1]
        r2 = row[2]
    
        idx = (rₚ .>= r1  .&& rₚ.<= r2)

        nums[i] = sum(idx)
        areas[i] = intersection_area(a,b,r1,r2)
        i+=1 # smarter way to do that?
    end
    area = π*(lims[:,2].^2 .- lims[:,1].^2)

    return nc./area
end

function radial_density_el(xy::Array{Float64,2}, nshells::Int64, a::Float64, b::Float64)
    r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
    rₚ = sqrt.(r)   

    thickness = b/nshells
    c = [b*i/nshells for i in 0:nshells-1]
    add(x) = x+thickness
    lims = hcat(c, map(add,c))

    nums = zeros(size(lims,1))
    areas = zeros(size(lims,1))
    i = 1
    for row in eachrow(lims) 
        r1 = row[1]
        r2 = row[2]
    
        idx = (rₚ .>= r1  .&& rₚ.<= r2)

        nums[i] = sum(idx)
        areas[i] = intersection_area(a,b,r1,r2)
        i+=1 # smarter way to do that?
    end

    return nums./areas
end

Random.seed!(666)
## MC checking
a,b = (1.,0.5)
N = 100
xy = (rand(N,2).-0.5).*repeat(2 .*[a b],N)

r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
rₚ = sqrt.(r)   
α =atan.(xy[:,2], xy[:,1]) 

rₑ = (a)*(b)./(sqrt.((((a)*sin.(α)).^2) .+ ((b)*cos.((α))).^2))

id = (rₚ .< (rₑ))
xy1 = [xy[id,1] xy[id,2]]

radial_density_el(xy1,10,a,b)



#=
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
=#