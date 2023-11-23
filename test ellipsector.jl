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
# println(ellipse_sector_area(1,0.5,0.4,0.2))
# println(ellipse_sector(1,0.5,0.6))



## MC checking
a,b = (1,0.5)
N = 10000
xy = (rand(N,2).-0.5).*repeat(2 .*[a b],N)

r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
rₚ = sqrt.(r)   
α =atan.(xy[:,2], xy[:,1]) 

rₑ = (a)*(b)./(sqrt.((((a)*sin.(α)).^2) .+ ((b)*cos.((α))).^2))

id = (rₚ .< (rₑ))
xy1 = [xy[id,1] xy[id,2]]

rₜ = sqrt.((xy1[:,1]).*(xy1[:,1]) + (xy1[:,2]).*(xy1[:,2]))

# idx = (rₜ .>= 0.6 .&& rₜ .< 0.7)
# xy2 = [xy1[idx,1] xy1[idx,2]]
# print(size(xy2,1))
# p = scatter(xy1[:,1], xy1[:,2], legend = false, aspect_ratio = :equal)
# p = scatter!(xy2[:,1], xy2[:,2], legend = false, aspect_ratio = :equal)

# p



rv = [[0.3, 0.4] [0.4, 0.6] [0.6, 0.7]]
for i in eachcol(rv)
    r1 = i[1]
    r2 = i[2]

    idx = (rₜ .>= r1  .&& rₜ.<= r2)

    xy2 = [xy1[idx,1] xy1[idx,2]]


    el_area = pi*a*b
    c = round(intersection_area(a,b,r1,r2)/el_area, digits = 5)
    d = round(size(xy2, 1)/size(xy1,1),digits = 5)

    println("$c,$d,$(size(xy2,1))")
end


# p = scatter(xy1[:,1], xy1[:,2], legend = false, aspect_ratio = :equal)
# p = scatter!(xy2[:,1], xy2[:,2], legend = false, aspect_ratio = :equal)
# p
