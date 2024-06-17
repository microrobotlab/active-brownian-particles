#TEsting
using Distances, LinearAlgebra

pwdist(x) =[a-b for a in x, b in x] 
function radial_directions(xy::Array{Float64,2})
    diff_x = pwdist(xy[:,1])./pairwise(Euclidean(), xy, xy, dims = 1)
    diff_y = pwdist(xy[:,2])./pairwise(Euclidean(), xy, xy, dims = 1)
    replace!(diff_x, NaN=>0.)
    replace!(diff_y, NaN=>0.)
    return diff_x, diff_y
end

#Function to calculate force vectors
function interactions(xy::Array{Float64,2}, R::Float64)
    ϵ=.1
    σ= 2*R
    k = .1 
    Np = size(xy,1)
    dists = zeros(Np,Np)

    dists .= pairwise(Euclidean(),xy,xy,dims=1)
    
    strength_param = 1.
    force = strength_param.*lennard_jones.(dists, σ, ϵ)
    replace!(force, NaN => 0.)

    dirs = radial_directions(xy)
    F_x = force.*dirs[1]
    F_y = force.*dirs[2]
    ΣFx = .-sum(F_x, dims = 1)
    ΣFy = .-sum(F_y, dims = 1)
    ΣF = vcat.(ΣFx, ΣFy)    
    return  reduce(vcat, transpose(ΣF))
end

lennard_jones(x, σ, ϵ) = 24*ϵ*(((2*σ^(12))./(x.^(13))).- (σ^(6)./(x.^(7))))

# Np = 10
# L = 10
# xy=(rand(10,2).-0.5).*repeat([L L], 10)

# interactions(xy,1.)
