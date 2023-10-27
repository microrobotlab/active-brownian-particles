using Distances
using LinearAlgebra

function diffusion_coeff(R::Float64, T::Float64=300.0, η::Float64=1e-3)
    # Boltzmann constant [J/K]
    kB = 1.38e-23
    # friction coefficient [Ns/m]
    γ = 6*pi*R*η
    # translational diffusion coefficient [m^2/s]
    DT = kB*T/γ
    # rotational diffusion coefficient [rad^2/s]
    DR = 6*DT/(8*R^2)
    return DT, DR, γ
end

# function radial_directions(xy::Array{Float64,2})
#     Np = size(xy,1)
#     dists = zeros(Np,Np)
#     dists .= pairwise(Euclidean(),xy,xy,dims=1)

#     diff_x = [a[1]-b[1] for a in eachrow(xy), b in eachrow(xy)]
#     diff_y = [a[2]-b[2] for a in eachrow(xy), b in eachrow(xy)]

    
#     diff_v = []
#     for i in 1:Np
#         for j in 1:Np
#             diff_v[i,j] = [diff_x[i,j] diff_y[i,j]]
#         end
#     end
#     # diff_v = [a-b for a in eachrow(xy), b in eachrow(xy)]
#     diff_v_norm = diff_v./dists

#     replace!(diff_v_norm, NaN => 0.)

#     return diff_v_norm
# end
function radial_directions(xy::Array{Float64,2})
    Np = size(xy,1)
    dists = zeros(Np,Np)
    dists .= pairwise(Euclidean(),xy,xy,dims=1)

    diff_v = [a-b for a in eachrow(xy), b in eachrow(xy)]
    diff_v_norm = diff_v./dists

    replace!(diff_v_norm, [NaN, NaN] => [0.,0.])

    return diff_v_norm
end

function attractive_interactions!(xy::Array{Float64,2}, R::Float64)

    ϵ=0.1
    σ= 2*R
    Np = size(xy,1)
    dists = zeros(Np,Np)
    #superpose = falses(Np,Np)
    #uptriang = falses(Np,Np)
    # uptriang = trues(Np,Np)
    # for i = 1:Np
        
    #     #uptriang[i,i+1:Np] .= true
    #     uptriang[i,i] = false
    # end

    dists .= pairwise(Euclidean(),xy,xy,dims=1)
    dists13= dists.^(13)

    dists7= dists.^(7)

    force= 24*ϵ*(((2*σ^(12))./dists13).- (σ^(6)./dists7))
    replace!(force, NaN => 0.)

    dirs = radial_directions(xy)
    F_v = force.*dirs    
    # force= force.* uptriang
    ΣF = .-sum(F_v, dims = 1)

    return  reduce(vcat, transpose(ΣF))

end

xy = [[1. 1.]; [3. 2.]; [3. 3.]]

a = attractive_interactions!(xy,2.)
# attractive_interactions!(xy, 2.)

 #the new step function should look like

# function step(abpe::ABPE, δt::Float64) where {ABPE <: ABPsEnsemble}
    
#     γ = diffusion_coeff(abpe.R)[3]

#     if size(position(abpe),2) == 2
#         δp = sqrt.(2*δt*abpe.DT)*randn(abpe.Np,2) .+ abpe.v*δt*[cos.(abpe.θ) sin.(abpe.θ)] .+ δt*attractive_interactions!(position(abpe),abpe.R)'/γ
#         δθ = sqrt(2*abpe.DR*δt)*randn(abpe.Np)


#     else
#         println("No step method available")
#     end

#     return (δp, δθ)
# end

