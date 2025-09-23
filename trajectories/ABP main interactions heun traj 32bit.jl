using .Threads
using CalculusWithJulia, Dates, Distributions, ForwardDiff, GeometryBasics, Random, Statistics, VoronoiCells
include(joinpath("..","sim_32bit","geo_toolbox_32bit.jl"))
include(joinpath("..", "force_functions.jl"))


# Define an "ABPsEnsemble" Type
abstract type ABPsEnsemble end

# Define a specific type for 2D ABPsEnsembles (CURRENTLY ASSUMING ALL PARTICLES ARE EQUAL)
struct ABPE2 <: ABPsEnsemble
    Np::Int64                      # number of particles®®
    L::Float32                      # size of observation space (μm)
	R::Float32  # Radius (μm)                                   --> Vector{Float32}(undef,Np)
    T::Float32 #temperature (K)
	v::Vector{Float32}  	# velocity (μm/s)                   --> Vector{Float32}(undef,Np)
    ω::Vector{Float32} #Angular velocity (rad/s)                --> Vector{Float32}(undef,Np)            
	DT::Float32 # translational diffusion coefficient (μm^2/s)  --> Vector{Float32}(undef,Np)
	DR::Float32 # rotational diffusion coefficient (rad^2/s)    --> Vector{Float32}(undef,Np)
	x::Vector{Float32}    # x position (μm)
	y::Vector{Float32}    # y position (μm)
	θ::Vector{Float32}    # orientation (rad)
end

#------------------------------------------------------------For square ---------------------------------------------------------------------------------------------------------

## Initialize ABP ensemble (CURRENTLY ONLY 2D) 
function initABPE(Np::Int64, L::Real, R::Real, T::Real, vd::Union{Real,Array{Float32,1},Distribution}, ωd::Union{Real,Array{Float32,1},Distribution}, int_func::Function, offcenter::Float32, int_params...; η::Float32=1f-3, xyθ = Float32[])
    # translational diffusion coefficient [m^2/s] & rotational diffusion coefficient [rad^2/s] - R [m]
    # Intial condition will be choosen as per the geometry under study
    (vd isa Float32) ? vd = [vd] : Nt
    (ωd isa Float32) ? ωd = [ωd] : Nt
    DT, DR = diffusion_coeff(1f-6R, T)
    if isempty(xyθ)
        xyθ = (rand(Float32,(Np,3)).-5.f-1).*repeat([L L Float32(2π)],Np)
    end
    xyθ[:,1:2], dists, superpose, uptriang = hardsphere_periodic(xyθ[:,1:2], R, L) #xyθ[:,1:2] gives x and y positions of intitial particles
    v = Float32.(rand(vd, Np))
    ω = Float32.(rand(ωd, Np))
    abpe = ABPE2( Np, L, R, T, v, ω, 1f12DT, DR, xyθ[:,1], xyθ[:,2], xyθ[:,3])
    return abpe, (dists, superpose, uptriang)
end

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Calculate diffusion coefficient and friction coefficient
function diffusion_coeff(R::Float32, T::Float32=3f2, η::Float32=1f-3)
    # Boltzmann constant [J/K]
    kB = 1.38f-23
    # friction coefficient [Ns/m]
    γ = 6*(pi*1f0)*R*η
    # translational diffusion coefficient [m^2/s]
    DT = kB*T/γ
    # rotational diffusion coefficient [rad^2/s]
    DR = 6*DT/(8*R^2)
    return DT, DR, γ
end
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
position(abpe::ABPE2) = [ abpe.x abpe.y ]
orientation(abpe::ABPE2) = abpe.θ
force(abpe::ABPE2) = [abpe.fx abpe.fy]
torque(abpe::ABPE2) = abpe.torque

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions to update particles for the next step

function update_heun(abpe::ABPE, matrices::Tuple{Matrix{Float32}, BitMatrix, BitMatrix}, δt::Float32, offcenter::Float32, range::Float32, int_func::Function, int_params...) where {ABPE <: ABPsEnsemble}

    δp_i = Array{Float32,2}(undef,abpe.Np,2)
    δθ_i = Array{Float32,1}(undef,abpe.Np)
    δp_f = Array{Float32,2}(undef,abpe.Np,2)
    δθ_f = Array{Float32,1}(undef,abpe.Np)
    γₜ = diffusion_coeff(1f-6*abpe.R, abpe.T)[3] #Output in international system units kg/s
    γᵣ = (8f-12γₜ / 6) * abpe.R^2                #Output in international system units
    oc_length = offcenter*abpe.R
    
    #intermediate step
    if (!isapprox(offcenter,0.f0))
        f_i, t_i = force_torque(position(abpe), orientation(abpe), abpe.L, oc_length, range, int_func, int_params...)
    else
        f_i = interactions_range(position(abpe), abpe.L, abpe.Np, int_func, int_params...)
        t_i = zeros(abpe.Np)
    end 

    det_part_i = (abpe.v.*[cos.(abpe.θ) sin.(abpe.θ)] .+ 1f-6f_i/γₜ, abpe.ω .+ 1f-18t_i/γᵣ)
    noise_part =  (sqrt.(2*δt*abpe.DT)*randn(abpe.Np,2), sqrt(2*abpe.DR*δt)*randn(abpe.Np))
    δp_i .= noise_part[1] .+ δt.*det_part_i[1]
    δθ_i .= noise_part[2] .+ δt.*det_part_i[2]

    pθ_i = (position(abpe), orientation(abpe)) .+ (δp_i, δθ_i)

    periodic_BC_array!(pθ_i[1], abpe.L, abpe.Np)

    #final step
    if (!isapprox(offcenter,0.f0))
        f_f, t_f = force_torque(pθ_i..., abpe.L, oc_length, range, int_func, int_params...)
    else
        f_f = interactions_range(pθ_i[1], abpe.L, abpe.Np, int_func, int_params...)
        t_f = zeros(abpe.Np)
    end

    det_part_f = (abpe.v.*[cos.(pθ_i[2]) sin.(pθ_i[2])] .+ 1f-6f_f/γₜ, abpe.ω .+ 1f-18t_f/γᵣ)
    δp_f .= noise_part[1] .+ δt.*(det_part_f[1] .+ det_part_i[1])/2
    δθ_f .= noise_part[2] .+ δt.*(det_part_f[2] + det_part_i[2])/2

    pθ = (position(abpe), orientation(abpe)) .+ (δp_f, δθ_f)

    periodic_BC_array!(pθ[1],abpe.L, abpe.Np)
    # println(maximum(abs.(pθ[1])))
    hardsphere_periodic!(pθ[1], matrices[1], matrices[2],abpe.R, abpe.L)
    periodic_BC_array!(pθ[1],abpe.L, abpe.Np) # Apply a second time for good measure

    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.T, abpe.v, abpe.ω, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2] )

    return new_abpe
end

function offcenter_nosuperpose!(abpe::ABPE2, δt::Float32, offcenter::Float32, int_func::Function, int_params...)
    xy = position(abpe)
    θ = orientation(abpe)
    R = abpe.R
    L = abpe.L
    γₜ = diffusion_coeff(1f-6*R, T)[3] #Output in international system units kg/s
    γᵣ = (8f-12γₜ / 6) * R^2                #Output in international system units
    oc = offcenter*R

    xy_chgcen = xy .+ oc* [cos.(θ) sin.(θ)]

    ocdists = pairwise(Euclidean(), xy_chgcen, xy_chgcen, dims=1)
    superpositions = sum(0.f0 .< ocdists .< 2R)/2
    j = 0
    while superpositions > 0
        # sup_per_particle = sum(ocdists .< 2R, dims=1)
        # println(sup_per_particle)
        # break
        # if j<=100
        #     t = force_torque_voronoi(xy_chgcen, θ, R, L, forward, offcenter, range, int_func, int_params...)[2]
        #     δθ = (δt)*1e-18t/γᵣ
        #     abpe.θ .+= δθ
        #     xy_chgcen = xy .+ (2*forward-1) .* [cos.(θ) sin.(θ)] .* R*offcenter
        #     superpositions = sum(0.0 .< pairwise(Euclidean(), xy_chgcen, xy_chgcen, dims=1) .< 2R)/2
        #     j +=1
        #     else
        ocdists = pairwise(Euclidean(), xy_chgcen, xy_chgcen, dims=1)
        sup_per_particle = sum(0.0 .< ocdists .< 2R, dims=1)
        id = vec(sup_per_particle .> 0)
        abpe.θ[id] .+= rand(length(abpe.θ[id.>0])).*2π
        xy_chgcen = xy .+ oc* [cos.(θ) sin.(θ)]
        superpositions = sum(0.0 .< pairwise(Euclidean(), xy_chgcen, xy_chgcen, dims=1) .< 2R)/2
        j +=1
        # end
    end
    # println(j)
    return nothing
end
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions for the hard sphere corrections

function hardsphere_periodic!(xy::Array{Float32,2}, periodicdists::Array{Float32,2}, superpose::BitArray{2}, R::Float32,L::Float32; tol::Float32=0.f0)
    superpositions = 1
    counter = 0
    Δxy = zeros(Float32, size(periodicdists)...,2)
    while superpositions > 0
        Threads.@threads for i in 1:2
            Δxy[:,:,i] .= pairwise(-,xy[:,i],xy[:,i])
            ix = abs.(Δxy[:,:,i]) .> abs.(abs.(Δxy[:,:,i]).-L)
            Δxy[ix,i] .= -sign.(Δxy[ix,i]).*abs.(abs.(Δxy[ix,i]).-L)
        end
        periodicdists .= sqrt.(Δxy[:,:,1].^2 .+ Δxy[:,:,2].^2)
        superpose .= 0.f0 .< periodicdists.<2R
        superpositions = sum(superpose)÷2
        # @show(superpositions)
        if superpositions > 0
            # println("correcting...")
            hardsphere_correction_periodic!(xy,Δxy,periodicdists,superpose,R,tol=1f-3)
            periodic_BC_array!(xy,L,size(xy,1))
        end
        counter += 1
        # @show(counter)
        if counter >= 100
            println("$superpositions superpositions remaining after 100 cycles")
            break
        end
    end
    return nothing
end

function hardsphere_correction_periodic!(xy::Array{Float32,2}, Δxy::Array{Float32,3}, dists::Array{Float32,2}, superpose::BitArray{2}, R::Float32; tol::Float32=1f-3)
    # Np = size(superpose,1)
    # Δxy = zeros(size(dists)...,2)
    Δpi(Δxi,d,s) = s==0 ? 0.f0 : Δxi*(((1+tol)*2R / d - 1) / 2 )
    Threads.@threads for i in 1:2
        # Δxy[:,:,i] .= pairwise(-,xy[:,i],xy[:,i])
        xy[:,i] .+= sum(Δpi.(Δxy[:,:,i],dists,superpose),dims=2)
    end
    return nothing
end

function hardsphere_periodic(xy::Array{Float32,2}, R::Float32, L::Float32; tol::Float32=1f-3) # called in initABPE
    Np = size(xy,1)
    periodicdists = zeros(Float32, Np,Np)
    superpose = falses(Np,Np)
    uptriang = triu(trues(Np,Np),1)
    hardsphere_periodic!(xy, periodicdists, superpose, R,L; tol=tol)
    return xy, periodicdists, superpose, uptriang
end


function hs_displacement(adjacency::Set{Tuple{Int64, Int64}}, xy::Array{Float32, 2}, Np::Int, R::Real; tol::Float32= 1f-8)
    displacement = zeros(Float32, Np, 2)
    for a in adjacency
        # Since tuples are ordered and the set is composed starting from points inside the simulation box, we are sure that a[1] <= Np
        dxy = xy[a[1],:] - xy[a[2],:]
        dist = sqrt(sum(abs2, dxy))
        if dist < 2R
            dir = dxy / dist
            displacement[a[1],:] += 0.5(1+tol)*(2R - dist) * dir
            if a[2] <= Np
                displacement[a[2],:] -= 0.5(1+tol)*(2R - dist) * dir
            end
        end
    end
    return displacement
end

function hs_voronoi!(xy::Array{Float32, 2}, L::Real, Np::Int, R::Real)
    for _ in 1:100
        #Step 1: Get the tessellation
        rect = Rectangle(Point2(-L/2, -L/2), Point2(L/2, L/2))

        tess = voronoicells(xy_to_points(xy), rect)

        # Step 2: Get the periodic projection of the tessellation

        xy_periodic_projection, proj_inds = find_boundary_points(xy, tess.Cells)
        xy_periodic_projection = vcat(xy, xy_periodic_projection)

        # Get the rectangle of the periodic projection
        minpoint = Point2(minimum(xy_periodic_projection[:,1]), minimum(xy_periodic_projection[:,2]))
        maxpoint = Point2(maximum(xy_periodic_projection[:,1]), maximum(xy_periodic_projection[:,2]))
        rect_periodic = Rectangle(minpoint, maxpoint)
        tess_periodic = voronoicells(xy_to_points(xy_periodic_projection), rect_periodic)

        # Step 3: Get the adjacency from the periodic tessellation
        adjacency = get_adjacency_from_points(tess_periodic.Cells, Np)

        #Step 4: Displace the particles according to the adjacency
        dis = hs_displacement(adjacency, xy_periodic_projection, Np, R)
        dis == zeros(Float32, Np, 2) && break
        xy.+= dis
        periodic_BC_array!(xy, L, Np)
    end
end

function hs_voronoi(xy::Array{Float32, 2}, L::Real, Np::Int, R::Real) #Called in initABPE
    # This function is a wrapper for hs_voronoi! to return the modified xy array
    hs_voronoi!(xy, L, Np, R)
    return xy
end

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function periodic_BC_array!(xy::AbstractArray{Float32,2}, L::Real, Np::Int)
    for i in 1:Np
        x, y = xy[i, 1], xy[i, 2]

        if abs(x) > L/2
            xy[i, 1] = x - sign(x) * L
        end

        if abs(y) > L/2
            xy[i, 2] = y - sign(y) * L
        end
    end
    return nothing
end

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Function to calculate force vectors

function interactions_range(xy::Array{Float64, 2}, L::Float64, l::Float64, Np::Int, int_func::Function, int_params...)
    # Preallocate the result array for efficiency
    ΣFtot = zeros(Float64, Np, 2)
    
    @threads for i in axes(xy,1)
        # Compute distances and apply periodic boundary conditions
        xy_shifted = xy .- [xy[i,1] xy[i,2]]
        periodic_BC_array!(xy_shifted, L, Np)
        xy_nonzero = @view xy_shifted[1:end .!=i, :]

        dists = d2(xy_nonzero)
        inside = vec(dists .<= l)
        # If no particles are inside the interaction range, force remains zero
        !any(inside) && continue
        xy_inside = @view xy_nonzero[inside, :]
        dirs = zeros(Float64, sum(inside), 2)

        # Compute distances and interaction forces
        dists_nonzero = @view dists[inside]
        forces = int_func.(dists_nonzero, int_params...)
        
        # Compute normalized direction vectors
        dirs .= xy_inside ./ dists_nonzero

        # Sum forces for each direction and assign to ΣFtot
        ΣFtot[i, :] = -sum(forces .* dirs, dims=1)'
    end
    return ΣFtot
end

#Function used to compute torques in aligning interactions
function force_torque(xy::Array{Float64,2}, θ::Array{Float64,1}, L::Real, oc::Float64, range::Real, int_func::Function, int_params...) #Forces are retuned in μN, torques in μN×μm
    xy_chgcen = xy .+ oc* [cos.(θ) sin.(θ)]
    forces = interactions_range(xy_chgcen, L, range, size(xy,1), int_func, int_params...)
    torques = oc*(forces[:,2] .* cos.(θ) .- forces[:,1] .* sin.(θ))
    return forces, torques
end
