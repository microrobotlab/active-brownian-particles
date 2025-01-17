using .Threads
using CalculusWithJulia, Dates, Distributions, ForwardDiff, Random, Statistics
include("geo_toolbox.jl")

# Define an "ABPsEnsemble" Type
abstract type ABPsEnsemble end

# Define a specific type for 2D ABPsEnsembles (CURRENTLY ASSUMING ALL PARTICLES ARE EQUAL)
struct ABPE2 <: ABPsEnsemble
    Np::Int64                      # number of particles®®
    L::Float64                      # size of observation space (μm)
	R::Float64  # Radius (μm)                                   --> Vector{Float64}(undef,Np)
    T::Float64 #temperature (K)
	v::Vector{Float64}  	# velocity (μm/s)                   --> Vector{Float64}(undef,Np)
    ω::Vector{Float64} #Angular velocity (rad/s)                --> Vector{Float64}(undef,Np)            
	DT::Float64 # translational diffusion coefficient (μm^2/s)  --> Vector{Float64}(undef,Np)
	DR::Float64 # rotational diffusion coefficient (rad^2/s)    --> Vector{Float64}(undef,Np)
	x::Vector{Float64}    # x position (μm)
	y::Vector{Float64}    # y position (μm)
	θ::Vector{Float64}    # orientation (rad)
end

#------------------------------------------------------------For square ---------------------------------------------------------------------------------------------------------

## Initialize ABP ensemble (CURRENTLY ONLY 2D) 
function initABPE(Np::Int64, L::Float64, R::Float64, T::Float64, vd::Union{Float64,Array{Float64,1},Distribution}, ωd::Union{Float64,Array{Float64,1},Distribution}, int_func::Function, forward::Bool, offcenter::Float64, range::Float64, int_params...; η::Float64=1e-3)
    # translational diffusion coefficient [m^2/s] & rotational diffusion coefficient [rad^2/s] - R [m]
    # Intial condition will be choosen as per the geometry under study
    (vd isa Float64) ? vd = [vd] : Nt
    (ωd isa Float64) ? ωd = [ωd] : Nt
    DT, DR = diffusion_coeff(1e-6R, T)
    xyθ = (rand(Np,3).-0.5).*repeat([L L 2π],Np)
    xyθ[:,1:2], dists, superpose, uptriang = hardsphere_periodic(xyθ[:,1:2], R, L) #xyθ[:,1:2] gives x and y positions of intitial particles
    v = rand(vd, Np)
    ω = rand(ωd,Np)
    abpe = ABPE2( Np, L, R, T, v, ω, 1e12DT, DR, xyθ[:,1], xyθ[:,2], xyθ[:,3])
    return abpe, (dists, superpose, uptriang)
end

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Calculate diffusion coefficient and friction coefficient
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
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions to simulate multiple spherical particles
function multiparticleE(Np::Integer, L::Float64, R::Float64, T::Float64,  v::Union{Float64,Array{Float64,1},Distribution}, ω::Union{Float64,Array{Float64,1},Distribution}, Nt::Int64, measevery::Int64, δt::Float64, int_func::Function, forward::Bool, offcenter::Float64, range::Float64, int_params...;)
    (Nt isa Int64) ? Nt : Nt=convert(Int64,Nt)

    ABPE_history = Vector{ABPE2}(undef,Nt÷(measevery)+1) # Nt is number of time steps

    ABPE, matrices = initABPE( Np, L, R, T, v, ω, int_func, forward, offcenter, range, int_params...,) # including initial hardsphere correction
    ABPE_history[1] = ABPE
    simulate!(ABPE_history, ABPE, matrices, Nt, measevery, δt, forward, offcenter, range, int_func, int_params...)

    return position.(ABPE_history), orientation.(ABPE_history), force.(ABPE_history), torque.(ABPE_history)
end

function simulate!(ABPE_history, ABPE, matrices, Nt, measevery, δt, forward, offcenter, range, int_func, int_params...)
    start = now()
    print_step = Nt÷100
    for nt in 1:Nt
        ABPE = update(ABPE,matrices,δt, forward, offcenter, range, int_func, int_params...)#updating information at every step
        if (nt) % measevery == 0
            ABPE_history[(nt)÷measevery+1] = ABPE
        end
        if nt % print_step == 0
            elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
            print("$((100*nt÷Nt))%... Step $nt, total elapsed time $(elapsed)\r")
        end
    end
    print("\n")
    return nothing
end

position(abpe::ABPE2) = [ abpe.x abpe.y ]
orientation(abpe::ABPE2) = abpe.θ
force(abpe::ABPE2) = [abpe.fx abpe.fy]
torque(abpe::ABPE2) = abpe.torque

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions to update particles for the next step

function update_heun(abpe::ABPE, matrices::Tuple{Matrix{Float64}, BitMatrix, BitMatrix}, δt::Float64, forward::Bool, offcenter::Float64, range::Float64, int_func::Function, int_params...) where {ABPE <: ABPsEnsemble}

    δp_i = Array{Float64,2}(undef,abpe.Np,2)
    δθ_i = Array{Float64,1}(undef,abpe.Np)
    δp_f = Array{Float64,2}(undef,abpe.Np,2)
    δθ_f = Array{Float64,1}(undef,abpe.Np)
    γₜ = diffusion_coeff(1e-6*abpe.R, abpe.T)[3] #Output in international system units kg/s
    γᵣ = (8e-12γₜ / 6) * abpe.R^2                #Output in international system units

    wca_sigma = 2abpe.R
    wca_range = wca_sigma*2^(1/6)
    #intermediate step
    if (!isapprox(offcenter,0.0))
        f_i, t_i = force_torque(position(abpe), orientation(abpe), abpe.R, abpe.L, forward, offcenter, range, int_func, int_params...)
        # f_i .+= interactions_range(position(abpe), abpe.R, abpe.L, wca_range, abpe.Np, weeks_chandler_andersen, wca_sigma, 1.)
    else
        f_i = interactions_range(position(abpe), abpe.R, abpe.L, range, abpe.Np, int_func, int_params...)
        # f_i .+= interactions_range(position(abpe), abpe.R, abpe.L, wca_range, abpe.Np, weeks_chandler_andersen, wca_sigma, 1.)
        t_i = zeros(abpe.Np)
    end 

    det_part_i = (abpe.v.*[cos.(abpe.θ) sin.(abpe.θ)] .+ 1e-6f_i/γₜ, abpe.ω .+ 1e-18t_i/γᵣ)
    noise_part =  (sqrt.(2*δt*abpe.DT)*randn(abpe.Np,2), sqrt(2*abpe.DR*δt)*randn(abpe.Np))
    δp_i .= noise_part[1] .+ δt.*det_part_i[1]
    δθ_i .= noise_part[2] .+ δt.*det_part_i[2]

    pθ_i = (position(abpe), orientation(abpe)) .+ (δp_i, δθ_i)

    periodic_BC_array!(pθ_i[1], abpe.L, abpe.R)

    #final step
    if (!isapprox(offcenter,0.0))
        f_f, t_f = force_torque(pθ_i..., abpe.R, abpe.L, forward, offcenter, range, int_func, int_params...)
        # f_f .+= interactions_range(pθ_i[1], abpe.R, abpe.L, wca_range, abpe.Np, weeks_chandler_andersen, wca_sigma, 1.)
    else
        f_f = interactions_range(pθ_i[1], abpe.R, abpe.L, range, abpe.Np, int_func, int_params...)
        # f_f .+= interactions_range(pθ_i[1], abpe.R, abpe.L, wca_range, abpe.Np, weeks_chandler_andersen, wca_sigma, 1.)
        t_f = zeros(abpe.Np)
    end 

    det_part_f = (abpe.v.*[cos.(pθ_i[2]) sin.(pθ_i[2])] .+ 1e-6f_f/γₜ, abpe.ω .+ 1e-18t_f/γᵣ)
    δp_f .= noise_part[1] .+ δt.*(det_part_f[1] .+ det_part_i[1])/2
    δθ_f .= noise_part[2] .+ δt.*(det_part_f[2] + det_part_i[2])/2

    pθ = (position(abpe), orientation(abpe)) .+ (δp_f, δθ_f)

    periodic_BC_array!(pθ[1],abpe.L, abpe.R)
    hardsphere_periodic!(pθ[1], matrices[1], matrices[2],abpe.R, abpe.L)

    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.T, abpe.v, abpe.ω, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2] )

    return new_abpe
end

function force_torque(xy::Array{Float64,2}, θ::Array{Float64,1}, R::Float64, L::Float64, forward::Bool, offcenter::Float64, range::Float64, int_func::Function, int_params...) #Forces are retuned in μN, torques in μN×μm
    xy_chgcen = xy .+ (2*forward-1) .* [cos.(θ) sin.(θ)] .* R*offcenter
    forces = interactions_range(xy_chgcen, R, L, range, size(xy,1), int_func, int_params...)
    torques = offcenter * R .* (forces[:,2] .* cos.(θ) .- forces[:,1] .* sin.(θ))
    return forces, torques
end

function offcenter_nosuperpose!(abpe::ABPE2, δt::Float64, forward::Bool, offcenter::Float64, range::Float64, int_func::Function, int_params...)
    xy = position(abpe)
    θ = orientation(abpe)
    R = abpe.R
    L = abpe.L
    γₜ = diffusion_coeff(1e-6*R, T)[3] #Output in international system units kg/s
    γᵣ = (8e-12γₜ / 6) * R^2                #Output in international system units

    xy_chgcen = xy .+ (2*forward-1) .* [cos.(θ) sin.(θ)] .* R*offcenter

    ocdists = pairwise(Euclidean(), xy_chgcen, xy_chgcen, dims=1)
    superpositions = sum(0.0 .< ocdists .< 2R)/2
    j = 0
    while superpositions > 0
        # sup_per_particle = sum(ocdists .< 2R, dims=1)
        # println(sup_per_particle)
        # break
        # if j<=100
        #     t = force_torque(xy_chgcen, θ, R, L, forward, offcenter, range, int_func, int_params...)[2]
        #     δθ = (δt)*1e-18t/γᵣ
        #     abpe.θ .+= δθ
        #     xy_chgcen = xy .+ (2*forward-1) .* [cos.(θ) sin.(θ)] .* R*offcenter
        #     superpositions = sum(0.0 .< pairwise(Euclidean(), xy_chgcen, xy_chgcen, dims=1) .< 2R)/2
        #     j +=1
        #     else
        ocdists = pairwise(Euclidean(), xy_chgcen, xy_chgcen, dims=1)
        sup_per_particle = sum(0.0 .< ocdists .< 2R, dims=1)
        id = vec(sup_per_particle .> 0)
        abpe.θ[id] .+= rand(abpe.θ[id.>0]).*2π
        xy_chgcen = xy .+ (2*forward-1) .* [cos.(θ) sin.(θ)] .* R*offcenter
        superpositions = sum(0.0 .< pairwise(Euclidean(), xy_chgcen, xy_chgcen, dims=1) .< 2R)/2
        j +=1
        # end
    end
    return nothing
end
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions for the hard sphere corrections

function hardsphere_periodic!(xy::Array{Float64,2}, periodicdists::Array{Float64,2}, superpose::BitArray{2}, R::Float64,L::Float64; tol::Float64=0.)
    superpositions = 1
    counter = 0
    Δxy = zeros(size(periodicdists)...,2)
    while superpositions > 0
        Threads.@threads for i in 1:2
            Δxy[:,:,i] .= pairwise(-,xy[:,i],xy[:,i])
            ix = abs.(Δxy[:,:,i]) .> abs.(abs.(Δxy[:,:,i]).-L)
            Δxy[ix,i] .= -sign.(Δxy[ix,i]).*abs.(abs.(Δxy[ix,i]).-L)
        end
        periodicdists .= sqrt.(Δxy[:,:,1].^2 .+ Δxy[:,:,2].^2)
        superpose .= 0. .< periodicdists.<2R
        superpositions = sum(superpose)÷2
        # @show(superpositions)
        if superpositions > 0
            # println("correcting...")
            hardsphere_correction_periodic!(xy,Δxy,periodicdists,superpose,R,tol=1e-3)
            periodic_BC_array!(xy,L,R)
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

function hardsphere_correction_periodic!(xy::Array{Float64,2}, Δxy::Array{Float64,3}, dists::Array{Float64,2}, superpose::BitArray{2}, R::Float64; tol::Float64=1e-3)
    # Np = size(superpose,1)
    # Δxy = zeros(size(dists)...,2)
    Δpi(Δxi,d,s) = s==0 ? 0.0 : Δxi*(((1+tol)*2R / d - 1) / 2 )
    Threads.@threads for i in 1:2
        # Δxy[:,:,i] .= pairwise(-,xy[:,i],xy[:,i])
        xy[:,i] .+= sum(Δpi.(Δxy[:,:,i],dists,superpose),dims=2)
    end
    return nothing
end

function hardsphere_periodic(xy::Array{Float64,2}, R::Float64, L::Float64; tol::Float64=1e-3) # called in initABPE
    Np = size(xy,1)
    periodicdists = zeros(Np,Np)
    superpose = falses(Np,Np)
    uptriang = triu(trues(Np,Np),1)
    hardsphere_periodic!(xy, periodicdists, superpose, R,L; tol=tol)
    return xy, periodicdists, superpose, uptriang
end

# Building
function hardsphere2_periodic!(xy::Vector{Tuple{Float64,Float64}}, R::Float64; tol::Float64=1e-3)
    superpositions = 1
    counter = 0
    periodicdists = zeros(length(xy),length(xy))
    Δxy = Matrix{Tuple{Float64,Float64}}(undef,size(periodicdists)...)
    Δp = Matrix{Tuple{Float64,Float64}}(undef,size(periodicdists)...)
    superpose = falses(size(periodicdists)...)

    displace(Δxy::Tuple{Float64,Float64},d::Float64) = Δxy.*(((1+tol)*2R / d - 1) / 2)
    correct(xy::Tuple{Float64,Float64},Δp_row) = xy.+reduce((x,y)->x.+y,Δp_row)
    periodic(d,L) = abs(d) > L-abs(d) ? -sign(d)*(L-abs(d)) : d

    while superpositions > 0
        pairwise!(.-,Δxy,xy,xy)
        Δxy .=  map(x->periodic.(x,L), Δxy)
        periodicdists = map(x->sqrt(x[1]^2+x[2]^2),Δxy)       
        superpose .= (0. .< periodicdists .< 2R*(1-tol))
        superpositions = sum(superpose)÷2
        if superpositions > 0
            fill!(Δp,(0.,0.))
            # Threads.@threads for s in findall(superpose)
            #     Δp[s] = displace(xy[s[1]].-xy[s[2]],dists[s])
            # end
            Δp[superpose] .= displace.(Δxy[superpose],periodicdists[superpose])
            xy = correct.(xy,eachrow(Δp)) 
        end
        counter += 1
        if counter >= 100
            println("$superpositions superpositions remaining after 100 cycles")
            break
        end
    end
    return nothing
end
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Function to calculate force vectors
function interactions(xy::Array{Float64,2}, R::Float64)
    ϵ=.1
    σ= 2R

    dists = pairwise(Euclidean(),xy,dims=1)
    
    strength_param = 1e0
    force = strength_param.*lennard_jones.(dists, σ, ϵ)
    replace!(force, NaN => 0.)

    dirs = radial_directions(xy)
    F_x = force.*dirs[1]
    F_y = force.*dirs[2]
    ΣFx = sum(F_x, dims = 1)
    ΣFy = sum(F_y, dims = 1)
    ΣF = vcat.(ΣFx, ΣFy)    
    return  reduce(vcat, transpose(ΣF))
end

function interactions_range(xy::Array{Float64, 2}, R::Float64, L::Float64, l::Float64, Np::Int, int_func::Function, int_params...)
    # Preallocate the result array for efficiency
    ΣFtot = Array{Float64}(undef, Np, 2)

    # Threshold for selecting particles within interaction range
    threshold = l
    
    @threads for i in axes(xy,1)
        # Compute distances and apply periodic boundary conditions
        xy_shifted = xy .- [xy[i,1] xy[i,2]]
        periodic_BC_array!(xy_shifted, L, R)

        # Filter to get particles within interaction range and not at the origin
        inside = all(reshape(abs.(xy_shifted) .<= threshold, Np, 2), dims=2)[:,1]
        xy_inside = xy_shifted[inside, :]

        # Remove central particle [0, 0] from interaction set
        xy_inside = xy_inside[xy_inside[:,1] .!= 0 .|| xy_inside[:,2] .!= 0, :]

        dirs = Array{Float64}(undef, size(xy_inside))

        if isempty(xy_inside)
            ΣFtot[i, :] .= 0.0
            continue
        end

        # Compute distances and interaction forces
        dists = d2(xy_inside)
        dists_nonzero = dists[dists .!= 0]
        forces = int_func.(dists_nonzero, int_params...)
        
        # Compute normalized direction vectors
        dirs .= xy_inside ./ dists_nonzero

        # Sum forces for each direction and assign to ΣFtot
        ΣFtot[i, :] .= .-sum(forces .* dirs, dims=1)'
    end

    return ΣFtot
end


lennard_jones(x, σ, ϵ) = 24*ϵ*(((2*σ^(12))/(x^(13)))- (σ^(6)/(x^(7))))
shifted_lennard_jones(x, σ, ϵ, shift) = 24*ϵ*(((2*σ^(12))/((x-shift)^(13)))- (σ^(6)/((x-shift)^(7))))


function weeks_chandler_andersen(x, σ::Float64, ϵ::Float64)
    rmin = σ*2^(1/6)
    if x > rmin
        return 0
    else
        return 24*ϵ*(((2*σ^(12))/(x^(13)))- (σ^(6)/(x^(7))))
    end
end

function purely_attractive_lennard_jones(x, σ::Float64, ϵ::Float64)
    rmin = σ*2^(1/6) + σ/2
    if x > rmin
        return 24*ϵ*(((2*σ^(12))/((x + σ/2)^(13)))- (σ^(6)/((x + σ/2)^(7))))
    else
        return 0
    end
end

function rlj_boundary(x::Float64, σ::Float64, ϵ::Float64)
    rmin = σ*2^(1/6) + σ/2
    if x > rmin
        return 0

    elseif x < σ
        return 390144*ϵ/σ
    else
        return 24*ϵ*(((2*σ^(12))/((x-σ/2)^(13)))- (σ^(6)/((x-σ/2)^(7)))) 
    end
end

function rlj_boundary(x::Float64, σ::Float64, ϵ::Float64)
    rmin = σ*2^(1/6) + σ/2
    if x > rmin
        return 0

    elseif x < σ
        return 390144*ϵ/σ
    else
        return 24*ϵ*(((2*σ^(12))/((x-σ/2)^(13)))- (σ^(6)/((x-σ/2)^(7)))) 
    end
end

function contact_lj(x, σ::Float64, ϵ::Float64)
    shift = 2^(1/6)-1
    return 24*ϵ*(((2*σ^(12))/((x+shift)^(13)))- (σ^(6)/((x+shift)^(7))))
end

coulomb(x,k) = k / x^2

excluded_volume(x::Float64,R::Float64, ϵ::Float64) = 12ϵ*((2R)^12)x^(-13)

spring(x,k, x0) = -k*(x-x0)
