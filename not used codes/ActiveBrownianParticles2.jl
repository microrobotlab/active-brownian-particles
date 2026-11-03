# Define an "ActiveBrownianParticle" Type
abstract type ActiveBrownianParticle end

# Define a specific type for 2D ABPs
struct ABP2 <: ActiveBrownianParticle
	x::Float64 		# x position (μm)
	y::Float64 		# y position (μm)
	θ::Float64  	# orientation (rad)
	R::Float64 		# Radius (μm)
	v::Float64 		# velocity (μm/s)
	DT::Float64 	# translational diffusion coefficient (μm^2/s)
	DR::Float64 	# rotational diffusion coefficient (rad^2/s)
end

# Define an "ABPsEnsemble" Type
abstract type ABPsEnsemble end

# Define a specific type for 2D ABPsEnsembles (CURRENTLY ASSUMING ALL PARTICLES ARE EQUAL)
struct ABPE2 <: ABPsEnsemble
    Np::Int64                      # number of particles
    L::Float64                      # size of observation space (μm)
	R::Float64  # Radius (μm)                                   --> Vector{Float64}(undef,Np)
	v::Float64 	# velocity (μm/s)                               --> Vector{Float64}(undef,Np)
	DT::Float64 # translational diffusion coefficient (μm^2/s)  --> Vector{Float64}(undef,Np)
	DR::Float64 # rotational diffusion coefficient (rad^2/s)    --> Vector{Float64}(undef,Np)
	x::Vector{Float64}    # x position (μm)
	y::Vector{Float64}    # y position (μm)
	θ::Vector{Float64}    # orientation (rad)
end

# Type methods
# ---

## Get position and orientation of the particle or ensemble (CURRENTLY ONLY 2D)
position(abp::ABP2) = ( abp.x, abp.y, abp.θ )

position(abpe::ABPE2) = [ abpe.x abpe.y ]
orientation(abpe::ABPE2) = abpe.θ

## Initialize ABP particle (CURRENTLY ONLY 2D)
function initABP(position::NTuple, R::Float64, v::Float64; T::Float64=300.0, η::Float64=1e-3)
    # translational diffusion coefficient [m^2/s] & rotational diffusion coefficient [rad^2/s] - R [m]
    DT, DR = diffusion_coeff(1e-6R)

    if length(position) == 3
        abp = ABP2( position..., R, v, 1e12DT, DR)
    else
        println("No init method available")
    end

    return abp
end

## Initialize ABP ensemble (CURRENTLY ONLY 2D)
function initABPE(Np::Int64, L::Float64, R::Float64, v::Float64; T::Float64=300.0, η::Float64=1e-3)
    # translational diffusion coefficient [m^2/s] & rotational diffusion coefficient [rad^2/s] - R [m]
    DT, DR = diffusion_coeff(1e-6R)

    # ONLY 2D!
    xyθ = (rand(Np,3).-0.5).*repeat([L L 2π],Np)
    xyθ[:,1:2], dists, superpose, uptriang = hardsphere(xyθ[:,1:2],R)
    abpe = ABPE2( Np, L, R, v, 1e12DT, DR, xyθ[:,1], xyθ[:,2], xyθ[:,3])

    return abpe, (dists, superpose, uptriang)
end

## Calculate the step position (and orientation) for the particle (CURRENTLY ONLY 2D)
function step(abp::ABP, δt::Float64) where {ABP <: ActiveBrownianParticle}
    if length(position(abp)) == 3
        δx = sqrt(2*abp.DT*δt)*randn() + abp.v*δt*cos(abp.θ)
        δy = sqrt(2*abp.DT*δt)*randn() + abp.v*δt*sin(abp.θ)
        δθ = sqrt(2*abp.DR*δt)*randn()
        δp = ( δx, δy, δθ )
    else
        println("No step method available")
    end
    return δp
end

## Update position and orientation of the particle (create new)
update(abp::ABP, step) where {ABP <: ActiveBrownianParticle} = ABP( (position(abp) .+ step)..., abp.R, abp.v, abp.DT, abp.DR )

# Functions
# ---

function diffusion_coeff(R::Float64, T::Float64=300.0, η::Float64=1e-3)
    # Boltzmann constant [J/K]
    kB = 1.38e-23
    # friction coefficient [Ns/m]
    γ = 6*pi*R*η
    # translational diffusion coefficient [m^2/s]
    DT = kB*T/γ
    # rotational diffusion coefficient [rad^2/s]
    DR = 6*DT/(8*R^2)
    return DT, DR
end;

## Create single particle trajectory
function trajectory(abp::ABP, N; δt::Float64=1e-3) where {ABP <: ActiveBrownianParticle}
	p = [position(abp)]
	t = 0:δt:δt*N
    
	for i in 1:N
		abp = update(abp, step(abp,δt))
		push!(p, position(abp))
	end
	return p, t
end

## Plot trajectory
function plot_trajectory(x,y)
    plot(x, y, c=:black, legend=false, aspect_ratio=:equal);
	scatter!([x[end]],[y[end]], legend=false);
    xlabel!("x [μm]");
    ylabel!("y [μm]");
end

## Animate single particle trajectory
## ---
##  works with:
# begin
# 	ticks_per_sec =10
# 	Δt = 1/ticks_per_sec
# 	tidx = findall(t .≈ round.(t/Δt).*Δt)
# 	time_vec = []
# end;
# @bind ticks Clock(Δt,true)
# begin
#     ticks
#     animate_trajectory(p,t,time_vec,Δt,tidx)
# end
## ---
function animate_trajectory(p, t, time_vec, tidx)
	push!(time_vec, time())
	if length(time_vec)>1
		framerate = round(1/(time_vec[end] - time_vec[end-1]),digits=1)
	else
		framerate = 0.0
	end
	nframe = length(time_vec)
	
	ti = tidx[1:nframe]
	x = first.(p)[ti]
	y = [ pi[2] for pi in p[ti] ]
	
	plot_trajectory(x, y)
	title!("t = $(round(t[ti[end]],digits=1)) s, $framerate fps")
end

function hardsphere_correction!(xy::Array{Float64,2}, dists::Array{Float64,2}, superpose::BitArray{2}, R::Float64; tol::Float64=1e-3)
    Np = size(superpose,1)
    for np1 in 1:Np
        if any(superpose[np1,:])
            np2 = findfirst(superpose[np1,:])
            Δp = (xy[np1,:] - xy[np2,:]) .* ( ( (1+tol)*2R / dists[np1,np2] - 1 ) / 2 )
            xy[np1,:] += Δp
            xy[np2,:] -= Δp
            dists[np2,np2+1:Np] = pairwise(Euclidean(), xy[np2:np2,:], xy[np2+1:Np,:], dims=1 )
            superpose[np2,np2+1:Np] = (dists[np2,np2+1:Np] .< 2R*(1-tol))
        end
    end
    return nothing
end

function hardsphere!(xy::Array{Float64,2}, dists::Array{Float64,2}, superpose::BitArray{2}, uptriang::BitArray{2}, R::Float64; tol::Float64=1e-3)
    superpositions = 1
    counter = 0
    # @time begin
    while superpositions > 0
        dists .= pairwise(Euclidean(),xy,xy,dims=1)
        superpose .= (dists .< 2R*(1-tol)).*uptriang
        # @show(findall(superpose))
        superpositions = sum(superpose)
        # @show(superpositions)
        if superpositions > 0
            hardsphere_correction!(xy,dists,superpose,R,tol=tol)
        end
        counter += 1
        # @show(counter)
        if counter >= 100
            println("$superpositions superpositions remaining after 100 cycles")
            break
        end
    end
    # end
    return nothing
end

function hardsphere(xy::Array{Float64,2}, R::Float64; tol::Float64=1e-3)
    Np = size(xy,1)
    dists = zeros(Np,Np)
    superpose = falses(Np,Np)
    uptriang = falses(Np,Np)
    for i = 1:Np-1
        uptriang[i,i+1:Np] .= true
    end
    hardsphere!(xy, dists, superpose, uptriang, R; tol=tol)
    return xy, dists, superpose, uptriang
end

# function periodic_BC_abp2(abp::ABP2,L::Float64)
#     xt = position(abp)
#     step = [0.0, 0.0, 0.0]
#     # Boundary conditions: horizontal edge
#     if abs(xt[1])>L/2
#         step[1] = -sign(xt[1])*L
#     end
#     # Boundary conditions: horizontal edge
#     if abs(xt[2])>L/2
#         step[2] = -sign(xt[2])*L
#     end
#     update(abp,Tuple(step))
#     return abp
# end