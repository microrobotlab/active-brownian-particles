using CalculusWithJulia, Dates, Distributions, ForwardDiff, Random, Statistics
include("geo_toolbox.jl")

# Define an "ABPsEnsemble" Type
abstract type ABPsEnsemble end

# Define a specific type for 2D ABPsEnsembles (CURRENTLY ASSUMING ALL PARTICLES ARE EQUAL)
struct ABPE2 <: ABPsEnsemble
    Np::Int64                      # number of particles®®
    L::Float64                      # size of observation space (μm)
	R::Float64  # Radius (μm)                                   --> Vector{Float64}(undef,Np)
	v::Vector{Float64}  	# velocity (μm/s)                   --> Vector{Float64}(undef,Np)
    ω::Vector{Float64} #Angular velocity (rad/s)                --> Vector{Float64}(undef,Np)            
	DT::Float64 # translational diffusion coefficient (μm^2/s)  --> Vector{Float64}(undef,Np)
	DR::Float64 # rotational diffusion coefficient (rad^2/s)    --> Vector{Float64}(undef,Np)
	x::Vector{Float64}    # x position (μm)
	y::Vector{Float64}    # y position (μm)
	θ::Vector{Float64}    # orientation (rad)
    fx::Vector{Float64}    # force x
    fy::Vector{Float64}    # force y
    torque::Vector{Float64}    # torque
end

#------------------------------------------------------------For square ---------------------------------------------------------------------------------------------------------

## Initialize ABP ensemble (CURRENTLY ONLY 2D) 
function initABPE(Np::Int64, L::Float64, R::Float64, vd::Union{Float64,Array{Float64,1},Distribution}, ωd::Union{Float64,Array{Float64,1},Distribution}, int_func::Function, forward::Bool, offcenter::Float64, range::Float64, int_params...; T::Float64=300.0, η::Float64=1e-3)
    # translational diffusion coefficient [m^2/s] & rotational diffusion coefficient [rad^2/s] - R [m]
    # Intial condition will be choosen as per the geometry under study
    (vd isa Float64) ? vd = [vd] : Nt
    (ωd isa Float64) ? ωd = [ωd] : Nt
    DT, DR = diffusion_coeff(1e-6R)
    xyθ = (rand(Np,3).-0.5).*repeat([L L 2π],Np)
    xyθ[:,1:2], dists, superpose, uptriang = hardsphere(xyθ[:,1:2],R) #xyθ[:,1:2] gives x and y positions of intitial particles
    v = rand(vd, Np)
    ω = rand(ωd,Np)
    force, torque = force_torque(xyθ[:,1:2], xyθ[:,3], R, L, forward, offcenter, range, int_func, int_params...)
    fx = force[:,1]
    fy = force[:,2]
    abpe = ABPE2( Np, L, R, v, ω, 1e12DT, DR, xyθ[:,1], xyθ[:,2], xyθ[:,3], fx, fy, torque)
    return abpe, (dists, superpose, uptriang)
end

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Generating inside ellipse

# function initABPE(Np::Int64, L::Float64, R::Float64, v::Float64; T::Float64=300.0, η::Float64=1e-3)
#     # translational diffusion coefficient [m^2/s] & rotational diffusion coefficient [rad^2/s] - R [m]
#     # Intial condition will be choosen as per the geometry under study
#     DT, DR = diffusion_coeff(1e-6R)


#     α = rand(Np).*2π

#     rₑ = (a-R)*(b-R)./(sqrt.((((a-R)*sin.(α)).^2) .+ ((b-R)*cos.((α))).^2))  # r value for boundary

#     r = rand(Np).*rₑ
#     θ = rand(Np).*2π

#     xyθ = [r.*cos.(α) r.*sin.(α) θ]

#     Np1= size(xyθ,1)    # number of particles inside the boundary while Np is total number of particles
#     #xyθ = (rand(Np,3).-0.0).*repeat([L L 2π],Np)
#     xyθ[:,1:2], dists, superpose, uptriang = hardsphere(xyθ[:,1:2],R) #xyθ[:,1:2] gives x and y positions of intitial particles
#     abpe = ABPE2( Np1, L, R, v, 1e12DT, DR, xyθ[:,1], xyθ[:,2], xyθ[:,3])

#     return abpe, (dists, superpose, uptriang)
# end

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
function multiparticleE(Np::Integer, L::Float64, R::Float64, v::Union{Float64,Array{Float64,1},Distribution}, ω::Union{Float64,Array{Float64,1},Distribution}, Nt::Int64, δt::Float64, int_func::Function, forward::Bool, offcenter::Float64, range::Float64, int_params...)
    (Nt isa Int64) ? Nt : Nt=convert(Int64,Nt)
    
    ABPE = Vector{ABPE2}(undef,Nt+1) # Nt is number of time steps
    ABPE[1], matrices = initABPE( Np, L, R, v, ω, int_func, forward, offcenter, range, int_params...) # including initial hardsphere correction

    simulate!(ABPE, matrices, Nt, δt, forward, offcenter, range, int_func, int_params...)

    return position.(ABPE), orientation.(ABPE), force.(ABPE), torque.(ABPE)
end

function simulate!(ABPE, matrices, Nt, δt, forward, offcenter, range, int_func, int_params...)
    start = now()
    print_step = Nt÷100
    for nt in 1:Nt
        ABPE[nt+1] = update(ABPE[nt],matrices,δt, forward, offcenter, range, int_func, int_params...)#updating information at every step
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
function update(abpe::ABPE, matrices::Tuple{Matrix{Float64}, BitMatrix, BitMatrix}, δt::Float64, forward::Bool, offcenter::Float64, range::Float64, int_func::Function, int_params...) where {ABPE <: ABPsEnsemble}

    pθ = ( position(abpe), orientation(abpe) ) .+ step(abpe,δt, force(abpe), abpe.torque)

    periodic_BC_array!(pθ[1],abpe.L, abpe.R)
    #circular_wall_condition!(pθ[1],L::Float64, R, step_mem::Array{Float64,2})
    hardsphere!(pθ[1], matrices[1], matrices[2], matrices[3], abpe.R)
    # @btime hardsphere!($p[:,1:2], $matrices[1], $matrices[2], $matrices[3], $params.R)

    new_force, new_torque = force_torque(pθ[1], pθ[2], abpe.R, abpe.L, forward, offcenter, range, int_func, int_params...)
    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.v, abpe.ω, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2], new_force[:,1], new_force[:,2], new_torque )

    return new_abpe
end

function step(abpe::ABPE, δt::Float64, force::Array{Float64,2}, torque::Array{Float64,1}) where {ABPE <: ABPsEnsemble}    
    γₜ = diffusion_coeff(1e-6*abpe.R)[3]
    γᵣ = γₜ*abpe.R*abpe.R*8/6   
    if size(position(abpe),2) == 2
        δp = sqrt.(2*δt*abpe.DT)*randn(abpe.Np,2) .+ δt.*abpe.v.*[cos.(abpe.θ) sin.(abpe.θ)] .+ δt*force/γₜ
        δθ = sqrt(2*abpe.DR*δt)*randn(abpe.Np) .+ δt.*abpe.ω #.+ δt*torque/γᵣ
    else
        println("No step method available")
    end

    return (δp, δθ)
end

function force_torque(xy::Array{Float64,2}, θ::Array{Float64,1}, R::Float64, L::Float64, forward::Bool, offcenter::Float64, range::Float64, int_func::Function, int_params...) #offcenter
    xy_chgcen = xy .+ (2*forward-1) .* [cos.(θ) sin.(θ)] .* R*offcenter
    forces = hcat(interactions_range(xy_chgcen, R, L, range, size(xy,1), int_func, int_params...), zeros(size(xy,1)))
    orientations = [cos.(θ) sin.(θ) zeros(size(xy,1))]
    torques = offcenter * R * (cross.(eachrow(orientations), eachrow(forces)))
    return forces[:,1:2], reduce(hcat, torques)'[:,3]
end

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions for the hard sphere corrections
function hardsphere_correction!(xy::Array{Float64,2}, dists::Array{Float64,2}, superpose::BitArray{2}, R::Float64; tol::Float64=1e-3)
    Np = size(superpose,1) ##gives me the number of superpose lines 
    for np1 in 1:Np
        if any(superpose[np1,:]) #if at least one value of the np1 row is true 
            np2 = findfirst(superpose[np1,:])
            Δp = (xy[np1,:] - xy[np2,:]) .* ( ( (1+tol)*2R / dists[np1,np2] - 1 ) / 2 )
            xy[np1,:] += Δp
            xy[np2,:] -= Δp
            dists[np2,np2+1:Np] = pairwise(Euclidean(), xy[np2:np2,:], xy[np2+1:Np,:], dims=1 )  # distances for the row wise pair operation
            superpose[np2,np2+1:Np] = (dists[np2,np2+1:Np] .< 2R*(1-tol))  #????
        end
    end
    return nothing
end

function hardsphere!(xy::Array{Float64,2}, dists::Array{Float64,2}, superpose::BitArray{2}, uptriang::BitArray{2}, R::Float64; tol::Float64=1e-3)
    superpositions = 1
    counter = 0
    # @time begin
    while superpositions > 0
        dists .= pairwise(Euclidean(),xy,dims=1)
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

function hardsphere(xy::Array{Float64,2}, R::Float64; tol::Float64=1e-3) # called in initABPE
    Np = size(xy,1)
    dists = zeros(Np,Np)
    superpose = falses(Np,Np)
    uptriang = triu(trues(Np,Np),1)
    hardsphere!(xy, dists, superpose, uptriang, R; tol=tol)
    return xy, dists, superpose, uptriang
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
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Functions for updating reflective boundary AND WALL UPDATE

function multiparticleE_wall(Np::Integer, L::Float64, R::Float64, v::Float64, Nt::Int64=2, δt::Float64=1e-3)
    (Nt isa Int64) ? Nt : Nt=convert(Int64,Nt)
    
    ABPE = Vector{ABPE2}(undef,Nt+1)
    ABPE[1], matrices = initABPE( Np, L, R, v ) # including initial hardsphere correction
    
    simulate_wall!(ABPE, matrices, Nt, δt)
    println("I am in multiwall update")
    return position.(ABPE), orientation.(ABPE)
end

function simulate_wall!(ABPE, matrices, Nt, δt)
    # PΘ = [ (position(abpe), orientation(abpe)) ]
    # pθ = PΘ[1]
    print_step = Nt÷100
    start = now()
    for nt in 1:Nt
        # start_step = now()
        ABPE[nt+1] = update_wall(ABPE[nt],matrices,δt)
        if nt % print_step == 0
            elapsed = Dates.canonicalize(now()-start)
            # per_step = Dates.canonicalize(now()-start_step)
            print("\r$((100*nt÷Nt))%... Step $nt, total elapsed time $(elapsed)")#, time per step $per_step")
        end
    end
    print("\n")
    return nothing
end

function update_wall(abpe::ABPE, matrices::Tuple{Matrix{Float64}, BitMatrix, BitMatrix}, δt::Float64) where {ABPE <: ABPsEnsemble}
    memory_step = step(abpe,δt)
  
    pθ = ( position(abpe), orientation(abpe) ) .+ memory_step

    wall_condition!(pθ[1],abpe.L, abpe.R, memory_step[1])
    #elliptical_wall_condition!(pθ[1],abpe.L, abpe.R, memory_step[1])
    # elliptical_wall_condition!(pθ[2],pθ[1],abpe.L, abpe.R, memory_step[1])

    hardsphere!(pθ[1], matrices[1], matrices[2], matrices[3], abpe.R)
    # @btime hardsphere!($p[:,1:2], $matrices[1], $matrices[2], $matrices[3], $params.R)
    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.v, abpe.ω, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2] )

    return new_abpe
end

function wall_condition!(xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2}) # this condition is for square reflective boundary
  
	# Boundary conditions: horizontal edge
	idx = abs.(xy[:,1]) .> (L/2 - R)
	if any(idx)
		xy[idx,1] .-= 2*sign.(xy[idx,1]).*(abs.(xy[idx,1]) .- (L/2 - R)) 
	end
	# Boundary conditions: vertical edge
	idy = abs.(xy[:,2]) .> (L/2 - R)
	if any(idy)
        xy[idy,2] .-= 2*sign.(xy[idy,2]).*(abs.(xy[idy,2]) .- (L/2 - R))
	end
    #println("I am in square wall")
	return nothing
end

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#=function circular_wall_condition!(xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2}) # this condition is for cicular reflective boundary
  # here the condition is calculated w.r.t to r value of the particle and have no edges here
  # this is first method used 
     
    r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])

    rr = sqrt.(r)

    rθ= atan.(xy[:,2], xy[:,1])   # angle which the particle make with the origin 
     
    id = rr.> (L/2 - R)        # checking the condition for r vector id = 1 when true, 0 when false
    
    Δr = zeros(length(r),1)
    #println("$id")
             
        Δr[id] = rr[id].- (L/2 - R)
        rr.-= 2*(Δr)
        #println(size(rr)) # gives size of an array
        xy[id,1] = rr[id].*(cos.(rθ[id]))
        xy[id,2] = rr[id].*(sin.(rθ[id]))


	
	#println("I am in circular wall")
	return nothing
end
=#
function circular_wall_condition1g!(xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2}) # this condition is for cicular reflective boundary
    # here the condition is calculated w.r.t to the normal at the intersection point of the radial distance with the wall
    # this is second method used 
       
      r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
  
      rr = sqrt.(r)
  
      rθ= atan.(xy[:,2], xy[:,1])   # angle which the particle make with the origin 
       
      id = rr.> (L/2 - R)        # checking the condition for r vector id = 1 when true, 0 when false
      
      correction = zeros(length(r),1)
      projection = zeros(length(r),1)
      inside = zeros(length(r),1)
      hat_normal = zeros(length(r),1)
      x = [[p[1],p[2]] for p in eachrow(xy)]
      
      function grad(x::Array{Float64},R::Float64)
       f(x) = x[1]^2 + x[2]^2 - (L*L)/4
       df=ForwardDiff.gradient(f, [x[1],x[2]])
       return df
    end
    normal = grad.(x,R)
     
    hat_normal = normal./norm.(normal)  # unit normal vector
   
    correction = (rr.-(L/2 - 2*R)).*[[cos.(θ),sin.(θ)] for θ in rθ]
    
    projection = dot.(correction,hat_normal)
                                  
    
    inside = (0.5*L.-1*projection).*[[cos.(θ),sin.(θ)] for θ in rθ]         # the vector which is inside the boundry now
      #println("$id")

      xpos= [p[1] for p in inside]

      ypos= [p[2] for p in inside]
          
          xy[id,1] = xpos[id,1]
          xy[id,2] = ypos[id,1]
  
        println("I am in circular wall 1g")
    
      return nothing
  end

function elliptical_wall_condition!(orientation::Array{Float64,1},xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2}) # this condition is for cicular reflective boundary
    # here the condition is calculated w.r.t to the normal at the intersection point of the radial distance with the wall
    # this is second method used 
    # orientation of particle will also change   
        a= L/2
        b= L/4

        a1= a-R
        b1= b-R
         
        #e = sqrt(1-(b/a)^2)
      r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
  
      rₚ = sqrt.(r)                # particle r co ordinate
  
      rθ= atan.(xy[:,2], xy[:,1])   # angle which the particle make with the origin 

      rₒ = ((a*b)./(sqrt.(((a*sin.(rθ)).^2) .+ (b*cos.(rθ)).^2)))

      rₑ = rₒ .- R
   
      #rₑ = b/sqrt.(1 .-((e*cos.(rθ)).^2))
      id = (rₚ .> (rₑ))
    
      x = [[p[1],p[2]] for p in eachrow(xy)]
      
      function grad(x::Array{Float64},a::Float64,b::Float64)
        
       f(x) = (x[1]^2)*b^2 + (x[2]^2)*a^2 - (a*b)^2
       df=ForwardDiff.gradient(f, [x[1],x[2]])
       return df
       end
    normal = grad.(x,a1,b1)     # gradient calculated from each particle position onto the ellipse boundary 
     
    hat_normal = normal./norm.(normal)  # unit normal vector

    
    correction_x = (rₚ .- rₑ) .*[(cos.(θ)) for θ in rθ]

    correction_y = (rₚ .- rₑ) .*[(sin.(θ)) for θ in rθ]
    
    c= [correction_x correction_y]

    correction = [[p[1],p[2]] for p in eachrow(c)]
  
    projection = dot.(correction,hat_normal)


    cᵥ = 2*(projection.+0*R).* hat_normal #vector of vector
    # to access this vector I am breaking it in x and y in the following lines

    cᵥx =  [p[1] for p in cᵥ]

    cᵥy =  [p[2] for p in cᵥ]
############################# calculations for orientation vector###############################
###########################3 conversion of vector into matrix was job of this part###########################
        xy[id,1] .-= cᵥx[id]
        xy[id,2] .-= cᵥy[id]

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

    for (i, xyc) in enumerate(eachrow(xy))
        # Compute distances and apply periodic boundary conditions
        xy_shifted = xy .- xyc'
        periodic_BC_array!(xy_shifted, L, R)

        # Filter to get particles within interaction range and not at the origin
        inside = all(reshape(abs.(xy_shifted) .<= threshold, Np, 2), dims=2)[:,1]
        xy_inside = xy_shifted[inside, :]

        # Remove central particle [0, 0] from interaction set
        xy_inside = xy_inside[xy_inside[:,1] .!= 0 .|| xy_inside[:,2] .!= 0, :]

        if isempty(xy_inside)
            ΣFtot[i, :] .= 0.0
            continue
        end

        # Compute distances and interaction forces
        dists = d2(xy_inside)
        dists_nonzero = dists[dists .!= 0]
        forces = int_func.(dists_nonzero, int_params...)
        
        # Compute normalized direction vectors
        dirs = xy_inside ./ dists_nonzero

        # Sum forces for each direction and assign to ΣFtot
        ΣFtot[i, :] .= .-sum(forces .* dirs, dims=1)'
    end

    return ΣFtot
end


lennard_jones(x, σ, ϵ) = 24*ϵ*(((2*σ^(12))./(x.^(13))).- (σ^(6)./(x.^(7))))
shifted_lennard_jones(x, σ, ϵ, shift) = 24*ϵ*(((2*σ^(12))./((x-shift).^(13))).- (σ^(6)./((x-shift).^(7))))


function rectified_lennard_jones(x::Array{Float64,2}, σ::Float64, ϵ::Float64)
    rmin = σ*2^(1/6)
    if x > rmin
        return 0
    else
        return 24*ϵ*(((2*σ^(12))./(x^(13))).- (σ^(6)./(x^(7))))
    end
end

function purely_attractive_lennard_jones(x, σ::Float64, ϵ::Float64)
    rmin = σ*2^(1/6) + σ/2
    if x > rmin
        return 24*ϵ*(((2*σ^(12))./((x + σ/2)^(13))).- (σ^(6)./((x + σ/2)^(7))))
    else
        return 0
    end
end

function rlj_boundary(x::Vector{Float64}, σ::Float64, ϵ::Float64)
    rmin = σ*2^(1/6) + σ/2
    if x > rmin
        return 0

    elseif x < σ
        return 390144*ϵ/σ
    else
        return 24*ϵ*(((2*σ^(12))./((x.-σ/2)^(13))).- (σ^(6)./((x.-σ/2)^(7)))) 
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
    return 24*ϵ*(((2*σ^(12))./((x+shift)^(13))).- (σ^(6)./((x+shift)^(7))))
end

coulomb(x,k) = k/(x*x)
