using CalculusWithJulia, ForwardDiff

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

## Get position and orientation of the particle or ensemble (CURRENTLY ONLY 2D)
position(abp::ABP2) = ( abp.x, abp.y, abp.θ )

## Initialize ABP particle (2D)
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

## Calculate the step position (and orientation) for the particle (2D)
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

## Update position and orientation of the particle 
update(abp::ABP, step) where {ABP <: ActiveBrownianParticle} = ABP( (position(abp) .+ step)..., abp.R, abp.v, abp.DT, abp.DR )

##Calculate diffusion coefficient
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
function trajectory(abp::ABP, N, δt::Float64) where {ABP <: ActiveBrownianParticle}
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
    plot(x, y, legend=false, aspect_ratio=:equal);
	scatter!([x[end]],[y[end]], legend=false);
    xlabel!("x [μm]");
    ylabel!("y [μm]");
end


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

#-------------------------------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE MSD FOR ONE TRAJECTORY
function MSDcalculate(x,y,N_Max,N)
    ltrack= N+1
    msd=zeros(N_Max+1)
    for tau in 1:N_Max
        for i in tau+1:ltrack
            msd[tau+1]+=((x[i]-x[i-tau])^2+(y[i]-y[i-tau])^2)/(ltrack-tau)
        end
    end
    return msd
end

#-----------------------------------------------------------------------------------------------------------------
# FUNCTION THAT GENERETAES n TRAJECTORIES, PLOTS FIRST FOUR OF THEM, CALCULATE AVERAGE MSD OF ALL OF THEM
function traj_and_MSD(x0, y0, R::Float64, v::Float64, num_traj::Int64, N, Delta_t::Float64, N_Max)
    graph = plot();
    matrMSD = fill(NaN, N_Max+1, num_traj)

    for i in 1:num_traj
        orientazione = rand()*2*pi
        orientazione = round(orientazione, digits=3)
        abp = initABP( (x0, y0, orientazione), R, v);

        p, t = trajectory( abp, N, Delta_t);

        x = [pi[1] for pi in p]
        y = [pi[2] for pi in p]

        if i <= 4 
            plot!(x,y, range=[-50,50],  title = "ActiveParticle (R=$R µm, v=$v µm/s)", aspect_ratio= :equal, legend=false)
            xlabel!("x [μm]")
            ylabel!("y [μm]")
        end

        matrMSD[1:N_Max+1, i] = MSDcalculate(x,y, N_Max, N)
    end

    return graph, matrMSD 

end

#--------------------------------------------------------------------------------------------
# Define an "ABPsEnsemble" Type
abstract type ABPsEnsemble end

# Define a specific type for 2D ABPsEnsembles (CURRENTLY ASSUMING ALL PARTICLES ARE EQUAL)
struct ABPE2 <: ABPsEnsemble
    Np::Int64                       # number of particles®®
    L::Float64                      # size of observation space (μm)
	R::Float64                      # Radius (μm)                                   --> Vector{Float64}(undef,Np)
	v::Float64 	                    # velocity (μm/s)                               --> Vector{Float64}(undef,Np)
	DT::Float64                     # translational diffusion coefficient (μm^2/s)  --> Vector{Float64}(undef,Np)
	DR::Float64                     # rotational diffusion coefficient (rad^2/s)    -->
	x::Vector{Float64}              # x position (μm)
	y::Vector{Float64}              # y position (μm)
	θ::Vector{Float64}              # orientation (rad)
end

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  INITIALIZE ABP ENSEMBLE (2D)
function initABPE(Np::Int64, L::Float64, R::Float64, v::Float64; T::Float64=300.0, η::Float64=1e-3)
    # translational diffusion coefficient [m^2/s] & rotational diffusion coefficient [rad^2/s] - R [m]
    # Intial condition will be choosen as per the geometry under study
    DT, DR = diffusion_coeff(1e-6R)
    k=0.5
    xyθ1 = (rand(Np,3).-k).*repeat([L L 2π],Np) #  3dim matrix with x, y and θ

    #calculations for ellipitical boundary
    r = (xyθ1[:,1]).*(xyθ1[:,1]) + (xyθ1[:,2]).*(xyθ1[:,2])
  
    rₚ = sqrt.(r)   
    α =atan.(xyθ1[:,2], xyθ1[:,1]) 
    a= L/2
    b= L/4
   
    rₑ = (a-R)*(b-R)./(sqrt.((((a-R)*sin.(α)).^2) .+ ((b-R)*cos.((α))).^2))  

    id = (rₚ .< (rₑ))
    xyθ = [xyθ1[id,1] xyθ1[id,2] xyθ1[id,3]]

    Np1= size(xyθ,1)                  # Number of particles inside the boundary while Np is total number of particles
 
    xyθ[:,1:2], dists, superpose, uptriang = hardsphere(xyθ[:,1:2],R) #xyθ[:,1:2] gives initial x and y positions of particles
    abpe = ABPE2( Np1, L, R, v, 1e12DT, DR, xyθ[:,1], xyθ[:,2], xyθ[:,3])

    return abpe, (dists, superpose, uptriang)
end

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  SIMULATE ENSEMBLE OF SPHERICAL PARTICLES WITH OPEN BOUNDARY
function multiparticleE(Np::Integer, L::Float64, R::Float64, v::Float64, Nt::Int64=2, δt::Float64=1e-3)
    (Nt isa Int64) ? Nt : Nt=convert(Int64,Nt)
    
    ABPE = Vector{ABPE2}(undef,Nt+1) # Nt is number of time steps
    ABPE[1], matrices = initABPE( Np, L, R, v ) 
    
    simulate!(ABPE, matrices, Nt, δt)

    return position.(ABPE), orientation.(ABPE)
end

function simulate!(ABPE, matrices, Nt, δt)
 
    for nt in 1:Nt
        ABPE[nt+1] = update(ABPE[nt],matrices,δt) 
        println("Step $nt")
    end
    return nothing
end

position(abpe::ABPE2) = [ abpe.x abpe.y ]
orientation(abpe::ABPE2) = abpe.θ

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# UPDATE PARTICLES FOR THE NEXT TIME STEP 
function update(abpe::ABPE, matrices::Tuple{Matrix{Float64}, BitMatrix, BitMatrix}, δt::Float64) where {ABPE <: ABPsEnsemble}
    pθ = ( position(abpe), orientation(abpe) ) .+ step(abpe,δt)

    periodic_BC_array!(pθ[1],abpe.L, abpe.R)
 
    hardsphere!(pθ[1], matrices[1], matrices[2], matrices[3], abpe.R)

    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.v, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2] )

    return new_abpe
end


function step(abpe::ABPE, δt::Float64) where {ABPE <: ABPsEnsemble}
    
    if size(position(abpe),2) == 2
        δp = sqrt.(2*δt*abpe.DT)*randn(abpe.Np,2) .+ abpe.v*δt*[cos.(abpe.θ) sin.(abpe.θ)]
        δθ = sqrt(2*abpe.DR*δt)*randn(abpe.Np)
    else
        println("No step method available")
    end
 
    return (δp, δθ)
end

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# HARD SPHERE CORRECTIONS 
function hardsphere_correction!(xy::Array{Float64,2}, dists::Array{Float64,2}, superpose::BitArray{2}, R::Float64; tol::Float64=1e-3)
    Np = size(superpose,1) 
    for np1 in 1:Np
        if any(superpose[np1,:]) #if at least one value of the np1 row is true 
            np2 = findfirst(superpose[np1,:])
            Δp = (xy[np1,:] - xy[np2,:]) .* ( ( (1+tol)*2R / dists[np1,np2] - 1 ) / 2 )
            xy[np1,:] += Δp
            xy[np2,:] -= Δp
            dists[np2,np2+1:Np] = pairwise(Euclidean(), xy[np2:np2,:], xy[np2+1:Np,:], dims=1 )  # distances for the row wise pair operation
            superpose[np2,np2+1:Np] = (dists[np2,np2+1:Np] .< 2R*(1-tol))  
        end
    end
    return nothing
end

function hardsphere!(xy::Array{Float64,2}, dists::Array{Float64,2}, superpose::BitArray{2}, uptriang::BitArray{2}, R::Float64; tol::Float64=1e-3)
    superpositions = 1
    counter = 0
    while superpositions > 0
        dists .= pairwise(Euclidean(),xy,xy,dims=1)
        superpose .= (dists .< 2R*(1-tol)).*uptriang
     
        superpositions = sum(superpose)
        
        if superpositions > 0
            hardsphere_correction!(xy,dists,superpose,R,tol=tol)
        end
        counter += 1
     
        if counter >= 100
            println("$superpositions superpositions remaining after 100 cycles")
            break
        end
    end
    
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
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OPEN BOUNDARY CONDITION
# When a particle crosses an edge it reappears on the opposite side

function periodic_BC_array!(xy::Array{Float64,2},L::Float64, R)  
	# Boundary conditions: horizontal edge
	idx = abs.(xy[:,1]) .> L/2 + R 
	if any(idx)
		xy[idx,1] .-= sign.(xy[idx,1]).*L  
	end
	# Boundary conditions: vertical edge
	idy = abs.(xy[:,2]) .> L/2 + R
	if any(idy)
		xy[idy,2] .-= sign.(xy[idy,2]).*L
	end
	return nothing
end
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SIMULATE ENSEMBLE OF SPHERICAL PARTICLES WITH CLOSED BOUNDARY/CONFINEMENT

function multiparticleE_wall(Np::Integer, L::Float64, R::Float64, v::Float64, Nt::Int64=2, δt::Float64=1e-3)
    (Nt isa Int64) ? Nt : Nt=convert(Int64,Nt)
    
    ABPE = Vector{ABPE2}(undef,Nt+1)
    ABPE[1], matrices = initABPE( Np, L, R, v ) 
    
    simulate_wall!(ABPE, matrices, Nt, δt)
    println("I am in multiwall update")
    return position.(ABPE), orientation.(ABPE)
end

function simulate_wall!(ABPE, matrices, Nt, δt)
       for nt in 1:Nt
        ABPE[nt+1] = update_wall(ABPE[nt],matrices,δt)
        
      
    end
    return nothing
end

function update_wall(abpe::ABPE, matrices::Tuple{Matrix{Float64}, BitMatrix, BitMatrix}, δt::Float64) where {ABPE <: ABPsEnsemble}
    memory_step = step(abpe,δt)
  
    pθ = ( position(abpe), orientation(abpe) ) .+ memory_step

    #wall_condition!(pθ[1],abpe.L, abpe.R, memory_step[1])    # Square boundary
  
    elliptical_wall_condition!(pθ[2],pθ[1],abpe.L, abpe.R, memory_step[1]) # Elliptical boundary

    hardsphere!(pθ[1], matrices[1], matrices[2], matrices[3], abpe.R)
 
    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.v, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2] )

    return new_abpe
end

######################################## SQUARE REFLECTIVE BOUNDARY#######################################################
function wall_condition!(xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2}) 
  
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
    println("I am in square wall")
	return nothing
end

############################## ELLIPTICAL REFLECTIVE BOUNDARY#################################################################
function elliptical_wall_condition!(orientation::Array{Float64,1},xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2}) 
        a= L/2
        b= L/4

        a1= a-R
        b1= b-R
         
        #e = sqrt(1-(b/a)^2)
      r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
  
      rₚ = sqrt.(r)                     # particle r co ordinate
  
      rθ= atan.(xy[:,2], xy[:,1])       # angle which the particle make with the origin 

      rₒ = ((a*b)./(sqrt.(((a*sin.(rθ)).^2) .+ (b*cos.(rθ)).^2)))

      rₑ = rₒ .- R
   
      id = (rₚ .> (rₑ))
     
      ###################### Calculations only on the particles outside the boundary #####################################
      if any(id)     
      positions = [[p[1],p[2]] for p in eachrow(xy)]

      function grad(x::Array{Float64},a::Float64,b::Float64)
        
       f(x) = (x[1]^2)*b^2 + (x[2]^2)*a^2 - (a*b)^2
       df=ForwardDiff.gradient(f, [x[1],x[2]])
       return df
       end
    normal = grad.(positions[id],a1,b1)     # gradient calculated from each particle position onto the ellipse boundary 
     
    hat_normal = normal./norm.(normal)       # unit normal vector

    
    correction_x = (rₚ[id] .- rₑ[id]) .*[(cos.(θ)) for θ in rθ[id]]

    correction_y = (rₚ[id] .- rₑ[id]) .*[(sin.(θ)) for θ in rθ[id]]
    
    c= [correction_x correction_y]

    correction = [[p[1],p[2]] for p in eachrow(c)]
  
    projection = dot.(correction,hat_normal)


    cᵥ = 2*(projection).* hat_normal 

    cᵥx =  [p[1] for p in cᵥ]

    cᵥy =  [p[2] for p in cᵥ]

########################### updating the posiitons ###########################
       
        xy[id,1] .-= cᵥx
        xy[id,2] .-= cᵥy
       end
      return nothing
  end
