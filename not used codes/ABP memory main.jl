using CalculusWithJulia, ForwardDiff
using Plots,Distances,NaNStatistics,CSV, DataFrames
using Dates
# Define an "ABPsEnsemble" Type
abstract type ABPsEnsemble end
include("ABP file.jl")
# Define a specific type for 2D ABPsEnsembles (CURRENTLY ASSUMING ALL PARTICLES ARE EQUAL)
struct ABPE2 <: ABPsEnsemble
    Np::Int64                      # number of particles®®
    L::Float64                      # size of observation space (μm)
	R::Float64  # Radius (μm)                                   --> Vector{Float64}(undef,Np)
	v::Float64 	# velocity (μm/s)                               --> Vector{Float64}(undef,Np)
	DT::Float64 # translational diffusion coefficient (μm^2/s)  --> Vector{Float64}(undef,Np)
	DR::Float64 # rotational diffusion coefficient (rad^2/s)    --> Vector{Float64}(undef,Np)
	x::Vector{Float64}    # x position (μm)
	y::Vector{Float64}    # y position (μm)
	θ::Vector{Float64}    # orientation (rad)
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

@time function hardsphere!(xy::Array{Float64,2}, dists::Array{Float64,2}, superpose::BitArray{2}, uptriang::BitArray{2}, R::Float64; tol::Float64=1e-3)
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

function hardsphere(xy::Array{Float64,2}, R::Float64; tol::Float64=1e-3) # called in initABPE
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
#-------------------

function initABPE(Np::Int64, L::Float64, R::Float64, v::Float64; T::Float64=300.0, η::Float64=1e-3)
    # translational diffusion coefficient [m^2/s] & rotational diffusion coefficient [rad^2/s] - R [m]
    # Intial condition will be choosen as per the geometry under study
    DT, DR = diffusion_coeff(1e-6R)
     
    # ONLY 2D!
    k=0.5
    xyθ1 = (rand(Np,3).-k).*repeat([L L 2π],Np) # 3 dim matrix with x, y and θ 
    r = (xyθ1[:,1]).*(xyθ1[:,1]) + (xyθ1[:,2]).*(xyθ1[:,2])
  
    rₚ = sqrt.(r)   
    α =atan.(xyθ1[:,2], xyθ1[:,1]) 
    a= L/2
        b= L/4
    #rₑ = (a*b)./(sqrt.(((a*sin.(xyθ1[:,3])).^2) .+ (b*cos.((xyθ1[:,3]))).^2))  # r value for boundary
    rₑ = (a-R)*(b-R)./(sqrt.((((a-R)*sin.(α)).^2) .+ ((b-R)*cos.((α))).^2))  # r value for boundary
    #rₑ = b/sqrt.(1 .-((e*cos.(rθ)).^2))
    id = (rₚ .< (rₑ))
    xyθ = [xyθ1[id,1] xyθ1[id,2] xyθ1[id,3]]

    Np1= size(xyθ,1)  
   #xyθ = (rand(Np,3).-0.0).*repeat([L L 2π],Np)
    xyθ[:,1:2], dists, superpose, uptriang = hardsphere(xyθ[:,1:2],R) #xyθ[:,1:2] gives x and y positions of intitial particles
    abpe = ABPE2( Np1, L, R, v, 1e12DT, DR, xyθ[:,1], xyθ[:,2], xyθ[:,3])

    return abpe, (dists, superpose, uptriang)
end

function diffusion_coeff(R::Float64, T::Float64=300.0, η::Float64=1e-2)
    # Boltzmann constant [J/K]
    kB = 1.38e-23
    # friction coefficient [Ns/m]
    γ = 6*pi*R*η
    # translational diffusion coefficient [m^2/s]
    DT = kB*T/γ
    # rotational diffusion coefficient [rad^2/s]
    DR = 6*DT/(8*R^2)
    return DT, DR
end

function multiparticleE_wall(Np::Integer, L::Float64, R::Float64, v::Float64, Nt::Int64=2, δt::Float64=1.0e-03)
    (Nt isa Int64) ? Nt : Nt=convert(Int64,Nt)
    #println("time steps $δt")
    Nt1=Int(Nt/10)
    ABPE = Vector{ABPE2}(undef,1)
    ABPE[1], matrices = initABPE( Np, L, R, v ) # including initial hardsphere correction
   # @show matrices
   current_value = deepcopy(ABPE[1])  # This will be used for updates at each step
 
   for i in 2:Nt
    #temp_index = Int(ceil(i/10))
  
    current_value = update_wall(current_value,matrices,δt)
        if mod(i,10) == 0
            # temp_index += 1
            push!(ABPE,current_value)
            # ABPE[temp_index] = ABPE[temp_index-1]
            # println( ABPE)
            # println("index= $temp_index step = $i")
            # if temp_index == 11
            #     break
            # end
           
    end
        
        # # ABPE[i+1] = update_wall(ABPE[i-count],matrices,δt)
        # # push!(ABPE,ABPE[i+1])
        # println("Step $i")
    end
    # simulate_wall!(ABPE, matrices, Nt, δt)
 
    return position.(ABPE), orientation.(ABPE)
 end

function simulate_wall!(ABPE, matrices, Nt, δt)
    # PΘ = [ (position(abpe), orientation(abpe)) ]
    # pθ = PΘ[1]
    
    for nt in 1:Nt
        ABPE[nt+1] = update_wall(ABPE[nt],matrices,δt)
        println("Step $nt")
      
    end
    return nothing
end
position(abpe::ABPE2) = [ abpe.x abpe.y ]
orientation(abpe::ABPE2) = abpe.θ

function update_wall(abpe::ABPE, matrices::Tuple{Matrix{Float64}, BitMatrix, BitMatrix}, δt::Float64) where {ABPE <: ABPsEnsemble}
    memory_step = step(abpe,δt)
  
    pθ = ( position(abpe), orientation(abpe) ) .+ memory_step

    #wall_condition!(pθ[1],abpe.L, abpe.R, memory_step[1])
    #elliptical_wall_condition!(pθ[1],abpe.L, abpe.R, memory_step[1])
    #elliptical_wall_condition!(pθ[2],pθ[1],abpe.L, abpe.R, memory_step[1])

    hardsphere!(pθ[1], matrices[1], matrices[2], matrices[3], abpe.R)
    # @btime hardsphere!($p[:,1:2], $matrices[1], $matrices[2], $matrices[3], $params.R)
    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.v, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2] )

    return new_abpe
end

function step(abpe::ABPE, δt::Float64) where {ABPE <: ABPsEnsemble}
    
    if size(position(abpe),2) == 2
        δp = sqrt.(2*δt*abpe.DT)*randn(abpe.Np,2) .+ abpe.v*δt*[cos.(abpe.θ) sin.(abpe.θ)] #.+ δt*δt*  attractive_interactions!(position(abpe),2.0)
        δθ = sqrt(2*abpe.DR*δt)*randn(abpe.Np)
    else
        println("No step method available")
    end
    return (δp, δθ)
end


L = 100.0 	# μm box length
R = 2.0	# μm particle radius
v = 50.0 	# μm/s particle velocity
a=L/2
b=L/4

pf_factor = (R^2)
DT, DR = diffusion_coeff(R).*[1e12, 1]
packing_fraction = 0.2

Np = round(Int,packing_fraction*a*b/(R^2))  #Np is the number of particles inside the ellipse
#π
Nt = 100# Nt is the number of steps 
δt = 1.0e-2 #L/(v*Nt) # δt is the time step

graph_wall = multiparticleE_wall(Np,L,R,v,Nt,δt) 

file_store(graph_wall,Nt,pathf)