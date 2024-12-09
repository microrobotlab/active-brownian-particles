# PURPOSE: Output of ABP main 
# all codes in repository are complied in this function
#VARIABLES: Destination folder path and filename
include("ABP main interactions opt.jl")
# include("ABP main.jl")
include("ABP_file.jl")
include("ABP radialdistribution.jl")
# include("ABP SD.jl")
include("ABP multifolder.jl")
# include("ABP radialdensity.jl")
include("ABP plot_animation.jl")
# include("ABP average.jl")
using CSV, DataFrames, Dates, Distances, Distributions, JLD2, Logging, NaNStatistics, Printf, Random
gr()

## USER INTERFACE
# destination folders selection
path="C:\\Users\\NiccoloP\\Documents\\thesis\\simulations" # destination directory path

## PARAMETERS SET
# Simulation parameters
Nt = Int(1e3)           # number of steps
Delta_t = 1e-3          # s step time
ICS=10                  # Number of intial conditons to be scanned 
animation_ds = 4     # Downsampling in animation
measevery = Int(1e1)           # Downsampling in file
animation = false
radialdensity = false

# Physical parameters
BC_type = :periodic    # :periodic or :wall
box_shape = :square    # shapes: :square, :circle, :ellipse
R = 2.0		           # μm particle radius
L = 200.0 	           # μm box length
packing_fraction = (pi*R^2/L^2)*500 # Largest pf for spherical beads π/4 = 0.7853981633974483
# Velocities can also be distributions e.g. v = Normal(0.,0.025)
v = [20.] 	            # μm/s particle s
ω = 0.        # s⁻¹ particle angular velocity
T = 300. # K temperature

# Interaction parameters
int_func = coulomb
forward = true
intrange = 5. # interaction range
offcenter = 1e-4
int_params = (1.) # σ and ϵ in the case of LJ 

#-------------------------------------------------------------------------------------------------------------------

## PRELIMINARY CALCULATIONS
DT, DR, γ = diffusion_coeff(1e-6R,T).*[1e12, 1, 1] # Translational and Rotational diffusion coefficients, drag coefficient
T_tot = Delta_t*Nt
actual_steps = (Nt÷measevery)+1

if box_shape == :square
    Np = round(Int,packing_fraction*L^2/(pi*R^2))
    density = Np/L^2
    a,b = L,L
elseif box_shape == :ellipse
    a=L/2 # Larger semiaxis
    b=L/4 # Smaller semiaxis
    area_el = π*a*b
    Np = round(Int,packing_fraction*area_el/(pi*R^2))
    density = Np/area_el
end

## NAMING AND FILE MANAGEMENT
datestamp=Dates.format(now(),"YYYYmmdd-HHMMSS")  # today's date
println(datestamp)

pathmain= joinpath(path, datestamp)

mainfolder= mkdir(pathmain)    # creates a folder named today's date

patht= joinpath(pathmain,"data\\")

mainfolder1= mkdir(patht)

folders=  multipledir(patht,ICS) 

# Info printing on shell and file
infos = @sprintf "Box shape: %s\nNumber of particles = %i\nNumber density = %s μm⁻²\nR=%.1f μm \nT = %.1f (K)\nv=%s (μm/s) \nω=%s (rad/s)\nCharacteristic lengths: (a=%.1f b=%.1f) μm\npf=%s\nIntegration step: dt=%.0e s \nSimulation downsampling: %i\nNumber of steps: Nt=%.1e\nTotal simulated time T_tot = %.2e s\n\nInteraction function: %s with parameters: %s\nRange: %.1f μm\nOffcenter: %.2f" box_shape Np density R T v ω a b packing_fraction Delta_t measevery  Nt T_tot int_func int_params intrange offcenter*(2*forward-1)

println(infos)

info_file_path = joinpath(mainfolder, "simulation_info.txt")
open(info_file_path, "w") do infile
    write(infile, infos)
end

info_dict = Dict(
    "Box shape" => box_shape,
    "Np" => Np,
    "numdensity" => density,
    "R" => R,
    "T" => T,
    "v" => v,
    "ω" => ω,
    "a" => a,
    "b" => b,
    "pf" => packing_fraction,
    "dt" => Delta_t,
    "measevery" => measevery,
    "Nt" => Nt,
    "T_tot" => T_tot,
    "int_func" => string(int_func),
    "int_params" => int_params,
    "intrange" => intrange,
    "offcenter" => offcenter,
)

JLD2.save(joinpath(mainfolder, "siminfo_dict.jld2"), info_dict)

if box_shape == :square
    for i=1:ICS
        pathf= joinpath(patht, "run$i\\")
        filename= "$datestamp"*"_run$i"
        @show pathf= pathf*filename

        # Simulation and file storage
        @info "$(now()) Started simulation #$i"
        if BC_type == :periodic
            history = multiparticleE(Np,L,R,T,v,ω,Nt,measevery,Delta_t, int_func, forward, offcenter, intrange, int_params...,) # has values of x and y position in each frame in history[1]
        end
        #---------------------------------------------------------------------------------------------------------------------
        # file store
        file_store_txt(history,actual_steps,pathf,downsampling = 1)
        #---------------------------------------------------------------------------------------------------------------------
        # animation
        if i==1
            animation_from_history(history,pathf,L,R,Np,Delta_t,actual_steps,measevery,animation_ds, show = false, record=true, final_format = "mkv", color_code_dir = true)
        end
        #---------------------------------------------------------------------------------------------------------------------
        # analysis
        if radialdensity
            physicalanalysis1(pathf, 500, L, R,".txt")
        end
    end
end