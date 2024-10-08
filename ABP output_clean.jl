# PURPOSE: Output of ABP main 
# all codes in repository are complied in this function
#VARIABLES: Destination folder path and filename
include("ABP main interactions.jl")
# include("ABP main.jl")
include("ABP file.jl")
include("ABP analysis.jl")
# include("ABP SD.jl")
include("ABP multifolder.jl")
include("ABP radialdensity.jl")
# include("ABP average.jl")
using CSV, DataFrames, Dates, Distances, Distributions, Logging, NaNStatistics, Plots, Printf, Random
gr()

## USER INTERFACE
# destination folders selection
path="C:\\Users\\nikko\\OneDrive\\Documents\\Uni\\magistrale\\tesi\\simulations\\" # destination directory path

## PARAMETERS SET
# Simulation parameters
Nt = 5000             # number of steps
Delta_t = 1e-3          # s step time
ICS=1                   # Number of intial conditons to be scanned 
animation_ds = 100     # Downsampling in animation
file_ds = 100           # Downsampling in file


# Physical parameters
BC_type = :periodic    # :periodic or :wall
box_shape = :square    # shapes: :square, :circle, :ellipse
L = 300.0 	           # μm box length
R = 2.0		           # μm particle radius
packing_fraction = 0.1 # Largest pf for spherical beads π/4 = 0.7853981633974483
# Velocities can also be distributions e.g. v = Normal(0.,0.025)
v = 10. 	           # μm/s particle velocity
ω = 1.                 # s⁻¹ particle angular velocity

#-------------------------------------------------------------------------------------------------------------------

## PRELIMINARY CALCULATIONS
DT, DR, γ = diffusion_coeff(R).*[1e12, 1, 1] # Translational and Rotational diffusion coefficients, drag coefficient
T_tot = Delta_t*Nt

if box_shape == :square
    Np = round(Int,packing_fraction*L^2/(pi*R^2))
    density = Np/L^2
    # @printf "Box shape: %s\nNumber of particles = %i\nNumber density = %.s μm⁻²" box_shape Np density
    a,b = L,L
elseif box_shape == :ellipse
    a=L/2 # Larger semiaxis
    b=L/4 # Smaller semiaxis
    area_el = π*a*b
    Np = round(Int,packing_fraction*area_el/(pi*R^2))
    density = Np/area_el
    # @printf "Box shape: %s\nNumber of particles = %i\nNumber density = %.s μm⁻²" box_shape Np density
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
infos = @sprintf "Box shape: %s\nNumber of particles = %i\nNumber density = %s μm⁻²\nR=%.1f μm \nv=%s (μm/s) \nω=%s (rad/s)\nCharacteristic lengths: (a=%.1f b=%.1f) μm\npf=%s\nIntegration step: dt=%.0e s\nNumber of steps: Nt=%.1e\nTotal simulated time T_tot = %.2e s" box_shape Np density R v ω a b packing_fraction Delta_t Nt T_tot

println(infos)

info_file_path = joinpath(patht, "simulation_info.txt")

open(info_file_path, "w") do infile
    write(infile, infos)
end

if box_shape == :square
    for i=1:ICS
        pathf= patht*"\\run$i\\"
        filename= "$datestamp"*"_run$i"
        pathf= pathf*filename

        # Simulation and file storage
        start_sim = now()
        @info "$start_sim Started simulation #$i"
        if BC_type == :periodic
            graph_wall = multiparticleE(Np,L,R,v,ω,Nt,Delta_t) # has values of x and y position in each frame in graph_wall[1]
            elapsed_time = Dates.canonicalize(now()-start_sim)
            println("multiparticleE complied: elapsed time $elapsed_time\n")
        end
        #---------------------------------------------------------------------------------------------------------------------
        # file store
        file_store(graph_wall,Nt,pathf,downsampling = file_ds)
                #---------------------------------------------------------------------------------------------------------------------
        # animation
        anim = @animate for i = 1:animation_ds:Nt
            scatter(graph_wall[1][i][:,1], graph_wall[1][i][:,2], aspect_ratio=:equal, lims=(-L/2, L/2),markersize=350R/L,marker =:circle,legend=false, title = "$(Np) particles, steps $i, ",)
            
            plot!([L/2], seriestype="vline", color=:black)  #square
            plot!([-L/2], seriestype="vline", color=:black)
            plot!([L/2], seriestype="hline", color=:black)
            plot!([-L/2], seriestype="hline", color=:black)
            # quiver!(graph_wall[1][i][:,1],graph_wall[1][i][:,2],quiver=(4*cos.(graph_wall[2][i,1]),4*sin.(graph_wall[2][i,1])), color=:red)
        end
   
        f1= pathf*".gif"
        gif(anim, f1)
    end
end
