# PURPOSE: Output of ABP main 
# all codes in repository are complied in this function
#VARIABLES: Destination folder path and filename
include("ABP main interactions heun.jl")
# include("ABP main.jl")
include("ABP_file.jl")
include("ABP radialdistribution.jl")
# include("ABP SD.jl")
include("ABP multifolder.jl")
# include("ABP radialdensity.jl")
include("ABP plot_animation.jl")
include("ABP orderparameters.jl")
# include("ABP average.jl")
using  CSV, DataFrames, Dates, Distributions, IterTools, JLD2, Logging, Printf

path = "D:\\nic_simulations\\tests"

## PARAMETERS SET
# Simulation parameters
Nt = Int(2e4)           # number of steps
δt = 5e-3       # s step time
ICS=1                  # Number of intial conditons to be scanned 
animation_ds = 1     # Downsampling in animation
measevery = Int(1)           # Downsampling in file
animation = true
radialdensity = false

# Physical parameters
BC_type = :periodic    # :periodic or :wall
box_shape = :square    # shapes: :square, :circle, :ellipse
R = 2.0		           # μm particle radius
L = 100.0 	           # μm box length
packing_fraction = (pi*R^2/L^2)*20 # Largest pf for spherical beads π/4 = 0.7853981633974483

# Velocities can also be distributions e.g. v = Normal(0.,0.025)
v =           # μm/s particle s
ω = [0.]      # s⁻¹ particle angular velocity
T = 300. # K temperature

# Interaction parameters
int_func = lennard_jones
forward = false
intrange = 20R # interaction range
offcenter = 0.25   
int_params = (2R, .1) # σ and ϵ in the case of LJ 

## PRELIMINARY CALCULATIONS
DT, DR, γ = diffusion_coeff(1e-6R,T).*[1e12, 1, 1] # Translational and Rotational diffusion coefficients, drag coefficient
T_tot = δt*Nt
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
   

for ω in ωs
    ## NAMING AND FILE MANAGEMENT
    datestamp=Dates.format(now(),"YYYYmmdd-HHMMSS")  # today's date
    println(datestamp)
    println("number of threads = $( Threads.nthreads() )")

    pathmain= joinpath(path, datestamp)

    mainfolder= mkdir(pathmain)    # creates a folder named today's date

    patht= joinpath(pathmain,"data\\")

    mainfolder1= mkdir(patht)

    folders=  multipledir(patht, ICS) 

    # Info printing on shell and file
    infos = @sprintf "Box shape: %s\nNumber of particles = %i\nNumber density = %s μm⁻²\nR=%.1f μm \nT = %.1f (K)\nv=%s (μm/s) \nω=%s (rad/s)\nCharacteristic lengths: (a=%.1f b=%.1f) μm\npf=%s\nIntegration step: dt=%.0e s \nSimulation downsampling: %i\nNumber of steps: Nt=%.1e\nTotal simulated time T_tot = %.2e s\n\nInteraction function: %s with parameters: %s\nRange: %.1f μm\nOffcenter: %s" box_shape Np density R T v ω a b packing_fraction δt measevery  Nt T_tot int_func int_params intrange offcenter*(2*forward-1)

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
        "dt" => δt,
        "measevery" => measevery,
        "Nt" => Nt,
        "T_tot" => T_tot,
        "int_func" => string(int_func),
        "int_params" => int_params,
        "intrange" => intrange,
        "offcenter" => offcenter*(2*forward-1),
    )


    JLD2.save(joinpath(mainfolder, "siminfo_dict.jld2"), info_dict)


    for i in 1:ICS
        pathf= joinpath(patht, "run$i\\")

        filename= "$datestamp"*"_run$i"
        pathf= pathf*filename

        datafname = pathf*".txt"
        # polarfname = pathf*"_polarization.txt"

        # Simulation and file storage
        open(datafname, "w") do infile
            writedlm(infile, ["N" "Time" "xpos" "ypos" "orientation"], ",")
        end

        start = now()
        @info "$(start) Started simulation #$i"


        ABPE, matrices = initABPE( Np, L, R, T, v, ω, int_func, forward, offcenter, intrange, int_params...,)
        if int_func == lennard_jones
            offcenter_nosuperpose!(ABPE, δt, forward, offcenter, 2R+1e-3, int_func, int_params...)
        end
        for nt in 0:Nt
            if nt % measevery == 0
                pnumber = collect(1:Np)
                time = fill(nt, Np)
                #creating Data
                data = DataFrame(
                    N= pnumber,
                    Time= time,
                    xpos= ABPE.x,
                    ypos= ABPE.y,
                    orientation=ABPE.θ,
                )  
                CSV.write(datafname, data, append = true)
                # open(polarfname, "a") do polfile
                #     write(polfile, "$(mean_polarization(ABPE.θ))\n")
                # end
            end
            ABPE =update_heun(ABPE,matrices,δt, forward, offcenter, intrange, int_func, int_params...)
            if nt % (Nt÷100) == 0
                elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
                print("$((100*nt÷Nt))%... Step $nt, total elapsed time $(elapsed)\r")
            end
        end
        @info "$(now()) Simulation and file writing finished"
        if animation
            try
                animation_from_file(pathf,L,R,δt,measevery,25, show = false, record=true, final_format = "mp4", color_code_dir = true)
            catch
                println("No vid")
            end
        end
    end
end

