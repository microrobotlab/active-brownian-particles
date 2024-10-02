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

N = 1000
Delta_t = 1e-3
t_tot= N*Delta_t

# tauMax = t_tot/10               #is the actual maximum delta t over which I can calculate the MSD
# N_Max = Int64(tauMax/Delta_t)   # is the maximum number of frames of the camera. For me that I am simulating and that's it it will be given by the delta_t_MAX / delta_t_MIN on which I can calculate the MSD

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# THIS IS THE CODE TO CALL MAIN FUNCTION
# We plot the set of particles considering the correction of hard spheres

L = 300.0 	# μm box length
R = 2.0		# μm particle radius
v = 10. 	# μm/s particle velocity
a=L/2
b=L/4
ICS=1      # number of intial conditons to be scanned 
#pf_factor = (R^2)/(a*b)
pf_factor = (R^2)
DT, DR, γ = diffusion_coeff(R).*[1e12, 1, 1]  #1
packing_fraction = 0.1 #Largest pf for spherical beads π/4 = 0.7853981633974483

area_el = π*a*b
Np = round(Int,packing_fraction*L^2/(pi*R^2))
density = Np/area_el
#π
Nt = N# Nt is the number of steps 
#println(" Number of particles: $Np") 
#-------------------------------------------------------------------------------------------------------------------

# destination folders selection
path="C:\\Users\\nikko\\OneDrive\\Documents\\Uni\\magistrale\\tesi\\simulations\\trq" # destination folder path

datestamp=Dates.format(now(),"YYYYmmdd-HHMMSS")  # todays date
println(datestamp)

mainfolder= mkdir(path*"$datestamp")    # creates a folder names todays'late

path1= path*"$datestamp\\"

mainname = @sprintf "R=%.1f v=%.1f a=%.1f b=%.1f pf=%.4f dt=%.0e Nt=%.1e" R v a b packing_fraction Delta_t Nt

mainfolder1= mkdir(path1*mainname)

patht= joinpath([path*"$datestamp", mainname])

folders=  multipledir(patht,ICS) 

println("$mainname; total time simulated $(Delta_t*Nt); Np = $Np")

for i=1:ICS
    pathf= patht*"\\run$i\\"
    filename= "$datestamp "*mainname*" run$i"
    pathf= pathf*filename
    # println("simulation=$i")
    @info "$(now()) Started simulation #$i"
    # graph = multiparticleE(Np,L,R,v,Nt);    # function to simulation particles with open boundary 

    # println("multiparticleE complied\n")
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------
    # THIS IS THE CODE TO CALL WALL FUNCTIONS IN THE MAIN FUNCTION for square wall condition

    # Same with the wall condition (particles bounce off the edge)  
    # start_sim = now()
    # graph_wall = multiparticleE_wall(Np,L,R,v,Nt,Delta_t) 


    # scatter(graph[1][:,1], graph[1][:,2], markersize=350R/L, legend=false, aspect_ratio=:equal, title = "$Np particles, step n°1")
    # plot!([L/2], seriestype="vline")
    # plot!([-L/2], seriestype="vline")

    # scatter(graph[20][:,1], graph[20][:,2], markersize=350R/L, legend=false, aspect_ratio=:equal, title = "$Np particles, step n°$Nt" )
    # plot!([L/2], seriestype="vline")
    # plot!([-L/2], seriestype="vline")

    #-------------------------------------------------------------------------------------------------------------------------------------------------------------
    # THIS IS THE CODE TO CALL WALL FUNCTIONS IN THE MAIN FUNCTION for elliptical wall condition

    start_sim = now()
    graph_wall = multiparticleE(Np,L,R,Exponential(12/log(50)),Normal(0.,0.025),Nt,Delta_t) # has values of x and y position in each frame in graph_wall[1]
    elapsed_time = Dates.canonicalize(now()-start_sim)
    #sim_time = @elapsed graph_wall = multiparticleE(Np,L,R,v,Nt) # has values of x and y position in each frame in graph_wall[1]

    println("multiparticleE complied: elapsed time $elapsed_time\n")
    #---------------------------------------------------------------------------------------------------------------------
    # file store
    file_store(graph_wall,Nt,pathf)
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------
    # analysis
    # inside_Np=stat_analysis1(a,b,R,pathf)
    # mean and standard deviation

    #analysis_SD= stat_analysis2(a,b,R,pathf)

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #radial density
    # radial_density_profile(L,15, pathf, density)

    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # making animation


    #------------------------------------------------------------------------------For square-------------------------------------------------------------------------

    anim = @animate for i = 1:1000:Nt
        scatter(graph_wall[1][i][:,1], graph_wall[1][i][:,2], aspect_ratio=:equal, lims=(-L/2, L/2),markersize=350R/L,marker =:circle,legend=false, title = "$(Np) particles, steps $i, ",)
        
        plot!([L/2], seriestype="vline", color=:black)  #square
        plot!([-L/2], seriestype="vline", color=:black)
        plot!([L/2], seriestype="hline", color=:black)
        plot!([-L/2], seriestype="hline", color=:black)
        # quiver!(graph_wall[1][i][:,1],graph_wall[1][i][:,2],quiver=(4*cos.(graph_wall[2][i,1]),4*sin.(graph_wall[2][i,1])), color=:red)
    end
    #marker_z=graph_wall[2][i,1], color=:rainbow, for 

    f1= pathf*".gif"
    gif(anim, f1)
end

#------------------------------------------------------------------------------for ellipse-------------------------------------------------------------------------
# anim = @animate for i = 1:100:Nt
#     scatter(graph_wall[1][i][:,1], graph_wall[1][i][:,2], aspect_ratio=:equal, lims=(-L/2, L/2),markersize=350R/L,marker =:circle,legend=false, title = "$(inside_Np) particles, steps $i, ellipse a=L/2, b= L/4")
#     plot!(L/2*cos.(-π:0.01:π), L/4*sin.(-π:0.01:π)) # ellipse
#     quiver!(graph_wall[1][i][:,1],graph_wall[1][i][:,2],quiver=(4*cos.(graph_wall[2][i,1]),4*sin.(graph_wall[2][i,1])), color=:red)
# end
#marker_z=graph_wall[2][i,1], color=:rainbow, for 

# f1= pathf*".gif"
# gif(anim, f1)


# end

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# AVERAGE OF THE MULTIPLE OUTPUT FILES DATA
#mainfolder="C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-111614\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\"
# (average(patht))   # passing path of the main folders which has all the runs
