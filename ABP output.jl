# PURPOSE: Output of ABP main 
# all codes in repository are complied in this function
#VARIABLES: Destination folder path and filename
@show Threads.nthreads()

include("ABP main.jl")
include("ABP_file.jl")
include("ABP analysis.jl")
include("ABP SD.jl")
include("ABP multifolder.jl")
# include("ABP average.jl")
include("generation.jl")
#include("ABP freq analysis.jl")
include("ABP gif.jl")
using Plots,Distances,NaNStatistics,CSV, DataFrames
using Dates
gr()

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# THIS IS THE CODE TO CALL MAIN FUNCTION
# We plot the set of particles considering the correction of hard spheres

L = 100.0 	# μm box length
R = 1.5	# μm particle radius
v = 5.0 	# μm/s particle velocity
a=L/2
b=L/8
ICS=1
   # number of intial conditons to be scanned 
#pf_factor = (R^2)/(a*b)
pf_factor = (R^2)
DT, DR = diffusion_coeff(R).*[1e12, 1]
packing_fraction = 0.1

Np = round(Int,packing_fraction*a*b/(R^2))  #Np is the number of particles inside the ellipse
#π
Nt = 20000000# Nt is the number of steps 
resample=1000
Nt_store= Int(Nt/resample)  # time steps at which data has to be stored, not the actual simulation time step
δt = 1.0e-3 #L/(v*Nt) # δt is the time step
#println(" Number of particles: $Np") 
#-------------------------------------------------------------------------------------------------------------------

# destination folders selection
path= raw"D:\j.sharma\P07\workstationMRL\2025\03.March\\" # destination folder path

datestamp=Dates.format(now(),"YYYYmmdd-HHMMSS")  # todays date

mainfolder= mkdir(path*"$datestamp")    # creates a folder names todays'late

path1= path*"$datestamp\\"

mainfolder1= mkdir(path1*"R=$R v=$v a=$a b=$b pf=$packing_fraction")

patht= path*"$datestamp\\R=$R v=$v a=$a b=$b pf=$packing_fraction\\"


folders=  multipledir(patht,ICS) 
@show start_sim= time()
for i=1:ICS
  start = time()
    pathf= patht*"\\run$i\\"
    filename= "\\$datestamp R=$R v=$v a=$a b=$b pf=$packing_fraction run$i"
    pathf= pathf*filename
    println("simulation=$i")
    #graph = multiparticleE(Np,L,R,v,Nt);    # function to simulation particles with open boundary 

     # println("multiparticleE complied\n")
  #-------------------------------------------------------------------------------------------------------------------------------------------------------------
# THIS IS THE CODE TO CALL WALL FUNCTIONS IN THE MAIN FUNCTION for square wall condition

# Same with the wall condition (particles bounce off the edge)  
#graph_wall = multiparticleE_wall(Np,L,R,v,Nt) 

#=
scatter(graph[1][:,1], graph[1][:,2], markersize=350R/L, legend=false, aspect_ratio=:equal, title = "$Np particles, step n°1")
plot!([L/2], seriestype="vline")
plot!([-L/2], seriestype="vline")

scatter(graph[20][:,1], graph[20][:,2], markersize=350R/L, legend=false, aspect_ratio=:equal, title = "$Np particles, step n°$Nt" )
plot!([L/2], seriestype="vline")
plot!([-L/2], seriestype="vline")
=#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
# THIS IS THE CODE TO CALL WALL FUNCTIONS IN THE MAIN FUNCTION for elliptical wall condition
graph_wall = multiparticleE_wall(Np,L,R,v,Nt,δt) # has values of x and y posiiton in each frame in graph_wall[1]
#  for i in 1:Nt_store
#     k=plot!(graph_wall[1][i][:,1], graph_wall[1][i][:,2], aspect_ratio=:equal, lims=(-L/2, L/2),markersize=350R/L,marker =:circle,legend=false, title = "$(Np) particles, steps $i, ")
#  end  #square
# display(k)
println("multiparticleE_wall compiled\n")
#---------------------------------------------------------------------------------------------------------------------
# file store  
file_store_csv(graph_wall,Nt_store,pathf)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
# analysis
#  inside_Np=stat_analysis1(a,b,R,pathf,δt,2)
# mean and standard deviation


#analysis_SD= stat_analysis2(a,b,R,pathf)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# making animation


#------------------------------------------------------------------------------For square-------------------------------------------------------------------------
#=
anim = @animate for i = 1:100:Nt
    scatter(graph_wall[1][i][:,1], graph_wall[1][i][:,2], aspect_ratio=:equal, lims=(-L/2, L/2),markersize=350R/L,marker =:circle,legend=false, title = "$(Np) particles, steps $i, ")
    
    plot!([L/2], seriestype="vline")  #square
    plot!([-L/2], seriestype="vline")
    quiver!(graph_wall[1][i][:,1],graph_wall[1][i][:,2],quiver=(4*cos.(graph_wall[2][i,1]),4*sin.(graph_wall[2][i,1])), color=:red)
end
#marker_z=graph_wall[2][i,1], color=:rainbow, for 

f1= pathf*".gif"
gif(anim, f1)
end
=#
#------------------------------------------------------------------------------for ellipse-------------------------------------------------------------------------

# anim = @animate for i = 1:Nt_store
#     scatter(graph_wall[1][i][:,1], graph_wall[1][i][:,2], aspect_ratio=:equal, lims=(-L/2, L/2),markersize=350R/L,marker =:circle,legend=false, title = "$Np particles, steps $(i*resample), ellipse a=L/2, b= L/4")
#     plot!(a*cos.(-π:0.01:π), b*sin.(-π:0.01:π)) # ellipse
#     quiver!(graph_wall[1][i][:,1],graph_wall[1][i][:,2],quiver=(4*cos.(graph_wall[2][i,1]),4*sin.(graph_wall[2][i,1])), color=:red)
# end
#marker_z=graph_wall[2][i,1], color=:rainbow, for 

# f1= pathf*".gif" # gif(anim, f1)
    finish = time()
    println("Time taken for simulation run$i: $(round((finish-start)/60.0, digits=3)) minutes")

end
@show time_end = time()
@show time_end-start_sim
# println("Total time taken for all runs: $(((time_end-start_sim), digits=3)) seconds")
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# AVERAGE OF THE MULTIPLE OUTPUT FILES DATA
# mainfolder="C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-205011\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\"
# (average(mainfolder))   # passing path of the main folders which has all the run Sant'Anna\P07Coding\2024\11.November\20241128-183816\R=2.0 v=15.0 a=50.0 b=25.0 pf=0.2

 mainfolder1= "D:\\j.sharma\\P07\\workstationMRL\\2025\\03.March\\20250324-152734\\R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1\\run1\\"
 
   filename="20250324-152734 R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1 run1\\"
  t= mainfolder1*filename
  
  f1000= joinpath(mainfolder1,filename*".csv") 

#    f1= "D:\\j.sharma\\P07\\workstationMRL\\20241104-121620\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.2\\run1\\20241104-121620 R=2.0 v=10.0 a=50.0 b=25.0 pf=0.2 run1_p.csv\\"
   df= CSV.read(t*".csv",DataFrame) 
#    FFT_analysis(t,δt)

#  inside_Np=stat_analysis1(a,b,R,t,δt,2) # 0 for pole, equator, 1 for only right left, 2 for entire

  gif(mainfolder1,t)

