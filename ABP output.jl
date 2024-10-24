# PURPOSE: Output of ABP main 
# all codes in repository are complied in this function
#VARIABLES: Destination folder path and filename

include("ABP main.jl")
include("ABP file.jl")
include("ABP analysis.jl")
include("ABP SD.jl")
include("ABP multifolder.jl")
# include("ABP average.jl")
include("generation.jl")
using Plots,Distances,NaNStatistics,CSV, DataFrames
using Dates
gr()

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# THIS IS THE CODE TO CALL MAIN FUNCTION
# We plot the set of particles considering the correction of hard spheres

L = 100.0 	# μm box length
R = 2.0	# μm particle radius
v = 10.0 	# μm/s particle velocity
a=L/2
b=L/4
ICS=1
   # number of intial conditons to be scanned 
#pf_factor = (R^2)/(a*b)
pf_factor = (R^2)
DT, DR = diffusion_coeff(R).*[1e12, 1]
packing_fraction = 0.2

Np = round(Int,packing_fraction*a*b/(R^2))  #Np is the number of particles inside the ellipse
#π
Nt = 1000000# Nt is the number of steps 
δt = 1.0e-4 #L/(v*Nt) # δt is the time step
#println(" Number of particles: $Np") 
#-------------------------------------------------------------------------------------------------------------------

# destination folders selection
path= raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\2024\10.October\ellipse\\" # destination folder path

datestamp=Dates.format(now(),"YYYYmmdd-HHMMSS")  # todays date

mainfolder= mkdir(path*"$datestamp")    # creates a folder names todays'late

path1= path*"$datestamp\\"

mainfolder1= mkdir(path1*"R=$R v=$v a=$a b=$b pf=$packing_fraction")

patht= path*"$datestamp\\R=$R v=$v a=$a b=$b pf=$packing_fraction\\"


folders=  multipledir(patht,ICS) 

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

println("multiparticleE_wall compiled\n")
#---------------------------------------------------------------------------------------------------------------------
# file store
file_store(graph_wall,Nt,pathf)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
# analysis
inside_Np=stat_analysis1(a,b,R,pathf,2)
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

anim = @animate for i = 1:100:Nt
    scatter(graph_wall[1][i][:,1], graph_wall[1][i][:,2], aspect_ratio=:equal, lims=(-L/2, L/2),markersize=350R/L,marker =:circle,legend=false, title = "$(inside_Np) particles, steps $i, ellipse a=L/2, b= L/4")
    plot!(a*cos.(-π:0.01:π), b*sin.(-π:0.01:π)) # ellipse
    quiver!(graph_wall[1][i][:,1],graph_wall[1][i][:,2],quiver=(4*cos.(graph_wall[2][i,1]),4*sin.(graph_wall[2][i,1])), color=:red)
end
#marker_z=graph_wall[2][i,1], color=:rainbow, for 

f1= pathf*".gif"
gif(anim, f1)

    finish = time()
    println("Time taken for simulation run$i: $(finish-start) seconds")

end
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# AVERAGE OF THE MULTIPLE OUTPUT FILES DATA
# mainfolder="C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-205011\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\"
# (average(mainfolder))   # passing path of the main folders which has all the runs

# mainfolder= raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07 Coding\2024\10.October\ellipse\20241016-151502\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\run1\\"
# filename="20241016-151502 R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1 run1"

# inside_Np=stat_analysis1(a,b,R,mainfolder*filename,0) # 0 for pole, ewuator, 1 for only right left, 2 for entire
