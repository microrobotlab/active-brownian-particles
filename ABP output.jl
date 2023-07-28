# PURPOSE: Output of ABP main.jl. All codes in repository are complied in this function
# INPUT: Destination folder path and ICS
# ICS is the number of Initial conditions scan given same parameters
# VARIABLES: N, L, R, V, a,b, packing_fraction, ICS
include("ABP main.jl")
include("ABP file.jl")
include("ABP analysis.jl")
include("ABP multifolder.jl")

using Plots,Distances,NaNStatistics,CSV, DataFrames, Dates

gr()

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# PARAMETERS SET
L = 100.0 	# μm box length
R = 2.0		# μm particle radius
v = 10.0 	# μm/s particle velocity
a = L/2     # μm major axis
b = L/4     # μm minor axis
Nt = 100000  # Nt is the number of steps 
ICS = 2     # Number of intitial condiitons 
pf_factor = (R^2)
DT, DR = diffusion_coeff(R).*[1e12, 1]
packing_fraction = 0.1
Np = round(Int,packing_fraction*L^2/(2R^2))  #Np is the number of particles 

#-------------------------------------------------------------------------------------------------------------------
# DESTINATION FOLDERS
path="C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\July\\clean\\" # destination folder path

datestamp=Dates.format(now(),"YYYYmmdd-HHMMSS")  # todays date

mainfolder= mkdir(path*"$datestamp")    

path1= path*"$datestamp\\"

mainfolder1= mkdir(path1*"R=$R v=$v a=$a b=$b")

patht= path*"$datestamp\\R=$R v=$v a=$a b=$b\\"

folders=  multipledir(patht,ICS) 

for i=1:ICS

    pathf= patht*"run$i\\"
    filename= "$datestamp R=$R v=$v a=$a b=$b run$i"
    pathf= pathf*filename
# -----------------------------------------------------------------------------------------------------------------------------------------------------

# CALL MAIN FUNCTION FOR OPEN BOUNDARY CONDITION
# PARTICLES INTERACT VIA HARD SPHERE CORRECTION
#=
graph = multiparticleE(Np,L,R,v,Nt);

#println("multiparticleE complied\n")
=#

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
# CALL MAIN FUNCTION FOR CLOSED BOUNDARY CONDITION

# PARTICLES INTERACT VIA HARD SPHERE CORRECTION AND BOUNCE OFF AFTER REACHING THE BOUNDARY

graph_wall = multiparticleE_wall(Np,L,R,v,Nt) # has values of x and y posiiton at each time step in graph_wall[1]

println("multiparticleE_wall complied\n")
#---------------------------------------------------------------------------------------------------------------------
# file store
file_store(graph_wall,Nt,pathf)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
# ANALYSIS AND PLOTING OF PARTICLES
inside_Np=stat_analysis1(a,b,R,pathf)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ANIMATION OF DYNAMICS
anim = @animate for i = 1:100:Nt
    scatter(graph_wall[1][i][:,1], graph_wall[1][i][:,2], aspect_ratio=:equal, lims=(-L/2, L/2),markersize=350R/L,marker =:circle,legend=false, title = "$(inside_Np) particles, steps $i, ellipse a=L/2, b= L/4")
    plot!(L/2*cos.(-π:0.01:π), L/4*sin.(-π:0.01:π))
    quiver!(graph_wall[1][i][:,1],graph_wall[1][i][:,2],quiver=(4*cos.(graph_wall[2][i,1]),4*sin.(graph_wall[2][i,1])), color=:red)
end
f1= pathf*".gif"
gif(anim, f1)

end



