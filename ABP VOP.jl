# this funcion will calculate the velocity and velocity polarization of particles , i.e., velocity order parameter
# Date: 29-05-2023
# Method: Absoulte velocity calculation at each instant
# Input: File from the output.jl, which has positions and orientation at each instant
# Output: Velocity of each particle at every instant and average velocity of the ensemble

using  Plots, LaTeXStrings, Statistics, CSV, DataFrames,CategoricalArrays
 

# destination folder path

path= raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07 Coding\2024\09.September\circle\20240926-105431\R=1.5 v=0.6 a=35.0 b=35.0 pf=0.01\run1"
filename="\\test file"   # filename base for all data

pathf= path*filename 
    f= pathf*".csv"
    f1= pathf*"_v.csv"
    f2= pathf*"_vop.png"
   
    df= CSV.read(f, DataFrame)
    #steps= df[!, :N]
    time= df[!,:Time]
    x= df[!,:xpos]
    y= df[!,:ypos]
    df[!,:Time] = categorical(df[!,:Time],compress=true) # it sorts out time step data 
    ## Group dataframe by values in categorical column
    gdf = groupby(df,:Time,sort=true) # only 1000 data groups because I have omitted 100 time steps means 1 s
    #xy = [x y]
   vel = [(sqrt.((diff(g[!,:xpos])).^2+(diff(g[!,:ypos]).^2))) for g in gdf]  # dr vector and velocity magnitude at each time second
   vel_x = [(diff(g[!, :xpos])) for g in gdf] 
   vel_y= [(diff(g[!,:ypos])) for g in gdf]

    
    vp= mean((vel_x./vel).+ (vel_y./vel))         # polarization vector 
    xgdf= [g[!,:xpos] for g in gdf]  # one data column

    ygdf= [g[!,:ypos] for g in gdf]  # one data column
       
  
   xy= ([g[!, [:xpos, :ypos]] for g in gdf]) # x and y position for g in gdf]
   
  xyc= vcat(xy...)
 
   xym=Matrix(xyc[:,[:xpos,:ypos]])
    distances = (pairwise(Euclidean(), xym[1:8,:], xym[1:8,:],dims =1)) # distance matrix

    up_dis= triu(distances, 1)
   plot(time_actual,first.(vel))

   plot(time_actual,getindex.(vel,3))

   scatter(time_actual[2:end],distances[5,200:299])
time_actual=[first(g[!,:Time]) for g in gdf]

j=[30.117 -5.32911]
y=[-30.4604 -5.96872]

rjy= sqrt((j[1]-y[1])^2+(j[2]-y[2])^2)
