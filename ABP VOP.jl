# this funcion will calculate the velocity and velocity polarization of particles , i.e., velocity order parameter
# Date: 29-05-2023
# Method: Absoulte velocity calculation at each instant
# Input: File from the output.jl, which has positions and orientation at each instant
# Output: Velocity of each particle at every instant and average velocity of the ensemble

using BenchmarkTools, Plots, LaTeXStrings, Statistics, CSV, DataFrames,CategoricalArrays
 

path="C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\July\\07\\v5.0\\"  # destination folder path

filename="data_ellipse"   # filename base for all data

pathf= path*filename 
    f= pathf*".csv"
    f1= pathf*"_v.csv"
    f2= pathf*"_vop.png"
   
    df= CSV.read(f, DataFrame)
    #steps= df[!, :N]
    time= df[!,:Time]
    x= df[!,:xpos]
    y= df[!,:ypos]
    df[!,:N] = categorical(df[!,:N],compress=true) # it sorts out time step data 
    ## Group dataframe by values in categorical column
    gdf = groupby(df,:N,sort=true) # only 1000 data groups because I have omitted 100 time steps means 1 s
   vel = [mean(sqrt.((diff(g[!,:xpos])).^2+(diff(g[!,:ypos]).^2))) for g in gdf]  # dr vector and velocity magnitude at each time second
   vel_x = [mean(diff(g[!, :xpos])) for g in gdf] 
   vel_y= [mean(diff(g[!,:ypos])) for g in gdf]


    vp= mean((vel_x./vel).+ (vel_y./vel))         # polarization vector 


   #plot(vel[100],vel[1000])
