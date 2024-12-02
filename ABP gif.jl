# This function is to make the gif of ABP simulation
# you need to pass the file and path of the gif

function makegif(mainfolder,pathf)


    fmain= pathf*".csv"
    df= CSV.read(fmain, DataFrame)
    df[!,:StepN] = categorical(df[!,:StepN],compress=true) # it sorts out time step data 
  
    ## Group dataframe by values in categorical column
    gdf = groupby(df,:StepN,sort=true) # only 1000 data groups because I have omitted 100 time steps means 1 s
  
    anim = @animate for i = 1:10:length(gdf)
        scatter(gdf[i][!,:xpos], gdf[i][!,:ypos], aspect_ratio=:equal, lims=(-L/2, L/2),markersize=350R/L,marker =:circle,legend=false, title = "$Np particles, steps $(i*resample)")
        plot!(a*cos.(-π:0.01:π), b*sin.(-π:0.01:π)) # ellipse
        quiver!(gdf[i][!,:xpos],gdf[i][!,:ypos],quiver=(4*cos.(gdf[i][!,:orientation]),4*sin.(gdf[i][!,:orientation])), color=:red)
    end
    f1= pathf*"skip10.gif"
    gif(anim, f1)
    
end