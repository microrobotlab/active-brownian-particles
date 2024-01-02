using Random
using  Plots, LaTeXStrings, Statistics, CSV, DataFrames,CategoricalArrays


function radial_density_sq(xy::Array{Float64,2}, nshells::Int64, L::Float64)
    thickness = .5*L/nshells
    lims = [i*thickness for i in 0:nshells-1]
    nums = zeros(nshells)
    for i in 1:nshells
        id = (abs.(xy[:,1]) .> lims[i]) .| (abs.(xy[:,2]) .> lims[i])
        nums[i] = sum(id)
    end
    area = (lims.+thickness).^2 .- lims.^2
    nc = [nums[i]-nums[i+1] for i in 1:length(nums)-1]
    nc = vcat(nc, nums[end])

    return nc./area, lims
end

function radial_density_cr(xy::Array{Float64,2}, nshells::Int64, R::Float64)
    thickness = 0.5R/nshells
    lims = [i*thickness for i in 0:nshells-1]

    r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
  
    rₚ = sqrt.(r)   
    # α =atan.(xy[:,2], xy[:,1]) 
    nums = zeros(nshells)

    for i in 1:nshells
        id = rₚ .> lims[i]
        nums[i] = sum(id)
    end
    area = π*((lims.+thickness).^2 .- lims.^2)
    nc = [nums[i]-nums[i+1] for i in 1:length(nums)-1]
    nc = vcat(nc, nums[end])

    return nc./area, lims
end

# function radial_density_el(xy::Array{Float64,2}, nshells::Int64, L::Float64)
#     r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
#     rₚ = sqrt.(r)   
#     α =atan.(xy[:,2], xy[:,1]) 

    
#     a = L/2 #These are the proportions used in the main code
#     b = L/4
#     rₑ = (a-R)*(b-R)./(sqrt.((((a-R)*sin.(α)).^2) .+ ((b-R)*cos.((α))).^2))
#     nums = zeros(nshells)

#     for i in 1:nshells
#         id = rₚ .> lims[i]
#         nums[i] = sum(id)
#     end
#     area = π*((lims.+thickness).^2 .- lims.^2)
#     nc = [nums[i]-nums[i+1] for i in 1:length(nums)-1]
#     nc = vcat(nc, nums[end])

#     return nc
# end

function radial_density_profile(R::Float64,nshells, pathf,ρtot)
    f= pathf*".csv"
    f1= pathf*"_d.csv"
    f2= pathf*"_dp.png"
    df= CSV.read(f, DataFrame)
    time= df[!,:Time]
    #x= df[!,:xpos]
    #y= df[!,:ypos]
    df[!,:Time] = categorical(df[!,:Time],compress=true) # it sorts out time step data 
    ## Group dataframe by values in categorical column
    gdf = groupby(df,:Time,sort=true) # only 1000 data groups because I have omitted 100 time steps means 1 s
    lastxy = gdf[end][!,[:xpos, :ypos]]
    xy = Matrix(lastxy)
    dens,l = radial_density_sq(xy,nshells,R)
    dens = dens/ρtot
    #image
    t1=plot();
    scatter!(t1, l,dens,legend = false)
    xlabel!("radius (μm)", xguidefont=font(16), xtickfont=font(11))
    ylabel!("ρ/Ρ", xguidefont=font(16), xtickfont=font(11))
    title!("Radial density profile")
    savefig(f2)

    #file wriiting
    touch(f1)

    efg = open(f1, "w")

    #creating DataFrame for number of particles at equators n1, and at poles n2
    data = DataFrame(Limits = l, Density= dens) 

    CSV.write(f1, data)


    println("I am out of ABP radialdensity")
    return dens
end