# PURPOSE: Output of ABP main 
# all codes in repository are complied in this function
#VARIABLES: Destination folder path and filename
include("ABP main interactions opt.jl")
# include("ABP main.jl")
include("ABP_file.jl")
include("ABP radialdistribution.jl")
# include("ABP SD.jl")
include("ABP multifolder.jl")
# include("ABP radialdensity.jl")
include("ABP plot_animation.jl")
# include("ABP average.jl")
using CSV, DataFrames, Dates, Distances, Distributions, JLD2, Logging, NaNStatistics, Printf, Random