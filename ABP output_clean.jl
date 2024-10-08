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

## PARAMETERS SET
# Simulation parameters
N = 100000       # number of steps
Delta_t = 1e-3   # s step time
ICS=1            # Number of intial conditons to be scanned 

# Geometrical parameters

L = 300.0 	           # μm box length
R = 2.0		           # μm particle radius
v = 10. 	           # μm/s particle velocity
packing_fraction = 0.1 # Largest pf for spherical beads π/4 = 0.7853981633974483



pf_factor = (R^2)
DT, DR, γ = diffusion_coeff(R).*[1e12, 1, 1]  #1

area_el = π*a*b
Np = round(Int,packing_fraction*L^2/(pi*R^2))
density = Np/area_el
#π
Nt = N# Nt is the number of steps 
#println(" Number of particles: $Np") 

t_tot= N*Delta_t
