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

## USER INTERFACE
# destination folders selection
path="C:\\Users\\nikko\\OneDrive\\Documents\\Uni\\magistrale\\tesi\\simulations\\" # destination folder path

## PARAMETERS SET
# Simulation parameters
Nt = 100000       # number of steps
Delta_t = 1e-3   # s step time
ICS=1            # Number of intial conditons to be scanned 

# Geometrical parameters
box_shape = :square    # shapes: :square, :circle, :ellipse
L = 300.0 	           # μm box length
R = 2.0		           # μm particle radius
v = 10. 	           # μm/s particle velocity
packing_fraction = 0.1 # Largest pf for spherical beads π/4 = 0.7853981633974483

#-------------------------------------------------------------------------------------------------------------------

## PRELIMINARY CALCULATIONS
DT, DR, γ = diffusion_coeff(R).*[1e12, 1, 1] # Translational and Rotational diffusion coefficients, drag coefficient

if box_shape == :square
    Np = round(Int,packing_fraction*L^2/(pi*R^2))
    density = Np/L^2
    @printf "Box shape: %s\nNumber of particles = %i\nNumber density = %.4f μm⁻²" box_shape Np density

elseif box_shape == :ellipse
    a=L/2 # Larger semiaxis
    b=L/4 # Smaller semiaxis
    area_el = π*a*b
    Np = round(Int,packing_fraction*area_el/(pi*R^2))
    density = Np/area_el
    @printf "Box shape: %s\nNumber of particles = %i\nNumber density = %.4f μm⁻²" box_shape Np density
end


