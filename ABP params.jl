using JLD2, Distributions, Statistics

include("ABP main interactions heun.jl")

simfolder = "..\\simulations\\tests"
datestamp = "20250611-172658"
patht = joinpath(simfolder, datestamp)
infodict = JLD2.load(joinpath(patht, "siminfo_dict.jld2"))

T = infodict["T"] # Temperature (K)
R = infodict["R"] # Particles' radius (μm)
v = infodict["v"] # Particles' velocity (μm/s or distribution)
σ = 2R

if isa(v, Real)
    vm = v
elseif isa(v, Distribution) || isa(v, Array)
    vm = mean(v)
end

diffusion_coeff(1e-6R,T,1e-3)

τ = 1/diffusion_coeff(1e-6R,T,1e-3)[2] # Persistence time 1/Dᵣ (s)
Pe = vm*τ/σ # Following definition in 10.1103/PhysRevLett.130.148202
