# this is a code for curvature calculations

using CalculusWithJulia, Plots, Unitful
L=100
a= L/2
b= L/4
path="C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2024\\"
f= path*"curvature1.png"
θ1= [ θ for θ in -π:0.01:π]  #in radian

angle= [atan(b/a)] # angle at which the pole area= equator area 

# curvature in  polar co ordinates
κ= [a*b./((a*a*sin.(θ)*sin.(θ))+(b*b*cos.(θ)*cos.(θ)))^(3/2) for θ in -π:0.01:π]

κ_angle= [a*b./((a*a*sin.(θ)*sin.(θ))+(b*b*cos.(θ)*cos.(θ)))^(3/2) for θ in angle]

kx= [a*cos.(θ) for θ in -π:0.01:π]

ky= [b*sin.(θ) for θ in -π:0.01:π]

k=plot(kx,ky, marker_z= κ,color = :lajolla,seriestype=:scatter,label="",aspect_ratio=:equal, xlabel="x[μm]", ylabel="y[μm]")
display(k)
savefig(k,f)

#lajolla