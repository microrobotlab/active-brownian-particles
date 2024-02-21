# this is a rough code for any calculations
# not to be taken seriously
using CalculusWithJulia, Plots
L=100
a= L/2
b= L/4

θ1= [ θ for θ in -π:0.01:π]  #in radian

angle= [atan(b/a)] # angle at which the pole area= equator area 

# curvature in  polar co ordinates
κ= [a*b./((a*a*sin.(θ)*sin.(θ))+(b*b*cos.(θ)*cos.(θ)))^(3/2) for θ in -π:0.01:π]

κ_angle= [a*b./((a*a*sin.(θ)*sin.(θ))+(b*b*cos.(θ)*cos.(θ)))^(3/2) for θ in angle]



j= plot(L/2*cos.(-π:0.01:π), L/4*sin.(-π:0.01:π),marker_z=κ , color=:rainbow,aspect_ratio=:equal,legend=true) 



#k=scatter(θ1, κ, marker_z= κ, color=:rainbow,legend=true)
display(j)

#display(scatter(θ1,κ))

#display(scatter!(angle,κ_angle))
#k= 8*a*b/((3*a*a+b*b)^(3/2))


k1= 8*8/(L*((13)^(3/2)))
