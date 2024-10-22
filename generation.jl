# Generating uniform distribution of particles inside an ellipse 

function genellipse(Nt,a,b,R)
    xy = (2 .*rand(Nt, 2).-1).*repeat([a b], Nt)

    r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
    rₚ = sqrt.(r)   
    θ =atan.(xy[:,2], xy[:,1]) 

    r_max= (a-R)*(b-R)./(sqrt.((((a-R)*sin.(θ)).^2) .+ ((b-R)*cos.((θ))).^2))
    id = (rₚ .< (r_max))
    xyacc = [xy[id,1] xy[id,2]]
    Nacc = size(xyacc,1)
    return xyacc, Nacc
end

function iterative_gen(Np, a, b, R)
    ##First generation
    xyacc, Nacc = genellipse(Np,a,b,R)
    while Nacc < Np
        numgen = Np-Nacc
        xyadd = genellipse(numgen,a,b,R)[1]
        xyacc = vcat(xyacc, xyadd)
        Nacc = size(xyacc,1)
    end
    return xyacc, Nacc
end

# xy = iterative_gen(100000, 100, 25, 2.0)[1]  # number of particles inside the boundary while Np is total number of particles
# xyθ = [xy 2π*rand(100000)]
# scatter(xyθ[:,1], xyθ[:,2], markersize=1, color="red")