include("geo_toolbox.jl")

## First way. Problem: center biased. Advantage: fast, simple and correct

function circgen(Np, a, b, R)
    α = rand(Np).*2π
    r_max= (a-R)*(b-R)./(sqrt.((((a-R)*sin.(α)).^2) .+ ((b-R)*cos.((α))).^2))
    r = rand(Np).*r_max
    x = r.*cos.(α)
    y = r.*sin.(α)
    return hcat(x,y)
end

## Second way. Presumably slower

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
    # println(Nacc)
    while Nacc < Np
        numgen = Np-Nacc
        xyadd = genellipse(numgen,a,b,R)[1]
        xyacc = vcat(xyacc, xyadd)
        Nacc = size(xyacc,1)
    end
    # println(Nacc)
    # x = xyacc[:,1]
    # y = xyacc[:,2]
    return xyacc, Nacc
end