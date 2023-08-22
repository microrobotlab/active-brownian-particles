


function attractive_interactions!(xy::Array{Float64,2}, R::Float64)

    ϵ=0.1
    σ= 2*R
    Np = size(xy,1)
    dists = zeros(Np,Np)
    #superpose = falses(Np,Np)
    #uptriang = falses(Np,Np)
    uptriang = trues(Np,Np)
    for i = 1:Np
        
        #uptriang[i,i+1:Np] .= true
        uptriang[i,i] = false
    end
    dists .= pairwise(Euclidean(),xy,xy,dims=1)
    dists13= dists.^(13)

    dists7= dists.^(7)

    force= 24*ϵ*(((2*σ^(13))./dists13).- (σ^(7)./dists7))

    force= force.* uptriang

    return  sum(force, dims=2)

end

xy= [[1.0 2.0 ] ;[1.0 1.0 ]; [1.5 2.0] ]
#xy = [[p[1],p[2]] for p in eachrow(x)]

AT= attractive_interactions!(xy, 2.0)

#AT1= sum(AT, dims=2)
