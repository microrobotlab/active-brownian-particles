using Distances, LinearAlgebra, Markdown

"""
    pwdist(x::Vector{Float64}[, y::Vector{Float64}])
Compute the pairwise difference matrix between 'x' and 'y'.

pwdist(x) is a short method for pwdist(x,x).

# Examples
```julia-repl
julia> pwdist([1.,2.],[3.,4.])
2×2 Matrix{Float64}:
 2.0  3.0
 1.0  2.0

julia> pwdist([1.,2.])
2×2 Matrix{Float64}:
  0.0  1.0
 -1.0  0.0
```
"""
pwdist(x::Vector{Float64}) =[b-a for a in x, b in x]

pwdist(x::Vector{Float64}, y::Vector{Float64}) = [b-a for a in x, b in y] 

"""
    radial_directions(xy::Array{Float64,2})

Compute the cos and sin matrix of the angles between any two positions of xy.

# Examples
```julia-repl
julia> xy = rand(3,2)
3×2 Matrix{Float64}:
 0.689357   0.36113
 0.0486193  0.780761
 0.923794   0.989613

julia> radial_directions(xy)[1]
3×3 Matrix{Float64}:
  0.0       -0.836558  0.349497
  0.836558   0.0       0.972687
 -0.349497  -0.972687  0.0

julia> radial_directions(xy)[2]
3×3 Matrix{Float64}:
  0.0        0.547878  0.936938
 -0.547878   0.0       0.232123
 -0.936938  -0.232123  0.0
```
"""
function radial_directions(xy::Array{Float64,2})
    dist = pairwise(Euclidean(), xy, dims = 1)
    dist[diagind(dist)].=eps()
    diff_x = pwdist(xy[:,1])./dist
    diff_y = pwdist(xy[:,2])./dist
    return diff_x, diff_y 
end

"""
    angdiff(ϕ1::Float64, ϕ2::Float64)
    angdiff(ϕ1::Array{Float64,1}, ϕ2::Array{Float64,1})

Compute the angle between two particle directions.

# Examples
```julia-repl
julia> rad2deg(angdiff(deg2rad(45.), deg2rad(30.)))
15.0
julia> rad2deg(angdiff(deg2rad(45.), deg2rad(315.)))
90.0
```
"""
function angdiff(ϕ1::Float64, ϕ2::Float64)
	Δϕ = abs(ϕ1-ϕ2) % 2π
	(Δϕ < 2π - Δϕ) ? diff = Δϕ : diff = 2π - Δϕ
	return diff
end

function angdiff(ϕ1::Array{Float64,1}, ϕ2::Array{Float64,1})
	Δϕ = pwdist(ϕ1, ϕ2) .% 2π
	pirotation!(Δϕ)
	Δϕ[Δϕ .> 2π .- Δϕ] .= 2π .- Δϕ[Δϕ .> 2π .- Δϕ]
	return Δϕ
end

"""
    pirotation(θ::Array{Float64})

Shifts all the angles from any domain to [0, 2π]

# Examples
```julia-repl
julia> rad2deg(pirotation(deg2rad(-45)))
315.0
```
"""
function pirotation(θ::Real)
	θr = θ%2π
	θr<0 ? θr+2π : θr
end

function pirotation(θ::Array{Float64})
	θr = θ.%2π
	θr[θr.<0].+=2π
	return θr
end

"""
    pirotation!(θ::Array{Float64})

Shifts all the angles from any domain to [0, 2π], reassigning.
```
"""
function pirotation!(θ::Float64)
    θ%=2π
	if θ < 0
        θ+=2π
    end
    return nothing
end

function pirotation!(θ::Array{Float64})
	θ.%=2π
	θ[θ.<0].+=2π
    return nothing
end
