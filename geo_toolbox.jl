using Distances
using GeometryBasics
using LinearAlgebra
using Markdown

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
function pirotation!(θ::Array{Float64})
	θ.%=2π
	θ[θ.<0].+=2π
    return nothing
end
"""
Calculate the 2-norm of a vector.
"""
d2(xy::Array{Float64,2}) = sqrt.(sum(abs2, xy, dims=2))
d2(xy) = sqrt.(sum(abs2, xy, dims=2))


"""
    ellipse_sectors(a::Float64,b::Float64,r::Float64,shell::Bool)

Compute the area of the interception between an ellipse and a circle with radius r>b as the sum of all the sectors given by the interceptions.
a larger semiaxis, b smaller semiaxis, r circle radius, shell: if function is used to compute shells area; if true, functions returns area of the ellipse when r = a.
```
"""
function ellipse_sectors(a::Float64,b::Float64,r::Float64,shell::Bool) # area of intersection between ellipse and circle with r >b
    r<=a ? r : throw(DomainError(r, "r must be smaller than the larger semiaxis or no intersection will exist."))
    if shell && isapprox(r,a)
        return pi*a*b
    end

    x = a*sqrt((r^2-b^2)/(a^2-b^2))
    y = b*sqrt((a^2-r^2)/(a^2-b^2))

    θ₁ = atan(y,x)
    θ₂ = atan(y,-x)
    Aₑ = 0.5a*b*atan(a*tan(θ₁)/b)
    Aᵢₙₜ = pi*a*b - 4*Aₑ
    Aₒ = 2*θ₁*r^2

    return Aᵢₙₜ+Aₒ
end

"""
    intersection_area_el(a::Float64,b::Float64,r::Float64,shell::Bool)

Compute the area of the shell between R₁ and R₂ inside the ellipse
```
"""
function intersection_area_el(a::Float64,b::Float64,R₁::Float64,R₂::Float64) #a larger semiaxis, b smaller semiaxis, R₁ smaller circle radius
    if R₂ <= b
        area = pi*(R₂^2-R₁^2)
    end

    if R₁ < b && R₂ > b
        area = ellipse_sector(a,b,R₂,true)-(pi*R₁*R₁)
    end

    if R₁ >= b
        area = ellipse_sector(a,b,R₂,true) - ellipse_sector(a,b,R₁,true)
    end

    return area
end

"""
    intersection_area_sq(r::Float64,L::Float64)

Compute the area of the intersection between a circle of radius r and a square of side L
```
"""
function intersection_area_sq(r::Float64, L::Float64)
    l = L/2

    if r<=l 
        A = pi*r^2
    else
        x = sqrt(r^2-l^2)
        y = l 
        θ = atan(y,x)-atan(x,y)

        A = 4*x*l + 2*θ*r^2
    end

    return A
end

"""
Still developing
"""
function get_angles(xyθ::Array{Float64,2})
	diffs = radial_directions(xyθ)
	angs_c2c = atan.(diffs[2], diffs[1])
	angs_c2c[angs_c2c .< 0].+=π
	return (xyθ[:,3].%2π).-angs_c2c
	#now you can combine an angle and its opposite through cat in 3d or in some other way depending on what you need
	#in [i,j] you have the angle that the i-th particle's direction forms with the direction between j and i
end

#Functions for Voronoi tessellation

"""
    xy_to_points(xy::Array{Float64,2})

Transforms a 2D array of coordinates into an array of Point2 objects.

# Examples
```julia-repl
julia> xy = [0.4 1.2; π ℯ]
2×2 Matrix{Float64}:
 0.4      1.2
 3.14159  2.71828

julia> xy_to_points(xy)
2-element Vector{Point{2, Float64}}:
 [0.4, 1.2]
 [3.141592653589793, 2.718281828459045]
```
"""
function xy_to_points(xy::Array{Float64,2})
    return [Point2(row) for row in eachrow(xy)]
end

"""
    points_to_xy(points::Vector{Point{2, Float64}})

Transforms a vector of 2D points into an 2D vector of coordinates.

# Examples
```julia-repl
julia> xy = [0.4 1.2; π ℯ]
2×2 Matrix{Float64}:
 0.4      1.2
 3.14159  2.71828

julia> xy_points = xy_to_points(xy)
2-element Vector{Point{2, Float64}}:
 [0.4, 1.2]
 [3.141592653589793, 2.718281828459045]

julia> points_to_xy(xy_points)
2×2 Matrix{Float64}:
 0.4      1.2
 3.14159  2.71828
```
"""
function points_to_xy(points::Vector{Point{2, Float64}})
    xy = zeros(Float64, length(points), 2)
    for (i, p) in enumerate(points)
        xy[i, 1] = p[1]
        xy[i, 2] = p[2]
    end
    return xy
end

"""
    round_point(p::Point{2,Float64}, tol=1e-8)

Transforms a 2D array of coordinates into an array of Point2 objects.

# Examples
```julia-repl
julia> xy = [0.4 1.2; π ℯ]
2×2 Matrix{Float64}:
 0.4      1.2
 3.14159  2.71828

julia> xyp = xy_to_points(xy)
2-element Vector{Point{2, Float64}}:
 [0.4, 1.2]
 [3.141592653589793, 2.718281828459045]

julia> round_point.(xyp)
2-element Vector{Point{2, Float64}}:
 [0.4, 1.2]
 [3.14159265, 2.71828183]
```
"""
function round_point(p::Point{2,Float64}, tol=1e-8)
    return round.(p, digits=ceil(Int, -log10(tol)))
end

"""
    get_adjacency_from_points(cells::Vector{Vector{Point{2,Float64}}}; tol=1e-8)

Returns a set of tuples representing the adjacency of Voronoi cells with particles at their centers.
"""
function get_adjacency_from_points(cells::Vector{Vector{Point{2,Float64}}}, Np::Int; tol=1e-8)
    # Step 1: Normalize and index all unique points
    point_to_index = Dict{Point{2,Float64}, Int}()
    index_counter = 1

    normalized_cells = Vector{Vector{Int}}(undef, length(cells))

    for (i, cell) in enumerate(cells)
        normalized_cell = Int[]
        for p in cell
            rp = Point(round_point(p, tol))
            if !haskey(point_to_index, rp)
                point_to_index[rp] = index_counter
                index_counter += 1
            end
            push!(normalized_cell, point_to_index[rp])
        end
        normalized_cells[i] = normalized_cell
    end

    # Step 2: Build edge -> cell map
    edge_map = Dict{Set{Int}, Vector{Int}}()

    for (cell_id, verts) in enumerate(normalized_cells)
        n = length(verts)
        for i in 1:n
            v1 = verts[i]
            v2 = verts[mod1(i+1, n)]  # wrap around
            edge = Set([v1, v2])
            edge_map[edge] = get(edge_map, edge, Int[])
            push!(edge_map[edge], cell_id)
        end
    end

    # Remove edges between two points if both points are outside of the simulation box
    # This is done by checking if both points are greater than the center length

    for (key, points) in collect(edge_map)
        if all(points .> Np)
            delete!(edge_map, key)
        end
    end

    # Step 3: Extract adjacency
    adjacency = Set{Tuple{Int, Int}}()
    for ids in values(edge_map)
        if length(ids) == 2
            i, j = sort(ids)
            push!(adjacency, (i, j))
        end
    end

    return adjacency
end

"""
    find_boundary_points(cells::Vector{Vector{Point{2,Float64}}})

Finds the boundary points of Voronoi cells that are outside the simulation box and returns their periodic projections.
"""
function find_boundary_points(xy::Array{Float64, 2}, cells)
    xy_periodic_projection = Array{Float64, 2}(undef, 0, 2)
    proj_inds = Array{Int, 1}(undef, 0)
    for (i, cell) in enumerate(points_to_xy.(cells))
        if any(abs.(cell[:,1]) .== L/2)
            xy_periodic_projection = vcat(xy_periodic_projection, [xy[i,1] - sign(xy[i,1])*L xy[i,2]])
            push!(proj_inds, i)
            end
        if any(abs.(cell[:,2]) .== L/2)
            xy_periodic_projection = vcat(xy_periodic_projection, [xy[i,1] xy[i,2] - sign(xy[i,2])*L])
            push!(proj_inds, i)
        end
    end
    return xy_periodic_projection, proj_inds
end
