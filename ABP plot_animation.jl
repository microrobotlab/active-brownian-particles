using GLMakie, GeometryBasics
using CSV, DataFrames, DelimitedFiles,LinearAlgebra, Markdown, Parquet

function animation_from_file(pathf::String, L::Float64, R::Float64, timestep::Float64, downsampling::Int, ext::String=".txt"; record::Bool=false)
    fname = pathf*ext
    timestep *= downsampling
    df = CSV.read(fname, DataFrame)    
    Np = maximum(df[!,:N])
    xpos = Array(df[!,:xpos])
    ypos = Array(df[!,:ypos])
    θ = Array(df[!,:orientation])

    xpos = reshape(xpos, Np,:)[:,1:downsampling:end]
    ypos = reshape(ypos, Np,:)[:,1:downsampling:end]
    θ = reshape(θ, Np,:)[:,1:downsampling:end]
    u,v = cos.(θ), sin.(θ)

    simstep = Observable(1)

    xs = @lift xpos[:,$simstep]
    ys = @lift ypos[:,$simstep]
    us = @lift u[:,$simstep]
    vs = @lift v[:,$simstep]

    fig = Figure()
    ax = Axis(fig[1,1], limits = (-L/2, L/2, -L/2, L/2), aspect = 1)
    mrk = decompose(Point2f,Circle(Point2f0(0), 2R))
    sc = scatter!(ax, xs, ys, marker = Polygon(mrk),markersize = 200/L)
    ar = arrows!(ax, xs, ys, us, vs, color = :black, lengthscale=R, arrowsize = 300R/L)

    timestamps = 1:maximum(df[!,:Time])÷downsampling

    slider = Slider(fig[2, 1], range=1:size(xpos, 2), startvalue=1)

    on(slider.value) do val
        simstep[] = round(Int, val)
    end
    fig
    if record
        record(fig, pathf*".gif", timestamps;
        framerate = Int(1/timestep)) do t
            simstep[] = t
        end
    end

end
