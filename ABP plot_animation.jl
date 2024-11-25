using GLMakie
using GeometryBasics: Point2f, Circle, Point2f0, Polygon
using CSV, DataFrames, Dates, Logging

function animation_from_file(pathf::String, L::Float64, R::Float64, timestep::Float64, measevery::Int, downsampling::Int=1; ext::String=".txt", record::Bool=false, final_format::String="gif")
    fname = pathf*ext
    timestep *= measevery*downsampling
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

    fig = GLMakie.Figure()
    ax = GLMakie.Axis(fig[1,1], limits = (-L/2, L/2, -L/2, L/2), aspect = 1)
    mrk = GLMakie.decompose(Point2f,Circle(Point2f0(0), 2R))
    sc = GLMakie.scatter!(ax, xs, ys, marker = Polygon(mrk),markersize = 200/L)
    ar = GLMakie.arrows!(ax, xs, ys, us, vs, color = :black, lengthscale=R, arrowsize = 300R/L)

    timestamps = 1:downsampling:(maximum(df[!,:Time])÷downsampling)

    slider = GLMakie.Slider(fig[2, 1], range=1:size(xpos, 2), startvalue=1)

    on(slider.value) do val
        simstep[] = round(Int, val)
    end
    
    if record
        GLMakie.record(fig, pathf*".$final_format", timestamps;
        framerate = Int(1/timestep)) do t
            simstep[] = t
        end
        @info "$(now()) Animation saved to $pathf.$final_format"
    end
    return fig

end

function animation_from_graph(graph_wall, pathf, L::Float64, R::Float64, Np::Int, timestep::Float64,Nt::Int,measevery::Int, downsampling::Int=1; record::Bool=false, final_format::String="gif")
    pos = vcat(graph_wall[1]...)
    orient = vcat(graph_wall[2]...)

    timestep *= measevery*downsampling

    xpos = reshape(pos[:,1], Np,:)[:,1:downsampling:end]
    ypos = reshape(pos[:,2], Np,:)[:,1:downsampling:end]
    θ = reshape(orient, Np,:)[:,1:downsampling:end]
    u,v = cos.(θ), sin.(θ)

    simstep = Observable(1)
    xs = @lift xpos[:,$simstep]
    ys = @lift ypos[:,$simstep]
    us = @lift u[:,$simstep]
    vs = @lift v[:,$simstep]

    fig = GLMakie.Figure()
    ax = GLMakie.Axis(fig[1,1], limits = (-L/2, L/2, -L/2, L/2), aspect = 1)
    mrk = GLMakie.decompose(Point2f,Circle(Point2f0(0), 2R))
    sc = GLMakie.scatter!(ax, xs, ys, marker = Polygon(mrk),markersize = 200/L)
    ar = GLMakie.arrows!(ax, xs, ys, us, vs, color = :black, lengthscale=R, arrowsize = 300R/L)

    timestamps = 1:downsampling:(Nt÷downsampling)

    slider = GLMakie.Slider(fig[2, 1], range=1:size(xpos, 2), startvalue=1)

    on(slider.value) do val
        simstep[] = round(Int, val)
    end
    
    if record
        GLMakie.record(fig, pathf*".$final_format", timestamps;
        framerate = Int(1/timestep)) do t
            simstep[] = t
        end
        @info "$(now()) Animation saved to $pathf.$final_format"
    end
    return fig

end
