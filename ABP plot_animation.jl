using GLMakie, CairoMakie
using GeometryBasics: Point2f, Circle, Point2f0, Polygon
using CSV, DataFrames, Dates, Logging
include("geo_toolbox.jl")
GLMakie.activate!(inline=false)

function animation_from_file(pathf::String, L::Float64, R::Float64, timestep::Float64, measevery::Int, output_framerate::Int=25; ext::String=".txt", show::Bool=true, record::Bool=false, final_format::String="gif", color_code_dir::Bool=false)
    GLMakie.activate!(inline=false)

    fname = pathf*ext
    timestep *= measevery
    if output_framerate > 1/timestep
        output_framerate = Int(1/timestep)
        downsampling = 1
    else
        downsampling = Int(1/(timestep*output_framerate))
    end
    df = CSV.read(fname, DataFrame)    
    Np = maximum(df[!,:N])
    xpos = Array(df[!,:xpos])
    ypos = Array(df[!,:ypos])
    θ = Array(df[!,:orientation])

    xpos = reshape(xpos, Np,:)[:,1:downsampling:end]
    ypos = reshape(ypos, Np,:)[:,1:downsampling:end]
    θ = reshape(θ, Np,:)[:,1:downsampling:end]
    pirotation!(θ)

    simstep = Observable(1)
    xs = @lift xpos[:,$simstep]
    ys = @lift ypos[:,$simstep]
    θs = @lift θ[:,$simstep]

    fig = GLMakie.Figure(size = (1080,1080))
    ax = GLMakie.Axis(fig[1,1], limits = (-L/2, L/2, -L/2, L/2), aspect = 1)
    mrk = GLMakie.decompose(Point2f,Circle(Point2f0(0), 1.))
    mrkdir = Polygon(Point2f[(0.6,0), (-0.4,0.4), (-0.4, -0.4)])
    if color_code_dir
        sc = GLMakie.scatter!(ax, xs, ys, marker = Polygon(mrk), color = θs, colormap = :cyclic_mygbm_30_95_c78_n256_s25, alpha = 0.5, markerspace = :data, markersize = R)
    else
        sc = GLMakie.scatter!(ax, xs, ys, marker = Polygon(mrk), color = :slategrey, markerspace = :data, markersize = R)
    end
    GLMakie.scatter!(ax, xs, ys, marker = mrkdir, rotation = θs,  color = :black, markerspace = :data, markersize = 3.)

    if record
        timestamps = 1:(size(xpos,2))
        GLMakie.record(fig, pathf*".$final_format", timestamps;
        framerate = output_framerate) do t
            simstep[] = t
        end
        @info "$(now()) Animation saved to $pathf.$final_format"
    end

    if show
        slider = GLMakie.Slider(fig[2, 1], range=1:size(xpos, 2), startvalue=1, color_active =:grey12, color_inactive = :grey60, color_active_dimmed = :grey30)

        on(slider.value) do val
            simstep[] = round(Int, val)
        end
        GLMakie.display(fig)
    end
    return nothing

end

function animation_from_history(history, pathf, L::Float64, R::Float64, Np::Int, timestep::Float64,Nt::Int,measevery::Int, downsampling::Int=1; show::Bool=true, record::Bool=false, final_format::String="gif", color_code_dir::Bool=false)
    GLMakie.activate!(inline=false)

    pos = vcat(history[1]...)
    orient = vcat(history[2]...)

    timestep *= measevery*downsampling

    xpos = reshape(pos[:,1], Np,:)[:,1:downsampling:end]
    ypos = reshape(pos[:,2], Np,:)[:,1:downsampling:end]
    θ = reshape(orient, Np,:)[:,1:downsampling:end]
    pirotation!(θ)

    simstep = Observable(1)
    xs = @lift xpos[:,$simstep]
    ys = @lift ypos[:,$simstep]
    θs = @lift θ[:,$simstep]

    fig = GLMakie.Figure(size = (600,600))
    ax = GLMakie.Axis(fig[1,1], limits = (-L/2, L/2, -L/2, L/2), aspect = 1)
    mrk = GLMakie.decompose(Point2f,Circle(Point2f0(0), 1))
    mrkdir = Polygon(Point2f[(0.6,0), (-0.4,0.4), (-0.4, -0.4)])
    if color_code_dir
        sc = GLMakie.scatter!(ax, xs, ys, marker = Polygon(mrk), color = θs, colormap = :cyclic_mygbm_30_95_c78_n256_s25, alpha = 0.5, markerspace = :data, markersize = R)
    else
        sc = GLMakie.scatter!(ax, xs, ys, marker = Polygon(mrk), color = :slategrey, markerspace = :data, markersize = R)
    end
    GLMakie.scatter!(ax, xs, ys, marker = mrkdir, rotation = θs,  color = :black, markerspace = :data, markersize = 3.)

    if record
        timestamps = 1:(size(xpos,2))
        GLMakie.record(fig, pathf*".$final_format", timestamps;
        framerate = Int(1/timestep)) do t
            simstep[] = t
        end
        @info "$(now()) Animation saved to $pathf.$final_format"
    end

    if show
        slider = GLMakie.Slider(fig[2, 1], range=1:size(xpos, 2), startvalue=1, color_active =:grey12, color_inactive = :grey60, color_active_dimmed = :grey30)

        on(slider.value) do val
            simstep[] = round(Int, val)
        end
        GLMakie.display(fig)
    end
    return nothing

end

function plot_one_timestep(df::DataFrame, R::Number, L::Number, timestep::Int; savepath = nothing, title= nothing)
    CairoMakie.activate!(inline=true)

    Np = maximum(df[!,:N])
    xpos = Array(df[!,:xpos])
    ypos = Array(df[!,:ypos])
    θ = Array(df[!,:orientation])

    xpos = reshape(xpos, Np,:)
    ypos = reshape(ypos, Np,:)
    θ = reshape(θ, Np,:)
    pirotation!(θ)

    fig = CairoMakie.Figure(size = (1080,1080))
    if title !== nothing
        ax = CairoMakie.Axis(fig[1,1], limits = (-L/2, L/2, -L/2, L/2), aspect = 1, title = title)
    else
        ax = CairoMakie.Axis(fig[1,1], limits = (-L/2, L/2, -L/2, L/2), aspect = 1)
    end

    mrk = CairoMakie.decompose(Point2f,Circle(Point2f0(0), 1))
    sc = CairoMakie.scatter!(ax, xpos[:,timestep], ypos[:,timestep], marker = Polygon(mrk), color = θ[:,timestep], colormap = :cyclic_mygbm_30_95_c78_n256_s25, alpha = 0.5, markerspace = :data, markersize = R)
    mrkdir = Polygon(Point2f[(0.6,0), (-0.4,0.4), (-0.4, -0.4)])
    CairoMakie.scatter!(ax, xpos[:,timestep], ypos[:,timestep], marker = mrkdir, rotation = θ[:,timestep],  color = :black, markerspace = :data, markersize = 3.)

    CairoMakie.display(fig)
    if savepath !== nothing
        CairoMakie.save(savepath, fig)
    end
end
