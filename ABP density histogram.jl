using CSV, DataFrames, Plots, StatsBase
pyplot()    # use pyplot backend

# -------------------------------------------------------

function ellipse_density(csv_path::String, a::Real, b::Real, outfile::String;
                         bins::Int=60)

    df = CSV.read(csv_path, DataFrame)

    # particles inside the ellipse
    inside = @. (df.xpos / a)^2 + (df.ypos / b)^2 <= 1.0
    x = df.xpos[inside]
    y = df.ypos[inside]

   # bin centres (perfectly symmetric)
   
    xcent = range(-a, a, length=bins)
    ycent = range(-b, b, length=bins)

    # bin width (safe way)
    dx = (2a) / (bins - 1)
    dy = (2b) / (bins - 1)

    xedges = range(-a - dx/2, a + dx/2, length=bins+1)
    yedges = range(-b - dy/2, b + dy/2, length=bins+1)
     # use StatsBase.fit explicitly to avoid name conflicts
    h = StatsBase.fit(Histogram, (x, y), (xedges, yedges))
    dens = h.weights ./ maximum(h.weights .+ 1e-12)
    # mask outside the ellipse
      mask = [(xc/a)^2 + (yc/b)^2 > 1 for xc in xcent, yc in ycent]
      dens[mask'] .= NaN
 # custom colour map: off-white → blue → red
    cmap = cgrad([:white, :blue, :deepskyblue, :orangered, :red], [0.0, 0.15, 0.45, 0.7, 1.0])
    p = heatmap(xcent, ycent, dens',
                c=cmap,
                aspect_ratio=:equal,
                xlims=(-a*1.25, a*1.25),
                ylims=(-b*1.25, b*1.25),
                xlabel="x (μm)", ylabel="y (μm)",
                #title="Time-averaged density",
                colorbar_title="relative density",
                framestyle=:box,
                colorbar_ticks=0:0.5:1)

    # ellipse boundary
    θ = range(0, 2π, length=400)
    plot!(p, a .* cos.(θ), b .* sin.(θ),
          lw=2.5, lc=:black, label="")

    savefig(p, outfile)
    println("Saved → $outfile")
    return p
end
# ---------- your paths (edit only this block) ----------
parent   = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2026\10.October\final paper data\R1.5 pf 0.3\aspect 0.55\20261105-164611\R=1.5 v=5.0 a=33.67003367003367 b=18.51851851851852 pf=0.3\run6\\"
filename = "20261105-164611 R=1.5 v=5.0 a=33.67003367003367 b=18.51851851851852 pf=0.3 run6"

input  = joinpath(parent, filename * ".csv")
output = joinpath(parent, filename * "_density.png")

a = 33.67003367003367
b = 18.51851851851852
# run
ellipse_density(input, a, b, output)