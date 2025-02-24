using DelimitedFiles, CairoMakie, JLD2, LaTeXStrings, Colors
include("stat_toolbox.jl")

Np = 250

datestamps = ["20250213-120012" "20250110-111152" "20250213-120105"]

fig = Figure(size = (600,800), fontsize = 20)
polax = Axis(fig[1,1], xlabel = "α", ylabel = "Polarization", xgridvisible = false, ygridvisible = false)
suscax = Axis(fig[2,1], xlabel = "α", ylabel = "Polarization Susceptivity", xgridvisible=false,ygridvisible=false)
colors = distinguishable_colors(length(datestamps), parse(Colorant, "firebrick2"))

for (i,datestamp) in enumerate(datestamps)
    folder = "D:\\NiccoloP\\simulations\\flocking2\\$(datestamp)\\data"
    info_dict = JLD2.load("D:\\NiccoloP\\simulations\\flocking2\\$(datestamp)\\siminfo_dict.jld2")
    v = info_dict["v"]

    filevec = readdir(folder)
    oc=Float64[]
    pol_list = Vector{Float64}[]
    for foldername in filevec
        from = findfirst("oc", foldername)[end] +1
        offc = foldername[from:end]
        push!(oc, parse(Float64, offc))
        polarfile = filter(name -> occursin("polarization", name), readdir(joinpath(folder, foldername)))
        println(polarfile)
        pv = readdlm(joinpath(folder, foldername, polarfile[1]))
        push!(pol_list, vec(pv))
    end
    sort_ind = sortperm(oc)
    oc = oc[sort_ind]
    pol_list = pol_list[sort_ind]

    therm = 20000
    blocksize = 2500
    if isapprox(v, 10)
        thermfig = Figure(fontsize = 20)
        ax = Axis(thermfig[1,1], xlabel = "Simulation step", ylabel = "Polarization", xgridvisible = false, ygridvisible = false)
        j = 7
        colstherm = distinguishable_colors(5, parse(Colorant, "turquoise2"))
        for pol in pol_list[7:end]
            lines!(ax, pol, alpha = 0.75, label = "α = $(oc[i])", color = colstherm[j-6])
            j+=1
        end
        vlines!(ax, 20000, ymin = 0, ymax=1., color = :black, linewidth = 2., linestyle = :dash)
        axislegend(ax, position = :rb)
        display(thermfig)
        save("D:\\NiccoloP\\simulations\\flocking2\\imgs\\thermalization_v10.pdf", thermfig)
    

        sizes = collect(1:10:10000)
        jkfig = Figure()
        axjk = Axis(jkfig[1,1], xlabel = "Block Size [steps]", ylabel = "Sample Standard Deviation", xgridvisible=false,ygridvisible=false)
        for (j,pol) in enumerate(pol_list)
            j !==7 && continue
            stds = [samp_std(jackknife(pol[therm:end], s)) for s in sizes]
            stdsusc = [samp_std(Np*(jackknife(pol[therm:end].^2,s).-(jackknife(pol[therm:end], s)).^2)) for s in sizes]
            lines!(axjk, sizes, stds, linewidth = 2., label = "Polarization")
            lines!(axjk, sizes, stdsusc, linewidth = 2., label = "Susceptivity")    
        end
        vlines!(axjk, 2500, color = :black, linewidth = 2., linestyle = :dash)
        axislegend(axjk, position = :lt)
        display(jkfig)
        save("D:\\NiccoloP\\simulations\\flocking2\\imgs\\blocksize_v10.pdf", jkfig)

    end

    pjack = zeros(size(pol_list))
    perrs = zeros(size(pol_list))
    χ = zeros(size(pol_list))
    χerrs = zeros(size(pol_list))
    for (i,pol) in enumerate(pol_list)
        jked = jackknife(pol[therm:end], blocksize)
        jked2 = jackknife(pol[therm:end].^2, blocksize)
        pjack[i] = mean(jked)
        perrs[i] = samp_std(jked)
        susc = Np*(jked2.-jked.^2)
        χ[i] = mean(susc)
        χerrs[i] = samp_std(susc)
    end
    scatterlines!(polax, oc, pjack, linewidth = 0.5, color = colors[i])
    errorbars!(polax, oc, pjack, perrs, whiskerwidth = 4, linewidth = 2, label = "v = $v", color = colors[i])
    scatterlines!(suscax, oc, χ, linewidth = 0.5, color = colors[i])
    errorbars!(suscax, oc, χ, χerrs, whiskerwidth = 4, linewidth = 2, label = "v = $v", color = colors[i])
end

axislegend(polax, position = :lt)
axislegend(suscax, position = :lt)
display(fig)
save("D:\\NiccoloP\\simulations\\flocking2\\imgs\\pol_susc.pdf", fig)
