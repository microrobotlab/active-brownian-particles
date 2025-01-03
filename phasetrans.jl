using DelimitedFiles, Plots
include("stat_toolbox.jl")

folder = joinpath("..", "simulations", "flocking2")

filevec = readdir(folder)
oc=Float64[]
pol_list = Vector{Float64}[]
for filename in filevec
    from = findfirst("oc", filename)[end] +1
    to = findfirst("_p", filename)[1] -1
    offc = filename[from:to]
    push!(oc, parse(Float64, offc))
    pv = readdlm(joinpath(folder, filename))
    push!(pol_list, vec(pv))
end
sort_ind = sortperm(oc)
oc = oc[sort_ind]
pol_list = pol_list[sort_ind]

sizes = collect(1:10:10000)
p = plot()
for pol in pol_list[1:5]
    stds = [samp_std(jackknife(pol, s)) for s in sizes]
    plot!(p, sizes, stds)
end
display(p)

pjack = zeros(size(pol_list))
perrs = zeros(size(pol_list))
χ = zeros(size(pol_list))
χerrs = zeros(size(pol_list))
therm = 20000
for (i,pol) in enumerate(pol_list)
    jked = jackknife(pol[therm:end], 2500)
    jked2 = jackknife(pol[therm:end].^2, 2500)
    pjack[i] = mean(jked)
    perrs[i] = samp_std(jked)
    susc = length(jked)*(jked2.-jked.^2)
    χ[i] = mean(susc)
    χerrs[i] = samp_std(susc)
end

display(scatter(oc, pjack, yerr = perrs))
display(scatter(oc, χ, yerr = χerrs))

# plot(pol_list[2])
