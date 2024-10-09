using CSV
using DataFrames
using FFTW
using FindPeaks1D
using JSON3
using Plots
using Statistics

rdf = CSV.read("..\\simulations\\20241009-094744\\data\\run1\\20241009-094744_run1_rdf.csv", DataFrame)
bindata = CSV.read("..\\simulations\\20241009-094744\\data\\run1\\20241009-094744_run1_rdfbin.csv", DataFrame)
rs = bindata[!,:Radius]

rdf.RadialDistributionFunction = [JSON3.read(x) for x in rdf.RadialDistributionFunction]
begin
    ngroups = 100
    lendf = size(rdf)[1]
    div = Int(lendf//ngroups)
    grouplabel = [(i÷div) for i in 0:lendf-1]
    rdf[!,:GroupLabel] = grouplabel

    rdf_avg = transform(groupby(rdf,:GroupLabel), :RadialDistributionFunction => ((x)->(mean(x),)) => :RDF_mean)
    transform!(rdf_avg, :RDF_mean => ByRow(x -> x[1]) => :RDF_mean)
end

begin
    p = plot()
    start,end_ = 0,20
    xticks!(p, start:1:end_)
    xlims!(p,start,end_)
    times = Int[]
    peaks = Float64[]
    for i in unique(div .* grouplabel .+1)[1:10:end]
        averaged_rdf = rdf_avg[i,:RDF_mean]
        pkindices, properties = findpeaks1d(averaged_rdf;
                                            height = 2., 
                                            prominence = 0.2,)
        push!(peaks, properties["peak_heights"][1])
        push!(times, rdf_avg[i,:Time])
        plot!(p, rs, averaged_rdf)
        scatter!(p, rs[pkindices], averaged_rdf[pkindices])
    end
    # p2 = plot()
    # scatter(p2, times, peaks)
end
display(p)
println(argmin(peaks))