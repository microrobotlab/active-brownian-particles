using CSV, DataFrames, FFTW, JSON3, LinearAlgebra, Plots, Random, Statistics

rdf = CSV.read("..\\simulations\\20241008-160813\\data\\run1\\20241008-160813_run1_rdf.csv", DataFrame)

rdf.RadialDistributionFunction = [JSON3.read(x) for x in rdf.RadialDistributionFunction]
begin
    lendf = size(rdf)[1]
    div = Int(lendf//groupsize)
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
    for i in unique(10 .* grouplabel .+1)
        plot!(p, rs, rdf_avg[i,:RDF_mean])
    end
end
display(p)