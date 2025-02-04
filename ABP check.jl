using CSV, DataFrames

path="D:\\nic_simulations\\lj_offcenter"

datestamp = "20250124-102030"

df_list = []
i = 0
for p in readdir(joinpath(path,datestamp,"data"), join = true)
    global i+=1
    i>1 && break
    if isdir(p)
        num = p[91:end]
        df = CSV.read(joinpath(p, datestamp*"_run$num.txt"), DataFrame)
        push!(df_list, df)
    end
end

L= 175.
R = 2.


println(any(reduce(vcat,[abs.(df.xpos) .> L/2.0+R for df in df_list])))
println(any(reduce(vcat,[abs.(df.ypos) .> L/2.0+R for df in df_list])))

println(maximum(maximum([abs.(df.xpos) for df in df_list])),"\n", maximum(maximum([abs.(df.ypos) for df in df_list])) )