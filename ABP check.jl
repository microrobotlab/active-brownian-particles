using CSV, DataFrames

path="C:\\Users\\nikko\\OneDrive\\Documents\\Uni\\magistrale\\tesi\\simulations\\"

datestamp = "20241125-174109"

df_list = []
for p in readdir(joinpath(path,datestamp,"data"), join = true)
    if isdir(p)
        num = p[91:end]
        df = CSV.read(joinpath(p, datestamp*"_run$num.txt"), DataFrame)
        push!(df_list, df)
    end
end

L= 25.
R = 2.


println(any(reduce(vcat,[abs.(df.xpos) .> L/2.0+R for df in df_list])))
println(any(reduce(vcat,[abs.(df.ypos) .> L/2.0+R for df in df_list])))

println(maximum(maximum([abs.(df.xpos) for df in df_list])),"\n", maximum(maximum([abs.(df.ypos) for df in df_list])) )