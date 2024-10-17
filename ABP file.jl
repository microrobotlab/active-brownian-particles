# FILE WRITING
# PURPOSE:Store the output of main code ABP output.jl
# METHOD: Dataframes are used 
# INPUT: array of postions and orientation
# OUTPUT: saved excel file (f1)

using CSV, DataFrames, DelimitedFiles

function file_store_csv(graph_wall,Nt,pathf;downsampling::Int=100)
    f1= pathf*".csv"               # destination file name

pnumber=[]
time=[]
x=[]
y=[]
θ=[]
for i = 1:downsampling:Nt+1
    for j = 1:length(graph_wall[1][i][:,1])
        push!(pnumber,j)
        push!(time, i)
    push!(x, graph_wall[1][i][j,1])
    push!(y, graph_wall[1][i][j,2])
    push!(θ, graph_wall[2][i,1][j])
    end
end

touch(f1)

    efg = open(f1, "w")
    
    #creating DataFrame
    data = DataFrame(
    N= pnumber,
    Time= time,
    xpos= x,
    ypos= y,
    orientation=θ) 

    CSV.write(f1, data)
    
    close(efg)
    @info "$(now()) Out of ABP file"
    return nothing
end

function file_store_txt(graph_wall,Nt,Np,pathf;downsampling::Int=100) # This is cuter
    f1= pathf*".txt"               # destination file name

    Nt +=1
    pnumber = Int.(repeat(1:Np, Nt÷downsampling))
    time = Int.(repeat(1:Nt, inner = [Np]))
    xy = vcat(graph_wall[1]...)
    x = xy[:,1]
    y = xy[:,2]
    θ = vcat(graph_wall[2]...)

    touch(f1)
    
    #creating DataFrame
    data = DataFrame(
    N = pnumber,
    Time = time,
    xpos = x,
    ypos = y,
    orientation =θ) 

    CSV.write(f1, data, delim = ",")
    
    @info "$(now()) Out of ABP file"
    return nothing
end

function file_store_txt(graph_wall,Nt,pathf;downsampling::Int=100) #... But this is faster 😢
    f1= pathf*".txt"              # destination file name

    pnumber=[]
    time=[]
    x=[]
    y=[]
    θ=[]
    for i = 1:downsampling:Nt+1
        for j = 1:length(graph_wall[1][i][:,1])
            push!(pnumber,j)
            push!(time, i)
        push!(x, graph_wall[1][i][j,1])
        push!(y, graph_wall[1][i][j,2])
        push!(θ, graph_wall[2][i,1][j])
        end
    end

    touch(f1)
    
    #creating DataFrame
    data = DataFrame(
    N= pnumber,
    Time= time,
    xpos= x,
    ypos= y,
    orientation=θ) 

    CSV.write(f1, data)
    
    @info "$(now()) Out of ABP file"
    return nothing
end