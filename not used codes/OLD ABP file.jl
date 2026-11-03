# FILE WRITING
# PURPOSE:Store the output of main code ABP output.jl
# METHOD: Dataframes are used 
# INPUT: array of postions and orientation
# OUTPUT: saved excel file (f1)

using CSV, DataFrames

function old_file_store(graph_wall,Nt,pathf)
    f1= pathf*".csv"               # destination file name

pnumber=[]
time=[]
x=[]
y=[]
θ=[]
for i = 1:1:Nt
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
  
    println("I am out of ABP file")
    return nothing
end