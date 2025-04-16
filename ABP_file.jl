# FILE WRITING
# PURPOSE:Store the output of main code ABP output.jl
# METHOD: Dataframes are used 
# INPUT: array of postions and orientation
# OUTPUT: saved excel file (f1)

using CSV, DataFrames, DelimitedFiles, Markdown

# """
#     file_store_csv(graph_wall,Nt,pathf;downsampling = 100)
#     Stores the content of graphw_wall in a .csv file
# """

function file_store_csv(graph_wall,Nt_store,pathf,resample)
    f1= pathf*".csv"               # destination file name

pnumber=[]
timestep=[]
x=[]
y=[]
θ=[]
for i = 1:Nt_store
    for j = 1:length(graph_wall[1][i][:,1])
        push!(pnumber,j)
        push!(timestep, i*resample)
    push!(x, graph_wall[1][i][j,1])
    push!(y, graph_wall[1][i][j,2])
    push!(θ, graph_wall[2][i,1][j])
    end
end

touch(f1)

    efg = open(f1, "w")
    
    #creating DataFrame
    data = DataFrame(
    Pn= pnumber,
    StepN= timestep,
    xpos= x,
    ypos= y,
    orientation=θ) 

    CSV.write(f1, data)
    
    close(efg)
    @info "$(now()) Out of ABP file"
    return nothing
end


# function file_store_txt(graph_wall,Nt,pathf;downsampling::Int=100)
#     f1= pathf*".txt"              # destination file name

#     pnumber=[]
#     time=[]
#     x=[]
#     y=[]
#     θ=[]
#     for i = 1:downsampling:Nt+1
#         for j = 1:length(graph_wall[1][i][:,1])
#             push!(pnumber,j)
#             push!(time, i)
#         push!(x, graph_wall[1][i][j,1])
#         push!(y, graph_wall[1][i][j,2])
#         push!(θ, graph_wall[2][i,1][j])
#         end
#     end

#     touch(f1)
    
#     #creating DataFrame
#     data = DataFrame(
#     N= pnumber,
#     Time= time,
#     xpos= x,
#     ypos= y,
#     orientation=θ) 

#     CSV.write(f1, data)
    
#     @info "$(now()) Out of ABP file"
#     return nothing
# end
