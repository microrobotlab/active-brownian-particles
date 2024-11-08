# FILE WRITING
# PURPOSE:Store the output of main code ABP output.jl
# METHOD: Dataframes are used 
# INPUT: array of postions and orientation
# OUTPUT: saved excel file (f1)

using CSV, DataFrames, DelimitedFiles, Markdown, Parquet

"""
    file_store_csv(graph_wall,Nt,pathf;downsampling = 100)
    Stores the content of graphw_wall in a .csv file
"""

function file_store_csv(graph_wall,Nt,pathf;downsampling::Int=100)
    f1= pathf*".csv"               # destination file name

    pnumber=[]
    time=[]
    x=[]
    y=[]
    θ=[]
    fx=[]
    fy=[]
    torque=[]
    for i = 1:downsampling:Nt
        for j = 1:length(graph_wall[1][i][:,1])
            push!(pnumber,j)
            push!(time, i)
        push!(x, graph_wall[1][i][j,1])
        push!(y, graph_wall[1][i][j,2])
        push!(θ, graph_wall[2][i,1][j])
        push!(fx, graph_wall[3][i][j,1])
        push!(fy, graph_wall[3][i][j,2])
        push!(torque, graph_wall[4][i,1][j])
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
    orientation=θ,
    fx=fx,
    fy=fy,
    torque=torque,
    )  

    CSV.write(f1, data)
    
    close(efg)
    @info "$(now()) Out of ABP file"
    return nothing
end

"""
    file_store_txt(graph_wall,Nt,pathf;downsampling = 100)
    Stores the content of graphw_wall in a .txt file.
"""

function file_store_txt(graph_wall,Nt,pathf;downsampling::Int=100)
    f1= pathf*".txt"              # destination file name

    pnumber=[]
    time=[]
    x=[]
    y=[]
    θ=[]
    fx=[]
    fy=[]
    torque=[]
    for i = 1:downsampling:Nt
        for j = 1:length(graph_wall[1][i][:,1])
            push!(pnumber,j)
            push!(time, i)
        push!(x, graph_wall[1][i][j,1])
        push!(y, graph_wall[1][i][j,2])
        push!(θ, graph_wall[2][i,1][j])
        push!(fx, graph_wall[3][i][j,1])
        push!(fy, graph_wall[3][i][j,2])
        push!(torque, graph_wall[4][i,1][j])
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
    orientation=θ,
    fx=fx,
    fy=fy,
    torque=torque,
    )  

    CSV.write(f1, data)
    
    @info "$(now()) Out of ABP file"
    return nothing
end

function file_store_parquet(graph_wall,Nt,pathf;downsampling::Int=100)
    f1= pathf*".parquet"               # destination file name

    pnumber=[]
    time=[]
    x=[]
    y=[]
    θ=[]
    fx=[]
    fy=[]
    torque=[]
    for i = 1:downsampling:Nt
        for j = 1:length(graph_wall[1][i][:,1])
            push!(pnumber,j)
            push!(time, i)
        push!(x, graph_wall[1][i][j,1])
        push!(y, graph_wall[1][i][j,2])
        push!(θ, graph_wall[2][i,1][j])
        push!(fx, graph_wall[3][i][j,1])
        push!(fy, graph_wall[3][i][j,2])
        push!(torque, graph_wall[4][i,1][j])
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
    orientation=θ,
    fx=fx,
    fy=fy,
    torque=torque,
    ) 

    Parquet.write_parquet(f1, data)
    
    close(efg)
    @info "$(now()) Out of ABP file"
    return nothing
end