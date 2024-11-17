# In this code, the time one particle is spending near the bounday is calculated
# and the total time spent in the boundary is also calculated
# Method: From the number of particles at the pole and eqautors, I am adding that number and gettin the time
# Time as output will be sun time. like this much time stend not a regration of time 

using CSV, DataFrames, Plots

data=DataFrame( run=Int64[], time_poles=Int64[], time_equators=Int64[])
for i in 1:90
runfolder= "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07Coding\\2024\\11.November\\20241112-142941\\R=2.0 v=5.0 a=50.0 b=25.0 pf=0.2\\run$i\\"
main_folder= joinpath(runfolder, "..\\")  # takes back to the previous folder
filename="20241112-142941 R=2.0 v=5.0 a=50.0 b=25.0 pf=0.2 run$(i)_p"
 
t = joinpath(runfolder, filename)
f= main_folder*"time_analysis.csv"
f1= t*".csv"
f2= main_folder*"time_analysis_run.png"
df= CSV.read(f1,DataFrame)

start_frame= 1  
end_frame= 100000
Neq1= df[start_frame:end_frame,:NeqL]
Neq2= df[start_frame:end_frame,:NeqR]
Npole1= df[start_frame:end_frame,:NpoleU]
Npole2= df[start_frame:end_frame,:NpoleD]

# total time spent in the boundary
total_time_poles = sum(Npole1) + sum(Npole2) 

tolal_time_equators = sum(Neq1) + sum(Neq2)

push!(data, (run=i, time_poles=total_time_poles, time_equators=tolal_time_equators))


CSV.write(f, data)
end
mainfolder= raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\2024\11.November\20241112-142941\R=2.0 v=5.0 a=50.0 b=25.0 pf=0.2\\"
f1= mainfolder*"time_analysis.csv"
df= CSV.read(f1,DataFrame)
plot(df.run,df.time_poles)
f2= mainfolder*"time_analysis_run.png"
# k= plot(df.run,df.time_poles,legend=false)
j= plot(df.run,(df.time_equators-df.time_poles)/100000,legend=false,marker=:circ, ms=2.25, mc=:red, ml=:thin)
m=plot(j,layout=(1,2),title="Poles Equators",xlabel="run", ylabel="T(eq)-T(poles)",legend=false,)
savefig(m,f2)
