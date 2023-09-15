# code can take files from multiple folders and can calculate their running average (select window_size) and FFT

using DataFrames, CSV, Statistics,  FFTW, Plots
#path1 = "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-205011\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\"

all_data = DataFrame( R=Float64[],v=Float64[],a=Float64[],b=Float64[],pf=Float64[],run=Int64[],timestep=Float64[],f=Float64[], FFT_real= Float64[], FFT_img= Float64[])

param_df = DataFrame(Parameter=String[], Value=String[])
 # creates an empaty data frame

for i in 1:10

path = "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-205011\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\run$i\\"

main_folder= joinpath(path, "..\\")  # takes back to the previous folder
filename= "20230824-205011 R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1 run$(i)_p"

pathf= path*filename
f=  pathf*".csv"
f1= path*"running_avg_FFT_run$i.png"
f2= main_folder*"FFT_data.csv"

df= CSV.read(f,DataFrame)            # e is equators particles , p is poles particles

################################################################### printing filename parameters#################################################################################
parts= split(filename, ' ') # splits the strings after every space

for part in parts
    if occursin("=",part)
        param_value = split(part, '=') # furter split the filename with = 
        if length(param_value) == 2
            param, value = param_value
            push!(param_df, [param, value])
        end
    end
   
end

########################################################################## Running average #############################################################################################################
window_size= 1

function running_avg!(df, window_size)
#max(1, i -window_size +1) , takes care of index not getting negative max(1,-5) =1 so even if i -window_size + 1 is neagtive for intial data, net index is positive
    df[!,:running_avg]= [mean(df[max(1,i-window_size +1):i, :p1]) for i in 1:nrow(df)] #it will add eunning avg in df
end
 time= df[!,:t]

running_avg!(df,window_size)
################################################################################ FFT of data ######################################################################################################

start_frame= 1000
end_frame= 10000
fs=1.0
freq= fftshift(fft(df[start_frame:end_frame,:running_avg]))
freqs = fftshift(fftfreq(length(df[start_frame:end_frame,:running_avg]), fs))

data= DataFrame( R=param_df[1,:Value],v= param_df[2,:Value],a= param_df[3,:Value],b= param_df[4,:Value],pf= param_df[5,:Value], run=i,timestep=df[start_frame:end_frame,:t],f=freqs, FFT_real= real.(freq), FFT_img= imag.(freq)) 
global all_data=vcat(all_data,data)

CSV.write(f2,all_data)
k= plot(freqs,real.(freq), xlimit=(-0.3,0.3), ylimit=(0.02,2000),seriestype=:stem)

display(k)
#savefig(k,f1)

end




