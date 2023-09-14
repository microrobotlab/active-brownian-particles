# code can take files from multiple folders and can calculate their running average (select window_size) and FFT

using DataFrames, CSV, Statistics,  FFTW, Plots

for i in 1:1

path = "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-205011\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\run$i\\"
path1= "..\\"
filename= "20230824-205011 R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1 run$(i)_p"
pathf= path*filename
f= pathf*".csv"
f1= path*"running_avg_FFT_run$i.png"
f2=path*"FFT_data.csv"

df= CSV.read(f,DataFrame) # e is equators particles , p is poles particles


window_size= 10

function running_avg!(df, window_size)
#start_index= max(1, i -window_size +1) , takes care of index not getting negative max(1,-5) =1 so even if i -window_size + 1 is neagtive for intial data, net index is positive
    df[!,:running_avg]= [mean(df[max(1,i-window_size +1):i, :p1]) for i in 1:nrow(df)] #it will add eunning avg in df
end
 time= df[!,:t]

running_avg!(df,window_size)

start_frame= 1000
end_frame= 10000
fs=1.0
freq= fftshift(fft(df[start_frame:end_frame,:running_avg]))
freqs = fftshift(fftfreq(length(df[start_frame:end_frame,:running_avg]), fs))

k= plot(freqs,abs.(freq), xlimit=(-0.3,0.3), ylimit=(0.02,1000),seriestype=:stem)

display(k)
savefig(k,f1)

end
#y=plot(df_avg[start_frame:end_frame,:t]./100.0,df_avg[start_frame:end_frame,:running_avg], seriestype= :scatter, xlimit=(start_frame,end_frame))



#savefig(y,f2)