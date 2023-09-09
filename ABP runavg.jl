using DataFrames, CSV, Statistics,  FFTW, Plots

path = "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-205011\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\"
f= path*"average100Copy.csv"
f1= path*"running_avg_FFT1.png"
f2=path*"running_avg_curve1.png"
df= CSV.read(f,DataFrame) # e is equators particles , p is poles particles


window_size= 1

function running_avg!(df, window_size)
#start_index= max(1, i -window_size +1) , takes care of index not getting negative max(1,-5) =1 so even if i -window_size + 1 is neagtive for intial data, net index is positive
    df[!,:running_avg]= [mean(df[max(1,i-window_size +1):i, :e]) for i in 1:nrow(df)] #it will add eunning avg in df
end
 time= df[!,:t]

running_avg!(df,window_size)

start_frame= 1000
end_frame= 10000
fs=1.0
freq= fftshift(fft(df[start_frame:end_frame,:running_avg]))
freqs = fftshift(fftfreq(length(df[start_frame:end_frame,:running_avg]), fs))

k= plot(freqs,abs.(freq), xlimit=(-0.01,0.01), ylimit=(0.02,1000),seriestype=:stem)

display(k)

#y=plot(df_avg[start_frame:end_frame,:t]./100.0,df_avg[start_frame:end_frame,:running_avg], seriestype= :scatter, xlimit=(start_frame,end_frame))

#savefig(k,f1)

#savefig(y,f2)