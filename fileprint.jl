using CSV, FileIO, DataFrames, Plots, LaTeXStrings, Statistics, PlotlyJS, FFTW

mainfolder = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07 Coding\2024\10.October\ellipse\20241017-105529\R=2.0 v=5.0 a=50.0 b=25.0 pf=0.1\run1"
filename= "\\20241017-105529 R=2.0 v=5.0 a=50.0 b=25.0 pf=0.1 run1_p.csv"

path= mainfolder*filename

f= mainfolder*"\\number of particles.png"
f1= mainfolder*"\\FFT_difference_equators_fs1000.png"
f2= mainfolder*"\\FFT_difference_poles_fs1000.png"
df= CSV.read(path,DataFrame)
t5=plot()
t6=plot()

start_frame= 1
end_frame= 100000
time= df[start_frame:end_frame,:t]
Neq1= df[start_frame:end_frame,:Neq1]
Neq2= df[start_frame:end_frame,:Neq2]
Npole1= df[start_frame:end_frame,:Npole1]
Npole2= df[start_frame:end_frame,:Npole2]
time1=time./100.0
yl= 30
#  title!(" Equators ")
t1= scatter(time./100.0, df[start_frame:end_frame,:Neq1],legend=false) 
xlabel!("Time (s)", xguidefont=font(14),ylimit=(0,30))#xticks=(0:4000:10000), xtickfont=font(11))
 plot!(ylabel=L"\mathrm{N_{eq(L)}}",xlimit=(0,100),yguidefont=font(20), ytickfont=font(11))
#  title!(" Poles ")
 t2= scatter(time./100.0, df[start_frame:end_frame,:Neq2],legend=false) 
xlabel!("Time (s)", xguidefont=font(14),ylimit=(0,30))#xticks=(0:4000:10000), xtickfont=font(11))
 plot!(ylabel=L"\mathrm{N_{eq(R)}}",xlimit=(0,100),yguidefont=font(20), ytickfont=font(11))
#  title!(" Poles ")
t3= scatter(time./100.0, df[start_frame:end_frame,:Npole1],legend=false) 
xlabel!("Time (s)", xguidefont=font(14),ylimit=(0,30))#xticks=(0:4000:10000), xtickfont=font(11))
 plot!(ylabel=L"\mathrm{N_{Pole(u)}}",xlimit=(0,100),yguidefont=font(20), ytickfont=font(11))

 t4= scatter(time./100.0, df[start_frame:end_frame,:Npole2],legend=false) 
xlabel!("Time (s)", xguidefont=font(14),ylimit=(0,30))#xticks=(0:4000:10000), xtickfont=font(11))
 plot!(ylabel=L"\mathrm{N_{Pole(d)}}",xlimit=(0,100),yguidefont=font(20), ytickfont=font(11))

 scatter!(t5,time1,Neq1.-Neq2, ylimit=(-yl,yl),mode="markers",markersize=0.5,seriescolor=:red,legend=false,ylabel=L"\mathrm{N_{Eq(L)}}-\mathrm{N_{Eq(R)}}") 

   scatter!(t6,time1,Npole1.-Npole2, ylimit=(-yl,yl),mode="markers",markersize=0.5,markercolor=:blue,legend=false,ylabel=L"\mathrm{N_{Pole(U)}}-\mathrm{N_{Pole(D)}}") 
p= plot(t1,t2,t3,t4, layout=(2,2))

q= plot(t5,t6,layout=(2,1))

display(p)
savefig(q,f)

ydata1 = Neq1.-Neq2
ydata2 = Npole1.-Npole2

fs=1000
freq= fftshift(fft(ydata1))
freqs = fftshift(fftfreq(length(ydata1), fs))

k= plot(freqs,real.(freq), xlimit=(0,500), ylimit=(0.02,20000),seriestype=:stem, xlabel="Frequency(Hz)", ylabel="Power",legend=false)

display(k)
savefig(k,f1)

# t1= plot(df[start_frame:end_frame,:f], df[start_frame:end_frame,:real],legend=false)  
# xlabel!("Frquency (HZ)", seriestype=:stem,xguidefont=font(20),xlimit=(0,0.5),ylimit=(0,200), xtickfont=font(11))
#  plot!(ylabel="Power",yguidefont=font(20), ytickfont=font(11))#