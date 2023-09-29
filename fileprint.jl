using CSV, FileIO, DataFrames, Plots, LaTeXStrings, Statistics 

mainfolder="C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-205011\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\"


filename= "average of FTT 100_data.csv"
path= mainfolder*filename

f= mainfolder*"poster_FFT.png"
df= CSV.read(path,DataFrame)


start_frame= 1
end_frame= 10000
#time= df[start_frame:end_frame,:t]
 t1= plot(df[start_frame:end_frame,:f], df[start_frame:end_frame,:real],legend=false)  
xlabel!("Frquency (HZ)", seriestype=:stem,xguidefont=font(20),xlimit=(0,0.5),ylimit=(0,200), xtickfont=font(11))
 plot!(ylabel="Power",yguidefont=font(20), ytickfont=font(11))
#  title!(" Equators ")

# t2= scatter(time./100.0, df[start_frame:end_frame,:p],legend=false) 
# xlabel!("Time (s)", xguidefont=font(14),ylimit=(15,30),xticks=(0:4000:10000), xtickfont=font(11))
# plot!(ylabel=L"\mathrm{N_{poles}}",xlimit=(0,10000),yguidefont=font(20), ytickfont=font(11))
# title!(" Poles ")

p= plot(t1)
display(p)
savefig(p,f)
