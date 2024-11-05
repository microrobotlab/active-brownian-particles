using CSV, FileIO, DataFrames, Plots, LaTeXStrings, Statistics, FFTW

mainfolder = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\2024\10.October\ellipse\20241022-185034\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.2\run10\\"
filename= "20241022-185034 R=2.0 v=10.0 a=50.0 b=25.0 pf=0.2 run10_p.csv"

path= mainfolder*filename

f= mainfolder*"\\number of particles.png"
# f1= mainfolder*"\\FFT_difference_equators_fs1000.png"
f2= mainfolder*"\\FFT_difference.png"
df= CSV.read(path,DataFrame)

start_frame= 1
end_frame= 1000000
dt=0.001
time_step= df[:,:t]
Neq1= df[:,:Neq1]
Neq2= df[:,:Neq2]
Npole1= df[:,:Npole1]
Npole2= df[:,:Npole2]
Waqt = time_step .* dt
########################################## FFT Analysis start #############################################################àà

# please choose the value of i, 1 for poles and 2 for equators
i=2
if i==1
    ydata = Npole1.-Npole2
    f1= mainfolder*"\\FFt_diff_poles.png"
elseif i==2
    ydata = Neq1.-Neq2
    f1= mainfolder*"\\FFT_diff_equators.png"
elseif i==3
    ydata = (Neq1.+Neq2).-(Npole1.+Npole2)
    f1= mainfolder*"\\FFT_diff_curvature.png"
   # f= mainfolder*"\\number of particles.png"
end
# tt= plot(Waqt,ydata,legend=false)
# display(tt) 
# savefig(tt,f)
ydata_corr= ydata[start_frame:end_frame].-mean(ydata[start_frame:end_frame])
fs=Int(1/dt) #sampling rate = sampling frequency = 1/dt if at every time step data is printed
freq= fftshift(fft(ydata_corr))
freqs = fftshift(fftfreq(length(ydata_corr), fs))
peaks, peak_values = findmax(real.(freq)) 
# Find all peaks in the real part of the frequency data
real_freq = (real.(freq))
all_peaks = findall(i -> (i > 1 && i < length(real_freq) && real_freq[i] > real_freq[i-1] && real_freq[i] > real_freq[i+1]), 1:length(real_freq))

# Sort peaks by their heights and select the top three
sorted_peaks = sort(all_peaks, by = x -> real_freq[x], rev = true)
top_peaks = sorted_peaks[1:min(10, length(sorted_peaks))]  # get up to the top 3 peaks
 top_frequencies = freqs[top_peaks]
 top_values = real_freq[top_peaks]

# Plot the frequency data with peaks highlighted
k=plot(freqs, real_freq, seriestype=:stem, xlim=(0.00, 0.02), ylim=(0.02, 1000000), xlabel="Frequency (Hz)", ylabel="Power", legend=false)

# Add the top peaks to the plot with annotations
scatter!(top_frequencies, top_values, color=:red, label="Peaks", markersize=3)
for (i, peak) in enumerate(top_frequencies)
    annotate!(peak, i,text("$(round(peak, digits=3))", :center, 8))
end

# Show or save the plot
title!("Peaks in Frequency= $(round(top_frequencies[1], digits=5)) $(round(top_frequencies[2], digits=3)) $(round(top_frequencies[3], digits=3))")
display(k)
savefig(k, f1)
Int(round(peaks))
# Plot with peaks highlighted
# k = plot!(freqs, real.(freq), xlimit=(0, 0.5), ylimit=(0.02, 10000), seriestype=:stem, xlabel="Frequency(Hz)", ylabel="Power", legend=false)
#scatter!(freqs[peak_values], peaks, color=:red, label="Peaks")

# kkk=peaks
# scatter!([abs(freqs[peak_values])], [peaks], color=:red, title="Most frequency= $(abs(freqs[peak_values]))", label="Peaks", markersize=5)
# # Save the plot
# savefig(k, f1)

# k= plot(freqs,real.(freq), xlimit=(-0.5,0.5), ylimit=(0.02,10000),seriestype=:stem, xlabel="Frequency(Hz)", ylabel="Power",legend=false)

# display(k)
# savefig(k,f1)
########################################## Plots start #############################################################
# yl= 30
# #  title!(" Equators ")
# t1= scatter(time./100.0, df[start_frame:end_frame,:Neq1],legend=false) 
# xlabel!("Time (s)", xguidefont=font(14),ylimit=(0,30))#xticks=(0:4000:10000), xtickfont=font(11))
#  plot!(ylabel=L"\mathrm{N_{eq(L)}}",xlimit=(0,100),yguidefont=font(20), ytickfont=font(11))
# #  title!(" Poles ")
#  t2= scatter(time./100.0, df[start_frame:end_frame,:Neq2],legend=false) 
# xlabel!("Time (s)", xguidefont=font(14),ylimit=(0,30))#xticks=(0:4000:10000), xtickfont=font(11))
#  plot!(ylabel=L"\mathrm{N_{eq(R)}}",xlimit=(0,100),yguidefont=font(20), ytickfont=font(11))
# #  title!(" Poles ")
# t3= scatter(time./100.0, df[start_frame:end_frame,:Npole1],legend=false) 
# xlabel!("Time (s)", xguidefont=font(14),ylimit=(0,30))#xticks=(0:4000:10000), xtickfont=font(11))
#  plot!(ylabel=L"\mathrm{N_{Pole(u)}}",xlimit=(0,100),yguidefont=font(20), ytickfont=font(11))

#  t4= scatter(time./100.0, df[start_frame:end_frame,:Npole2],legend=false) 
# xlabel!("Time (s)", xguidefont=font(14),ylimit=(0,30))#xticks=(0:4000:10000), xtickfont=font(11))
#  plot!(ylabel=L"\mathrm{N_{Pole(d)}}",xlimit=(0,100),yguidefont=font(20), ytickfont=font(11))

#  scatter!(t5,time1,Neq1.-Neq2, ylimit=(-yl,yl),mode="markers",markersize=0.5,seriescolor=:red,legend=false,ylabel=L"\mathrm{N_{Eq(L)}}-\mathrm{N_{Eq(R)}}") 

#    scatter!(t6,time1,Npole1.-Npole2, ylimit=(-yl,yl),mode="markers",markersize=0.5,markercolor=:blue,legend=false,ylabel=L"\mathrm{N_{Pole(U)}}-\mathrm{N_{Pole(D)}}") 
# p= plot(t1,t2,t3,t4, layout=(2,2))

# q= plot(t5,t6,layout=(2,1))

# display(p)
# savefig(q,f)

# t1= plot(df[start_frame:end_frame,:f], df[start_frame:end_frame,:real],legend=false)  
# xlabel!("Frquency (HZ)", seriestype=:stem,xguidefont=font(20),xlimit=(0,0.5),ylimit=(0,200), xtickfont=font(11))
#  plot!(ylabel="Power",yguidefont=font(20), ytickfont=font(11))#