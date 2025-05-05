using CSV, FileIO, DataFrames, Plots, LaTeXStrings, Statistics, FFTW, Peaks

function FFT_analysis(mainfolder,dt,resample,task)
    
f1= mainfolder*"_parameter_p.csv"

f= mainfolder*"\\number of particles.png"
# f1= mainfolder*"\\FFT_difference_equators_fs1000.png"
f2= mainfolder*"\\FFT_difference.png"
f3=mainfolder*"_FFT.csv"
df= CSV.read(f1,DataFrame)

start_frame= 1000
end_frame= 20000
time_step= df[:,:t]
Neq1= df[:,:NeqL]
Neq2= df[:,:NeqR]
Npole1= df[:,:NpoleU]
Npole2= df[:,:NpoleD]
Waqt = time_step .* dt
########################################## FFT Analysis start #############################################################àà

# please choose the value of i, 1 for poles and 2 for equators
i=task
if i==1
    ydata = Npole1.-Npole2
    f1= mainfolder*"_perimeter_FFt_diff_poles.png"
elseif i==2
    ydata = Neq1.-Neq2
    f1= mainfolder*"_perimeter_FFT_diff_equators.png"
elseif i==3
    ydata = (Neq1.+Neq2).-(Npole1.+Npole2)
    f1= mainfolder*"_perimeter_FFT_diff_curvature.png"
   # f= mainfolder*"\\number of particles.png"
end

ydata_corr= ydata[start_frame:end_frame].-mean(ydata[start_frame:end_frame])
@show fs=Int(1/(dt*resample)) #sampling rate = sampling frequency = 1/dt if at every time step data is printed
#freq= fftshift(fft(ydata_corr))
freq=(fft(ydata_corr))
freqs = (fftfreq(length(ydata_corr), fs))
real_freq = (abs.(freq))
indices, heights = findmaxima(real_freq)
h_ranks = sortperm(heights, rev=true)
top_indices = indices[h_ranks[1:4]]
top_frequencies = freqs[top_indices]
#=
peaks, peak_values = findmax(abs.(freq)) 
# Find all peaks in the real part of the frequency data
real_freq = (abs.(freq))
all_peaks = findall(i -> (i > 1 && i < length(real_freq) && real_freq[i] > real_freq[i-1] && real_freq[i] > real_freq[i+1]), 1:length(real_freq))

# Sort peaks by their heights and select the top three
sorted_peaks = sort(all_peaks, by = x -> real_freq[x], rev = true)
top_peaks = sorted_peaks[1:min(4, length(sorted_peaks))]  # get up to the top 3 peaks
top_frequencies = freqs[top_peaks]
top_values = real_freq[top_peaks]
=#

# Plot the frequency data with peaks highlighted
k=plot(freqs*1000, real_freq, seriestype=:stem, xlim=(0, 2.0), ylim=(0.02, 55000), xlabel="Frequency (mHz)", ylabel="Power", legend=false)
# save the data in a csv file
freq_data= DataFrame(f=freqs, real= real_freq)
CSV.write(f3,freq_data)
#=
# Add the top peaks to the plot with annotations
scatter!(top_frequencies, top_values, color=:red, label="Peaks", markersize=3)
for (i, peak) in enumerate(top_frequencies)
    annotate!(peak, i,text("$(round(peak, digits=3))", :center, 8))
end
=#
# Show or save the plot
title!("Peaks in Freq(mHZ)= $(round(top_frequencies[1]*1000, digits=3)) $(round(top_frequencies[2]*1000, digits=3)) $(round(top_frequencies[3]*1000, digits=3))")
display(k)
savefig(k, f1)
#Int(round(peaks))
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

end