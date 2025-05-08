exp_folder="C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P10 Microfabrication\\Experiments\\2024\\12.December\\05\\exp1\\VID001\\data_ellipse21_V15\\"
exp_file= "data.csv"
exp_file= joinpath(exp_folder,exp_file)
df= CSV.read(exp_file,DataFrame)
ydata= df[:,:pdiff]
ydata_corr= ydata.-mean(ydata)
freq=(fft(ydata_corr))
freqs = (fftfreq(length(ydata_corr), 12))
real_freq = (abs.(freq))
findmaxima(real_freq)
indices, heights = findmaxima(real_freq)
k =plot(freqs*1000, real_freq, seriestype=:stem, xlim=(0, 100.0), ylim=(0.02, 500), xlabel="Frequency (mHz)", ylabel="Power", legend=false)
display(k)
savefig(k,joinpath(exp_folder,"FFT_analysis1.png"))
#FFT_analysis_window(exp_folder,exp_folder,δt,resample,6)