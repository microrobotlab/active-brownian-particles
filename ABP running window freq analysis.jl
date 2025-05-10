#FFT analysis of the number of particles in the system using welch_pgram.jl method for digital signal processing in julia
using CSV, FileIO, DataFrames, Plots,DSP,Polynomials

function FFT_analysis_window(mainfolder,run_folder,dt,resample,task)
    
# f1= mainfolder*"_parameter_p.csv"
f1= mainfolder*"data.csv"
f= run_folder*"\\FFT_analysis.png"

df= CSV.read(f1,DataFrame)

start_frame= 1
end_frame= 4300
# time_step= df[:,:t]
# Neq1= df[:,:NeqL]
# Neq2= df[:,:NeqR]
# Npole1= df[:,:NpoleU]
# Npole2= df[:,:NpoleD]
ydata= df[:,:pdiff]

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

    t_vec = 1:length(ydata)  # Time vector for the data points
    p = fit(t_vec, ydata, 1)  # Linear fit
    trend = p.(t_vec)
    time_series_detrended = ydata .- trend
    time_series_centered = time_series_detrended .- mean(time_series_detrended)


time_series = time_series_centered[start_frame:end_frame] .-mean(ydata[start_frame:end_frame])

# = 20000
# fs=Int(1/(dt*resample)) #sampling rate = sampling frequency = 1/dt if at every time step data is printed
#assert length(time_series) == N

fs=12
# Parameters for Welch's method
nperseg = 1024 # Number of data points in each segment. You can increase this value to get more frequency resolution
noverlap = nperseg ÷ 2

# Compute PSD
p = welch_pgram(time_series, nperseg, noverlap; fs=fs)

# Extract frequencies and PSD
 freqs = p.freq
 psd = p.power
freqs_mHz = freqs * 1000

# Ignore the zero frequency because it corresponds to the not changing number of particles 
freqs_mHz_no_zero = freqs_mHz[2:end]
psd_no_zero = psd[2:end]

# Find the peak frequency excluding zero frequency
peak_idx_no_zero = argmax(psd_no_zero)
peak_freq_mHz = freqs_mHz_no_zero[peak_idx_no_zero]
peak_period = 1 / (peak_freq_mHz / 1000)
# Plot the PSD (optional)
idx = findall(f -> 0.01 <= f <= 40.01, freqs_mHz)# sort of high pass filter
pt_linear= plot(freqs_mHz[idx], psd[idx], xlabel="Frequency (mHz)", ylabel="Power", label="PSD", title="Power Spectral Density")
vline!([peak_freq_mHz], label="Peak at $peak_freq_mHz mHz", linestyle=:dash)
savefig(pt_linear,joinpath(run_folder, "FFT_analysis.png"))
# pt_log= plot(freqs_mHz[idx], psd[idx], xscale=:log10,yscale=:log10,xlabel="log(Frequency)(mHz)", ylabel="Power", label="PSD", title="Power Spectral Density")
# vline!([peak_freq_mHz], label="Peak at $peak_freq_mHz mHz", linestyle=:dash)
# savefig(pt_log,joinpath(run_folder, "FFT_analysis_log.png"))
return peak_freq_mHz, peak_period, freqs_mHz, psd

end