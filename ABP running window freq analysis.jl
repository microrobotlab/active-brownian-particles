#FFT analysis of the number of particles in the system using welch_pgram.jl method for digital signal processing in julia
using CSV, FileIO, DataFrames, Plots,DSP

function FFT_analysis_window(mainfolder,dt,resample,task)
    
f1= mainfolder*"_parameter_p.csv"

f= mainfolder*"\\number of particles.png"

df= CSV.read(f1,DataFrame)

start_frame= 1000
end_frame= 20000
time_step= df[:,:t]
Neq1= df[:,:NeqL]
Neq2= df[:,:NeqR]
Npole1= df[:,:NpoleU]
Npole2= df[:,:NpoleD]

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



time_series = ydata[start_frame:end_frame] .-mean(ydata[start_frame:end_frame])

# = 20000
fs=Int(1/(dt*resample)) #sampling rate = sampling frequency = 1/dt if at every time step data is printed
#assert length(time_series) == N


# Parameters for Welch's method
nperseg = 8192  # Number of data points in each segment. You can increase this value to get more frequency resolution
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
idx = findall(f -> 0 <= f <= 20.01, freqs_mHz)
plot(freqs_mHz[idx], psd[idx], xlabel="Frequency (mHz)", ylabel="Power", label="PSD", title="Power Spectral Density")
vline!([peak_freq_mHz], label="Peak at $peak_freq_mHz mHz", linestyle=:dash)

return peak_freq_mHz, peak_period, freqs_mHz, psd

end