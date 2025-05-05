using DSP
using Plots
using CSV
using DataFrames

# Load the time series data for one run (adjust as needed)
# Example: Loading from a CSV file
mainfolder= "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07Coding\\workstationMRL\\2025\\04.April\\20250415-120846\\R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1\\" # destination folder path
mainfolder1= mainfolder*"run1\\"
filename="20250415-120846 R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1 run1_parameter_p" # dont put \\ after the filename
t= mainfolder1 * filename * ".csv"

df = CSV.read(t, DataFrame)  # Replace with your file path
start_frame= 1000
end_frame= 20000
time_step= df[:,:t]
Neq1= df[:,:NeqL]
Neq2= df[:,:NeqR]
Npole1= df[:,:NpoleU]
Npole2= df[:,:NpoleD]
time_series =(Neq1.+Neq2).-(Npole1.+Npole2)
time_series = time_series .-mean(time_series)


# Alternatively, if you have the data as a vector:
# time_series = ...  # Your vector of 20,000 samples

# Verify data length
N = 20000
fs = 1.0  # Sampling rate (Hz)
@assert length(time_series) == N
println("Time series length: $N samples")

# Parameters for Welch's method
nperseg = 4096  #you can increase this value to get more frequency resolution
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

println("Characteristic frequency for this run: $peak_freq_mHz mHz")
println("Corresponding period: $peak_period seconds")

# Plot the PSD (optional)
idx = findall(f -> 0 <= f <= 20.01, freqs_mHz)
plot(freqs_mHz[idx], psd[idx], xlabel="Frequency (mHz)", ylabel="Power", label="PSD", title="Power Spectral Density")
vline!([peak_freq_mHz], label="Peak at $peak_freq_mHz mHz", linestyle=:dash)