# to check FFt for constant functions
using CSV, FileIO, DataFrames, Plots, LaTeXStrings, Statistics, FFTW
t = []
v = []
for i in 1:1000
 push!(t,i)
 push!(v,randn())
end

df = DataFrame(t=t,v=v)
df.v
plot(df.t,df.v, seriestype=:line)
freq = fftshift(fft(df.v))
freqs = fftshift(fftfreq(length(df.v), 1/maximum(df.t)))
plot(freqs,real.(fft(df.v)), seriestype=:stem, title="FFT of v(t)")

t0 = 0              # Start time 
fs = 44100          # Sampling rate (Hz)
tmax = 0.1          # End time       

t = t0:1/fs:tmax;   
signal = sin.(2π * 60 .* t)

F = fftshift(fft(signal))
freqs = fftshift(fftfreq(length(t), fs))

# plots 
time_domain = plot(t, signal, title = "Signal", label='f',legend=:top)
freq_domain = plot(freqs, abs.(F), title = "Spectrum", xlim=(-100, +100), xticks=-100:20:100, label="abs.(F)",legend=:top) 
plot(time_domain, freq_domain, layout = (2,1))