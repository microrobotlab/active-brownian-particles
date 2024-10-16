# for fitting the exponential curve in the average number of particles at equators and at poles.
# input can be any single file or average file
# FFT file can be saved and fitting could be done

using LsqFit, Plots, DataFrames, CSV, FFTW

path = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07 Coding\2024\10.October\ellipse\20241015-105629\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\run1"
filename= "\\20241015-105629 R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1 run1_p.csv"
f= path*filename
f1= path*"fitting_best_equators.png"
f2= path*"\\FFT_difference_equators.png"
df= CSV.read(f,DataFrame)

# define model
#y(t) = a * exp(b * t)
#odel(y, p) = p[1] * exp.(p[2] * y) #.+ p[3]  # exponential fitting

model(t, p) = p[1] .+ p[2] * (1 .- exp.(-t/p[3]))

#model(y, p) = p[1] .+ (p[2] * y)         # linear fitting
# linear fitting

#model(t, p) = p[1] + (p[2] * t)  

start_frame= 1
end_frame= 1000

tdata = df[start_frame:end_frame,:t]./100.0
Neq1= df[start_frame:end_frame,:Neq1]
Neq2= df[start_frame:end_frame,:Neq2]
Npole1= df[start_frame:end_frame,:Npole1]
Npole2= df[start_frame:end_frame,:Npole2]
ydata1 = Neq1.-Neq2
ydata2 = Npole1.-Npole2

#=
p0 = [23.0, 3.0, 100.0]
fit = curve_fit(model, tdata, ydata, p0)
param = fit.param
yfit= model(tdata, param)

p=plot(tdata,ydata, seriestype=:scatter, label="Data", xlabel= "Time (s)", ylabel="N_eq")
q= plot!(tdata,yfit,label="Fitted ($(param[1])+ $(param[2]) * [1- exp(-t/$(param[3])])", title= "p0 = $p0")

display(q)
savefig(q,f1)
=#
#######################################################################################FFT##############################################################


fs=1
freq= fftshift(fft(ydata1))
freqs = fftshift(fftfreq(length(ydata1), fs))

k= plot(freqs,real.(freq), xlimit=(0,0.6), ylimit=(0.02,200),seriestype=:stem, xlabel="Frequency(Hz)", ylabel="Power",legend=false)

display(k)
savefig(k,f2)


