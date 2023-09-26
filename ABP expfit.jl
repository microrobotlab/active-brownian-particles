# for fitting the exponential curve in the average number of particles at equators and at poles.
# input will be average data generated from ABP average.jl 

using LsqFit, Plots, DataFrames, CSV, FFTW

path = "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-205011\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\"
f= path*"average100.csv"
f1= path*"fitting_best_equators.png"
df= CSV.read(f,DataFrame)

# define model
#y(t) = a * exp(b * t)
#odel(y, p) = p[1] * exp.(p[2] * y) #.+ p[3]  # exponential fitting

model(t, p) = p[1] .+ p[2] * (1 .- exp.(-t/p[3]))

#model(y, p) = p[1] .+ (p[2] * y)         # linear fitting
# linear fitting

#model(t, p) = p[1] + (p[2] * t)  

tdata = df[1:10000,:t]./100.0
neq= df[1:10000,:e]
np= df[1:10000,:p]
ydata= neq

p0 = [23.0, 3.0, 100.0]
fs=1


fit = curve_fit(model, tdata, ydata, p0)
param = fit.param
yfit= model(tdata, param)

p=plot(tdata,ydata, seriestype=:scatter, label="Data", xlabel= "Time (s)", ylabel="N_eq")
q= plot!(tdata,yfit,label="Fitted ($(param[1])+ $(param[2]) * [1- exp(-t/$(param[3])])", title= "p0 = $p0")

display(q)
savefig(q,f1)
#######################################################################################

#FFT
#=
freq= fftshift(fft(neq))
freqs = fftshift(fftfreq(length(tdata), fs))

k= plot(freqs,abs.(freq), xlimit=(-0.02,0.02), seriestype=:scatter)



display(k)


=#