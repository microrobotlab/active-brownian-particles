# for fitting the exponential curve in the average number of particles at equators and at poles.
# input will be average data generated from ABP average.jl 
# model funcion is

# m(t,p)=p1*exp(−p2*t)
# testing the terminal
# writing commit from terminal
# checking pushing and committing together 

using LsqFit, Plots, DataFrames, CSV, FFTW

path = "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-205011\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\"
f= path*"average100.csv"
f1= path*"fitting_best.png"
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
ydata= np

p0 = [23.0, 3.0, 100.0]
fs=1


fit = curve_fit(model, tdata, ydata, p0)
param = fit.param
yfit= model(tdata, param)

p=plot(tdata,ydata, seriestype=:scatter, label="Data")
q= plot!(tdata,yfit,label="Fitted $((param[1]))exp($(param[2])t)", title= "p0 = $p0")
# k= plot!(jdata,ydata)#label="Fitted $(param[1])t+ $(param[1])")
display(q)
#savefig(k,f1)

#jdata = 0.26 * exp.(0.3 * neq) ##+ 0.1*randn(length(neq))
#tdata= log10.(tdata)
#ydata = model(tdata, [1.0 2.0]) + 0.01*randn(length(tdata))
#ydata= log10.(ydata)
#plot(tdata,ydata, seriestype=:scatter, label="Data")
#ydata=log10.(neq)
#######################################################################################

#FFT
#=
freq= fftshift(fft(neq))
freqs = fftshift(fftfreq(length(tdata), fs))

k= plot(freqs,abs.(freq), xlimit=(-0.02,0.02), seriestype=:scatter)



display(k)


=#