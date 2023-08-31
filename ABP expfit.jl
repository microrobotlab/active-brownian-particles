# for fitting the exponential curve in the average number of particles at equators and at poles.
# input will be average data generated from ABP average.jl 
# model funcion is

# m(t,p)=p1*exp(−p2*t)

using LsqFit, Plots, DataFrames, CSV, CurveFit, FFTW

path = "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-205011\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\"
f= path*"average100.csv"
f1= path*"fitting.png"
df= CSV.read(f,DataFrame)

# define model
#y(t) = a * exp(b * t)
#model(y, p) = p[1] * exp.(p[2] * y) #.+ p[3]  # exponential fitting

model(y, p) = p[1] .+ (p[2] * y)         # linear fitting
# linear fitting

#model(t, p) = p[1] + (p[2] * t)  

tdata = df[9250:10000,:t]./100.0
neq= df[9250:10000,:e]
np= df[9250:10000,:p]
ydata= neq
jdata = 0.0026 * exp.(0.471 * neq) ##+ 0.1*randn(length(neq))
#tdata= log10.(tdata)
#ydata = model(tdata, [1.0 2.0]) + 0.01*randn(length(tdata))
#ydata= log10.(ydata)
plot(tdata,ydata, seriestype=:scatter, label="Data")
#ydata=log10.(neq)
p0 = [25.0, 0.01]
fs=1.0

freq= fftshift(fft(neq))
freqs = fftshift(fftfreq(length(tdata), fs))

k= plot(freqs,abs.(freq), xlimit=(-0.02,0.02), seriestype=:scatter)



display(k)


#=
fit = curve_fit(model, ydata, tdata, p0)
param = fit.param
xfit= model(ydata, param)

p=plot(tdata,ydata, seriestype=:scatter, label="Data")
#q= plot!(xfit,ydata,label="Fitted $((param[1]))exp($(param[2])t)")
q= plot!(xfit,ydata,label="Fitted $(param[1])t+ $(param[1])")
display(q)
#savefig(q,f1)
=#