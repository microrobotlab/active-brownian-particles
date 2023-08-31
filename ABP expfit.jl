# for fitting the exponential curve in the average number of particles at equators and at poles.
# input will be average data generated from ABP average.jl 
# model funcion is

# m(t,p)=p1*exp(−p2*t)

using LsqFit, Plots, DataFrames, CSV

path = "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-205011\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\"
f= path*"average100.csv"
df= CSV.read(f,DataFrame)
model(t, p) = p[1] * exp.(-p[2] * t)

tdata = df[!,:t]./100.0
neq= df[!,:e]
np= df[!,:p]
#ydata = model(tdata, [1.0 2.0]) + 0.01*randn(length(tdata))
ydata=neq
p0 = [0.5, 10.0]

fit = curve_fit(model, tdata, ydata, p0)
param = fit.param
yfit= model(tdata, param)



p=plot(tdata,ydata, seriestype=:scatter, label="Data")
q= plot!(tdata,yfit,label="Fitted $((param[1]))exp(-$(param[2])t)")
display(q)