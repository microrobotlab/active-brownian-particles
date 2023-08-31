# for fitting the exponential curve in the average number of particles at equators and at poles.
# input will be average data generated from ABP average.jl 
# model funcion is

# m(t,p)=p1*exp(−p2*t)

using LsqFit, Plots, DataFrames, CSV, CurveFit

path = "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\08.Aug\\ellipse\\20230824-205011\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\"
f= path*"average100.csv"
df= CSV.read(f,DataFrame)

# define model
y(t) = a * exp(b * t)
model(t, p) = p[1] * exp.(p[2] * t)   # exponential fitting

# linear fitting

#model(t, p) = p[1] + (p[2] * t)  

tdata = 0:0.1:2.0
ydata = 3.0 * exp.(-2.0 * tdata) + 0.1*randn(length(tdata))


#tdata = df[100:1000,:t]./100.0
#neq= df[100:1000,:e]
#np= df[100:1000,:p]
#tdata= log10.(tdata)
#ydata = model(tdata, [1.0 2.0]) + 0.01*randn(length(tdata))
#ydata= log10.(ydata)
#ydata= neq
#ydata=log10.(neq)
p0 = [2.0, 1.5]

fit = curve_fit(model, tdata, ydata, p0)
param = fit.param
yfit= model(tdata, param)

p=plot(tdata,ydata, seriestype=:scatter, label="Data")
q= plot!(tdata,yfit,label="Fitted $((param[1]))exp($(param[2])t)")
display(q)