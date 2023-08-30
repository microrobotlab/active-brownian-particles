# for fitting the exponential curve in the average number of particles at equators and at poles.
# input will be average data generated from ABP average.jl 
# model funcion is

# m(t,p)=p1*exp(−p2*t)

using LsqFit, Plots


model(t, p) = p[1] * exp.(-p[2] * t)

tdata = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
ydata = model(tdata, [1.0 2.0]) + 0.01*randn(length(tdata))

p0 = [0.5, 0.5]

fit = curve_fit(model, tdata, ydata, p0)
param = fit.param

p=plot(tdata,ydata)
q= plot(tdata,fit)
display(q)