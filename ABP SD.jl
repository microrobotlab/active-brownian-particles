
# PURPOSE: Calculate the mean and standard deviation of the packing fraction
# INPUT: Data file genereated by ABP analysis.jl

using Plots, LaTeXStrings, Statistics, CSV, DataFrames

gr()
pathf = "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\July\\parameter scan\\20230729-181916\\R=2.0 v=1.0 a=50.0 b=25.0 pf=0.1\\run2\\"
function stat_analysis2(a,b,R,pathf)
  f1= pathf*"20230729-181916 R=2.0 v=1.0 a=50.0 b=25.0 pf=0.1 run2_p.csv"
  df= CSV.read(f1, DataFrame)
  neq= df[!,:p1]
  npole= df[!,:p2]
  eθ = atan(b/a)   # angle at which area will be same
pf_factor = R*R
Aeq= a*b*(atan(a*tan(eθ)/b))  # equator area

Ap= a*b*(atan(b/(tan(eθ)*a)))
  meq= mean(neq)
  sdeq= stdm(neq,meq)    # standard deviation of equator
  mpole= mean(npole)
  sdpole= stdm(npole,mpole)
  mpfe = meq*(0.5*π/Aeq)*pf_factor 
  mpfp = mpole*(0.5*π/Ap)*pf_factor 
  println("i am in ABP SD")
  sdpfe= sdeq*(0.5*π/Aeq)*pf_factor  # standard deviation of pf at equators
  sdpfp= sdpole*(0.5*π/Ap)*pf_factor 

#plot( [1:2], [mpfe,mpfp], yerror = [sdpfe,sdpfp], xtick=0:3, ytick= 0:0.2:0.5, xlimit= (0,3), ylimit= (0,0.5), grid= false, legend = false, ytickfont=font(17), xtickfont=font(17), linewidth=2, linecolor= "black")

#savefig(pathf*"_error.png")
return mpfe, sdpfe,mpfp, sdpfp
#return meq, sdeq, mpole,sdpole

end

L = 100.0 	# μm box length
R = 2.0		# μm particle radius

a=L/2
b=L/4
out= stat_analysis2(a,b,R,pathf)