# Here the analysis for multiple files is done, the data is read from multiple folders 
# Analysis is done on each file and 
# One can alsso plot average of all data from multiple files
# INPUT: mainfolder path of the main folder which has all the run folders
# OUTPUT: Analysis of each file 

function multianalysis(mainfolder)

    for i in 1:10
    mainfolder1= mainfolder*"run$i\\"
    filename="20250415-120846 R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1 run$i" # dont put \\ after the filename
    t= mainfolder1 * filename 
   
  #f1000= joinpath(mainfolder1,filename * ".csv") 
 
 # #    f1= "D:\\j.sharma\\P07\\workstationMRL\\20241104-121620\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.2\\run1\\20241104-121620 R=2.0 v=10.0 a=50.0 b=25.0 pf=0.2 run1_p.csv\\"
# df= CSV.read(f1000, DataFrame) 
  #inside_Np=stat_analysis_perimeter(a,b,R,t,δt,0) # 0 for pole, equator, 1 for only right left, 2 for entire
  FFT_analysis(t,δt,resample,3) # FFT analysis of the data, number at the end is for different type of analyis
    end
    
end