# Here the analysis for multiple files is done, the data is read from multiple folders 
# Analysis is done on each file and 
# One can alsso plot average of all data from multiple files
# INPUT: mainfolder path of the main folder which has all the run folders
# OUTPUT: Analysis of each file 

function multianalysis(mainfolder,ICS)
  # Initialize storage for peak frequencies and PSDs
  peak_frequencies = Float64[]
  psds = []
    for i in 1:ICS # number of runs
      println("analysis for run $i")
    mainfolder1= mainfolder*"run$i\\"
    filename="20250415-120846 R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1 run$i" # dont put \\ after the filename
    t= mainfolder1 * filename 
   
  #f1000= joinpath(mainfolder1,filename * ".csv") 
 
 # #    f1= "D:\\j.sharma\\P07\\workstationMRL\\20241104-121620\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.2\\run1\\20241104-121620 R=2.0 v=10.0 a=50.0 b=25.0 pf=0.2 run1_p.csv\\"
# df= CSV.read(f1000, DataFrame) 
   #inside_Np=stat_analysis_perimeter(a,b,R,t,δt,2) # 0 for pole, equator, 1 for only right left, 2 for entire
 #FFT_analysis(t,δt,resample,3) # FFT analysis of the data, number at the end is for different type of analyis

 ########################## FFT analysis of the data using Welch's method ##########################
 peak_freq_mHz, peak_period, freqs_mHz, psd = FFT_analysis_window(t,δt,resample,3)

 push!(peak_frequencies, peak_freq_mHz)
 push!(psds, psd)

end
   # Save peak frequencies to a CSV file
   output_file = joinpath(mainfolder, "peak_frequencies.csv")
   peak_df = DataFrame(Run = 1:ICS, Peak_Frequency_mHz = peak_frequencies)
   CSV.write(output_file, peak_df)
   plot(1:ICS, peak_frequencies, title="Peak Frequencies for $ICS runs", xlabel="Run", ylabel="Frequency (mHz)", legend=false)
   savefig(mainfolder*"peak_frequencies.png")
  
    
end