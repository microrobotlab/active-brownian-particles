parent_folder  = raw"C:\Users\j.sharma\OneDrive - Scuola Superiore Sant'Anna\P07Coding\workstationMRL\2025\03.March\20250325-165047\R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1\run1\\"
filename = "20250325-165047 R=1.5 v=5.0 a=50.0 b=12.5 pf=0.1 run1" # dont put \\ after the filename
input = joinpath(parent_folder, filename * ".csv")
output = joinpath(parent_folder, filename * "_angle_distance.csv")
#output = raw"C:\Users\j.sharma\...\20250325-165047_angle_distance.csv"

export_angle_distance(input, output, 50.0, 12.5; plot_step=500000)