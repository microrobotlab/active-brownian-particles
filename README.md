## PURPOSE: 
Simulate the dynamics of active brownian particles moving in 2D (in Julia). The particles are interacting via hard sphere correction. 
Code can simulate both open (periodic) and closed hard boundary condition.
When confined, particles are reflected from the boundary.

## AUTHORS: 
Jyoti Sharma <sup>1</sup>, Lapo Corti<sup>2</sup>, Stefano Palagi<sup>1</sup>

## AFFILIATION:
 1) The BioRobotics Institute, Scuola Superiore Sant’Anna , Pontedera (Pisa), Italy
 2) University of Pisa, Pisa, Italy.

## RELEASE DATE 
28-July-2023

## PREREQUISITES:  
Julia
## INSTRUCTIONS TO RUN THE CODE:
1) Open ABP output.jl file
2) Set destination folder path
3) Set the parameters
4) Run the code
5) After compliation, there will be a folder named "date/time" of run and inside subfolder named "simulation parameters" and inside it there will be multiple folders, named run1, run2---- runICS, where ICS is the  number of initial conditions scan for a given set of parameters
6) Output will be a gif, two csv files and one plot
7) Gif file shows the evolution of the particles
7) File named XXXX.csv has information of the position and orientation of all simulated particles
8) File named XXXX_p.csv has information of number of particles at equators or poles of the ellipse
9) Plot named XXXX.png shows the temporal evolution of packing fraction of particles at equators or poles of the ellipse


## REFERENCES
1) Numerical Simulations of Active Brownian Particles, *Agnese Callegari and Giovanni Volpe*, Book series (SOBIMA)

DOI: https://link.springer.com/chapter/10.1007/978-3-030-23370-9_7

2) Simulation of the active Brownian motion of a microswimmer, *Giorgio Volpe, Sylvain Gigan and Giovanni Volpe*, American Journal of Physics, 82, 659 (2014)

DOI: https://doi.org/10.1119/1.4870398