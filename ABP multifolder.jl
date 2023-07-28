# PURPOSE: CREATE MULTIPLE FOLDERS
# INPUT: PATH AND NUMBER OF FOLDERS
using FilePaths

function multipledir(path, ICS)
    rundirs = []
    for i=1:ICS
        rundir = joinpath(path, "run$i")
        if !isdir(rundir)
            mkdir(rundir)
            push!(rundirs, rundir) # Save the path of the created directory.
        else
            println("Directory already exists: ", rundir)
        end
    end
    if isempty(rundirs)
        error("All directories already exist, time for a cleanup!")
    end
    return rundirs # Return the list of created directories.
end
