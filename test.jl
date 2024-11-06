# this code is for testing the index change in the julia and to save space in the code
# discussion with Stefano on 5-11-2024


function test(Nt)
    Nt1=Int(Nt/10)
    ABPE = Vector{Vector{Float64}}(undef, Nt1 + 2)
    ABPE[1] = [1.0]
      # Initialize current_value in the same format as ABPE[1]
      current_value = deepcopy(ABPE[1])  # This will be used for updates at each step
      abpe_index = 2  # Position in ABPE to store every 10th step
    for i = 1:Nt
        ttt= Int(ceil(i/10))
       
        current_value = current_value .+δt
        # println(kk)
    #   ABPE[i+1] = ABPE[i].+δt
    # #     # ABPE[1+1] = ABPE[1].+δt
      if mod(i,10) == 0
            println("Step $i")
           ABPE[abpe_index] = deepcopy(current_value)
           abpe_index += 1
    println("Step $i")
        
   
       
     end
     end
    return ABPE
end

Np = round(Int,packing_fraction*a*b/(R^2))  #Np is the number of particles inside the ellipse
#π
Nt = 1000# Nt is the number of steps 
δt = 1.0

test(Nt)


######################################################
