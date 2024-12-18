using Statistics

samp_std(x::Vector{Float64}) = std(x, corrected = false)*sqrt(length(x)-1)
susceptibility(x::Vector{Float64}) = length(x)*(mean(y->y^2, x) - (mean(x))^2)

function blocking(x::Array, blocksize::Int)
    length = size(x,1) ÷ blocksize 
    x_cut = x[1:length*blocksize]
    x_b = reshape(x_cut, (blocksize, length))

    return mean(x_b, dims = 1), length
end

function jackknife(x::Array, blocksize::Int)
    x_blocked, length = blocking(x, blocksize)
    jk(a) = let s = sum(a); [s-v for v in a] end
    x_j = vec(jk(x_blocked)/(length-1))
    return x_j
end